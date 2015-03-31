import scala.sys.process._
import scala.collection.JavaConversions._
import scala.collection.immutable.StringOps

import org.broadinstitute.sting.queue.QScript
import org.broadinstitute.sting.queue.util.QScriptUtils
import org.broadinstitute.sting.queue.extensions.gatk._
import org.broadinstitute.sting.queue.extensions.picard._
import org.broadinstitute.sting.queue.function.JavaCommandLineFunction

import net.sf.samtools.SAMFileReader
import net.sf.samtools.SAMFileHeader.SortOrder

//import org.apache.commons.io.FilenameUtils
//import java.io.File

class RefineBams extends QScript {
  qscript =>

  @Input(doc="Input BAM file(s)", fullName="input", shortName="i", required=true)
  var inputBams: List[File] = _

  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=true)
  var reference: File = _

  @Input(doc="Output directory", fullName="out_dir", shortName="O", required=true)
  var outPath: File = _

  @Input(doc="The path to the binary of bwa", fullName="path_to_bwa", shortName="bwa", required=true)
  var bwaPath: File = _

  @Input(doc="The path to the binary of bwa", fullName="path_to_samtools", shortName="samtools", required=true)
  var samtoolsPath: File = _

  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=true)
  var dbSNP: Seq[File] = Seq()

  @Input(doc="extra VCF files to use as reference indels for Indel Realignment", fullName="extra_indels", shortName="indels", required=true)
  var indels: Seq[File] = Seq()

  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = 4

  @Argument(doc="How many threads to use", fullName="num_threads", shortName="nt", required=false)
  var numThreads: Int = 4

  @Argument(doc="Number of threads BWA should use", fullName="bwa_threads", shortName="bt", required=false)
  var bwaThreads: Int = 4

  @Argument(doc="The path to create log files", fullName="queue_log_dir", shortName="logDir", required=false)
  val queueLogDir: String = ".qlog/"  // Gracefully hide Queue's output

  def withParallelism[A](n : Int)(block : => A) : A = {
    import collection.parallel.ForkJoinTasks.defaultForkJoinPool._
    val defaultParLevel = getParallelism
    setParallelism(n)
    val ret = block
    setParallelism(defaultParLevel)
    ret
  }

  def script() {
    qscript.outPath.mkdir()
    qscript.queueLogDir.mkdir()
    for (bam <- qscript.inputBams) {
    //withParallelism(4) {
      //qscript.inputBams.toList.par.foreach { bam =>
        // Extract sampld ID and run ID from file and directory names
        val basename        = bam.getName()
        var runId           = bam.getParent().getName()
        val sampleIdPattern = "^(TCGA-\\S{2}-\\S{4}-\\S{3})\\S+?\\_(\\S+?)(\\_s\\_\\S+|)\\.\\S+\\.bam".r
        val sampleIdPattern(sampleId, sampleRunId, lane) = basename
        runId = if (runId.toCharArray()(0).isDigit) runId else sampleRunId
        val newStem         = sampleId + "___" + runId
        val readGroup       = new ReadGroup(id = "0", lb = sampleId, pl = "Illumina", pu = runId, sm = sampleId, cn = "Harvard_GCC", ds = "WGS")
        // Revert original bam to name-sorted bam
        val revertedBam     = qscript.outPath + "/" + swapExt(newStem, ".bam", ".rv.bam")
        add(revertBam(bam, revertedBam, true))
        val sai1              = qscript.outPath + "/" + swapExt(revertedBam, ".rv.bam", ".1.sai")
        val sai2              = qscript.outPath + "/" + swapExt(revertedBam, ".rv.bam", ".2.sai")
        val newBam            = qscript.outPath + "/" + swapExt(revertedBam, ".rv.bam", ".new.bam")
        val sortedBam         = qscript.outPath + "/" + swapExt(revertedBam, ".rv.bam", ".new.st.bam")
        val dupMarkedBam      = qscript.outPath + "/" + swapExt(revertedBam, ".rv.bam", ".new.st.dm.bam")
        val metricsFile       = qscript.outPath + "/" + swapExt(revertedBam, ".rv.bam", ".new.st.dm.metrics")
        val realignInterval   = qscript.outPath + "/" + swapExt(revertedBam, ".rv.bam", ".new.st.dm.ra.intervals")
        val realignedBam      = qscript.outPath + "/" + swapExt(revertedBam, ".rv.bam", ".new.st.dm.ra.bam")
        val recalData         = qscript.outPath + "/" + swapExt(revertedBam, ".rv.bam", ".new.st.dm.ra.rc.data")
        val recalibratedBam   = qscript.outPath + "/" + swapExt(revertedBam, ".rv.bam", ".new.st.dm.ra.rc.bam")
        //val reducedBam        = qscript.outPath + "/" + swapExt(revertedBam, ".rv.bam", ".new.st.dm.ra.rc.rd.bam")
        add(bwaAlnPe(revertedBam, sai1, 1),
            bwaAlnPe(revertedBam, sai2, 2),
            bwaSamPe(revertedBam, sai1, sai2, newBam, readGroup),
            sortSam(newBam, sortedBam, SortOrder.coordinate),
            markDup(sortedBam, dupMarkedBam, metricsFile),
            targetIntervals(Seq(dupMarkedBam), realignInterval),
            realignIndel(Seq(dupMarkedBam), realignInterval, realignedBam),
            recalBase(Seq(realignedBam), recalData),
            recalBam(Seq(realignedBam), recalData, recalibratedBam))
            //reduceBam(Seq(recalibratedBam), reducedBam))
      //}
    }
  }

  class ReadGroup (val id: String,
                   val lb: String,
                   val pl: String,
                   val pu: String,
                   val sm: String,
                   val cn: String,
                   val ds: String) {}

  trait ExternalCommonArgs extends CommandLineFunction {
    this.memoryLimit = 4
    this.isIntermediate = true
  }

  trait SAMargs extends PicardBamFunction with ExternalCommonArgs {
      this.maxRecordsInRam = 100000
  }

  trait CommandLineGATKArgs extends CommandLineGATK with ExternalCommonArgs {
    this.reference_sequence = qscript.reference
  }

  case class revertBam (inBam: File, outBam: File, removeAlignmentInfo: Boolean) extends RevertSam with ExternalCommonArgs {
    this.input                      :+= inBam
    this.output                     = outBam
    this.removeAlignmentInformation = removeAlignmentInfo;
    this.sortOrder                  = if (removeAlignmentInfo) { SortOrder.queryname } else { SortOrder.coordinate }
    this.analysisName               = queueLogDir + outBam + ".reverted"
    this.jobName                    = queueLogDir + outBam + ".reverted"
    this.jobQueue                   = "long"
  }

  case class addReadGroup (inBam: File, outBam: File, readGroup: ReadGroup) extends AddOrReplaceReadGroups with ExternalCommonArgs {
    this.input :+= inBam
    this.output = outBam
    this.RGID = readGroup.id
    this.RGCN = readGroup.cn
    this.RGDS = readGroup.ds
    this.RGLB = readGroup.lb
    this.RGPL = readGroup.pl
    this.RGPU = readGroup.pu
    this.RGSM = readGroup.sm
    this.analysisName   = queueLogDir + outBam + ".rg"
    this.jobName        = queueLogDir + outBam + ".rg"
  }

  case class markDup (inBam: File, outBam: File, metricsFile: File) extends MarkDuplicates with ExternalCommonArgs {
    this.input        :+= inBam
    this.output       = outBam
    this.metrics      = metricsFile
    this.memoryLimit  = 4
    this.analysisName = queueLogDir + outBam + ".dm"
    this.jobName      = queueLogDir + outBam + ".dm"
    this.jobQueue     = "long"
  }

  case class sortSam (inSam: File, outBam: File, sortOrderP: SortOrder) extends SortSam with ExternalCommonArgs {
    this.input        :+= inSam
    this.output       = outBam
    this.sortOrder    = sortOrderP
    this.analysisName = queueLogDir + outBam + ".st"
    this.jobName      = queueLogDir + outBam + ".st"
    this.jobQueue     = "long"
  }

  //case class mergeBams (inBams: Seq[File], outBam: File) extends MergeSamFiles with ExternalCommonArgs {
    //this.input        = inBams
    //this.output       = outBam
    //this.sortOrder    = SortOrder.coordinate
    //this.analysisName = queueLogDir + "/" + outBam + ".mg"
    //this.jobName      = queueLogDir + "/" + outBam + ".mg"
    //this.jobQueue     = "long"
  //}

  //case class countReads (inBam: File, readCount: File) extends CommandLineFunction with ExternalCommonArgs {
    //@Input(doc="bam file to be sliced", required = true) var ibam = inBam
    //@Output(doc="file containing total number of reads in an input bam", required = true) var rcFile = readCount
    //def commandLine = "samtools view -c " + ibam + " > " + rcFile
    //this.analysisName   = queueLogDir + "/" + obam + ".readCount"
    //this.jobName        = queueLogDir + "/" + obam + ".readCount"
  //}

  //case class sliceBam (inBam: File, startPos: Int, stopPos: Int, outBam: File) extends CommandLineFunction with ExternalCommonArgs {
    //@Input(doc="bam file to be sliced", required = true) var ibam = inBam
    //@Argument(doc="strat position to be sliced", required = true) var posStart = startPos
    //@Argument(doc="end position to be sliced", required = true) var posStop = stopPos
    //@Output(doc="output bam file", required = true) var obam = outBam
    //val nextStartPos = stopPos + 1
    //val osam = qscript.outPath + "/" + swapExt(outBam, ".bam", ".sam")
    //def commandLine = "samtools view -H " + ibam + " > " + osam + " && samtools view " + ibam + " | sed -n '" + posStart + "," + posStop + "p;" + nextStartPos + "q' >> " + osam + " && samtools view -bSh " + osam + " > " + obam + " && rm " + osam
    //this.analysisName   = queueLogDir + "/" + obam + ".sliced"
    //this.jobName        = queueLogDir + "/" + obam + ".sliced"
  //}

  case class bwaAlnPe (inBam: File, outSai: File, index: Int) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Output(doc="output sai file for mating pair") var sai = outSai
    def commandLine   = bwaPath + " aln -t " + qscript.bwaThreads + " " + reference + " -b" + index + " " + bam + " > " + sai
    this.analysisName = queueLogDir + outSai + ".bwaAlnPe"
    this.jobName      = queueLogDir + outSai + ".bwaAlnPe"
    if (qscript.bwaThreads > 1) {
      this.jobQueue       = "mcore"
      this.nCoresRequest  = qscript.bwaThreads
    } else {
      this.jobQueue =      "long"
    }
  }

  case class bwaSamPe (inBam: File, inSai1: File, inSai2:File, outBam: File, readGroup: ReadGroup) extends CommandLineFunction with ExternalCommonArgs {
    @Input(doc="bam file to be aligned") var bam = inBam
    @Input(doc="bwa alignment index file for 1st mating pair") var sai1 = inSai1
    @Input(doc="bwa alignment index file for 2nd mating pair") var sai2 = inSai2
    @Output(doc="output aligned bam file") var alignedBam = outBam
    //@RG     ID:0    PL:Illumina     PUReducedBams   LB:TCGA-DK-A1A3-01A___120921_SN208_0431_BD1D5FACXX.mf.1_10000000.bam    SM:TCGA-DK-A1A3-01A___120921_SN208_0431_BD1D5FACXX.mf.1_10000000.bam
    //def commandLine = bwaPath + " sampe " + reference + " " + sai1 + " " + sai2 + " " + bam + " " + bam + " > " + alignedBam
    def commandLine   = bwaPath + " sampe " + reference + " " + sai1 + " " + sai2 + " " + bam + " " + bam + 
      " -r '@RG\\tID:" + readGroup.id + "\\tPL:" + readGroup.pl + "\\tPU:" + readGroup.pu + "\\tLB:" + readGroup.lb + "\\tSM:" + readGroup.sm + "' | " + samtoolsPath + " view -bSho " + alignedBam + " -"
    this.memoryLimit  = 6
    this.analysisName = queueLogDir + outBam + ".bwaSamPe"
    this.jobName      = queueLogDir + outBam + ".bwaSamPe"
    this.jobQueue     = "long"
  }

  case class targetIntervals (inBams: Seq[File], outIntervals: File) extends RealignerTargetCreator with CommandLineGATKArgs {
    this.input_file       = inBams
    this.out              = outIntervals
    this.mismatchFraction = 0.0
    this.known            ++= qscript.dbSNP
    this.known            ++= qscript.indels
    this.analysisName     = queueLogDir + outIntervals + ".ti"
    this.jobName          = queueLogDir + outIntervals + ".ti"
    this.num_threads      = qscript.numThreads
    if (qscript.numThreads > 1) {
      this.jobQueue       = "mcore"
      this.nCoresRequest  = qscript.numThreads
    } else {
      this.jobQueue       = "long"
    }
    //this.fix_misencoded_quality_scores = true
  }

  case class realignIndel (inBams: Seq[File], realignIntervals: File, outBam: File) extends IndelRealigner with CommandLineGATKArgs {
    this.input_file       = inBams
    this.targetIntervals  = realignIntervals
    this.out              = outBam
    this.known            ++= qscript.dbSNP
    this.known            ++= qscript.indels
    this.compress         = 0
    this.scatterCount     = qscript.nContigs
    this.analysisName     = queueLogDir + outBam + ".ra"
    this.jobName          = queueLogDir + outBam + ".ra"
    //this.fix_misencoded_quality_scores = true
  }

  case class recalBase (inBams: Seq[File], recalTable: File) extends BaseRecalibrator with CommandLineGATKArgs {
    this.input_file   = inBams
    this.out          = recalTable
    this.knownSites   ++= qscript.dbSNP
    this.knownSites   ++= qscript.indels
    this.scatterCount = qscript.nContigs
    if (qscript.numThreads > 1) {
      this.jobQueue       = "mcore"
      this.nCoresRequest  = qscript.numThreads
    }
    this.analysisName = queueLogDir + recalTable + ".rc"
    this.jobName      = queueLogDir + recalTable + ".rc"
  }

  case class recalBam (inBams: Seq[File], recalTable: File, outBam: File) extends PrintReads with CommandLineGATKArgs {
    this.input_file         = inBams
    this.out                = outBam
    this.BQSR               = recalTable
    if (qscript.numThreads > 1) {
      this.jobQueue         = "mcore"
      this.nCoresRequest    = qscript.numThreads
      this.num_cpu_threads_per_data_thread = qscript.numThreads
    } else {
      this.jobQueue         = "long"
    }
    this.analysisName       = queueLogDir + outBam + ".pr"
    this.jobName            = queueLogDir + outBam + ".pr"
  }

  case class reduceBam (inBams: Seq[File], outBam: File) extends ReduceReads with CommandLineGATKArgs {
    this.input_file         = inBams
    this.out                = outBam
    this.scatterCount       = qscript.nContigs
    this.analysisName       = queueLogDir + outBam + ".rd"
    this.jobName            = queueLogDir + outBam + ".rd"
    this.isIntermediate     = false
  }

}

