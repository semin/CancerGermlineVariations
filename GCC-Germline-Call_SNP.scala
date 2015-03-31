import scala.sys.process._
import scala.collection.JavaConversions._
import scala.collection.immutable.StringOps

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.util.QScriptUtils
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.picard._
import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction

import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model

import htsjdk.samtools.SAMFileReader
import htsjdk.samtools.SAMFileHeader.SortOrder

class CallSnpsIndels extends QScript {
  qscript =>

  @Input(doc="Input BAM file of list of BAM file(s)", fullName="input", shortName="i", required=true)
  var input: File = _

  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=true)
  var reference: File = _

  @Input(doc="Output directory", fullName="out_dir", shortName="O", required=true)
  var outPath: File = _

  @Input(doc="dbsnp ROD to use (must be in VCF format)", fullName="dbsnp", shortName="D", required=true)
  var dbSNP: File = _

  @Argument(doc="Output stem", fullName="out_stem", shortName="st", required=false)
  var outStem: String = "Harvard_GCC_WGS"

  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = 50

  @Argument(doc="How many threads to use", fullName="num_threads", shortName="nt", required=false)
  var numThreads: Int = 4

  def script() {
    qscript.outPath.mkdir()
    qSettings.logDirectory.mkdir()

    val bams = QScriptUtils.createSeqFromFile(input)
    val rawVcf = qscript.outPath + "/" + outStem + ".vcf.gz"
    val ug = new UnifiedGenotyper
    ug.scatterCount = nContigs
    ug.reference_sequence = reference
    ug.input_file = bams
    ug.out = rawVcf
    ug.genotype_likelihoods_model = Model.BOTH
    ug.stand_call_conf = 10
    ug.stand_emit_conf = 4
    ug.read_filter :+= "BadCigar"
    ug.isIntermediate = false
    ug.analysisName = rawVcf.getName() + ".ug"
    ug.jobName = rawVcf.getName() + ".ug"
    //ug.analysisName = qSettings.logDirectory + "/" + rawVcf.getName() + ".ug"
    //ug.jobName = qSettings.logDirectory + "/" + rawVcf.getName() + ".ug"
    //ug.fix_misencoded_quality_scores  = true
    ug.allow_potentially_misencoded_quality_scores = true
    ug.memoryLimit = 10
    ug.isIntermediate = false
    add(ug)

    // Annotate dbSNP information
    val dbVcf = qscript.outPath + "/" + outStem + ".dbsnp.vcf.gz"
    val va = new VariantAnnotator
    va.reference_sequence = reference
    va.scatterCount = qscript.nContigs
    va.variant = rawVcf
    va.intervals = Seq(rawVcf)
    va.out = dbVcf
    va.dbsnp = dbSNP
    add(va)
  }

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
}

