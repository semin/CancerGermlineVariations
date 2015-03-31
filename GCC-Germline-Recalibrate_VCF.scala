import scala.sys.process._
import scala.collection.JavaConversions._
import scala.collection.immutable.StringOps

import org.broadinstitute.gatk.queue.QScript
import org.broadinstitute.gatk.queue.util.QScriptUtils
import org.broadinstitute.gatk.queue.extensions.gatk._
import org.broadinstitute.gatk.queue.extensions.picard._
import org.broadinstitute.gatk.queue.function.JavaCommandLineFunction

import org.broadinstitute.gatk.tools.walkers.genotyper.GenotypeLikelihoodsCalculationModel.Model
import org.broadinstitute.gatk.tools.walkers.variantrecalibration.VariantRecalibratorArgumentCollection.Mode

import htsjdk.samtools.SAMFileReader
import htsjdk.samtools.SAMFileHeader.SortOrder
import htsjdk.variant.variantcontext.VariantContext.Type

class RecalibrateVcf extends QScript {
  qscript =>

  @Input(doc="Input VCF file of list of VCF file(s)", fullName="input", shortName="i", required=true)
  var rawVcf: File = _

  @Input(doc="Reference fasta file", fullName="reference", shortName="R", required=true)
  var reference: File = _

  @Input(doc="Output directory", fullName="out_dir", shortName="O", required=true)
  var outPath: File = _

  @Input(doc="dbsnp to use (must be in VCF format)", fullName="dbsnp", shortName="db", required=true)
  var dbSNP: File = _

  @Input(doc="Hapmap VCF to use", fullName="hapmap", shortName="hm", required=true)
  var hapmap: File = _

  @Input(doc="Omni genotypying array VCF to use", fullName="omni", shortName="om", required=true)
  var omni: File = _

  @Input(doc="1K genomes high confidnence SNPs", fullName="genomes1k", shortName="kg", required=true)
  var genomes1k: File = _

  @Input(doc="Mills INDELs", fullName="mills", shortName="ml", required=true)
  var mills: File = _

  @Argument(doc="How many ways to scatter/gather", fullName="scatter_gather", shortName="sg", required=false)
  var nContigs: Int = 50

  @Argument(doc="How many threads to use", fullName="num_threads", shortName="nt", required=false)
  var numThreads: Int = 4

  def script() {
    qscript.outPath.mkdir()
    qSettings.logDirectory.mkdir()

    // SNP
    val rawSnp = qscript.outPath + "/" + swapExt(rawVcf.getName(), ".vcf.gz", ".snp.vcf.gz")
    val svSnp  = new SelectVariants
    svSnp.reference_sequence = reference
    svSnp.variant = rawVcf
    svSnp.out = rawSnp
    svSnp.selectTypeToInclude = Seq(Type.SNP, Type.MNP)
    add(svSnp)

    val recalSnp = qscript.outPath + "/" + swapExt(rawSnp, ".vcf.gz", ".recal.vcf.gz")
    val tranchesSnp = qscript.outPath + "/" + swapExt(rawSnp, ".vcf.gz", ".recal.tranches")
    val rscriptSnp = qscript.outPath + "/" + swapExt(rawSnp, ".vcf.gz", ".recal.plots.R")
    val vrSnp = new VariantRecalibrator
    vrSnp.reference_sequence = reference
    vrSnp.input = Seq(new File(rawSnp))
    vrSnp.resource = Seq( new TaggedFile(hapmap, "known=false,training=true,truth=true,prior=15.0"),
                          new TaggedFile(omni, "known=false,training=true,truth=true,prior=12.0"),
                          new TaggedFile(genomes1k, "known=false,training=true,truth=false,prior=10.0"),
                          new TaggedFile(dbSNP, "known=true,training=false,truth=false,prior=2.0"))
    vrSnp.use_annotation = Seq("DP", "QD", "FS", "MQRankSum", "ReadPosRankSum")
    vrSnp.mode = Mode.SNP
    vrSnp.minNumBadVariants = 3000
    vrSnp.recal_file = recalSnp
    vrSnp.tranches_file = tranchesSnp
    vrSnp.rscript_file = rscriptSnp
    vrSnp.memoryLimit = 80
    vrSnp.residentLimit = 80
    vrSnp.analysisName = rawVcf.getName() + ".vrSnp"
    vrSnp.jobName = rawVcf.getName() + ".vrSnp"
    //vrSnp.jobNativeArgs = Seq("-M 70000000")
    add(vrSnp)

    val vqsrSnp = qscript.outPath + "/" + swapExt(rawSnp.getName(), ".vcf.gz", ".vqsr.vcf.gz")
    val arSnp = new ApplyRecalibration
    arSnp.reference_sequence = reference
    arSnp.input = Seq(new File(rawSnp))
    arSnp.out = vqsrSnp
    arSnp.mode = Mode.SNP
    arSnp.ts_filter_level = 99.0
    arSnp.recal_file = recalSnp
    arSnp.tranches_file = tranchesSnp
    arSnp.memoryLimit = 80
    arSnp.residentLimit = 80
    arSnp.analysisName = rawVcf.getName() + ".arSnp"
    arSnp.jobName = rawVcf.getName() + ".arSnp"
    //arSnp.jobNativeArgs = Seq("-M 70000000")
    add(arSnp)

    // Indel
    val rawIndel = qscript.outPath + "/" + swapExt(rawVcf.getName(), ".vcf.gz", ".indel.vcf.gz")
    val svIndel  = new SelectVariants
    svIndel.reference_sequence = reference
    svIndel.variant = rawVcf
    svIndel.out = rawIndel
    svIndel.selectTypeToInclude = Seq(Type.INDEL)
    add(svIndel)

    val recalIndel = qscript.outPath + "/" + swapExt(rawIndel, ".vcf.gz", ".recal.vcf.gz")
    val tranchesIndel = qscript.outPath + "/" + swapExt(rawIndel, ".vcf.gz", ".recal.tranches")
    val rscriptIndel = qscript.outPath + "/" + swapExt(rawIndel, ".vcf.gz", ".recal.plots.R")
    val vrIndel = new VariantRecalibrator
    vrIndel.reference_sequence = reference
    vrIndel.input = Seq(new File(rawIndel))
    vrIndel.resource = Seq(new TaggedFile(mills, "known=true,training=true,truth=true,prior=12.0"))
    vrIndel.use_annotation = Seq("DP", "FS", "MQRankSum", "ReadPosRankSum")
    vrIndel.mode = Mode.INDEL
    vrIndel.minNumBadVariants = 3000
    vrIndel.recal_file = recalIndel
    vrIndel.tranches_file = tranchesIndel
    vrIndel.rscript_file = rscriptIndel
    vrIndel.maxGaussians = 4
    vrIndel.memoryLimit = 20 
    vrIndel.residentLimit = 25
    vrIndel.analysisName = rawVcf.getName() + ".vrIndel"
    vrIndel.jobName = rawVcf.getName() + ".vrIndel"
    add(vrIndel)

    val vqsrIndel = qscript.outPath + "/" + swapExt(rawIndel.getName(), ".vcf.gz", ".vqsr.vcf.gz")
    val arIndel = new ApplyRecalibration
    arIndel.reference_sequence = reference
    arIndel.input = Seq(new File(rawIndel))
    arIndel.out = vqsrIndel
    arIndel.mode = Mode.INDEL
    arIndel.ts_filter_level = 99.0
    arIndel.recal_file = recalIndel
    arIndel.tranches_file = tranchesIndel
    
    arIndel.residentLimit = 15
    arIndel.analysisName = rawVcf.getName() + ".arIndel"
    arIndel.jobName = rawVcf.getName() + ".arIndel"
    add(arIndel)
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

