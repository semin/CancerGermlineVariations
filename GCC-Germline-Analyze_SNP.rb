#!/usr/bin/env ruby

require 'open3'
require 'logger'
require "rinruby"
require 'sqlite3'
require "parallel"
require 'pathname'
require "google_drive"

$logger         = Logger.new(STDOUT)
$logger.level   = Logger::DEBUG
$javaBin        = Pathname.new("/opt/java/jdk7/bin/java")
$gccBamDir1     = Pathname.new("/groups/kucherlapati/GCC/LevelII")
$gccBamDir2     = Pathname.new("/files/CBMI/parklab/tcga/hg18BAMs")
$homeDir        = Pathname.new("/home/sl279")
$baseDir        = Pathname.new("/groups/kucherlapati/GCC/Germline")
$scriptDir      = $baseDir + "Scripts"
$gdacDir        = $baseDir + "GDAC"
$dbsnpDir       = $baseDir + "DBSNP"
$vcfDir         = $baseDir + "VCF.new";$vcfDir.mkpath
$canVcfDir      = $vcfDir + "CANCERS";$canVcfDir.mkpath
$panVcfDir      = $vcfDir + "PANCAN";$panVcfDir.mkpath
$installDir     = $homeDir + "BiO/Install"
$bwaBin         = $installDir + "bwa-0.7.5a/bwa"
$samtoolsBin    = $installDir + "samtools-0.1.19/samtools"
$queueBin       = $installDir + "Sting/dist/Queue.jar"
$gatkBin        = $installDir + "GATK2/GenomeAnalysisTK.jar"
$gatkBundleDir  = $installDir + "GATK-bundle"
$picardBinDir   = $installDir + "picard-tools"
#$refSeq         = $gatkBundleDir + "human_g1k_v37_decoy.fasta"
$refSeq         = $gccBamDir1 + "hg19Ref/Homo_sapiens_assembly19.fasta"
$dbsnp          = $gatkBundleDir + "dbsnp_137.b37.vcf"
$hapmap         = $gatkBundleDir + "hapmap_3.3.b37.vcf"
$g1kIndel       = $gatkBundleDir + "1000G_phase1.indels.b37.vcf"
$millsIndel     = $gatkBundleDir + "Mills_and_1000G_gold_standard.indels.b37.vcf"
$g1kOmni        = $gatkBundleDir + "1000G_omni2.5.b37.vcf"
$annovarDir     = Pathname.new "/groups/park/semin/BiO/Install/annovar"
$annovarBin     = $annovarDir + "annotate_variation.pl"
$convertBin     = $annovarDir + "convert2annovar.pl"
$annovarDbDir   = $annovarDir + "humandb"
$chrs           = (1..22).to_a << "X" << "Y" << "MT"
#$chrs           = (9..22).to_a << "X" << "Y" << "MT"
$chrchrs        = $chrs.map { |c| "chr#{c}" }
#$tmpDir         = Pathname.new("/hms/scratch1/sl279/tmp"); $tmpDir.mkpath
$tmpDir         = Pathname.new("/groups/park/semin/BiO/Temp"); $tmpDir.mkpath
$cancerTypes    = %w[BLCA BRCA CRC ESCA HNSC LGG LUAD PRAD SKCM STAD THCA UCEC]
$chromLenb37    = $gatkBundleDir + "hg19.genome"

class String
  def underscore
    self.gsub(/::/, '/').
      gsub(/([A-Z]+)([A-Z][a-z])/,'\1_\2').
      gsub(/([a-z\d])([A-Z])/,'\1_\2').
      tr("-", "_").
      downcase
  end
end

def submit(cmd)
  #puts cmd
  out = `#{cmd}`.chomp
  if out =~ /Job\s{1}<(\d+)>/
    puts "Job #{$1} submitted."
    return Integer($1)
  else
    abort out
  end
end

def sqlite3_import_tsv(dbFile, tsvFile, tableName)
  status = Open3::popen3("sqlite3 #{dbFile}") do |stdin, stdout, stderr, wait_thr|
    puts "Importing #{tsvFile} into #{dbFile}:#{tableName}..."
    pid = wait_thr.pid
    stdin.puts ".separator '\t'"
    stdin.puts ".import #{tsvFile} #{tableName}"
    stdin.puts ".exit"
    stdin.close
    stdout.close
    stderr.close
  end
  return status
end

def split_vqsr_vcfs_into_chromosomes
  vqsrVcfs = Pathname.glob($canVcfDir + "*/*.vqsr.vcf.gz")
  vqsrVcfs.each do |vqsrVcf|
    next if vqsrVcf.dirname.basename.to_s != "CRC"
    varType = vqsrVcf.basename.to_s.match(/indel/) ? "indel" : "snp"
    chrDir = vqsrVcf.dirname + varType + "chrs"; chrDir.mkpath
    $chrs.each do |chr|
      chrVcf = chrDir + vqsrVcf.basename(".gz").sub_ext(".#{chr}.vcf.gz")
      lsfOut = chrVcf.sub_ext(".gz.lsfout")
      splitCmd = <<-CMD
        bsub \\
          -q park_12h \\
          -g /gcc/germ/split \\
          -o #{lsfOut} \\
          "#{$javaBin} -Xmx5G -jar #{$gatkBin} \\
            -R #{$refSeq} \\
            -T SelectVariants \\
            --variant #{vqsrVcf} \\
            -o #{chrVcf} \\
            -L #{chr} &&
          tabix -p vcf #{chrVcf}"
      CMD
      submit splitCmd
    end
  end
end

def merge_chr_vqsr_vcfs
  $chrs.each do |chr|
    vqsrChrVcfs = Pathname.glob($canVcfDir + "*" + "chrs" + "*.vqsr.snp.#{chr}.vcf").sort
    if vqsrChrVcfs.size != $cancerTypes.size
      $logger.error "Found only #{vqsrChrVcfs.size} vcfs!"
      exit
    end
    panChrVcfDir  = $panVcfDir + "chrs"; panChrVcfDir.mkpath
    outVcf        = panChrVcfDir + "Harvard_GCC_WGS-PANCAN12-Normal.vqsr.snp.#{chr}.vcf"
    lsfOut        = outVcf.sub_ext(".vcf.lsfout")
    combineCmd = <<-CMD
      bsub \\
        -q i2b2_1d \\
        -g /gcc/vcf/comb \\
        -o #{lsfOut} \\
        #{$javaBin} -Xmx5G -jar #{$gatkBinDir}/GenomeAnalysisTK.jar \\
          -R #{$refSeq} \\
          -T CombineVariants \\
          #{vqsrChrVcfs.map { |v| "--variant:#{v.basename.to_s.split('-')[1]} #{v}" }.join(" ")} \\
          -o #{outVcf} \\
          --genotypemergeoption REQUIRE_UNIQUE
    CMD
    puts combineCmd
  end
end

def split_vqsr_vcfs_into_pieces
  # chunk genomic regions
  #chunkSize = 1000000  #(e.g. CRC)
  chunkSize = 500000 #(e.g. LUAD)
  chromLen  = {}
  $chromLenb37.each_line do |line|
    chrom, len = line.chomp.split(/\s+/)
    len = len.to_i
    next if len == 0
    chrom = "MT" if chrom == "M"
    chromLen[chrom] = len
  end
  vqsrVcfs = Pathname.glob($canVcfDir + "UCEC" + "*.vqsr.vcf.gz")
  vqsrVcfs.each do |vqsrVcf|
    varType = vqsrVcf.basename.to_s.match(/indel/) ? "indel" : "snp"
    #next if varType == "indel"
    chrDir = vqsrVcf.dirname + varType + "chrs"; chrDir.mkpath
    $chrs.each do |chr|
      slices      = (1..chromLen[chr.to_s]).each_slice(chunkSize)
      num_slices  = slices.to_a.size
      #slices.each_with_index do |slice, si|
      Parallel.each_with_index(slices, :in_threads => 4) do |slice, si|
        roi     = chr == "MT" ? "#{chr}:1-16569" : "#{chr}:#{slice[0]}-#{slice[-1]}"
        roiNum  = "#{(si+1).to_s.rjust(num_slices.to_s.size, '0')}_of_#{num_slices}"
        roiDir  = chrDir + chr.to_s + roiNum; roiDir.mkpath
        roiVcf  = roiDir + vqsrVcf.basename(".gz").sub_ext(".#{roiNum}.vcf.gz")
        lsfOut  = roiVcf.sub_ext(".gz.lsfout")
        if roiVcf.exist?
          if lsfOut.successful?
            next
          else
            roiVcf.delete
            lsfOut.delete
          end
        else
          lsfOut.delete if lsfOut.exist?
        end
        cmd = <<-CMD
          bsub \\
            -q mini \\
            -g /gcc/germ/split \\
            -o #{lsfOut} \\
            "#{$javaBin} -Xmx5G -jar #{$gatkBin} \\
              -R #{$refSeq} \\
              -T SelectVariants \\
              --variant #{vqsrVcf} \\
              -o #{roiVcf} \\
              -L #{roi} && \\
            tabix -p vcf #{roiVcf}"
        CMD
        submit cmd
      end
    end
  end
end

def split_vqsr_chr_vcfs_into_even_sized_chunks
  vqsrVcfs = Pathname.glob($canVcfDir + "UCEC" + "*.vqsr.vcf.gz")
  vqsrVcfs.each do |vqsrVcf|
    varType = vqsrVcf.basename.to_s.match(/indel/) ? "indel" : "snp"
    #next if varType == "indel"
    chrDir = vqsrVcf.dirname + varType + "chrs"; chrDir.mkpath
    $chrs.each do |chr|
      slices      = (1..chromLen[chr.to_s]).each_slice(chunkSize)
      num_slices  = slices.to_a.size
      #slices.each_with_index do |slice, si|
      Parallel.each_with_index(slices, :in_threads => 4) do |slice, si|
        roi     = chr == "MT" ? "#{chr}:1-16569" : "#{chr}:#{slice[0]}-#{slice[-1]}"
        roiNum  = "#{(si+1).to_s.rjust(num_slices.to_s.size, '0')}_of_#{num_slices}"
        roiDir  = chrDir + chr.to_s + roiNum; roiDir.mkpath
        roiVcf  = roiDir + vqsrVcf.basename(".gz").sub_ext(".#{roiNum}.vcf.gz")
        lsfOut  = roiVcf.sub_ext(".gz.lsfout")
        if roiVcf.exist?
          if lsfOut.successful?
            next
          else
            roiVcf.delete
            lsfOut.delete
          end
        else
          lsfOut.delete if lsfOut.exist?
        end
        cmd = <<-CMD
          bsub \\
            -q mini \\
            -g /gcc/germ/split \\
            -o #{lsfOut} \\
            "#{$javaBin} -Xmx5G -jar #{$gatkBin} \\
              -R #{$refSeq} \\
              -T SelectVariants \\
              --variant #{vqsrVcf} \\
              -o #{roiVcf} \\
              -L #{roi} && \\
            tabix -p vcf #{roiVcf}"
        CMD
        submit cmd
      end
    end
  end
end


def filter_vqsr_vcfs
  #Harvard_GCC_WGS-LUAD-Normal.dbsnp.indel.vqsr.001_of_499.vcf.gz
  vcfs = Pathname.glob($canVcfDir + "UCEC/{snp,indel}/chrs/*/*/*.vcf.gz").select { |x| x.to_s.match(/\d+\_of\_\d+\.vcf\.gz$/) }.sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    #vcfOut = vcf.dirname + vcf.basename(".gz").sub_ext(".af01.mss20.vcf.gz")
    vcfOut = vcf.dirname + vcf.basename(".gz").sub_ext(".mss20.vcf.gz")
    puts vcfOut
    lsfOut = vcfOut.sub_ext(".gz.lsfout")
    if vcfOut.exist?
      if lsfOut.successful?
        next
      else
        vcfOut.delete
        lsfOut.delete
      end
    else
      lsfOut.delete if lsfOut.exist?
    end
    cmd =<<-CMD
      bsub \\
        -q mini \\
        -g /gcc/germ/filter \\
        -o #{lsfOut} \\
        "zcat #{vcf} | ruby -ne '
          if (\\$_.start_with?(\\"#\\"))
            print \\$_;
          else
            cols              = \\$_.chomp.split(\\"\\\\t\\")
            next if cols[6]   != \\"PASS\\"
            alleleFrcs        = cols[7].match(/AF=(\\S+?);/)[1].split(\\",\\").map { |x| x.to_f }
            sampleSize        = cols.size - 9
            missingSamples    = cols[9..-1].select { |x| x.start_with?(\\"./.\\") }
            missingSampleFrc  = missingSamples.size / sampleSize.to_f
            print(\\$_) if (missingSampleFrc < 0.2)
          end' | bgzip -c > #{vcfOut} &&
          tabix -p vcf #{vcfOut}"
    CMD
            #print(\\$_) if (missingSampleFrc < 0.2 && alleleFrcs[0] > 0.01)
    submit cmd
  end
end

def annotate_vcfs_with_dbsnp
  vcfs = Pathname.glob($canVcfDir + "*/*.vqsr.{snp,indel}.mss20.vcf").sort
  vcfs.each_with_index do |vcf, vi|
    newVcf = vcf.sub_ext(".dbsnp.vcf")
    lsfOut = newVcf.sub_ext(".vcf.lsfout")
    cmd =<<-CMD
      bsub \\
        -q park_1d \\
        -g /gcc/germ/va \\
        -R "rusage[mem=20000]" -M 20000000 \\
        -o #{lsfOut} \\
        java -jar -Xmx18G #{$gatkBinDir}/GenomeAnalysisTK.jar \\
          -T VariantAnnotator \\
          -R #{$refSeq} \\
          -o #{newVcf} \\
          -D #{$dbsnp} \\
          --variant #{vcf} \\
          -L #{vcf}
    CMD
    submit cmd
  end
end

def annotate_vcfs_with_dbsnp_rsids
  dbsnp = $dbsnpDir + "human_9606_b142_GRCh37p13/VCF/All.rsid.exp.vcf.gz"
  vcfs = Pathname.glob($canVcfDir + "UCEC/{snp,indel}/chrs/*/*/*.mss20.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    dbsnp_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".dbsnp.vcf.gz")
    lsfout = dbsnp_vcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -g /gcc/germ/dbsnp \\
      -q mini \\
      -o #{lsfout} \\
      "zcat #{vcf} | vcf-annotate -a #{dbsnp} -d key=INFO,ID=DBSNP,Number=A,Type=String,Description=\\"dbSNP142 RSID\\" -c CHROM,POS,INFO/DBSNP,REF,ALT | bgzip -c > #{dbsnp_vcf} && tabix -p vcf #{dbsnp_vcf}"
    CMD
    submit cmd
  end
end

def annotate_vcfs_with_1000genomes_afs
  ###INFO=<ID=AF,Number=A,Type=Float,Description="Estimated allele frequency in the range (0,1)">
  ###INFO=<ID=EAS_AF,Number=A,Type=Float,Description="Allele frequency in the EAS populations calculated from AC and AN, in the range (0,1)">
  ###INFO=<ID=EUR_AF,Number=A,Type=Float,Description="Allele frequency in the EUR populations calculated from AC and AN, in the range (0,1)">
  ###INFO=<ID=AFR_AF,Number=A,Type=Float,Description="Allele frequency in the AFR populations calculated from AC and AN, in the range (0,1)">
  ###INFO=<ID=AMR_AF,Number=A,Type=Float,Description="Allele frequency in the AMR populations calculated from AC and AN, in the range (0,1)">
  ###INFO=<ID=SAS_AF,Number=A,Type=Float,Description="Allele frequency in the SAS populations calculated from AC and AN, in the range (0,1)">

  #tgaf = Pathname.new("/groups/kucherlapati/GCC/Germline/1000Genomes/ALL.wgs.phase1_release_v3.20101123.snps_indels_sv.sites.vcf.afs.tsv.gz")
  #pops = %w[ALL ASN AMR AFR EUR]

  tgaf = $baseDir + "1000Genomes" + "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5.20130502.sites.afs.tsv.gz"
  pops = %w[ALL EAS EUR AFR AMR SAS]

  vcfs = Pathname.glob($canVcfDir + "UCEC/{snp,indel}/chrs/*/*/*.mss20.dbsnp.vcf.gz").sort
  #vcfs = Pathname.glob($canVcfDir + "UCEC/snp/chrs/*/*/*.mss20.dbsnp.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    tgALL_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".tgALL.vcf.gz")
    tgEAS_vcf = tgALL_vcf.dirname + tgALL_vcf.basename(".gz").sub_ext(".tgEAS.vcf.gz")
    tgEUR_vcf = tgEAS_vcf.dirname + tgEAS_vcf.basename(".gz").sub_ext(".tgEUR.vcf.gz")
    tgAFR_vcf = tgEUR_vcf.dirname + tgEUR_vcf.basename(".gz").sub_ext(".tgAFR.vcf.gz")
    tgAMR_vcf = tgAFR_vcf.dirname + tgAFR_vcf.basename(".gz").sub_ext(".tgAMR.vcf.gz")
    tgSAS_vcf = tgAMR_vcf.dirname + tgAMR_vcf.basename(".gz").sub_ext(".tgSAS.vcf.gz")
    lsfout    = tgSAS_vcf.sub_ext(".gz.lsfout")
    cmd       = <<-CMD
    bsub \\
      -g /gcc/germ/tg \\
      -q mini \\
      -o #{lsfout} "
      zcat #{vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_ALL_AF,Number=A,Type=Float,Description=\\"Global Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,INFO/TG_ALL_AF,-,-,-,-,- | bgzip -c > #{tgALL_vcf} &&
      tabix -p vcf #{tgALL_vcf} &&
      zcat #{tgALL_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_EAS_AF,Number=A,Type=Float,Description=\\"EAS Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,INFO/TG_EAS_AF,-,-,-,- | bgzip -c > #{tgEAS_vcf} &&
      tabix -p vcf #{tgEAS_vcf} &&
      zcat #{tgEAS_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_EUR_AF,Number=A,Type=Float,Description=\\"EUR Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,INFO/TG_EUR_AF,-,-,- | bgzip -c > #{tgEUR_vcf} &&
      tabix -p vcf #{tgEUR_vcf} &&
      zcat #{tgEUR_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_AFR_AF,Number=A,Type=Float,Description=\\"AFR Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,-,INFO/TG_AFR_AF,-,- | bgzip -c > #{tgAFR_vcf} &&
      tabix -p vcf #{tgAFR_vcf} &&
      zcat #{tgAFR_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_AMR_AF,Number=A,Type=Float,Description=\\"AMR Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,-,-,INFO/TG_AMR_AF,- | bgzip -c > #{tgAMR_vcf} &&
      tabix -p vcf #{tgAMR_vcf} &&
      zcat #{tgAMR_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_SAS_AF,Number=A,Type=Float,Description=\\"SAS Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,-,-,-,INFO/TG_SAS_AF | bgzip -c > #{tgSAS_vcf} &&
      tabix -p vcf #{tgSAS_vcf} &&
      rm -rf #{tgALL_vcf} #{tgEAS_vcf} #{tgEUR_vcf} #{tgAFR_vcf} #{tgAMR_vcf} &&
      rm -rf #{tgALL_vcf}.tbi #{tgEAS_vcf}.tbi #{tgEUR_vcf}.tbi #{tgAFR_vcf}.tbi #{tgAMR_vcf}.tbi
      "
    CMD
    submit cmd
  end
end

def annotate_snp_vcfs_with_cadd_scores
  cadd = Pathname.new("/groups/kucherlapati/GCC/Germline/CADD/whole_genome_SNVs.tsv.gz")
  vcfs = Pathname.glob($canVcfDir + "UCEC/snp/chrs/*/*/*.tgSAS.vcf.gz").sort
  #vcfs.each_with_index do |vcf, vi|
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    cadd_raw_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".cadd_raw.vcf.gz")
    cadd_phred_vcf = cadd_raw_vcf.dirname + cadd_raw_vcf.basename(".gz").sub_ext(".cadd_phred.vcf.gz")
    lsfout = cadd_phred_vcf.sub_ext(".gz.lsfout")
    if cadd_phred_vcf.exist? && cadd_raw_vcf.exist?
      if lsfout.successful?
        next
      else
        cadd_raw_vcf.delete
        cadd_phred_vcf.delete
        lsfout.delete
      end
    else
      lsfout.delete if lsfout.exist?
      cadd_raw_vcf.delete if cadd_raw_vcf.exist?
      cadd_phred_vcf.delete if cadd_phred_vcf.exist?
    end
    cmd = <<-CMD
    bsub \\
      -g /gcc/germ/cadd \\
      -q mini \\
      -o #{lsfout} \\
      "zcat #{vcf} | vcf-annotate -a #{cadd} -d key=INFO,ID=CADD_RAW,Number=A,Type=Float,Description=\\"CADD Raw Score\\" -c CHROM,POS,REF,ALT,INFO/CADD_RAW,- | bgzip -c > #{cadd_raw_vcf} &&
      tabix -p vcf #{cadd_raw_vcf} &&
      zcat #{cadd_raw_vcf} | vcf-annotate -a #{cadd} -d key=INFO,ID=CADD_PHRED,Number=A,Type=Float,Description=\\"CADD Phred Score\\" -c CHROM,POS,REF,ALT,-,INFO/CADD_PHRED | bgzip -c > #{cadd_phred_vcf}
      tabix -p vcf #{cadd_phred_vcf}"
    CMD
    submit cmd
  end
end

def get_not_hispanic_white_patient_indices_of(cancerType, vcf)
  # read ethinicity info from clinical data
  #patient.bcrpatientbarcode       tcga-a6-2671    tcga-a6-2672    tcga-a6-2674    tcga-a6-2675    tcga-a6-2676    tcga-a6-2677    tcga-a6-2678    tcga-a6-2679    tcga-a6-2680    tcga-a6-2681    tcga-a6-2682    tcga-a6-2683    tcga-a6-2684    tcga-a6-2685    tcga-a6-2686    tcga-a6-3807    tcg
  #patient.ethnicity       not hispanic or latino  not hispanic or latino  not hispanic or latino  not hispanic or latino  not hispanic or latino  not hispanic or latino  not hispanic or latino  not hispanic or latino  not hispanic or latino  not hispanic or latino  not hispanic or latino  not
  #patient.race    white   white   white   white   white   white   white   white   black or african american       white   white   white   white   black or african american       white   white   white   white   white   white   white   white   black or african american       white   black or af
  #clinFile = Pathname.glob($gdacDir + "stddata__2014_10_17" + cancerType + "20141017" + "*" + "#{cancerType}.clin.stage2File.txt")[0]
  clinFile = Pathname.glob($gdacDir + "stddata__2014_10_17" + cancerType + "20141017" + "*" + "#{cancerType}.clin.merged.picked.txt")[0]
  unless clinFile.exist?
    abort "Cannot fild clinical information for #{cancerType}!"
  end
  barcodes = nil
  ethinicities = nil
  races = nil
  clinFile.each_line do |line|
    if line.start_with?("Hybridization REF")
      barcodes = line.chomp.split("\t")[1..-1]
    elsif line.start_with?("ethnicity")
      ethinicities = line.chomp.split("\t")[1..-1]
    elsif line.start_with?("race")
      races = line.chomp.split("\t")[1..-1]
    end
  end
  notHispanicIndices = ethinicities.each_index.select { |i| ethinicities[i] == "not hispanic or latino" }
  whiteIndices = races.each_index.select { |i| races[i] == "white" }
  notHispanicWhiteIndices = notHispanicIndices & whiteIndices
  notHispanicWhiteBarcodes = barcodes.values_at(*notHispanicWhiteIndices).map(&:upcase)

  # Read TCGA bardcodes in vcf
  patientIds = nil
  notHispanicWhitePatientIndices = nil
  gzippedVcfFile = open(vcf)
  gzippedVcfFh = Zlib::GzipReader.new(gzippedVcfFile)
  gzippedVcfFh.each_line do |line|
    if line.start_with?("#CHROM")
      cols = line.chomp.split("\t")
      patientIds = cols[9..-1].map { |x| x.split("-")[0..2].join("-") }
      notHispanicWhitePatientIndices = patientIds.each_index.select { |i| notHispanicWhiteBarcodes.include?(patientIds[i]) }
      return(notHispanicWhitePatientIndices)
    end
  end
end

def get_1k_genomes_phase1_population_counts
  samplesInTgVcf = %w[HG00096 HG00097 HG00099 HG00100 HG00101 HG00102 HG00103 HG00104 HG00106 HG00108 HG00109 HG00110 HG00111 HG00112 HG00113 HG00114 HG00116 HG00117 HG00118 HG00119 HG00120 HG00121 HG00122 HG00123 HG00124 HG00125 HG00126 HG00127 HG00128 HG00129 HG00130 HG00131 HG00133 HG00134 HG00135 HG00136 HG00137 HG00138 HG00139 HG00140 HG00141 HG00142 HG00143 HG00146 HG00148 HG00149 HG00150 HG00151 HG00152 HG00154 HG00155 HG00156 HG00158 HG00159 HG00160 HG00171 HG00173 HG00174 HG00176 HG00177 HG00178 HG00179 HG00180 HG00182 HG00183 HG00185 HG00186 HG00187 HG00188 HG00189 HG00190 HG00231 HG00232 HG00233 HG00234 HG00235 HG00236 HG00237 HG00238 HG00239 HG00240 HG00242 HG00243 HG00244 HG00245 HG00246 HG00247 HG00249 HG00250 HG00251 HG00252 HG00253 HG00254 HG00255 HG00256 HG00257 HG00258 HG00259 HG00260 HG00261 HG00262 HG00263 HG00264 HG00265 HG00266 HG00267 HG00268 HG00269 HG00270 HG00271 HG00272 HG00273 HG00274 HG00275 HG00276 HG00277 HG00278 HG00280 HG00281 HG00282 HG00284 HG00285 HG00306 HG00309 HG00310 HG00311 HG00312 HG00313 HG00315 HG00318 HG00319 HG00320 HG00321 HG00323 HG00324 HG00325 HG00326 HG00327 HG00328 HG00329 HG00330 HG00331 HG00332 HG00334 HG00335 HG00336 HG00337 HG00338 HG00339 HG00341 HG00342 HG00343 HG00344 HG00345 HG00346 HG00349 HG00350 HG00351 HG00353 HG00355 HG00356 HG00357 HG00358 HG00359 HG00360 HG00361 HG00362 HG00364 HG00366 HG00367 HG00369 HG00372 HG00373 HG00375 HG00376 HG00377 HG00378 HG00381 HG00382 HG00383 HG00384 HG00403 HG00404 HG00406 HG00407 HG00418 HG00419 HG00421 HG00422 HG00427 HG00428 HG00436 HG00437 HG00442 HG00443 HG00445 HG00446 HG00448 HG00449 HG00451 HG00452 HG00457 HG00458 HG00463 HG00464 HG00472 HG00473 HG00475 HG00476 HG00478 HG00479 HG00500 HG00501 HG00512 HG00513 HG00524 HG00525 HG00530 HG00531 HG00533 HG00534 HG00536 HG00537 HG00542 HG00543 HG00553 HG00554 HG00556 HG00557 HG00559 HG00560 HG00565 HG00566 HG00577 HG00578 HG00580 HG00581 HG00583 HG00584 HG00589 HG00590 HG00592 HG00593 HG00595 HG00596 HG00607 HG00608 HG00610 HG00611 HG00613 HG00614 HG00619 HG00620 HG00625 HG00626 HG00628 HG00629 HG00634 HG00635 HG00637 HG00638 HG00640 HG00641 HG00650 HG00651 HG00653 HG00654 HG00656 HG00657 HG00662 HG00663 HG00671 HG00672 HG00683 HG00684 HG00689 HG00690 HG00692 HG00693 HG00698 HG00699 HG00701 HG00702 HG00704 HG00705 HG00707 HG00708 HG00731 HG00732 HG00734 HG00736 HG00737 HG00740 HG01047 HG01048 HG01051 HG01052 HG01055 HG01060 HG01061 HG01066 HG01067 HG01069 HG01070 HG01072 HG01073 HG01075 HG01079 HG01080 HG01082 HG01083 HG01085 HG01095 HG01097 HG01098 HG01101 HG01102 HG01104 HG01105 HG01107 HG01108 HG01112 HG01113 HG01124 HG01125 HG01133 HG01134 HG01136 HG01137 HG01140 HG01148 HG01149 HG01167 HG01168 HG01170 HG01171 HG01173 HG01174 HG01176 HG01183 HG01187 HG01188 HG01190 HG01191 HG01197 HG01198 HG01204 HG01250 HG01251 HG01257 HG01259 HG01271 HG01272 HG01274 HG01275 HG01277 HG01278 HG01334 HG01342 HG01344 HG01345 HG01350 HG01351 HG01353 HG01354 HG01356 HG01357 HG01359 HG01360 HG01365 HG01366 HG01374 HG01375 HG01377 HG01378 HG01383 HG01384 HG01389 HG01390 HG01437 HG01440 HG01441 HG01455 HG01456 HG01461 HG01462 HG01465 HG01488 HG01489 HG01491 HG01492 HG01494 HG01495 HG01497 HG01498 HG01515 HG01516 HG01518 HG01519 HG01521 HG01522 HG01550 HG01551 HG01617 HG01618 HG01619 HG01620 HG01623 HG01624 HG01625 HG01626 NA06984 NA06986 NA06989 NA06994 NA07000 NA07037 NA07048 NA07051 NA07056 NA07347 NA07357 NA10847 NA10851 NA11829 NA11830 NA11831 NA11843 NA11892 NA11893 NA11894 NA11919 NA11920 NA11930 NA11931 NA11932 NA11933 NA11992 NA11993 NA11994 NA11995 NA12003 NA12004 NA12006 NA12043 NA12044 NA12045 NA12046 NA12058 NA12144 NA12154 NA12155 NA12249 NA12272 NA12273 NA12275 NA12282 NA12283 NA12286 NA12287 NA12340 NA12341 NA12342 NA12347 NA12348 NA12383 NA12399 NA12400 NA12413 NA12489 NA12546 NA12716 NA12717 NA12718 NA12748 NA12749 NA12750 NA12751 NA12761 NA12763 NA12775 NA12777 NA12778 NA12812 NA12814 NA12815 NA12827 NA12829 NA12830 NA12842 NA12843 NA12872 NA12873 NA12874 NA12889 NA12890 NA18486 NA18487 NA18489 NA18498 NA18499 NA18501 NA18502 NA18504 NA18505 NA18507 NA18508 NA18510 NA18511 NA18516 NA18517 NA18519 NA18520 NA18522 NA18523 NA18525 NA18526 NA18527 NA18528 NA18530 NA18532 NA18534 NA18535 NA18536 NA18537 NA18538 NA18539 NA18541 NA18542 NA18543 NA18544 NA18545 NA18546 NA18547 NA18548 NA18549 NA18550 NA18552 NA18553 NA18555 NA18557 NA18558 NA18559 NA18560 NA18561 NA18562 NA18563 NA18564 NA18565 NA18566 NA18567 NA18570 NA18571 NA18572 NA18573 NA18574 NA18576 NA18577 NA18579 NA18582 NA18592 NA18593 NA18595 NA18596 NA18597 NA18599 NA18602 NA18603 NA18605 NA18606 NA18608 NA18609 NA18610 NA18611 NA18612 NA18613 NA18614 NA18615 NA18616 NA18617 NA18618 NA18619 NA18620 NA18621 NA18622 NA18623 NA18624 NA18626 NA18627 NA18628 NA18630 NA18631 NA18632 NA18633 NA18634 NA18635 NA18636 NA18637 NA18638 NA18639 NA18640 NA18641 NA18642 NA18643 NA18645 NA18647 NA18740 NA18745 NA18747 NA18748 NA18749 NA18757 NA18853 NA18856 NA18858 NA18861 NA18867 NA18868 NA18870 NA18871 NA18873 NA18874 NA18907 NA18908 NA18909 NA18910 NA18912 NA18916 NA18917 NA18923 NA18924 NA18933 NA18934 NA18939 NA18940 NA18941 NA18942 NA18943 NA18944 NA18945 NA18946 NA18947 NA18948 NA18949 NA18950 NA18951 NA18952 NA18953 NA18954 NA18956 NA18957 NA18959 NA18960 NA18961 NA18962 NA18963 NA18964 NA18965 NA18966 NA18968 NA18971 NA18973 NA18974 NA18975 NA18976 NA18977 NA18978 NA18980 NA18981 NA18982 NA18983 NA18984 NA18985 NA18986 NA18987 NA18988 NA18989 NA18990 NA18992 NA18994 NA18995 NA18998 NA18999 NA19000 NA19002 NA19003 NA19004 NA19005 NA19007 NA19009 NA19010 NA19012 NA19020 NA19028 NA19035 NA19036 NA19038 NA19041 NA19044 NA19046 NA19054 NA19055 NA19056 NA19057 NA19058 NA19059 NA19060 NA19062 NA19063 NA19064 NA19065 NA19066 NA19067 NA19068 NA19070 NA19072 NA19074 NA19075 NA19076 NA19077 NA19078 NA19079 NA19080 NA19081 NA19082 NA19083 NA19084 NA19085 NA19087 NA19088 NA19093 NA19095 NA19096 NA19098 NA19099 NA19102 NA19107 NA19108 NA19113 NA19114 NA19116 NA19117 NA19118 NA19119 NA19121 NA19129 NA19130 NA19131 NA19137 NA19138 NA19146 NA19147 NA19149 NA19150 NA19152 NA19160 NA19171 NA19172 NA19175 NA19185 NA19189 NA19190 NA19197 NA19198 NA19200 NA19204 NA19207 NA19209 NA19213 NA19222 NA19223 NA19225 NA19235 NA19236 NA19247 NA19248 NA19256 NA19257 NA19307 NA19308 NA19309 NA19310 NA19311 NA19312 NA19313 NA19315 NA19316 NA19317 NA19318 NA19319 NA19321 NA19324 NA19327 NA19328 NA19331 NA19332 NA19334 NA19338 NA19346 NA19347 NA19350 NA19351 NA19352 NA19355 NA19359 NA19360 NA19371 NA19372 NA19373 NA19374 NA19375 NA19376 NA19377 NA19379 NA19380 NA19381 NA19382 NA19383 NA19384 NA19385 NA19390 NA19391 NA19393 NA19394 NA19395 NA19396 NA19397 NA19398 NA19399 NA19401 NA19403 NA19404 NA19428 NA19429 NA19430 NA19431 NA19434 NA19435 NA19436 NA19437 NA19438 NA19439 NA19440 NA19443 NA19444 NA19445 NA19446 NA19448 NA19449 NA19451 NA19452 NA19453 NA19455 NA19456 NA19457 NA19461 NA19462 NA19463 NA19466 NA19467 NA19468 NA19469 NA19470 NA19471 NA19472 NA19473 NA19474 NA19625 NA19648 NA19651 NA19652 NA19654 NA19655 NA19657 NA19660 NA19661 NA19663 NA19664 NA19672 NA19675 NA19676 NA19678 NA19679 NA19681 NA19682 NA19684 NA19685 NA19700 NA19701 NA19703 NA19704 NA19707 NA19711 NA19712 NA19713 NA19716 NA19717 NA19719 NA19720 NA19722 NA19723 NA19725 NA19726 NA19728 NA19729 NA19731 NA19732 NA19734 NA19735 NA19737 NA19738 NA19740 NA19741 NA19746 NA19747 NA19749 NA19750 NA19752 NA19753 NA19755 NA19756 NA19758 NA19759 NA19761 NA19762 NA19764 NA19770 NA19771 NA19773 NA19774 NA19776 NA19777 NA19779 NA19780 NA19782 NA19783 NA19785 NA19786 NA19788 NA19789 NA19794 NA19795 NA19818 NA19819 NA19834 NA19835 NA19900 NA19901 NA19904 NA19908 NA19909 NA19914 NA19916 NA19917 NA19920 NA19921 NA19922 NA19923 NA19982 NA19984 NA19985 NA20126 NA20127 NA20276 NA20278 NA20281 NA20282 NA20287 NA20289 NA20291 NA20294 NA20296 NA20298 NA20299 NA20314 NA20317 NA20322 NA20332 NA20334 NA20336 NA20339 NA20340 NA20341 NA20342 NA20344 NA20346 NA20348 NA20351 NA20356 NA20357 NA20359 NA20363 NA20412 NA20414 NA20502 NA20503 NA20504 NA20505 NA20506 NA20507 NA20508 NA20509 NA20510 NA20512 NA20513 NA20515 NA20516 NA20517 NA20518 NA20519 NA20520 NA20521 NA20522 NA20524 NA20525 NA20527 NA20528 NA20529 NA20530 NA20531 NA20532 NA20533 NA20534 NA20535 NA20536 NA20537 NA20538 NA20539 NA20540 NA20541 NA20542 NA20543 NA20544 NA20581 NA20582 NA20585 NA20586 NA20588 NA20589 NA20752 NA20753 NA20754 NA20755 NA20756 NA20757 NA20758 NA20759 NA20760 NA20761 NA20765 NA20766 NA20768 NA20769 NA20770 NA20771 NA20772 NA20773 NA20774 NA20775 NA20778 NA20783 NA20785 NA20786 NA20787 NA20790 NA20792 NA20795 NA20796 NA20797 NA20798 NA20799 NA20800 NA20801 NA20802 NA20803 NA20804 NA20805 NA20806 NA20807 NA20808 NA20809 NA20810 NA20811 NA20812 NA20813 NA20814 NA20815 NA20816 NA20818 NA20819 NA20826 NA20828]

  subpopToPop = {}
  subpopToPopFile = Pathname.new "/groups/park/semin/BiO/Research/GCC/Germline/DGV/Population_Code.tsv"
  subpopToPopFile.each_line do |line|
    next if line.empty? || line.start_with?("#")
    subpop, pop = line.chomp.split("\t")
    subpopToPop[subpop] = pop
  end

  sampleToPop = {}
  tgSampleToPopFile = Pathname.new "/groups/park/semin/BiO/Research/GCC/Germline/DGV/1000_Genomes_Sample_Population_Matching_Table.tsv"
  tgSampleToPopFile.each_line do |line|
    next if line.empty? || line.start_with?("#")
    sample, subpop = line.chomp.split("\t")
    pop = subpopToPop[subpop]
    sampleToPop[sample] = pop
  end

  allCnt, asnCnt, amrCnt, afrCnt, eurCnt = 0, 0, 0, 0, 0
  samplesInTgVcf.each do |sample|
    allCnt += 1
    case sampleToPop[sample]
    when "ASN"
      asnCnt += 1
    when "AMR"
      amrCnt += 1
    when "AFR"
      afrCnt += 1
    when "EUR"
      eurCnt += 1
    end
  end
  #puts "ALL:#{allCnt},ASN:#{asnCnt},AMR:#{amrCnt},AFR:#{afrCnt},EUR:#{eurCnt}"
  return({:ALL => allCnt, :ASN => asnCnt, :AMR => amrCnt, :AFR => afrCnt, :EUR => eurCnt})
end

def get_1k_genomes_phase3_population_counts
  allCnt, easCnt, eurCnt, afrCnt, amrCnt, sasCnt = 0, 0, 0, 0, 0, 0
  popTableFile = $baseDir + "1000Genomes/integrated_call_samples_v3.20130502.ALL.panel"
  popTableFile.readlines[1..-1].each do |line|
    next if line.empty? || line.start_with?("#")
    sample, pop, super_pop, gender = line.chomp.split("\t")
    allCnt += 1
    case super_pop
    when "EAS"
      easCnt += 1
    when "EUR"
      eurCnt += 1
    when "AFR"
      afrCnt += 1
    when "AMR"
      amrCnt += 1
    when "SAS"
      sasCnt += 1
    end
  end
  #puts "ALL:#{allCnt}, EAS:#{easCnt}, EUR:#{eurCnt}, AFR:#{afrCnt}, AMR:#{amrCnt}, SAS:#{sasCnt}"
  return({:ALL => allCnt, :EAS => easCnt, :EUR => eurCnt, :AFR => afrCnt, :AMR => amrCnt, :SAS => sasCnt})
end


class Pathname
  def successful?
    lsfOutTxt = Pathname.new(self).read
    if lsfOutTxt =~ /Successful/m
      return true
    else
      return false
    end
  end
end

def convert_snp_vcf_to_tsv
  #popCnts = get_1k_genomes_phase1_population_counts
  popCnts = get_1k_genomes_phase3_population_counts
  vcfs = Pathname.glob($canVcfDir + "UCEC/snp/chrs/*/*/*cadd_phred.vcf.gz").sort
  #vcfs.each_with_index do |vcf, vi|
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    cancerType = vcf.basename.to_s.split("-")[1]
    cancerType2 = cancerType == "CRC" ? "COADREAD" : cancerType
    tsv = vcf.sub_ext(".tsv")
    lsfOut = tsv.sub_ext(".tsv.lsfout")
    if tsv.exist?
      if lsfOut.successful?
        next
      else
        tsv.delete
        lsfOut.delete
      end
    end
    notHispanicWhitePatientIndices = get_not_hispanic_white_patient_indices_of(cancerType2, vcf)
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TCGA-A2-A0EU-10A-01D-A060
    #6       202777  rs12206347      C       T       594.64  PASS    AC=19;AF=0.086;AN=222;BaseQRankSum=1.051;DB;DP=409;Dels=0.00;FS=4.675;HaplotypeScore=0.2332;InbreedingCoeff=-0.0044;MLEAC=19;MLEAF=0.086;MQ=58.08;MQ0=0;MQRankSum=0.312;POSITIVE_TRAIN_SITE;QD=9.44;ReadPosRankSum=-0.677;VQSLOD=25.27;culprit=FS;DBSNP=12206347;TG_ALL_AF=0.03;TG_AMR_AF=0.03;TG_AFR_AF=0.01;TG_EUR_AF=0.07;CADD_RAW=-0.518325;CADD_PHRED=1.634        GT:AD:DP:GQ:PL  0/1:6,4:10:81:81,0,123
        #-q long -W 700:0 \\
    cmd =<<-CMD
      bsub \\
        -g /gcc/germ/tsv \\
        -q short -W 720 \\
        -R "rusage[mem=20000]" -M 20000000 \\
        -o #{lsfOut} \\
        "zcat #{vcf} |
          ruby -ne 'BEGIN { require \\"rinruby\\"; \\$rr = RinRuby.new(echo = false, interactive = false); }
            if \\$_.start_with?(\\"#\\")
              next
            else
              cols    = \\$_.chomp.split(\\"\\t\\")
              alts    = cols[4].split(\\",\\")
              nalts   = Array.new(3, \\"\\"); nalts[0...alts.size] = alts
              nrsids  = Array.new(3, \\"\\");
              if \\$_ =~ /DBSNP=(\\d+)/
                rsids = \\$1.split(\\",\\") if \\$_ =~ /DBSNP=(\\d+)/
                nrsids[0...rsids.size] = rsids
              end
              ntg_all_afs   = Array.new(3, \\"\\")
              if \\$_ =~ /TG_ALL_AF=(\\S+?)(;|\\s)/
                tg_all_afs = \\$1.split(\\",\\") if \\$_ =~ /TG_ALL_AF=(\\S+?)(;|\\s)/
                ntg_all_afs[0...tg_all_afs.size] = tg_all_afs
              end
              ntg_eas_afs   = Array.new(3, \\"\\")
              if \\$_ =~ /TG_EAS_AF=(\\S+?)(;|\\s)/
                tg_eas_afs = \\$1.split(\\",\\") if \\$_ =~ /TG_EAS_AF=(\\S+?)(;|\\s)/
                ntg_eas_afs[0...tg_eas_afs.size] = tg_eas_afs
              end
              ntg_eur_afs   = Array.new(3, \\"\\")
              if \\$_ =~ /TG_EUR_AF=(\\S+?)(;|\\s)/
                tg_eur_afs = \\$1.split(\\",\\") if \\$_ =~ /TG_EUR_AF=(\\S+?)(;|\\s)/
                ntg_eur_afs[0...tg_eur_afs.size] = tg_eur_afs
              end
              ntg_afr_afs   = Array.new(3, \\"\\")
              if \\$_ =~ /TG_AFR_AF=(\\S+?)(;|\\s)/
                tg_afr_afs = \\$1.split(\\",\\") if \\$_ =~ /TG_AFR_AF=(\\S+?)(;|\\s)/
                ntg_afr_afs[0...tg_afr_afs.size] = tg_afr_afs
              end
              ntg_amr_afs   = Array.new(3, \\"\\")
              if \\$_ =~ /TG_AMR_AF=(\\S+?)(;|\\s)/
                tg_amr_afs = \\$1.split(\\",\\") if \\$_ =~ /TG_AMR_AF=(\\S+?)(;|\\s)/
                ntg_amr_afs[0...tg_amr_afs.size] = tg_amr_afs
              end
              ntg_sas_afs   = Array.new(3, \\"\\")
              if \\$_ =~ /TG_SAS_AF=(\\S+?)(;|\\s)/
                tg_sas_afs = \\$1.split(\\",\\") if \\$_ =~ /TG_SAS_AF=(\\S+?)(;|\\s)/
                ntg_sas_afs[0...tg_sas_afs.size] = tg_sas_afs
              end
              ncadd_raws    = Array.new(3, \\"\\")
              if \\$_ =~ /CADD_RAW=(\\S+?)(;|\\s)/
                cadd_raws = \\$1.split(\\",\\")
                ncadd_raws[0...cadd_raws.size] = cadd_raws
              end
              ncadd_phreds    = Array.new(3, \\"\\")
              if \\$_ =~ /CADD_PHRED=(\\S+?)(;|\\s)/
                cadd_phreds = \\$1.split(\\",\\")
                ncadd_phreds[0...cadd_phreds.size] = cadd_phreds
              end
              an            = \\$1.split(\\",\\") if \\$_ =~ /AN=(\\S+?)(;|\\s)/
              acs           = \\$1.split(\\",\\") if \\$_ =~ /AC=(\\S+?)(;|\\s)/
              nacs          = Array.new(3, \\"\\"); nacs[0...acs.size] = acs
              afs           = \\$1.split(\\",\\") if \\$_ =~ /AF=(\\S+?)(;|\\s)/
              nafs          = Array.new(3, \\"\\"); nafs[0...afs.size] = afs
              all_as        = cols[9..-1].map { |x| x.split(\\":\\")[0].split(\\"/\\") }.flatten
              all_an        = all_as.select { |a| a != \\".\\" }.size
              all_ac1       = all_as.select { |a| a == \\"1\\" }.size
              all_ac2       = all_as.select { |a| a == \\"2\\" }.size
              all_ac3       = all_as.select { |a| a == \\"3\\" }.size
              all_af1       = all_ac1 / all_an.to_f
              all_af2       = all_ac2 / all_an.to_f
              all_af3       = all_ac3 / all_an.to_f
              eur_as        = cols[9..-1].values_at(#{notHispanicWhitePatientIndices.join(",")}).map { |x| x.split(\\":\\")[0].split(\\"/\\") }.flatten
              eur_an        = eur_as.select { |a| a != \\".\\" }.size
              eur_ac1       = eur_as.select { |a| a == \\"1\\" }.size
              eur_ac2       = eur_as.select { |a| a == \\"2\\" }.size
              eur_ac3       = eur_as.select { |a| a == \\"3\\" }.size
              eur_af1       = eur_ac1 / eur_an.to_f
              eur_af2       = eur_ac2 / eur_an.to_f
              eur_af3       = eur_ac3 / eur_an.to_f
              \\$rr.eval <<-EOF
                chisqTestPvalAllAc1   = \\"\\"
                chisqTestPvalAllAc2   = \\"\\"
                chisqTestPvalAllAc3   = \\"\\"
                fisherTestPvalAllAc1  = \\"\\"
                fisherTestPvalAllAc2  = \\"\\"
                fisherTestPvalAllAc3  = \\"\\"
                chisqTestPvalEurAc1   = \\"\\"
                chisqTestPvalEurAc2   = \\"\\"
                chisqTestPvalEurAc3   = \\"\\"
                fisherTestPvalEurAc1  = \\"\\"
                fisherTestPvalEurAc2  = \\"\\"
                fisherTestPvalEurAc3  = \\"\\"
              EOF
              if (ntg_all_afs[0] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableAllAc1       = as.table(rbind(c(\#{all_ac1}, \#{all_an - all_ac1}),
                                                        c(as.integer(2 * #{popCnts[:ALL]} * \#{ntg_all_afs[0]}),
                                                          as.integer(2 * #{popCnts[:ALL]} - 2 * #{popCnts[:ALL]} * \#{ntg_all_afs[0]}))))
                  chisqTestAllAc1       = chisq.test(contTableAllAc1, simulate.p.value = TRUE)
                  chisqTestPvalAllAc1   = chisqTestAllAc1\\$p.value
                  fisherTestAllAc1      = fisher.test(contTableAllAc1)
                  fisherTestPvalAllAc1  = fisherTestAllAc1\\$p.value
                EOF
              end
              if (ntg_all_afs[1] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableAllAc2       = as.table(rbind(c(\#{all_ac2}, \#{all_an - all_ac2}),
                                                        c(as.integer(2 * #{popCnts[:ALL]} * \#{ntg_all_afs[1]}),
                                                          as.integer(2 * #{popCnts[:ALL]} - 2 * #{popCnts[:ALL]} * \#{ntg_all_afs[1]}))))
                  chisqTestAllAc2       = chisq.test(contTableAllAc2, simulate.p.value = TRUE)
                  chisqTestPvalAllAc2   = chisqTestAllAc2\\$p.value
                  fisherTestAllAc2      = fisher.test(contTableAllAc2)
                  fisherTestPvalAllAc2  = fisherTestAllAc2\\$p.value
                EOF
              end
              if (ntg_all_afs[2] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableAllAc3       = as.table(rbind(c(\#{all_ac3}, \#{all_an - all_ac3}),
                                                        c(as.integer(2 * #{popCnts[:ALL]} * \#{ntg_all_afs[2]}),
                                                          as.integer(2 * #{popCnts[:ALL]} - 2 * #{popCnts[:ALL]} * \#{ntg_all_afs[2]}))))
                  chisqTestAllAc3       = chisq.test(contTableAllAc3, simulate.p.value = TRUE)
                  chisqTestPvalAllAc3   = chisqTestAllAc3\\$p.value
                  fisherTestAllAc3      = fisher.test(contTableAllAc3)
                  fisherTestPvalAllAc3  = fisherTestAllAc3\\$p.value
                EOF
              end
              if (ntg_eur_afs[0] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableEurAc1       = as.table(rbind(c(\#{eur_ac1}, \#{eur_an - eur_ac1}),
                                                        c(as.integer(2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[0]}),
                                                          as.integer(2 * #{popCnts[:EUR]} - 2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[0]}))))
                  chisqTestEurAc1       = chisq.test(contTableEurAc1, simulate.p.value = TRUE)
                  chisqTestPvalEurAc1   = chisqTestEurAc1\\$p.value
                  fisherTestEurAc1      = fisher.test(contTableEurAc1)
                  fisherTestPvalEurAc1  = fisherTestEurAc1\\$p.value
                EOF
              end
              if (ntg_eur_afs[1] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableEurAc2       = as.table(rbind(c(\#{eur_ac2}, \#{eur_an - eur_ac2}),
                                                        c(as.integer(2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[1]}),
                                                          as.integer(2 * #{popCnts[:EUR]} - 2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[1]}))))
                  chisqTestEurAc2       = chisq.test(contTableEurAc2, simulate.p.value = TRUE)
                  chisqTestPvalEurAc2   = chisqTestEurAc2\\$p.value
                  fisherTestEurAc2      = fisher.test(contTableEurAc2)
                  fisherTestPvalEurAc2  = fisherTestEurAc2\\$p.value
                EOF
              end
              if (ntg_eur_afs[2] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableEurAc3       = as.table(rbind(c(\#{eur_ac3}, \#{eur_an - eur_ac3}),
                                                        c(as.integer(2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[2]}),
                                                          as.integer(2 * #{popCnts[:EUR]} - 2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[2]}))))
                  chisqTestEurAc3       = chisq.test(contTableEurAc3, simulate.p.value = TRUE)
                  chisqTestPvalEurAc3   = chisqTestEurAc3\\$p.value
                  fisherTestEurAc3      = fisher.test(contTableEurAc3)
                  fisherTestPvalEurAc3  = fisherTestEurAc3\\$p.value
                EOF
              end
              puts([cols[0..5],
                    nalts, nrsids, 
                    ntg_all_afs, ntg_eas_afs, ntg_eur_afs, ntg_afr_afs, ntg_amr_afs, ntg_sas_afs,
                    ncadd_raws, ncadd_phreds,
                    an, nacs, nafs,
                    all_an, eur_an,
                    all_ac1, all_ac2, all_ac3,
                    all_af1, all_af2, all_af3,
                    eur_ac1, eur_ac2, eur_ac3,
                    eur_af1, eur_af2, eur_af3,
                    \\$rr.chisqTestPvalAllAc1,
                    \\$rr.chisqTestPvalAllAc2,
                    \\$rr.chisqTestPvalAllAc3,
                    \\$rr.fisherTestPvalAllAc1,
                    \\$rr.fisherTestPvalAllAc2,
                    \\$rr.fisherTestPvalAllAc3,
                    \\$rr.chisqTestPvalEurAc1,
                    \\$rr.chisqTestPvalEurAc2,
                    \\$rr.chisqTestPvalEurAc3,
                    \\$rr.fisherTestPvalEurAc1,
                    \\$rr.fisherTestPvalEurAc2,
                    \\$rr.fisherTestPvalEurAc3,
                    cols[9..-1].map { |x| x.split(\\":\\")[0] }
                    ].flatten.join(\\"\\t\\"))
            end
          ' > #{tsv}"
    CMD
    submit cmd
  end
end

def convert_indel_vcf_to_tsv
  #popCnts = get_1k_genomes_phase1_population_counts
  popCnts = get_1k_genomes_phase3_population_counts
  vcfs = Pathname.glob($canVcfDir + "UCEC/indel/chrs/*/*/*tgSAS.vcf.gz").sort
  #vcfs.each_with_index do |vcf, vi|
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    cancerType = vcf.basename.to_s.split("-")[1]
    cancerType2 = cancerType == "CRC" ? "COADREAD" : cancerType
    tsv = vcf.sub_ext(".tsv")
    lsfOut = tsv.sub_ext(".tsv.lsfout")
    if tsv.exist?
      if lsfOut.successful?
        next
      else
        tsv.delete
        lsfOut.delete
      end
    end
    notHispanicWhitePatientIndices = get_not_hispanic_white_patient_indices_of(cancerType2, vcf)
    #CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  TCGA-A2-A0EU-10A-01D-A060
    #3       66519   rs201594994     G       GA      400.97  PASS    AC=15;AF=0.069;AN=216;BaseQRankSum=-0.070;DB;DP=397;FS=1.765;InbreedingCoeff=-0.0322;MLEAC=15;MLEAF=0.069;MQ=55.08;MQ0=0;MQRankSum=0.086;QD=8.02;RPA=4,5;RU=A;ReadPosRankSum=0.356;STR;VQSLOD=5.90;culprit=DP;DBSNP=201594994;TG_ALL_AF=0.03;TG_ASN_AF=0.0017;TG_AMR_AF=0.04;TG_AFR_AF=0.03;TG_EUR_AF=0.05      GT:AD:DP:GQ:PL  0/0:3,0:3:9:0,9,88
    cmd =<<-CMD
      bsub \\
        -q short -W 12:0 \\
        -g /gcc/germ/tsv \\
        -R "rusage[mem=20000]" -M 20000000 \\
        -o #{lsfOut} \\
        "zcat #{vcf} |
          ruby -ne 'BEGIN { require \\"rinruby\\"; \\$rr = RinRuby.new(echo = false, interactive = false); }
            if \\$_.start_with?(\\"#\\")
              next
            else
              cols    = \\$_.chomp.split(\\"\\t\\")
              alts    = cols[4].split(\\",\\")
              nalts   = Array.new(6, \\"\\"); nalts[0...alts.size] = alts
              nrsids  = Array.new(6, \\"\\");
              if \\$_ =~ /DBSNP=(\\d+)/
                rsids = \\$1.split(\\",\\") if \\$_ =~ /DBSNP=(\\d+)/
                nrsids[0...rsids.size] = rsids
              end
              ntg_all_afs   = Array.new(6, \\"\\")
              if \\$_ =~ /TG_ALL_AF=(\\S+?)(;|\\s)/
                tg_all_afs = \\$1.split(\\",\\") if \\$_ =~ /TG_ALL_AF=(\\S+?)(;|\\s)/
                ntg_all_afs[0...tg_all_afs.size] = tg_all_afs
              end
              ntg_eas_afs   = Array.new(6, \\"\\")
              if \\$_ =~ /TG_EAS_AF=(\\S+?)(;|\\s)/
                tg_eas_afs = \\$1.split(\\",\\") if \\$_ =~ /TG_EAS_AF=(\\S+?)(;|\\s)/
                ntg_eas_afs[0...tg_eas_afs.size] = tg_eas_afs
              end
              ntg_eur_afs   = Array.new(6, \\"\\")
              if \\$_ =~ /TG_EUR_AF=(\\S+?)(;|\\s)/
                tg_eur_afs = \\$1.split(\\",\\") if \\$_ =~ /TG_EUR_AF=(\\S+?)(;|\\s)/
                ntg_eur_afs[0...tg_eur_afs.size] = tg_eur_afs
              end
              ntg_afr_afs   = Array.new(6, \\"\\")
              if \\$_ =~ /TG_AFR_AF=(\\S+?)(;|\\s)/
                tg_afr_afs = \\$1.split(\\",\\") if \\$_ =~ /TG_AFR_AF=(\\S+?)(;|\\s)/
                ntg_afr_afs[0...tg_afr_afs.size] = tg_afr_afs
              end
              ntg_amr_afs   = Array.new(6, \\"\\")
              if \\$_ =~ /TG_AMR_AF=(\\S+?)(;|\\s)/
                tg_amr_afs = \\$1.split(\\",\\") if \\$_ =~ /TG_AMR_AF=(\\S+?)(;|\\s)/
                ntg_amr_afs[0...tg_amr_afs.size] = tg_amr_afs
              end
              ntg_sas_afs   = Array.new(6, \\"\\")
              if \\$_ =~ /TG_SAS_AF=(\\S+?)(;|\\s)/
                tg_sas_afs = \\$1.split(\\",\\") if \\$_ =~ /TG_SAS_AF=(\\S+?)(;|\\s)/
                ntg_sas_afs[0...tg_sas_afs.size] = tg_sas_afs
              end
              an            = \\$1.split(\\",\\") if \\$_ =~ /AN=(\\S+?)(;|\\s)/
              acs           = \\$1.split(\\",\\") if \\$_ =~ /AC=(\\S+?)(;|\\s)/
              nacs          = Array.new(6, \\"\\"); nacs[0...acs.size] = acs
              afs           = \\$1.split(\\",\\") if \\$_ =~ /AF=(\\S+?)(;|\\s)/
              nafs          = Array.new(6, \\"\\"); nafs[0...afs.size] = afs
              all_as        = cols[9..-1].map { |x| x.split(\\":\\")[0].split(\\"/\\") }.flatten
              all_an        = all_as.select { |a| a != \\".\\" }.size
              all_ac1       = all_as.select { |a| a == \\"1\\" }.size
              all_ac2       = all_as.select { |a| a == \\"2\\" }.size
              all_ac3       = all_as.select { |a| a == \\"3\\" }.size
              all_ac4       = all_as.select { |a| a == \\"4\\" }.size
              all_ac5       = all_as.select { |a| a == \\"5\\" }.size
              all_ac6       = all_as.select { |a| a == \\"6\\" }.size
              all_af1       = all_ac1 / all_an.to_f
              all_af2       = all_ac2 / all_an.to_f
              all_af3       = all_ac3 / all_an.to_f
              all_af4       = all_ac4 / all_an.to_f
              all_af5       = all_ac5 / all_an.to_f
              all_af6       = all_ac6 / all_an.to_f
              eur_as        = cols[9..-1].values_at(#{notHispanicWhitePatientIndices.join(",")}).map { |x| x.split(\\":\\")[0].split(\\"/\\") }.flatten
              eur_an        = eur_as.select { |a| a != \\".\\" }.size
              eur_ac1       = eur_as.select { |a| a == \\"1\\" }.size
              eur_ac2       = eur_as.select { |a| a == \\"2\\" }.size
              eur_ac3       = eur_as.select { |a| a == \\"3\\" }.size
              eur_ac4       = eur_as.select { |a| a == \\"4\\" }.size
              eur_ac5       = eur_as.select { |a| a == \\"5\\" }.size
              eur_ac6       = eur_as.select { |a| a == \\"6\\" }.size
              eur_af1       = eur_ac1 / eur_an.to_f
              eur_af2       = eur_ac2 / eur_an.to_f
              eur_af3       = eur_ac3 / eur_an.to_f
              eur_af4       = eur_ac4 / eur_an.to_f
              eur_af5       = eur_ac5 / eur_an.to_f
              eur_af6       = eur_ac6 / eur_an.to_f
              \\$rr.eval <<-EOF
                chisqTestPvalAllAc1   = \\"\\"
                chisqTestPvalAllAc2   = \\"\\"
                chisqTestPvalAllAc3   = \\"\\"
                chisqTestPvalAllAc4   = \\"\\"
                chisqTestPvalAllAc5   = \\"\\"
                chisqTestPvalAllAc6   = \\"\\"
                fisherTestPvalAllAc1  = \\"\\"
                fisherTestPvalAllAc2  = \\"\\"
                fisherTestPvalAllAc3  = \\"\\"
                fisherTestPvalAllAc4  = \\"\\"
                fisherTestPvalAllAc5  = \\"\\"
                fisherTestPvalAllAc6  = \\"\\"
                chisqTestPvalEurAc1   = \\"\\"
                chisqTestPvalEurAc2   = \\"\\"
                chisqTestPvalEurAc3   = \\"\\"
                chisqTestPvalEurAc4   = \\"\\"
                chisqTestPvalEurAc5   = \\"\\"
                chisqTestPvalEurAc6   = \\"\\"
                fisherTestPvalEurAc1  = \\"\\"
                fisherTestPvalEurAc2  = \\"\\"
                fisherTestPvalEurAc3  = \\"\\"
                fisherTestPvalEurAc4  = \\"\\"
                fisherTestPvalEurAc5  = \\"\\"
                fisherTestPvalEurAc6  = \\"\\"
              EOF
              if (ntg_all_afs[0] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableAllAc1       = as.table(rbind(c(\#{all_ac1}, \#{all_an - all_ac1}),
                                                        c(as.integer(2 * #{popCnts[:ALL]} * \#{ntg_all_afs[0]}),
                                                          as.integer(2 * #{popCnts[:ALL]} - 2 * #{popCnts[:ALL]} * \#{ntg_all_afs[0]}))))
                  chisqTestAllAc1       = chisq.test(contTableAllAc1, simulate.p.value = TRUE)
                  chisqTestPvalAllAc1   = chisqTestAllAc1\\$p.value
                  fisherTestAllAc1      = fisher.test(contTableAllAc1)
                  fisherTestPvalAllAc1  = fisherTestAllAc1\\$p.value
                EOF
              end
              if (ntg_all_afs[1] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableAllAc2       = as.table(rbind(c(\#{all_ac2}, \#{all_an - all_ac2}),
                                                        c(as.integer(2 * #{popCnts[:ALL]} * \#{ntg_all_afs[1]}),
                                                          as.integer(2 * #{popCnts[:ALL]} - 2 * #{popCnts[:ALL]} * \#{ntg_all_afs[1]}))))
                  chisqTestAllAc2       = chisq.test(contTableAllAc2, simulate.p.value = TRUE)
                  chisqTestPvalAllAc2   = chisqTestAllAc2\\$p.value
                  fisherTestAllAc2      = fisher.test(contTableAllAc2)
                  fisherTestPvalAllAc2  = fisherTestAllAc2\\$p.value
                EOF
              end
              if (ntg_all_afs[2] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableAllAc3       = as.table(rbind(c(\#{all_ac3}, \#{all_an - all_ac3}),
                                                        c(as.integer(2 * #{popCnts[:ALL]} * \#{ntg_all_afs[2]}),
                                                          as.integer(2 * #{popCnts[:ALL]} - 2 * #{popCnts[:ALL]} * \#{ntg_all_afs[2]}))))
                  chisqTestAllAc3       = chisq.test(contTableAllAc3, simulate.p.value = TRUE)
                  chisqTestPvalAllAc3   = chisqTestAllAc3\\$p.value
                  fisherTestAllAc3      = fisher.test(contTableAllAc3)
                  fisherTestPvalAllAc3  = fisherTestAllAc3\\$p.value
                EOF
              end
              if (ntg_all_afs[3] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableAllAc4       = as.table(rbind(c(\#{all_ac4}, \#{all_an - all_ac4}),
                                                        c(as.integer(2 * #{popCnts[:ALL]} * \#{ntg_all_afs[3]}),
                                                          as.integer(2 * #{popCnts[:ALL]} - 2 * #{popCnts[:ALL]} * \#{ntg_all_afs[3]}))))
                  chisqTestAllAc4       = chisq.test(contTableAllAc4, simulate.p.value = TRUE)
                  chisqTestPvalAllAc4   = chisqTestAllAc4\\$p.value
                  fisherTestAllAc4      = fisher.test(contTableAllAc4)
                  fisherTestPvalAllAc4  = fisherTestAllAc4\\$p.value
                EOF
              end
              if (ntg_all_afs[4] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableAllAc5       = as.table(rbind(c(\#{all_ac5}, \#{all_an - all_ac5}),
                                                        c(as.integer(2 * #{popCnts[:ALL]} * \#{ntg_all_afs[4]}),
                                                          as.integer(2 * #{popCnts[:ALL]} - 2 * #{popCnts[:ALL]} * \#{ntg_all_afs[4]}))))
                  chisqTestAllAc5       = chisq.test(contTableAllAc5, simulate.p.value = TRUE)
                  chisqTestPvalAllAc5   = chisqTestAllAc5\\$p.value
                  fisherTestAllAc5      = fisher.test(contTableAllAc5)
                  fisherTestPvalAllAc5  = fisherTestAllAc5\\$p.value
                EOF
              end
              if (ntg_all_afs[5] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableAllAc6       = as.table(rbind(c(\#{all_ac6}, \#{all_an - all_ac6}),
                                                        c(as.integer(2 * #{popCnts[:ALL]} * \#{ntg_all_afs[5]}),
                                                          as.integer(2 * #{popCnts[:ALL]} - 2 * #{popCnts[:ALL]} * \#{ntg_all_afs[5]}))))
                  chisqTestAllAc6       = chisq.test(contTableAllAc6, simulate.p.value = TRUE)
                  chisqTestPvalAllAc6   = chisqTestAllAc6\\$p.value
                  fisherTestAllAc6      = fisher.test(contTableAllAc6)
                  fisherTestPvalAllAc6  = fisherTestAllAc6\\$p.value
                EOF
              end
              if (ntg_eur_afs[0] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableEurAc1       = as.table(rbind(c(\#{eur_ac1}, \#{eur_an - eur_ac1}),
                                                        c(as.integer(2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[0]}),
                                                          as.integer(2 * #{popCnts[:EUR]} - 2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[0]}))))
                  chisqTestEurAc1       = chisq.test(contTableEurAc1, simulate.p.value = TRUE)
                  chisqTestPvalEurAc1   = chisqTestEurAc1\\$p.value
                  fisherTestEurAc1      = fisher.test(contTableEurAc1)
                  fisherTestPvalEurAc1  = fisherTestEurAc1\\$p.value
                EOF
              end
              if (ntg_eur_afs[1] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableEurAc2       = as.table(rbind(c(\#{eur_ac2}, \#{eur_an - eur_ac2}),
                                                        c(as.integer(2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[1]}),
                                                          as.integer(2 * #{popCnts[:EUR]} - 2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[1]}))))
                  chisqTestEurAc2       = chisq.test(contTableEurAc2, simulate.p.value = TRUE)
                  chisqTestPvalEurAc2   = chisqTestEurAc2\\$p.value
                  fisherTestEurAc2      = fisher.test(contTableEurAc2)
                  fisherTestPvalEurAc2  = fisherTestEurAc2\\$p.value
                EOF
              end
              if (ntg_eur_afs[2] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableEurAc3       = as.table(rbind(c(\#{eur_ac3}, \#{eur_an - eur_ac3}),
                                                        c(as.integer(2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[2]}),
                                                          as.integer(2 * #{popCnts[:EUR]} - 2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[2]}))))
                  chisqTestEurAc3       = chisq.test(contTableEurAc3, simulate.p.value = TRUE)
                  chisqTestPvalEurAc3   = chisqTestEurAc3\\$p.value
                  fisherTestEurAc3      = fisher.test(contTableEurAc3)
                  fisherTestPvalEurAc3  = fisherTestEurAc3\\$p.value
                EOF
              end
              if (ntg_eur_afs[3] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableEurAc4       = as.table(rbind(c(\#{eur_ac4}, \#{eur_an - eur_ac4}),
                                                        c(as.integer(2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[3]}),
                                                          as.integer(2 * #{popCnts[:EUR]} - 2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[3]}))))
                  chisqTestEurAc4       = chisq.test(contTableEurAc4, simulate.p.value = TRUE)
                  chisqTestPvalEurAc4   = chisqTestEurAc4\\$p.value
                  fisherTestEurAc4      = fisher.test(contTableEurAc4)
                  fisherTestPvalEurAc4  = fisherTestEurAc4\\$p.value
                EOF
              end
              if (ntg_eur_afs[4] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableEurAc5       = as.table(rbind(c(\#{eur_ac5}, \#{eur_an - eur_ac5}),
                                                        c(as.integer(2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[4]}),
                                                          as.integer(2 * #{popCnts[:EUR]} - 2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[4]}))))
                  chisqTestEurAc5       = chisq.test(contTableEurAc5, simulate.p.value = TRUE)
                  chisqTestPvalEurAc5   = chisqTestEurAc5\\$p.value
                  fisherTestEurAc5      = fisher.test(contTableEurAc5)
                  fisherTestPvalEurAc5  = fisherTestEurAc5\\$p.value
                EOF
              end
              if (ntg_eur_afs[5] != \\"\\")
                \\$rr.eval <<-EOF
                  contTableEurAc6       = as.table(rbind(c(\#{eur_ac6}, \#{eur_an - eur_ac6}),
                                                        c(as.integer(2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[5]}),
                                                          as.integer(2 * #{popCnts[:EUR]} - 2 * #{popCnts[:EUR]} * \#{ntg_eur_afs[5]}))))
                  chisqTestEurAc6       = chisq.test(contTableEurAc6, simulate.p.value = TRUE)
                  chisqTestPvalEurAc6   = chisqTestEurAc6\\$p.value
                  fisherTestEurAc6      = fisher.test(contTableEurAc6)
                  fisherTestPvalEurAc6  = fisherTestEurAc6\\$p.value
                EOF
              end
              puts([cols[0..5],
                    nalts, nrsids,
                    ntg_all_afs, ntg_eas_afs, ntg_eur_afs, ntg_afr_afs, ntg_amr_afs, ntg_sas_afs,
                    an, nacs, nafs,
                    all_an, eur_an,
                    all_ac1, all_ac2, all_ac3, all_ac4, all_ac5, all_ac6,
                    all_af1, all_af2, all_af3, all_af4, all_af5, all_af6,
                    eur_ac1, eur_ac2, eur_ac3, eur_ac4, eur_ac5, eur_ac6,
                    eur_af1, eur_af2, eur_af3, eur_af4, eur_af5, eur_af6,
                    \\$rr.chisqTestPvalAllAc1,
                    \\$rr.chisqTestPvalAllAc2,
                    \\$rr.chisqTestPvalAllAc3,
                    \\$rr.chisqTestPvalAllAc4,
                    \\$rr.chisqTestPvalAllAc5,
                    \\$rr.chisqTestPvalAllAc6,
                    \\$rr.fisherTestPvalAllAc1,
                    \\$rr.fisherTestPvalAllAc2,
                    \\$rr.fisherTestPvalAllAc3,
                    \\$rr.fisherTestPvalAllAc4,
                    \\$rr.fisherTestPvalAllAc5,
                    \\$rr.fisherTestPvalAllAc6,
                    \\$rr.chisqTestPvalEurAc1,
                    \\$rr.chisqTestPvalEurAc2,
                    \\$rr.chisqTestPvalEurAc3,
                    \\$rr.chisqTestPvalEurAc4,
                    \\$rr.chisqTestPvalEurAc5,
                    \\$rr.chisqTestPvalEurAc6,
                    \\$rr.fisherTestPvalEurAc1,
                    \\$rr.fisherTestPvalEurAc2,
                    \\$rr.fisherTestPvalEurAc3,
                    \\$rr.fisherTestPvalEurAc4,
                    \\$rr.fisherTestPvalEurAc5,
                    \\$rr.fisherTestPvalEurAc6,
                    cols[9..-1].map { |x| x.split(\\":\\")[0] }
                    ].flatten.join(\\"\\t\\"))
            end
          ' > #{tsv}"
    CMD
    submit cmd
  end
end

def create_annovar_input
  vcfs = Pathname.glob($canVcfDir + "LUAD/{snp,indel}/chrs/*/*/*mss20.vcf.gz").sort
  #vcfs.each_with_index do |vcf, vi|
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    annovarIn = vcf.sub_ext(".av")
    lsfOut = annovarIn.sub_ext(".av.lsfout")
    #1       10108   10108   C       T       1       10108   .       C       T       29.58   PASS
    cmd =<<-CMD
      bsub \\
        -g /gcc/germ/annovar \\
        -q mini \\
        -g /gcc/germ/annovar \\
        -o #{lsfOut} \\
        "zcat #{vcf} | perl #{$convertBin} --format vcf4old --includeinfo stdin | awk -F $'\\\\t' 'BEGIN {OFS = FS} { print \\$1, \\$2, \\$3, \\$4, \\$5, \\$7, \\$9, \\$10 }' > #{annovarIn}"
    CMD
    submit cmd
  end
end

def run_annovar_gene
  annovarIns  = Pathname.glob($canVcfDir + "LUAD/{snp,indel}/*.av").sort
  annovarIns.each_with_index do |annovarIn, vi|
    lsfOut = annovarIn.sub_ext(".av.refgene.lsfout")
    cmd =<<-CMD
        bsub \\
          -q park_12h \\
          -g /gcc/germ/annovar \\
          -R "rusage[mem=10000]" -M 10000000 \\
          -o #{lsfOut} \\
          perl #{$annovarBin} \\
            --batchsize 1000000000 \\
            --memfree 1000000 \\
            --memtotal 1000000 \\
            --comment \\
            --geneanno \\
            --dbtype refgene \\
            --buildver hg19 \\
            #{annovarIn} \\
            #{$annovarDbDir}
    CMD
    submit cmd
  end
end

def run_annovar_filter
  annovarDbs  = %w[ snp137
                    snp138
                    1000g2012apr_all
                    1000g2012apr_eur
                    1000g2012apr_afr
                    1000g2012apr_amr
                    1000g2012apr_asn
                    esp6500si_all
                    esp6500si_ea
                    esp6500si_aa
                    cosmic68
                    ljb23_sift
                    ljb23_pp2hdiv
                    ljb23_pp2hvar
                    ljb23_ma
                    ljb23_mt
                    ]
  annovarIns  = Pathname.glob($canVcfDir + "LUAD/{snp,indel}/*.av").sort
  annovarIns.each_with_index do |annovarIn, vi|
    annovarDbs.each do |db|
      lsfOut = annovarIn.sub_ext(".av.filter.#{db}.lsfout")
          #-q short -W 12:0 \\
      cmd =<<-CMD
        bsub \\
          -q park_12h \\
          -g /gcc/germ/annovar \\
          -R "rusage[mem=10000]" -M 10000000 \\
          -o #{lsfOut} \\
          perl #{$annovarBin} \\
            --batchsize 1000000000 \\
            --memfree 1000000 \\
            --memtotal 1000000 \\
            --comment \\
            --filter \\
            --dbtype #{db} \\
            --buildver hg19 \\
            #{annovarIn} \\
            #{$annovarDbDir}
      CMD
      submit cmd
    end
  end
end

def run_annovar_region
  #annovarDbs  = %w[ 
                    #gwascatalog
                    #tfbsConsSites
                    #wgRna
                    #targetScanS
                    #dgvMerged
                    #wgEncodeRegDnaseClusteredV2
                    #wgEncodeRegTfbsClusteredV3
                    #wgEncodeBroadHmmGm12878HMM
                    #wgEncodeBroadHmmH1hescHMM
                    #wgEncodeBroadHmmHepg2HMM
                    #wgEncodeBroadHmmHmecHMM
                    #wgEncodeBroadHmmHsmmHMM
                    #wgEncodeBroadHmmHuvecHMM
                    #wgEncodeBroadHmmK562HMM
                    #wgEncodeBroadHmmNhekHMM
                    #wgEncodeBroadHmmNhlfHMM
                    #phastConsElements46way
                    #phastConsElements46wayPlacental
                    #phastConsElements46wayPrimates
                    #rmsk_mod
                    #simpleRepeat
                    #nestedRepeats
                    #microsat
                    #genomicSuperDups
                    #]
  annovarDbs  = %w[ 
                    genomicSuperDups
                    ]
  annovarIns  = Pathname.glob($canVcfDir + "CRC/{snp,indel}/*.av").sort
  annovarIns.each_with_index do |annovarIn, vi|
    annovarDbs.each do |db|
      lsfOut = annovarIn.sub_ext(".av.region.#{db}.lsfout")
      cmd =<<-CMD
        bsub \\
          -q short -W 12:0 \\
          -g /gcc/germ/annovar \\
          -R "rusage[mem=10000]" -M 10000000 \\
          -o #{lsfOut} \\
          perl #{$annovarBin} \\
            --batchsize 1000000000 \\
            --memfree 100000 \\
            --memtotal 100000 \\
            --comment \\
            --regionanno \\
            --dbtype #{db} \\
            --buildver hg19 \\
            #{annovarIn} \\
            #{$annovarDbDir}
      CMD
      submit cmd
    end
  end
end

def run_annovar_dhs_gene
  annovarDbs = %w[combinedDhsPeaks
                  distalDhsToPromoterDhs
                  dhsToGeneExpression
                  promoterDhsMasterKnown
                  promoterDhsMasterNovel]
  annovarIns = Pathname.glob($canVcfDir + "LUAD/{snp,indel}/*.av").sort
  annovarIns.each_with_index do |annovarIn, vi|
    annovarDbs.each do |db|
      lsfOut = annovarIn.sub_ext(".av.region.#{db}.lsfout")
      cmd =<<-CMD
        bsub \\
          -q short -W 12:0 \\
          -g /gcc/germ/annovar \\
          -R "rusage[mem=10000]" -M 10000000 \\
          -o #{lsfOut} \\
          perl #{$annovarBin} \\
            --batchsize 1000000000 \\
            --memfree 1000000 \\
            --memtotal 1000000 \\
            --comment \\
            --otherinfo \\
            --regionanno \\
            --dbtype #{db} \\
            --buildver hg19 \\
            #{annovarIn} \\
            #{$annovarDbDir}
      CMD
      submit cmd
    end
  end
end

def extract_cadd_scores
  vcfs = Pathname.glob($canVcfDir + "LUAD/*cadd.vcf.gz")
  vcfs.each_with_index do |vcf, vi|
    next if vcf.dirname.basename.to_s != "CRC"
    #CADD        0.42    1       100000012       100000012       G       T       100000012       G       T
    #1       58814   rs114420996     G       A       313.97  PASS    AC=15;AF=0.068;AN=222;BaseQRankSum=4.195;DB;DP=465;Dels=0.00;FS=0.000;HaplotypeScore=0.2083;InbreedingCoeff=0.0996;MLEAC=13;MLEAF=0.059;MQ=21.43;MQ0=16;MQRankSum=2.055;QD=6.16;ReadPosRankSum=-0.446;VQSLOD=22.86;culprit=FS;CADD=-0.906873,0.424
    #Harvard_GCC_WGS-CRC-Normal.dbsnp.snp.vqsr.mss20.af1.annovar.hg19_snp137_filtered
    caddTsv = vcf.dirname + vcf.basename(".gz").sub_ext(".scores.tsv")
    lsfOut = caddTsv.sub_ext(".tsv.lsfout")
    cmd = <<-CMD
      bsub \\
        -g /gcc/germ/cadd \\
        -q park_12h \\
        -o #{lsfOut} \\
        "zcat #{vcf} | ruby -ne 'next if \\$_.start_with?(\\"#\\"); cols = \\$_.chomp.split(\\"\\t\\"); cadds = \\$1 if \\$_ =~ /;CADD=(\\S+)/; puts([\\"CADD\\", cadds.split(\\",\\")[0..1], cols[0], cols[1], cols[1], cols[3], cols[4].split(\\",\\")[0], cols[1], cols[3], cols[4]].flatten.join(\\"\\t\\"))' > #{caddTsv}"
    CMD
    puts cmd
  end
end

def import_annovar_results
  tsvs = Pathname.glob($canVcfDir + "LUAD/{snp,indel}/*.tsv").sort
  #tsvs = Pathname.glob($canVcfDir + "CRC/indel/*.tsv").sort
  tsvs.each_with_index do |tsv, vi|
    varType         = tsv.basename.to_s.match(/indel/) ? "indel" : "snp"
    mutTable        = "#{varType}s"
    cancerType      = tsv.dirname.dirname.basename.to_s
    vcf             = tsv.dirname.dirname + tsv.basename.sub_ext(".vqsr.vcf.gz")
    gzippedVcfFile  = open(vcf)
    gzippedVcfFh    = Zlib::GzipReader.new(gzippedVcfFile)
    colNames        = []
    gzippedVcfFh.each_line do |line|
      if line.start_with?("#CHROM")
        colNames = line.strip.split(/\s+/)
        break
      end
    end
    samplesNames = colNames[9..-1].map { |n| n.split("-")[0..3].join("_") }
    dbFile = (tsv.dirname.dirname + "#{cancerType}.db").to_s
    db = SQLite3::Database.new dbFile

    db.execute("DROP TABLE IF EXISTS #{mutTable};")
    if varType == "snp"
      db.execute("CREATE TABLE #{mutTable} (
                  chrom   VARCHAR,
                  pos     INT,
                  id      VARCHAR,
                  ref     VARCHAR,
                  alt     VARCHAR,
                  qual    REAL,
                  alt1    VARCHAR, alt2    VARCHAR, alt3    VARCHAR,
                  rsid1   INT, rsid2   INT, rsid3   INT,
                  tg_all_af1  REAL, tg_all_af2  REAL, tg_all_af3  REAL,
                  tg_asn_af1  REAL, tg_asn_af2  REAL, tg_asn_af3  REAL,
                  tg_amr_af1  REAL, tg_amr_af2  REAL, tg_amr_af3  REAL,
                  tg_afr_af1  REAL, tg_afr_af2  REAL, tg_afr_af3  REAL,
                  tg_eur_af1  REAL, tg_eur_af2  REAL, tg_eur_af3  REAL,
                  cadd_raw1   REAL, cadd_raw2   REAL, cadd_raw3   REAL,
                  cadd_phred1 REAL, cadd_phred2 REAL, cadd_phred3 REAL,
                  an          INT,
                  ac1 INT,  ac2 INT,  ac3 INT,
                  af1 REAL, af2 REAL, af3 REAL,
                  all_an      INT,
                  eur_an      INT,
                  all_ac1     INT,  all_ac2     INT,  all_ac3     INT,
                  all_af1     REAL, all_af2     REAL, all_af3     REAL,
                  eur_ac1     INT,  eur_ac2     INT,  eur_ac3     INT,
                  eur_af1     REAL, eur_af2     REAL, eur_af3     REAL,
                  chisq_pval_all_ac1  REAL,  chisq_pval_all_ac2 REAL,  chisq_pval_all_ac3 REAL,
                  fisher_pval_all_ac1 REAL, fisher_pval_all_ac2 REAL, fisher_pval_all_ac3 REAL,
                  chisq_pval_eur_ac1  REAL,  chisq_pval_eur_ac2 REAL,  chisq_pval_eur_ac3 REAL,
                  fisher_pval_eur_ac1 REAL, fisher_pval_eur_ac2 REAL, fisher_pval_eur_ac3 REAL,
                  #{ samplesNames.map { |x| "#{x} VARCHAR" }.join(", ") });")
    else
      db.execute("CREATE TABLE #{mutTable} (
                  chrom   VARCHAR,
                  pos     INT,
                  id      VARCHAR,
                  ref     VARCHAR,
                  alt     VARCHAR,
                  qual    REAL,
                  alt1    VARCHAR, alt2    VARCHAR, alt3    VARCHAR,
                  alt4    VARCHAR, alt5    VARCHAR, alt6    VARCHAR,
                  rsid1   INT, rsid2   INT, rsid3   INT,
                  rsid4   INT, rsid5   INT, rsid6   INT,
                  tg_all_af1  REAL, tg_all_af2  REAL, tg_all_af3  REAL,
                  tg_all_af4  REAL, tg_all_af5  REAL, tg_all_af6  REAL,
                  tg_asn_af1  REAL, tg_asn_af2  REAL, tg_asn_af3  REAL,
                  tg_asn_af4  REAL, tg_asn_af5  REAL, tg_asn_af6  REAL,
                  tg_amr_af1  REAL, tg_amr_af2  REAL, tg_amr_af3  REAL,
                  tg_amr_af4  REAL, tg_amr_af5  REAL, tg_amr_af6  REAL,
                  tg_afr_af1  REAL, tg_afr_af2  REAL, tg_afr_af3  REAL,
                  tg_afr_af4  REAL, tg_afr_af5  REAL, tg_afr_af6  REAL,
                  tg_eur_af1  REAL, tg_eur_af2  REAL, tg_eur_af3  REAL,
                  tg_eur_af4  REAL, tg_eur_af5  REAL, tg_eur_af6  REAL,
                  an          INT,
                  ac1 INT,  ac2 INT,  ac3 INT,
                  ac4 INT,  ac5 INT,  ac6 INT,
                  af1 REAL, af2 REAL, af3 REAL,
                  af4 REAL, af5 REAL, af6 REAL,
                  all_an      INT,
                  eur_an      INT,
                  all_ac1     INT,  all_ac2     INT,  all_ac3     INT,
                  all_ac4     INT,  all_ac5     INT,  all_ac6     INT,
                  all_af1     REAL, all_af2     REAL, all_af3     REAL,
                  all_af4     REAL, all_af5     REAL, all_af6     REAL,
                  eur_ac1     INT,  eur_ac2     INT,  eur_ac3     INT,
                  eur_ac4     INT,  eur_ac5     INT,  eur_ac6     INT,
                  eur_af1     REAL, eur_af2     REAL, eur_af3     REAL,
                  eur_af4     REAL, eur_af5     REAL, eur_af6     REAL,
                  chisq_pval_all_ac1  REAL,  chisq_pval_all_ac2 REAL,  chisq_pval_all_ac3 REAL,
                  chisq_pval_all_ac4  REAL,  chisq_pval_all_ac5 REAL,  chisq_pval_all_ac6 REAL,
                  fisher_pval_all_ac1 REAL, fisher_pval_all_ac2 REAL, fisher_pval_all_ac3 REAL,
                  fisher_pval_all_ac4 REAL, fisher_pval_all_ac5 REAL, fisher_pval_all_ac6 REAL,
                  chisq_pval_eur_ac1  REAL,  chisq_pval_eur_ac2 REAL,  chisq_pval_eur_ac3 REAL,
                  chisq_pval_eur_ac4  REAL,  chisq_pval_eur_ac5 REAL,  chisq_pval_eur_ac6 REAL,
                  fisher_pval_eur_ac1 REAL, fisher_pval_eur_ac2 REAL, fisher_pval_eur_ac3 REAL,
                  fisher_pval_eur_ac4 REAL, fisher_pval_eur_ac5 REAL, fisher_pval_eur_ac6 REAL,
                  #{ samplesNames.map { |x| "#{x} VARCHAR" }.join(", ") });")
    end
    sqlite3_import_tsv(dbFile, tsv, "#{mutTable}")

    #Gene annotation
    varTable  = "#{varType}_variant_functions"
    varFile   = tsv.sub_ext(".av.variant_function")
    #intergenic      NONE(dist=NONE),DDX11L1(dist=1766)      1       10108   10108   C       T       10108   C       T
    db.execute("DROP TABLE IF EXISTS #{varTable};")
    db.execute("CREATE TABLE #{varTable} (
                region  VARCHAR,
                feature VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, varFile, varTable)

    # Exonic variations
    exoTable  = "#{varType}_exonic_variant_functions"
    exoFile   = tsv.sub_ext(".av.exonic_variant_function")
    #line676 nonsynonymous SNV       SAMD11:NM_152486:exon10:c.T1027C:p.W343R,       1       877831  877831  T       C       877831  T       C
    db.execute("DROP TABLE IF EXISTS #{exoTable};")
    db.execute("CREATE TABLE #{exoTable} (
                line_num        VARCHAR,
                exonic_vartype  VARCHAR,
                exonic_change   VARCHAR,
                vchrom          VARCHAR,
                vstart          INT,
                vstop           INT,
                vref            VARCHAR,
                valt            VARCHAR,
                pos             INT,
                ref             VARCHAR,
                alt             VARCHAR);")
    sqlite3_import_tsv(dbFile, exoFile, exoTable)

    #dbSNP
    snp137File  = tsv.sub_ext(".av.hg19_snp137_dropped")
    snp137Table = "#{varType}_dbsnp137"
    #snp137  rs200279319     1       10231   10231   C       A       10231   C       A
    db.execute("DROP TABLE IF EXISTS #{snp137Table};")
    db.execute("CREATE TABLE #{snp137Table} (
                db        VARCHAR,
                rsid      VARCHAR,
                vchrom    VARCHAR,
                vstart    INT,
                vstop     INT,
                vref      VARCHAR,
                valt      VARCHAR,
                pos       INT,
                ref       VARCHAR,
                alt       VARCHAR);")
    sqlite3_import_tsv(dbFile, snp137File, snp137Table)

    snp138File  = tsv.sub_ext(".av.hg19_snp138_dropped")
    snp138Table = "#{varType}_dbsnp138"
    #snp138  rs200279319     1       10231   10231   C       A       10231   C       A
    db.execute("DROP TABLE IF EXISTS #{snp138Table};")
    db.execute("CREATE TABLE #{snp138Table} (
                db        VARCHAR,
                rsid      VARCHAR,
                vchrom    VARCHAR,
                vstart    INT,
                vstop     INT,
                vref      VARCHAR,
                valt      VARCHAR,
                pos       INT,
                ref       VARCHAR,
                alt       VARCHAR);")
    sqlite3_import_tsv(dbFile, snp138File, snp138Table)

    ## 1K Genomes
    populations = %w[ALL EUR AFR AMR ASN]
    populations.each do |population|
      tgFile  = tsv.sub_ext(".av.hg19_#{population}.sites.2012_04_dropped")
      tgTable = "#{varType}_tg_#{population.downcase}_afs"
      #1000g2012apr_all        0.74    1       1000156 1000156 C       T       1000156 C       T
      db.execute("DROP TABLE IF EXISTS #{tgTable};")
      db.execute("CREATE TABLE #{tgTable} (
                  db          VARCHAR,
                  af          REAL,
                  vchrom      VARCHAR,
                  vstart      INT,
                  vstop       INT,
                  vref        VARCHAR,
                  valt        VARCHAR,
                  pos         INT,
                  ref         VARCHAR,
                  alt         VARCHAR);")
      sqlite3_import_tsv(dbFile, tgFile, tgTable)
    end

    ## NHLBI-ESP project
    populations = %w[all ea aa]
    populations.each do |population|
      espFile  = tsv.sub_ext(".av.hg19_esp6500si_#{population}_dropped")
      espTable = "#{varType}_esp_#{population}_afs"
      #esp6500si_aa    0.389015        1       100111956       100111956       A       G       100111956       A       G
      db.execute("DROP TABLE IF EXISTS #{espTable};")
      db.execute("CREATE TABLE #{espTable} (
                  db      VARCHAR,
                  af      REAL,
                  vchrom  VARCHAR,
                  vstart  INT,
                  vstop   INT,
                  vref    VARCHAR,
                  valt    VARCHAR,
                  pos     INT,
                  ref     VARCHAR,
                  alt     VARCHAR);")
      sqlite3_import_tsv(dbFile, espFile, espTable)
    end

    ## COSMIC
    cosmicFile  = tsv.sub_ext(".av.hg19_cosmic68_dropped")
    cosmicTable = "#{varType}_cosmic68"
    #cosmic68        ID=COSM1600588,COSM1600589;OCCURENCE=1(liver)   1       115226991       115226991       A       -       115226990       CA      C
    db.execute("DROP TABLE IF EXISTS #{cosmicTable};")
    db.execute("CREATE TABLE #{cosmicTable} (
                db      VARCHAR,
                id      VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, cosmicFile, cosmicTable)

    ## SIFT
    siftFile  = tsv.sub_ext(".av.hg19_ljb23_sift_dropped")
    siftTable = "#{varType}_ljb23_sift"
    #ljb23_sift      0.71    1       10190575        10190575        G       C       10190575        G       C
    db.execute("DROP TABLE IF EXISTS #{siftTable};")
    db.execute("CREATE TABLE #{siftTable} (
                db      VARCHAR,
                score   REAL,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, siftFile, siftTable)

    ## PP2HDIV
    pp2hdivFile  = tsv.sub_ext(".av.hg19_ljb23_pp2hdiv_dropped")
    pp2hdivTable = "#{varType}_ljb23_pp2hdiv"
    db.execute("DROP TABLE IF EXISTS #{pp2hdivTable};")
    db.execute("CREATE TABLE #{pp2hdivTable} (
                db      VARCHAR,
                score   REAL,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, pp2hdivFile, pp2hdivTable)

    ## PP2HVAR
    pp2hvarFile  = tsv.sub_ext(".av.hg19_ljb23_pp2hvar_dropped")
    pp2hvarTable = "#{varType}_ljb23_pp2hvar"
    #ljb2_pp2hvar    0.001   1       100174670       100174670       C       T       100174670       C       T
    db.execute("DROP TABLE IF EXISTS #{pp2hvarTable};")
    db.execute("CREATE TABLE #{pp2hvarTable} (
                db      VARCHAR,
                score   REAL,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, pp2hvarFile, pp2hvarTable)

    ## MutationAssessor
    maFile  = tsv.sub_ext(".av.hg19_ljb23_ma_dropped")
    maTable = "#{varType}_ljb23_ma"
    #ljb2_ma 1.32    1       100195237       100195237       G       A       100195237       G       A
    db.execute("DROP TABLE IF EXISTS #{maTable};")
    db.execute("CREATE TABLE #{maTable} (
                db      VARCHAR,
                score   VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, maFile, maTable)

    ## MutationTester
    mtFile  = tsv.sub_ext(".av.hg19_ljb23_mt_dropped")
    mtTable = "#{varType}_ljb23_mt"
    #ljb2_mt 0.185158        1       100195237       100195237       G       A       100195237       G       A
    db.execute("DROP TABLE IF EXISTS #{mtTable};")
    db.execute("CREATE TABLE #{mtTable} (
                db      VARCHAR,
                score   VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, mtFile, mtTable)

    ## gwascatalog
    gwasFile  = tsv.sub_ext(".av.hg19_gwasCatalog")
    gwasTable = "#{varType}_gwas_catalogs"
    #gwascatalog     Name=Body mass index    1       1005806 1005806 C       T       1005806 C       T
    db.execute("DROP TABLE IF EXISTS #{gwasTable};")
    db.execute("CREATE TABLE #{gwasTable} (
                db        VARCHAR,
                phenotype VARCHAR,
                vchrom    VARCHAR,
                vstart    INT,
                vstop     INT,
                vref      VARCHAR,
                valt      VARCHAR,
                pos       INT,
                ref       VARCHAR,
                alt       VARCHAR);")
    sqlite3_import_tsv(dbFile, gwasFile, gwasTable)

    ## PhastConsElements
    phastCons = %w[phastConsElements46way phastConsElements46wayPlacental phastConsElements46wayPrimates]
    phastCons.each do |phastCon|
      phatConFile  = tsv.sub_ext(".av.hg19_#{phastCon}")
      phatConTable = "#{varType}_#{phastCon.underscore}"
      #phastConsElements46way  Score=369;Name=lod=42   1       63671   63671   G       A       63671   G       A
      db.execute("DROP TABLE IF EXISTS #{phatConTable};")
      db.execute("CREATE TABLE #{phatConTable} (
                  db      VARCHAR,
                  score   VARCHAR,
                  vchrom  VARCHAR,
                  vstart  INT,
                  vstop   INT,
                  vref    VARCHAR,
                  valt    VARCHAR,
                  pos     INT,
                  ref     VARCHAR,
                  alt     VARCHAR);")
      sqlite3_import_tsv(dbFile, phatConFile, phatConTable)
    end

    ## TFBS
    tfbsFile  = tsv.sub_ext(".av.hg19_tfbsConsSites")
    tfbsTable = "#{varType}_tfbs_cons_sites"
    #tfbsConsSites   Score=896;Name=V$ZID_01 1       900505  900505  G       C       900505  G       C
    db.execute("DROP TABLE IF EXISTS #{tfbsTable};")
    db.execute("CREATE TABLE #{tfbsTable} (
                db      VARCHAR,
                tfs     VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, tfbsFile, tfbsTable)

    ## wgRna
    wgRnaFile  = tsv.sub_ext(".av.hg19_wgRna")
    wgRnaTable = "#{varType}_wg_rna"
    #wgRna   Name=hsa-mir-200b       1       1102490 1102490 C       G       1102490 C       G
    db.execute("DROP TABLE IF EXISTS #{wgRnaTable};")
    db.execute("CREATE TABLE #{wgRnaTable} (
                db      VARCHAR,
                rna     VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, wgRnaFile, wgRnaTable)

    ## targetScanS
    targetScanSFile  = tsv.sub_ext(".av.hg19_targetScanS")
    targetScanSTable = "#{varType}_target_scans"
    #targetScanS     Score=50;Name=AGRN:miR-224      1       990839  990839  C       T       990839  C       T
    db.execute("DROP TABLE IF EXISTS #{targetScanSTable};")
    db.execute("CREATE TABLE #{targetScanSTable} (
                db      VARCHAR,
                target  VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, targetScanSFile, targetScanSTable)

    ## DGV
    dgvFile  = tsv.sub_ext(".av.hg19_dgvMerged")
    dgvTable = "#{varType}_dgvs"
    #dgvMerged       Name=nsv7879,dgv2n71    1       10108   10108   C       T       10108   C       T
    db.execute("DROP TABLE IF EXISTS #{dgvTable};")
    db.execute("CREATE TABLE #{dgvTable} (
                db      VARCHAR,
                var     VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, dgvFile, dgvTable)

    ## DNase cluster
    dnaseFile  = tsv.sub_ext(".av.hg19_wgEncodeRegDnaseClusteredV2")
    dnaseTable = "#{varType}_wg_encode_reg_dnase_clustered_v2"
    #wgEncodeRegDnaseClusteredV2     Name=38 1       10108   10108   C       T       10108   C       T
    db.execute("DROP TABLE IF EXISTS #{dnaseTable};")
    db.execute("CREATE TABLE #{dnaseTable} (
                db                VARCHAR,
                num_of_celllines  VARCHAR,
                vchrom            VARCHAR,
                vstart            INT,
                vstop             INT,
                vref              VARCHAR,
                valt              VARCHAR,
                pos               INT,
                ref               VARCHAR,
                alt               VARCHAR);")
    sqlite3_import_tsv(dbFile, dnaseFile, dnaseTable)

    ## TFBS cluster
    tfbsClstFile  = tsv.sub_ext(".av.hg19_wgEncodeRegTfbsClusteredV3")
    tfbsClstTable = "#{varType}_wg_encode_reg_tfbs_clustered_v3"
    #wgEncodeRegTfbsClusteredV3      Name=JUND,USF1,BHLHE40,FOSL2    1       534316  534316  G       C       534316  G       C
    db.execute("DROP TABLE IF EXISTS #{tfbsClstTable};")
    db.execute("CREATE TABLE #{tfbsClstTable} (
                db      VARCHAR,
                tfs     VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, tfbsClstFile, tfbsClstTable)

    ## BROAD HMM
    cellLines = %w[Gm12878 H1hesc Hepg2 Hmec Hsmm Huvec K562 Nhek Nhlf]
    cellLines.each do |cellLine|
      hmmFile  = tsv.sub_ext(".av.hg19_wgEncodeBroadHmm#{cellLine}HMM")
      hmmTable = "#{varType}_wg_encode_broad_hmm_#{cellLine.downcase}"
      #wgEncodeBroadHmmGm12878HMM      Name=15_Repetitive/CNV  1       10108   10108   C       T       10108   C       T
      db.execute("DROP TABLE IF EXISTS #{hmmTable};")
      db.execute("CREATE TABLE #{hmmTable} (
                  db      VARCHAR,
                  state   VARCHAR,
                  vchrom  VARCHAR,
                  vstart  INT,
                  vstop   INT,
                  vref    VARCHAR,
                  valt    VARCHAR,
                  pos     INT,
                  ref     VARCHAR,
                  alt     VARCHAR);")
      sqlite3_import_tsv(dbFile, hmmFile, hmmTable)
    end

    ## RepeatMasker
    rmskFile  = tsv.sub_ext(".av.hg19_rmsk_mod")
    rmskTable = "#{varType}_rmsk"
    #rmsk_mod        Name=Simple_repeat      1       10108   10108   C       T       10108   C       T
    db.execute("DROP TABLE IF EXISTS #{rmskTable};")
    db.execute("CREATE TABLE #{rmskTable} (
                db      VARCHAR,
                type    VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, rmskFile, rmskTable)

    ## Simple Repeat
    simpleRepeatFile  = tsv.sub_ext(".av.hg19_simpleRepeat")
    simpleRepeatTable = "#{varType}_simple_repeats"
    #simpleRepeat    Name=trf        1       10108   10108   C       T       10108   C       T
    db.execute("DROP TABLE IF EXISTS #{simpleRepeatTable};")
    db.execute("CREATE TABLE #{simpleRepeatTable} (
                db      VARCHAR,
                type    VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, simpleRepeatFile, simpleRepeatTable)

    ## Nested Repeat
    nestedRepeatFile  = tsv.sub_ext(".av.hg19_nestedRepeats")
    nestedRepeatTable = "#{varType}_nested_repeats"
    #nestedRepeats   Name=L2 1       54490   54490   G       A       54490   G       A
    db.execute("DROP TABLE IF EXISTS #{nestedRepeatTable};")
    db.execute("CREATE TABLE #{nestedRepeatTable} (
                db      VARCHAR,
                type    VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, nestedRepeatFile, nestedRepeatTable)

    ## Microsatelite
    microsatFile  = tsv.sub_ext(".av.hg19_microsat")
    microsatTable = "#{varType}_microsatelites"
    #microsat        Name=27xTG      1       1705164 1705164 T       A       1705164 T       A
    db.execute("DROP TABLE IF EXISTS #{microsatTable};")
    db.execute("CREATE TABLE #{microsatTable} (
                db      VARCHAR,
                type    VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, microsatFile, microsatTable)

    ## Combined DHS Peaks
    dhsFile  = tsv.sub_ext(".av.hg19_combinedDhsPeaks")
    dhsTable = "#{varType}_combined_dhs_peaks"
    #combinedDhsPeaks        Name=MCV-37     1       10231   10231   C       A       10231   C       A
    db.execute("DROP TABLE IF EXISTS #{dhsTable};")
    db.execute("CREATE TABLE #{dhsTable} (
                db                VARCHAR,
                num_of_celllines  VARCHAR,
                vchrom            VARCHAR,
                vstart            INT,
                vstop             INT,
                vref              VARCHAR,
                valt              VARCHAR,
                pos               INT,
                ref               VARCHAR,
                alt               VARCHAR);")
    sqlite3_import_tsv(dbFile, dhsFile, dhsTable)


    ## DHS - DHS
    dhsDhsFile  = tsv.sub_ext(".annovar.hg19_distalDhsToPromoterDhs")
    dhsDhsTable = "#{varType}_distal_dhs_to_promoter_dhs"
    #distalDhsToPromoterDhs  Name=RP5-857K21.1,RP5-857K21.3  1       747575  747575  A       G       747575  A       G
    db.execute("DROP TABLE IF EXISTS #{dhsDhsTable};")
    db.execute("CREATE TABLE #{dhsDhsTable} (
                db      VARCHAR,
                genes   VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, dhsDhsFile, dhsDhsTable)

    ## DHS - RNAseq
    dhsRnaFile  = tsv.sub_ext(".av.hg19_dhsToGeneExpression")
    dhsRnaTable = "#{varType}_dhs_to_gene_expression"
    #dhsToGeneExpression     Name=RP11-206L10.2      1       710706  710706  G       A       710706  G       A
    db.execute("DROP TABLE IF EXISTS #{dhsRnaTable};")
    db.execute("CREATE TABLE #{dhsRnaTable} (
                db      VARCHAR,
                genes   VARCHAR,
                vchrom  VARCHAR,
                vstart  INT,
                vstop   INT,
                vref    VARCHAR,
                valt    VARCHAR,
                pos     INT,
                ref     VARCHAR,
                alt     VARCHAR);")
    sqlite3_import_tsv(dbFile, dhsRnaFile, dhsRnaTable)

    ## DHS Promoter Known
    promoKnownFile  = tsv.sub_ext(".av.hg19_promoterDhsMasterKnown")
    promoKnownTable = "#{varType}_promoter_dhs_master_known"
    #promoterDhsMasterKnown  Name=MCV-106    1       569933  569933  G       A       569933  G       A
    db.execute("DROP TABLE IF EXISTS #{promoKnownTable};")
    db.execute("CREATE TABLE #{promoKnownTable} (
                db                VARCHAR,
                num_of_celllines  VARCHAR,
                vchrom            VARCHAR,
                vstart            INT,
                vstop             INT,
                vref              VARCHAR,
                valt              VARCHAR,
                pos               INT,
                ref               VARCHAR,
                alt               VARCHAR);")
    sqlite3_import_tsv(dbFile, promoKnownFile, promoKnownTable)

    ## DHS Promoter Novel
    promoNovelFile  = tsv.sub_ext(".av.hg19_promoterDhsMasterNovel")
    promoNovelTable = "#{varType}_promoter_dhs_master_novel"
    #promoterDhsMasterNovel  Name=MCV-106    1       569933  569933  G       A       569933  G       A
    db.execute("DROP TABLE IF EXISTS #{promoNovelTable};")
    db.execute("CREATE TABLE #{promoNovelTable} (
                db                VARCHAR,
                num_of_celllines  VARCHAR,
                vchrom            VARCHAR,
                vstart            INT,
                vstop             INT,
                vref              VARCHAR,
                valt              VARCHAR,
                pos               INT,
                ref               VARCHAR,
                alt               VARCHAR);")
    sqlite3_import_tsv(dbFile, promoNovelFile, promoNovelTable)
  end
end

def summarize_annovar_results
  dbFiles = Pathname.glob($canVcfDir + "LUAD/*.db").sort
  dbFiles.each do |dbFile|
    db = SQLite3::Database.new(dbFile.to_s)
    %w[snp indel].each do |varType|
      sampleNames = []
      tableInfo = `sqlite3 #{dbFile} 'PRAGMA table_info(#{varType}s);'`
      tableInfo.split("\n").each do |line|
        next if line.empty?
        vals = line.split("\t")
        sampleNames << vals[1] if vals[1].start_with?("TCGA")
      end
      db.execute("DROP TABLE IF EXISTS #{varType}_annotations;")
      sql =<<-SQL
CREATE TABLE #{varType}_annotations AS 
  SELECT
    m.chrom, m.pos, m.id, m.qual, m.ref, m.alt,
      #{varType == "snp" ? "m.alt1, m.alt2, m.alt3, 
                            m.rsid1, m.rsid2, m.rsid3,
                            m.all_an, m.eur_an,
                            m.all_ac1, m.all_ac2, m.all_ac3,
                            m.all_af1, m.all_af2, m.all_af3,
                            m.eur_ac1, m.eur_ac2, m.eur_ac3,
                            m.eur_af1, m.eur_af2, m.eur_af3,
                            m.tg_all_af1, m.tg_all_af2, m.tg_all_af3,
                            m.tg_asn_af1, m.tg_asn_af2, m.tg_asn_af3,
                            m.tg_amr_af1, m.tg_amr_af2, m.tg_amr_af3,
                            m.tg_afr_af1, m.tg_afr_af2, m.tg_afr_af3,
                            m.tg_eur_af1, m.tg_eur_af2, m.tg_eur_af3,
                            m.chisq_pval_all_ac1,  m.chisq_pval_all_ac2,  m.chisq_pval_all_ac3,
                            m.fisher_pval_all_ac1, m.fisher_pval_all_ac2, m.fisher_pval_all_ac3,
                            m.chisq_pval_eur_ac1,  m.chisq_pval_eur_ac2,  m.chisq_pval_eur_ac3,
                            m.fisher_pval_eur_ac1, m.fisher_pval_eur_ac2, m.fisher_pval_eur_ac3,
                            m.cadd_raw1, m.cadd_raw2, m.cadd_raw3,
                            m.cadd_phred1, m.cadd_phred2, m.cadd_phred3," :
                           "m.alt1, m.alt2, m.alt3, m.alt4, m.alt5, m.alt6, 
                            m.rsid1, m.rsid2, m.rsid3, m.rsid4, m.rsid5, m.rsid6, 
                            m.all_an, m.eur_an,
                            m.all_ac1, m.all_ac2, m.all_ac3, m.all_ac4, m.all_ac5, m.all_ac6,
                            m.all_af1, m.all_af2, m.all_af3, m.all_af4, m.all_af5, m.all_af6,
                            m.eur_ac1, m.eur_ac2, m.eur_ac3, m.eur_ac4, m.eur_ac5, m.eur_ac6,
                            m.eur_af1, m.eur_af2, m.eur_af3, m.eur_af4, m.eur_af5, m.eur_af6,
                            m.tg_all_af1, m.tg_all_af2, m.tg_all_af3, m.tg_all_af4, m.tg_all_af5, m.tg_all_af6,
                            m.tg_asn_af1, m.tg_asn_af2, m.tg_asn_af3, m.tg_asn_af4, m.tg_asn_af5, m.tg_asn_af6,
                            m.tg_amr_af1, m.tg_amr_af2, m.tg_amr_af3, m.tg_amr_af4, m.tg_amr_af5, m.tg_amr_af6,
                            m.tg_afr_af1, m.tg_afr_af2, m.tg_afr_af3, m.tg_afr_af4, m.tg_afr_af5, m.tg_afr_af6,
                            m.tg_eur_af1, m.tg_eur_af2, m.tg_eur_af3, m.tg_eur_af4, m.tg_eur_af5, m.tg_eur_af6,
                            m.chisq_pval_all_ac1,  m.chisq_pval_all_ac2,  m.chisq_pval_all_ac3,
                            m.chisq_pval_all_ac4,  m.chisq_pval_all_ac5,  m.chisq_pval_all_ac6,
                            m.fisher_pval_all_ac1, m.fisher_pval_all_ac2, m.fisher_pval_all_ac3,
                            m.fisher_pval_all_ac4, m.fisher_pval_all_ac5, m.fisher_pval_all_ac6,
                            m.chisq_pval_eur_ac1,  m.chisq_pval_eur_ac2,  m.chisq_pval_eur_ac3,
                            m.chisq_pval_eur_ac4,  m.chisq_pval_eur_ac5,  m.chisq_pval_eur_ac6,
                            m.fisher_pval_eur_ac1, m.fisher_pval_eur_ac2, m.fisher_pval_eur_ac3,
                            m.fisher_pval_eur_ac4, m.fisher_pval_eur_ac5, m.fisher_pval_eur_ac6,"}
    v.region AS region, v.feature AS feature,
    e.exonic_vartype AS exonic_vartype, e.exonic_change AS exonic_change,
    s1.rsid AS dbsnp137, s2.rsid AS dbsnp138,
    t1.af AS av_tg_all_af, t2.af AS av_tg_eur_af, t3.af AS av_tg_afr_af, t4.af AS av_tg_amr_af, t5.af AS av_tg_asn_af,
    e1.af AS esp6500_all_af, e2.af AS esp6500_ea_af, e3.af AS esp6500_aa_af,
    cs.id AS cosmic,
    sf.score AS sift,
    pd.score AS pp2hdiv,
    pv.score AS pp2hvar,
    ma.score AS ma,
    mt.score AS mt,
    g.phenotype AS gwas_catalog,
    c1.score AS phast_cons_elements46way,
    c2.score AS phast_cons_elements46way_placental,
    c3.score AS phast_cons_elements46way_primates,
    dc.num_of_celllines AS wg_encode_reg_dnase_clustered_v2,
    t.tfs AS tfbs_cons_sites,
    tc.tfs AS wg_encode_reg_tfbs_clustered_v3,
    h1.state AS wg_encode_broad_hmm_gm12878,
    h2.state AS wg_encode_broad_hmm_h1hesc,
    h3.state AS wg_encode_broad_hmm_hepg2,
    h4.state AS wg_encode_broad_hmm_hmec,
    h5.state AS wg_encode_broad_hmm_hsmm,
    h6.state AS wg_encode_broad_hmm_huvec,
    h7.state AS wg_encode_broad_hmm_k562,
    h8.state AS wg_encode_broad_hmm_nhek,
    h9.state AS wg_encode_broad_hmm_nhlf,
    r.type AS repeat_masker,
    s.type AS simple_repeat,
    n.type AS nested_repeat,
    wr.rna AS wg_rna,
    ts.target AS mirna_target,
    d.var AS dgv,
    ms.type AS microsatelite,
    d1.num_of_celllines AS combined_dhs_peak,
    d2.genes AS distal_dhs_to_promoter_dhs,
    d3.genes AS dhs_to_gene_expression,
    p1.num_of_celllines AS promoter_dhs_master_known,
    p2.num_of_celllines AS promoter_dhs_master_novel,
    #{sampleNames.map { |s| "quote(m.#{s}) AS #{s}" }.join(", ")}
    FROM #{varType}s AS m
    LEFT JOIN #{varType}_variant_functions AS v on m.chrom = v.vchrom and m.pos = v.pos and m.ref = v.ref and m.alt1 = v.alt
    LEFT JOIN #{varType}_exonic_variant_functions AS e on m.chrom = e.vchrom and m.pos = e.pos and m.ref = e.ref and m.alt1 = e.alt
    LEFT JOIN #{varType}_dbsnp137 AS s1 on m.chrom = s1.vchrom and m.pos = s1.pos and m.ref = s1.ref and m.alt1 = s1.alt
    LEFT JOIN #{varType}_dbsnp138 AS s2 on m.chrom = s2.vchrom and m.pos = s2.pos and m.ref = s2.ref and m.alt1 = s2.alt
    LEFT JOIN #{varType}_tg_all_afs AS t1 on m.chrom = t1.vchrom and m.pos = t1.pos and m.ref = t1.ref and m.alt1 = t1.alt
    LEFT JOIN #{varType}_tg_eur_afs AS t2 on m.chrom = t2.vchrom and m.pos = t2.pos and m.ref = t2.ref and m.alt1 = t2.alt
    LEFT JOIN #{varType}_tg_afr_afs AS t3 on m.chrom = t3.vchrom and m.pos = t3.pos and m.ref = t3.ref and m.alt1 = t3.alt
    LEFT JOIN #{varType}_tg_amr_afs AS t4 on m.chrom = t4.vchrom and m.pos = t4.pos and m.ref = t4.ref and m.alt1 = t4.alt
    LEFT JOIN #{varType}_tg_asn_afs AS t5 on m.chrom = t5.vchrom and m.pos = t5.pos and m.ref = t5.ref and m.alt1 = t5.alt
    LEFT JOIN #{varType}_esp_all_afs AS e1 on m.chrom = e1.vchrom and m.pos = e1.pos and m.ref = e1.ref and m.alt1 = e1.alt
    LEFT JOIN #{varType}_esp_ea_afs AS e2 on m.chrom = e2.vchrom and m.pos = e2.pos and m.ref = e2.ref and m.alt1 = e2.alt
    LEFT JOIN #{varType}_esp_aa_afs AS e3 on m.chrom = e3.vchrom and m.pos = e3.pos and m.ref = e3.ref and m.alt1 = e3.alt
    LEFT JOIN #{varType}_cosmic68 AS cs on m.chrom = cs.vchrom and m.pos = cs.pos and m.ref = cs.ref and m.alt1 = cs.alt
    LEFT JOIN #{varType}_ljb23_sift AS sf on m.chrom = sf.vchrom and m.pos = sf.pos and m.ref = sf.ref and m.alt1 = sf.alt
    LEFT JOIN #{varType}_ljb23_pp2hdiv AS pd on m.chrom = pd.vchrom and m.pos = pd.pos and m.ref = pd.ref and m.alt1 = pd.alt
    LEFT JOIN #{varType}_ljb23_pp2hvar AS pv on m.chrom = pv.vchrom and m.pos = pv.pos and m.ref = pv.ref and m.alt1 = pv.alt
    LEFT JOIN #{varType}_ljb23_ma AS ma on m.chrom = ma.vchrom and m.pos = ma.pos and m.ref = ma.ref and m.alt1 = ma.alt
    LEFT JOIN #{varType}_ljb23_mt AS mt on m.chrom = mt.vchrom and m.pos = mt.pos and m.ref = mt.ref and m.alt1 = mt.alt
    LEFT JOIN #{varType}_gwas_catalogs AS g on m.chrom = g.vchrom and m.pos = g.pos and m.ref = g.ref and m.alt1 = g.alt
    LEFT JOIN #{varType}_phast_cons_elements46way AS c1 on m.chrom = c1.vchrom and m.pos = c1.pos and m.ref = c1.ref and m.alt1 = c1.alt
    LEFT JOIN #{varType}_phast_cons_elements46way_placental AS c2 on m.chrom = c2.vchrom and m.pos = c2.pos and m.ref = c2.ref and m.alt1 = c2.alt
    LEFT JOIN #{varType}_phast_cons_elements46way_primates AS c3 on m.chrom = c3.vchrom and m.pos = c3.pos and m.ref = c3.ref and m.alt1 = c3.alt
    LEFT JOIN #{varType}_tfbs_cons_sites AS t on m.chrom = t.vchrom and m.pos = t.pos and m.ref = t.ref and m.alt1 = t.alt
    LEFT JOIN #{varType}_wg_rna AS wr on m.chrom = wr.vchrom and m.pos = wr.pos and m.ref = wr.ref and m.alt1 = wr.alt
    LEFT JOIN #{varType}_target_scans AS ts on m.chrom = ts.vchrom and m.pos = ts.pos and m.ref = ts.ref and m.alt1 = ts.alt
    LEFT JOIN #{varType}_dgvs AS d on m.chrom = d.vchrom and m.pos = d.pos and m.ref = d.ref and m.alt1 = d.alt
    LEFT JOIN #{varType}_wg_encode_reg_dnase_clustered_v2 AS dc on m.chrom = dc.vchrom and m.pos = dc.pos and m.ref = dc.ref and m.alt1 = dc.alt
    LEFT JOIN #{varType}_wg_encode_reg_tfbs_clustered_v3 AS tc on m.chrom = tc.vchrom and m.pos = tc.pos and m.ref = tc.ref and m.alt1 = tc.alt
    LEFT JOIN #{varType}_wg_encode_broad_hmm_gm12878 AS h1 on m.chrom = h1.vchrom and m.pos = h1.pos and m.ref = h1.ref and m.alt1 = h1.alt
    LEFT JOIN #{varType}_wg_encode_broad_hmm_h1hesc AS h2 on m.chrom = h2.vchrom and m.pos = h2.pos and m.ref = h2.ref and m.alt1 = h2.alt
    LEFT JOIN #{varType}_wg_encode_broad_hmm_hepg2 AS h3 on m.chrom = h3.vchrom and m.pos = h3.pos and m.ref = h3.ref and m.alt1 = h3.alt
    LEFT JOIN #{varType}_wg_encode_broad_hmm_hmec AS h4 on m.chrom = h4.vchrom and m.pos = h4.pos and m.ref = h4.ref and m.alt1 = h4.alt
    LEFT JOIN #{varType}_wg_encode_broad_hmm_hsmm AS h5 on m.chrom = h5.vchrom and m.pos = h5.pos and m.ref = h5.ref and m.alt1 = h5.alt
    LEFT JOIN #{varType}_wg_encode_broad_hmm_huvec AS h6 on m.chrom = h6.vchrom and m.pos = h6.pos and m.ref = h6.ref and m.alt1 = h6.alt
    LEFT JOIN #{varType}_wg_encode_broad_hmm_k562 AS h7 on m.chrom = h7.vchrom and m.pos = h7.pos and m.ref = h7.ref and m.alt1 = h7.alt
    LEFT JOIN #{varType}_wg_encode_broad_hmm_nhek AS h8 on m.chrom = h8.vchrom and m.pos = h8.pos and m.ref = h8.ref and m.alt1 = h8.alt
    LEFT JOIN #{varType}_wg_encode_broad_hmm_nhlf AS h9 on m.chrom = h9.vchrom and m.pos = h9.pos and m.ref = h9.ref and m.alt1 = h9.alt
    LEFT JOIN #{varType}_rmsk AS r on m.chrom = r.vchrom and m.pos = r.pos and m.ref = r.ref and m.alt1 = r.alt
    LEFT JOIN #{varType}_simple_repeats AS s on m.chrom = s.vchrom and m.pos = s.pos and m.ref = s.ref and m.alt1 = s.alt
    LEFT JOIN #{varType}_nested_repeats AS n on m.chrom = n.vchrom and m.pos = n.pos and m.ref = n.ref and m.alt1 = n.alt
    LEFT JOIN #{varType}_microsatelites AS ms on m.chrom = ms.vchrom and m.pos = ms.pos and m.ref = ms.ref and m.alt1 = ms.alt
    LEFT JOIN #{varType}_combined_dhs_peaks AS d1 on m.chrom = d1.vchrom and m.pos = d1.pos and m.ref = d1.ref and m.alt1 = d1.alt
    LEFT JOIN #{varType}_distal_dhs_to_promoter_dhs AS d2 on m.chrom = d2.vchrom and m.pos = d2.pos and m.ref = d2.ref and m.alt1 = d2.alt
    LEFT JOIN #{varType}_dhs_to_gene_expression AS d3 on m.chrom = d3.vchrom and m.pos = d3.pos and m.ref = d3.ref and m.alt1 = d3.alt
    LEFT JOIN #{varType}_promoter_dhs_master_known AS p1 on m.chrom = p1.vchrom and m.pos = p1.pos and m.ref = p1.ref and m.alt1 = p1.alt
    LEFT JOIN #{varType}_promoter_dhs_master_novel AS p2 on m.chrom = p2.vchrom and m.pos = p2.pos and m.ref = p2.ref and m.alt1 = p2.alt;
      SQL
      puts sql
      db.execute(sql)
    end
    db.close
  end
end


##split_vqsr_vcfs_into_chromosomes
##merge_chr_vqsr_vcfs

##split_vqsr_vcfs_into_pieces
#split_vqsr_chr_vcfs_into_even_sized_chunks
#filter_vqsr_vcfs
#annotate_vcfs_with_dbsnp_rsids
#annotate_vcfs_with_1000genomes_afs
#annotate_snp_vcfs_with_cadd_scores
#convert_snp_vcf_to_tsv
#convert_indel_vcf_to_tsv
#create_annovar_input
#run_annovar_gene
#run_annovar_filter
#run_annovar_region
#run_annovar_dhs_gene
##extract_cadd_scores
#import_annovar_results
#summarize_annovar_results
