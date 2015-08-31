#!/usr/bin/env ruby

require 'logger'
require 'pathname'
require 'parallel'
#require "google_drive_v0"
#require "google/api_client"

$logger = Logger.new(STDOUT)
$logger.level = Logger::DEBUG

$java7_bin        = Pathname.new("/opt/java/jdk1.7.0_71/bin/java")
$gcc_bam_dir1     = Pathname.new("/groups/kucherlapati/GCC/LevelII")
$home_dir         = Pathname.new("/home/sl279")
$base_dir         = Pathname.new("/groups/kucherlapati/GCC/Germline")
$bam_dir          = $base_dir + "BAM"
$script_dir       = $base_dir + "Scripts"
$vcf_dir          = $base_dir + "VCF"
$gvcf_dir          = $base_dir + "GVCF"
$install_dir      = $home_dir + "BiO/Install"
$bwa_bin          = $install_dir + "bwa-0.7.5a/bwa"
$samtools_bin     = $install_dir + "samtools-0.1.19/samtools"
$gatk2_bin        = $install_dir + "GenomeAnalysisTK-2.8-1-g932cd3a/GenomeAnalysisTK.jar"
$gatk3_bin        = $install_dir + "GATK3/GenomeAnalysisTK.jar"
$queue3_bin       = $install_dir + "Queue-3.3-0/Queue.jar"
$gatk_bundle_dir  = $install_dir + "GATK-bundle/2.8/b37"
#$refseq_broad     = Pathname.new "/groups/kucherlapati/GCC/LevelII/hg19Ref/Homo_sapiens_assembly19.fasta"
$ref_dict          = Pathname.new "/groups/kucherlapati/GCC/LevelII/hg19Ref/Homo_sapiens_assembly19.dict"
$dbsnp_dir         = $base_dir + "DBSNP"
$refseq           = $gatk_bundle_dir + "human_g1k_v37_decoy.fasta"
$dbsnp            = $gatk_bundle_dir + "dbsnp_138.b37.vcf"
$hapmap           = $gatk_bundle_dir + "hapmap_3.3.b37.vcf"
$g1kIndel         = $gatk_bundle_dir + "1000G_phase1.indels.b37.vcf"
$g1kSnp           = $gatk_bundle_dir + "1000G_phase1.snps.high_confidence.b37.vcf"
$millsIndel       = $gatk_bundle_dir + "Mills_and_1000G_gold_standard.indels.b37.vcf"
$g1kOmni          = $gatk_bundle_dir + "1000G_omni2.5.b37.vcf"
$chromLenb37      = $gatk_bundle_dir + "hg19.genome"
$chrs             = (1..22).to_a << "X" << "Y"
$chrchrs          = $chrs.map { |c| "chr#{c}" }
$tmp_dir           = $home_dir + "BiO/Temp/GATK"; $tmp_dir.mkpath
$cancer_types     = %w[BLCA BRCA CRC ESCA HNSC LGG LUAD PRAD SKCM STAD THCA UCEC UVM]
#$normalWskey   = "0AqyF8I3ZSraYdFl3dmZReEhEOU1iOXdFVWs2TFFmcGc"
#$normalWskey   = "0AqyF8I3ZSraYdEN4LTcxdVZtS21PQUtsOUFaeFJDR2c"
$normalWskey   = "1ypZKmoAeE7OVbrp_Wfxasn1yEVpDPQBQcZD5z_zfOJs"

def submit(cmd)
  out = `#{cmd}`.chomp
  if out =~ /Job\s{1}<(\d+)>/
    puts "Job #{$1} submitted."
    return Integer($1)
  else
    abort out
  end
end

def call_snps_using_queue_using_old_bams_normal
  #client = Google::APIClient.new
  #auth = client.authorization
  #auth.client_id = "seminlee@gmail.com"
  #auth.client_secret = "xgmuriuxalzmrgqe"
  ##auth.scope = "https://www.googleapis.com/auth/drive " + "https://spreadsheets.google.com/feeds/"
  #auth.scope = "https://docs.google.com/spreadsheet"
  #auth.redirect_uri = "urn:ietf:wg:oauth:2.0:oob"
  #print("1. Open this page:\n%s\n\n" % auth.authorization_uri)
  #print("2. Enter the authorization code shown in the page: ")
  #auth.code = $stdin.gets.chomp
  #auth.fetch_access_token!
  #access_token = auth.access_token
  #gdSession     = GoogleDrive.login_with_oauth(access_token)
  gdSession     = GoogleDriveV0.login("seminlee@gmail.com", "xgmuriuxalzmrgqe")
  foundPatient  = {}
  readySamples  = {}
  totalCount    = 0
  $cancer_types.each do |cancerType|
    gdSsheet = gdSession.spreadsheet_by_key($normalWskey)
    gdWsheet = gdSsheet.worksheet_by_title(cancerType)
    readySamples[cancerType] = []
    if gdWsheet.respond_to?("rows")
      gdWsheet.rows.each do |row|
        next unless row[0].start_with?("TCGA")
        nsid, nrid, segSize, cnv, snp = row
        if snp.to_i == 1
          pid = nsid.split('-')[0..2].join('-')
          if foundPatient.has_key?(pid)
            $logger.error "#{pid} has beend already found in #{cancerType}!"
            exit
          else
            foundPatient[pid] = true
          end
          shortNsid = nsid.split('-')[0..3].join('-')
          key = [shortNsid, nrid].join('___')
          readySamples[cancerType] << key
        end
      end
    end
    totalCount += readySamples[cancerType].size
    $logger.info "[#{cancerType}] #{readySamples[cancerType].size} analysis-ready bams found"
  end
  $logger.info "[PanCan#{$cancer_types.size}] #{totalCount} analysis-ready bams found in total"

  # Create a list of bam files to be fed into
  oriQueueScript = $script_dir + "GCC-Germline-Call_SNP.scala"
  panQueueScript = $panVcfDir + oriQueueScript.basename.sub_ext("-PANCAN-Normal.scala")
  FileUtils.cp(oriQueueScript, panQueueScript)
  panBamListFile = $panVcfDir + "BamFileListForGATK-PanCan-Normal.list"
  panBamListFile.open('w') do |pfile|
    readySamples.keys.each do |cancerType|
      outDir = $vcf_dir + cancerType
      canQueueScript = outDir + oriQueueScript.basename.sub_ext("-#{cancerType}-Normal.scala")
      canQueueScript.dirname.mkpath
      FileUtils.cp(oriQueueScript, canQueueScript)
      canBamListFile = outDir + "BamFileListForGATK-#{cancerType}-Normal.list"
      canBamListFile.open('w') do |cfile|
        canNomalSamples = readySamples[cancerType]
        canNomalSamples.each do |canNormalSample|
          sid, rid    = canNormalSample.split("___")
          nrid        = rid.gsub(/\_L\d+$/, "")
          globPattern = nil
          bams        = if cancerType == "CRC"
                          globPattern = "/groups/kucherlapati/GCC/LevelII/#{cancerType}/#{sid}*.bam"
                          Pathname.glob(globPattern)
                        elsif cancerType == "CESC"
                          globPattern = "/home/sl279/BiO/Research/GCC/LevelII/#{cancerType}/#{nrid}/#{sid}*.bam"
                          Pathname.glob(globPattern)
                        else
                          globPattern = "/groups/kucherlapati/GCC/LevelII/#{cancerType}/#{nrid}/#{sid}*.bam"
                          Pathname.glob(globPattern)
                        end
          bam         = if bams.size == 0
                          $logger.warn "[#{cancerType}] Cannot find #{canNormalSample}: #{globPattern} !!!"
                          exit
                        elsif bams.size > 1
                          if (sid == "TCGA-CG-4466-10A" || sid == "TCGA-CG-4477-10A" || sid == "TCGA-A5-A0GE-10A" || sid == "TCGA-A5-A0GV-10A")
                            bams[0]
                          else
                            $logger.warn "[#{cancerType}] Found #{bams.size} bams: #{globPattern} !!!\n#{bams.map { |b| b.to_s }.join("\n")}"
                            exit
                          end
                        else
                          bams[0]
                        end
          pfile.puts bam
          cfile.puts bam
        end
      end

      lsfout = canQueueScript.sub_ext(".scala.lsfout")
      logDir = outDir + ".qlog-normal"
          #-q long -W 200:0 \\
      cmd = <<-CMD
        bsub \\
          -q i2b2_7d \\
          -g /gcc/ug/#{cancerType} \\
          -R "rusage[mem=25000]" -M 25000000 \\
          -o #{lsfout} \\
          #{$java7_bin} -Xmx23G -jar #{$queue3_bin} \\
            -S #{canQueueScript} \\
            -i #{canBamListFile} \\
            -O #{$vcf_dir + cancerType} \\
            -st Harvard_GCC_WGS-#{cancerType}-Normal \\
            -logDir #{logDir} \\
            -R #{$refseq_broad} \\
            -D #{$dbsnp} \\
            -tempDir #{$tmp_dir} \\
            -bsub \\
            -jobQueue short \\
            -sg 2000 \\
            -run \\
            -l INFO \\
            -startFromScratch 
      CMD
      puts cmd
    end
  end
end

def call_snps_using_queue_using_old_bams_normal
end


def recalibrate_vcf_using_queue
  canType = "CESC"
  oriRecalScala = $script_dir + "GCC-Germline-Recalibrate_VCF.scala"
  raw_vcfs = Pathname.glob($vcf_dir + "#{canType}/Harvard_GCC_WGS-*-Normal.dbsnp.vcf.gz")
  raw_vcfs.each do |raw_vcf|
    cancerDir = raw_vcf.dirname
    cancerType = raw_vcf.dirname.basename.to_s
    cancerRecalScala = cancerDir + oriRecalScala.basename.sub_ext("-#{cancerType}.scala")
    lsfout = cancerRecalScala.sub_ext(".scala.lsfout")
    FileUtils.cp(oriRecalScala, cancerRecalScala)
        #-q long -W 200:0 \\
    cmd = <<-CMD
      bsub \\
        -g /gcc/germ/#{cancerType.downcase} \\
        -q i2b2_7d \\
        -R "rusage[mem=25000]" -M 25000000 \\
        -o #{lsfout} \\
        #{$java7_bin} -jar -Xmx23G #{$queue3_bin} \\
          -S #{cancerRecalScala} \\
          -i #{raw_vcf} \\
          -O #{cancerDir} \\
          -logDir #{cancerDir} \\
          -R #{$refseq_broad} \\
          -db #{$dbsnp} \\
          -hm #{$hapmap} \\
          -om #{$g1kOmni} \\
          -kg #{$g1kSnp} \\
          -ml #{$millsIndel} \\
          -tempDir #{$tmp_dir} \\
          -bsub \\
          -jobQueue i2b2_1d \\
          -sg 2000 \\
          -run \\
          -l INFO \\
          -startFromScratch 
    CMD
    puts cmd
  end
end

def split_cancer_vcf_by_chr
  cancer_type = "BLCA"
  %w[indel snp].each do |var_type|
    vqsr_vcf_dir = $vcf_dir + "CANCER" + cancer_type
    vqsr_vcf = vqsr_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.#{var_type}.vqsr.vcf.gz"
    out_vcf_dir = vqsr_vcf_dir + var_type; out_vcf_dir.mkpath
    $chrs.each do |chr|
      out_vcf_chr_dir = out_vcf_dir + "chrs" + chr.to_s; out_vcf_chr_dir.mkpath
      vqsr_chr_vcf = out_vcf_chr_dir + vqsr_vcf.basename(".gz").sub_ext(".pass.an100.chr#{chr}.vcf.gz")
      lsfout = vqsr_chr_vcf.sub_ext(".gz.lsfout")
          #-q short -W 12:0 \\
      cmd = <<-CMD
        bsub \\
          -g /germ/sv \\
          -q i2b2_12h \\
          -R "rusage[mem=5500]" -M 5500000 \\
          -o #{lsfout} \\
          #{$java7_bin} -XX:+UseSerialGC -Xmx5G -jar #{$gatk3_bin} \\
            -T SelectVariants \
            -R #{$refseq} \\
            --variant #{vqsr_vcf} \\
            -o #{vqsr_chr_vcf} \\
            -L #{chr} \\
            --excludeFiltered \\
            -select "AN >= 100"
      CMD
      submit cmd
    end
  end
end

def split_pancan_vcf_by_chr
  %w[indel snp].each do |var_type|
    vqsr_vcf_dir = $vcf_dir + "PANCAN"
    vqsr_vcf = vqsr_vcf_dir + "Harvard_GCC-WGS-PANCAN.#{var_type}.vqsr.vcf.gz"
    out_vcf_dir = vqsr_vcf_dir + var_type; out_vcf_dir.mkpath
    $chrs.each do |chr|
      out_vcf_chr_dir = out_vcf_dir + "chrs" + chr.to_s; out_vcf_chr_dir.mkpath
      vqsr_chr_vcf = out_vcf_chr_dir + vqsr_vcf.basename(".gz").sub_ext(".pass.an100.chr#{chr}.vcf.gz")
      lsfout = vqsr_chr_vcf.sub_ext(".gz.lsfout")
      cmd = <<-CMD
        bsub \\
          -g /germ/sv \\
          -q i2b2_12h \\
          -R "rusage[mem=5500]" -M 5500000 \\
          -o #{lsfout} \\
          #{$java7_bin} -XX:+UseSerialGC -Xmx5G -jar #{$gatk3_bin} \\
            -T SelectVariants \\
            -R #{$refseq} \\
            --variant #{vqsr_vcf} \\
            -o #{vqsr_chr_vcf} \\
            -L #{chr} \\
            --excludeFiltered \\
            -select "AN >= 100"
      CMD
      submit cmd
    end
  end
end

def filter_cancer_chr_vcfs
  cancer_type = "BLCA"
  vqsr_chr_vcfs = Pathname.glob("#{$vcf_dir}/CANCER/#{cancer_type}/*/chrs/*/*vqsr.pass.an100.chr*.vcf.gz")
  vqsr_chr_vcfs.each do |vqsr_chr_vcf|
    vqsr_chr_filt_vcf = vqsr_chr_vcf.dirname + vqsr_chr_vcf.basename(".gz").sub_ext(".maf05.vcf.gz")
    lsfout = vqsr_chr_filt_vcf.sub_ext(".gz.lsfout")
        #-q short -W 12:0 \\
    cmd = <<-CMD
      bsub \\
        -g /germ/filter \\
        -q i2b2_12h \\
        -R "rusage[mem=5500]" -M 5500000 \\
        -o #{lsfout} "\\
        vcftools \\
          --gzvcf #{vqsr_chr_vcf} \\
          --recode \\
          --recode-INFO-all \\
          --remove-filtered-all \\
          --maf 0.05 \\
          --stdout | bgzip -c > #{vqsr_chr_filt_vcf} && tabix -p vcf #{vqsr_chr_filt_vcf}"
    CMD
    submit cmd
  end
end

def filter_pancan_chr_vcfs
  vqsr_chr_vcfs = Pathname.glob("#{$vcf_dir}/PANCAN/*/chrs/*/*vqsr.pass.an100.chr*.vcf.gz")
  vqsr_chr_vcfs.each do |vqsr_chr_vcf|
    vqsr_chr_filt_vcf = vqsr_chr_vcf.dirname + vqsr_chr_vcf.basename(".gz").sub_ext(".maf05.vcf.gz")
    lsfout = vqsr_chr_filt_vcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
      bsub \\
        -g /germ/filter \\
        -q short -W 12:0 \\
        -R "rusage[mem=5500]" -M 5500000 \\
        -o #{lsfout} "\\
        vcftools \\
          --gzvcf #{vqsr_chr_vcf} \\
          --recode \\
          --recode-INFO-all \\
          --remove-filtered-all \\
          --maf 0.05 \\
          --stdout | bgzip -c > #{vqsr_chr_filt_vcf} && tabix -p vcf #{vqsr_chr_filt_vcf}"
    CMD
    submit cmd
  end
end


def split_vqsr_cancer_chr_vcfs_into_even_sized_chunks
  chunk_size = 20000
  cancer_type = "BLCA"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  vqsr_chr_vcfs = Pathname.glob("#{cancer_vcf_dir}/{snp,indel}/chrs/*/*.pass.an100.chr*.vcf.gz").sort
  vqsr_chr_vcfs.each do |vqsr_chr_vcf|
    chunk_dir = vqsr_chr_vcf.dirname + "chunks"; chunk_dir.mkpath
    chunk_stem = chunk_dir + vqsr_chr_vcf.basename.sub_ext(".chunk")
    lsfout = chunk_stem.sub_ext(".chunks.lsfout")
    cmd = <<-CMD
      bsub \\
        -g /germ/split \\
        -q short -W 12:0 \\
        -o #{lsfout} \\
        'zgrep -vP "^#" #{vqsr_chr_vcf} | split -l #{chunk_size} -a 5 -d - #{chunk_stem}'
    CMD
    submit cmd
  end
end

def add_cancer_vcf_header_to_even_sized_chunks
  oriVcf = $vcf_dir + "LUAD/Harvard_GCC_WGS-LUAD-Normal.dbsnp.vcf.gz"
  oriVcfFh = Zlib::GzipReader.open(oriVcf)
  vcfHeader = oriVcf.sub_ext(".header")
  unless vcfHeader.exist?
    vcfHeader.open('w') do |file|
      oriVcfFh.each_line do |line|
        if line.start_with?("#")
          file.print(line)
        else
          exit
        end
      end
    end
  end
  vcfs = Pathname.glob($vcf_dir + "LUAD/{snp,indel}/chrs/*/chunks/*chunk[0-0][0-9][0-9][0-9][0-9]").sort
  vcfs.each do |vcf|
    headeredVcf = Pathname.new("#{vcf}.headered")
    cmd = "bsub -g /gcc/germ/header -q mini 'cat #{vcfHeader} #{vcf} > #{headeredVcf}'"
    submit cmd
  end
end

def remove_raw_even_sized_chunks
  vcfs = Pathname.glob($vcf_dir + "LUAD/{snp,indel}/chrs/*/chunks/*chunk[0-0][0-9][0-9][0-9][0-9]").sort
  vcfs.each { |vcf| vcf.delete }
end

def bgzip_headered_even_sized_cancer_vcf_chunks
  vcfs = Pathname.glob($vcf_dir + "LUAD/{snp,indel}/chrs/*/chunks/*chunk[0-0][0-9][0-9][0-9][0-9].headered").sort
  vcfs.each do |vcf|
    lsfout = Pathname.new("#{vcf}.gz.lsfout")
    cmd = <<-CMD
      bsub -g /gcc/germ/bgzip -q mini -o #{lsfout} 'bgzip #{vcf} && tabix -p vcf #{vcf}.gz'
    CMD
    submit cmd
  end
end

def annotate_cancer_chr_vcfs_with_dbsnp_rsids
  cancer_type = "BLCA"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  #dbsnp = $dbsnp_dir + "human_9606_b142_GRCh37p13/VCF/All.rsid.exp.vcf.gz"
  dbsnp = $dbsnp_dir + "human_9606_b144_GRCh37p13/VCF/00-All.rsid.exp.vcf.gz"
  vcfs = Pathname.glob("#{cancer_vcf_dir}/{snp,indel}/chrs/*/*.maf05.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    dbsnp_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".rsid.vcf.gz")
    lsfout = dbsnp_vcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -g /germ/rsid \\
      -q i2b2_12h \\
      -o #{lsfout} \\
      "zcat #{vcf} | vcf-annotate -a #{dbsnp} -d key=INFO,ID=DBSNP,Number=A,Type=String,Description=\\"dbSNP144 RSID\\" -c CHROM,POS,INFO/DBSNP,REF,ALT | bgzip -c > #{dbsnp_vcf} && tabix -p vcf #{dbsnp_vcf}"
    CMD
    submit cmd
  end
end

def annotate_pancan_chr_vcfs_with_dbsnp_rsids
  pancan_vcf_dir = $vcf_dir + "PANCAN"
  dbsnp = $dbsnp_dir + "human_9606_b144_GRCh37p13/VCF/00-All.rsid.exp.vcf.gz"
  vcfs = Pathname.glob("#{pancan_vcf_dir}/{snp,indel}/chrs/*/*.maf05.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    dbsnp_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".rsid.vcf.gz")
    lsfout = dbsnp_vcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -g /germ/rsid \\
      -q short -W 12:0 \\
      -o #{lsfout} \\
      "zcat #{vcf} | vcf-annotate -a #{dbsnp} -d key=INFO,ID=DBSNP,Number=A,Type=String,Description=\\"dbSNP144 RSID\\" -c CHROM,POS,INFO/DBSNP,REF,ALT | bgzip -c > #{dbsnp_vcf} && tabix -p vcf #{dbsnp_vcf}"
    CMD
    submit cmd
  end
end

def annotate_cancer_chr_vcfs_with_1000genomes_afs
  cancer_type = "BLCA"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  tgaf = $base_dir + "1000Genomes" + "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.afs.exp.tsv.gz"
  vcfs = Pathname.glob("#{cancer_vcf_dir}/{snp,indel}/chrs/*/*.maf05.rsid.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    tg_all_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".tg_all.vcf.gz")
    tg_eas_vcf = tg_all_vcf.dirname + tg_all_vcf.basename(".gz").sub_ext(".tg_eas.vcf.gz")
    tg_eur_vcf = tg_all_vcf.dirname + tg_eas_vcf.basename(".gz").sub_ext(".tg_eur.vcf.gz")
    tg_afr_vcf = tg_all_vcf.dirname + tg_eur_vcf.basename(".gz").sub_ext(".tg_afr.vcf.gz")
    tg_amr_vcf = tg_all_vcf.dirname + tg_afr_vcf.basename(".gz").sub_ext(".tg_amr.vcf.gz")
    tg_sas_vcf = tg_all_vcf.dirname + tg_amr_vcf.basename(".gz").sub_ext(".tg_sas.vcf.gz")
    final_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".tg.vcf.gz")
    lsfout = final_vcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -g /germ/tg \\
      -q i2b2_12h \\
      -o #{lsfout} "
      zcat #{vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_ALL_AF,Number=A,Type=Float,Description=\\"Global Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,INFO/TG_ALL_AF,-,-,-,-,- | bgzip -c > #{tg_all_vcf} &&
      tabix -p vcf #{tg_all_vcf} &&
      zcat #{tg_all_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_EAS_AF,Number=A,Type=Float,Description=\\"EAS Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,INFO/TG_EAS_AF,-,-,-,- | bgzip -c > #{tg_eas_vcf} &&
      tabix -p vcf #{tg_eas_vcf} &&
      zcat #{tg_eas_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_EUR_AF,Number=A,Type=Float,Description=\\"EUR Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,INFO/TG_EUR_AF,-,-,- | bgzip -c > #{tg_eur_vcf} &&
      tabix -p vcf #{tg_eur_vcf} &&
      zcat #{tg_eur_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_AFR_AF,Number=A,Type=Float,Description=\\"AFR Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,-,INFO/TG_AFR_AF,-,- | bgzip -c > #{tg_afr_vcf} &&
      tabix -p vcf #{tg_afr_vcf} &&
      zcat #{tg_afr_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_AMR_AF,Number=A,Type=Float,Description=\\"AMR Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,-,-,INFO/TG_AMR_AF,- | bgzip -c > #{tg_amr_vcf} &&
      tabix -p vcf #{tg_amr_vcf} &&
      zcat #{tg_amr_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_SAS_AF,Number=A,Type=Float,Description=\\"SAS Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,-,-,-,INFO/TG_SAS_AF | bgzip -c > #{tg_sas_vcf} &&
      rm -rf #{tg_all_vcf} #{tg_eas_vcf} #{tg_eur_vcf} #{tg_afr_vcf} #{tg_amr_vcf} &&
      rm -rf #{tg_all_vcf}.tbi #{tg_eas_vcf}.tbi #{tg_eur_vcf}.tbi #{tg_afr_vcf}.tbi #{tg_amr_vcf}.tbi &&
      mv #{tg_sas_vcf} #{final_vcf} &&
      tabix -p vcf #{final_vcf}
      "
    CMD
    submit cmd
  end
end

def annotate_pancan_chr_vcfs_with_1000genomes_afs
  pancan_vcf_dir = $vcf_dir + "PANCAN"
  tgaf = $base_dir + "1000Genomes" + "ALL.wgs.phase3_shapeit2_mvncall_integrated_v5a.20130502.sites.afs.exp.tsv.gz"
  vcfs = Pathname.glob("#{pancan_vcf_dir}/{snp,indel}/chrs/*/*.maf05.rsid.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    tg_all_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".tg_all.vcf.gz")
    tg_eas_vcf = tg_all_vcf.dirname + tg_all_vcf.basename(".gz").sub_ext(".tg_eas.vcf.gz")
    tg_eur_vcf = tg_all_vcf.dirname + tg_eas_vcf.basename(".gz").sub_ext(".tg_eur.vcf.gz")
    tg_afr_vcf = tg_all_vcf.dirname + tg_eur_vcf.basename(".gz").sub_ext(".tg_afr.vcf.gz")
    tg_amr_vcf = tg_all_vcf.dirname + tg_afr_vcf.basename(".gz").sub_ext(".tg_amr.vcf.gz")
    tg_sas_vcf = tg_all_vcf.dirname + tg_amr_vcf.basename(".gz").sub_ext(".tg_sas.vcf.gz")
    final_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".tg.vcf.gz")
    lsfout = final_vcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -g /germ/tg \\
      -q short -W 12:0 \\
      -o #{lsfout} "
      zcat #{vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_ALL_AF,Number=A,Type=Float,Description=\\"Global Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,INFO/TG_ALL_AF,-,-,-,-,- | bgzip -c > #{tg_all_vcf} &&
      tabix -p vcf #{tg_all_vcf} &&
      zcat #{tg_all_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_EAS_AF,Number=A,Type=Float,Description=\\"EAS Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,INFO/TG_EAS_AF,-,-,-,- | bgzip -c > #{tg_eas_vcf} &&
      tabix -p vcf #{tg_eas_vcf} &&
      zcat #{tg_eas_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_EUR_AF,Number=A,Type=Float,Description=\\"EUR Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,INFO/TG_EUR_AF,-,-,- | bgzip -c > #{tg_eur_vcf} &&
      tabix -p vcf #{tg_eur_vcf} &&
      zcat #{tg_eur_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_AFR_AF,Number=A,Type=Float,Description=\\"AFR Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,-,INFO/TG_AFR_AF,-,- | bgzip -c > #{tg_afr_vcf} &&
      tabix -p vcf #{tg_afr_vcf} &&
      zcat #{tg_afr_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_AMR_AF,Number=A,Type=Float,Description=\\"AMR Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,-,-,INFO/TG_AMR_AF,- | bgzip -c > #{tg_amr_vcf} &&
      tabix -p vcf #{tg_amr_vcf} &&
      zcat #{tg_amr_vcf} | vcf-annotate -a #{tgaf} -d key=INFO,ID=TG_SAS_AF,Number=A,Type=Float,Description=\\"SAS Allele Frequency from 1000 Genomes Project\\" -c CHROM,POS,REF,ALT,-,-,-,-,-,INFO/TG_SAS_AF | bgzip -c > #{tg_sas_vcf} &&
      rm -rf #{tg_all_vcf} #{tg_eas_vcf} #{tg_eur_vcf} #{tg_afr_vcf} #{tg_amr_vcf} &&
      rm -rf #{tg_all_vcf}.tbi #{tg_eas_vcf}.tbi #{tg_eur_vcf}.tbi #{tg_afr_vcf}.tbi #{tg_amr_vcf}.tbi &&
      mv #{tg_sas_vcf} #{final_vcf} &&
      tabix -p vcf #{final_vcf}
      "
    CMD
    submit cmd
  end
end

def annotate_cancer_chr_vcfs_with_esp6500si_afs
  cancer_type = "BLCA"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  esp_af_file = $base_dir + "NHLBI" + "ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.afs.sorted.txt.gz"
  vcfs = Pathname.glob(cancer_vcf_dir + "{snp,indel}/chrs/*/*.tg.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    esp_ta_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".esp_ta.vcf.gz")
    esp_ea_vcf = esp_ta_vcf.dirname + esp_ta_vcf.basename(".gz").sub_ext(".esp_ea.vcf.gz")
    esp_aa_vcf = esp_ta_vcf.dirname + esp_ta_vcf.basename(".gz").sub_ext(".esp_aa.vcf.gz")
    final_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".esp.vcf.gz")
    lsfout = final_vcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -g /gcc/germ/esp \\
      -q i2b2_12h \\
      -o #{lsfout} "
      zcat #{vcf} | vcf-annotate -a #{esp_af_file} -d key=INFO,ID=TA_AF,Number=A,Type=Float,Description=\\"Total American Allele Frequency from NHLBI ESP6500SI-V2-SSA137\\" -c CHROM,POS,REF,ALT,-,-,INFO/TA_AF | bgzip -c > #{esp_ta_vcf} &&
      tabix -p vcf #{esp_ta_vcf} &&
      zcat #{esp_ta_vcf} | vcf-annotate -a #{esp_af_file} -d key=INFO,ID=EA_AF,Number=A,Type=Float,Description=\\"European American Allele Frequency from NHLBI ESP6500SI-V2-SSA137\\" -c CHROM,POS,REF,ALT,INFO/EA_AF,-,- | bgzip -c > #{esp_ea_vcf} &&
      tabix -p vcf #{esp_ea_vcf} &&
      zcat #{esp_ea_vcf} | vcf-annotate -a #{esp_af_file} -d key=INFO,ID=AA_AF,Number=A,Type=Float,Description=\\"African American Allele Frequency from NHLBI ESP6500SI-V2-SSA137\\" -c CHROM,POS,REF,ALT,-,INFO/AA_AF,- | bgzip -c > #{esp_aa_vcf} &&
      rm -rf #{esp_ta_vcf} #{esp_ea_vcf} &&
      rm -rf #{esp_ta_vcf}.tbi #{esp_ea_vcf}.tbi &&
      mv #{esp_aa_vcf} #{final_vcf} &&
      tabix -p vcf #{final_vcf}
      "
    CMD
    submit cmd
  end
end

def annotate_pancan_chr_vcfs_with_esp6500si_afs
  pancan_vcf_dir = $vcf_dir + "PANCAN"
  esp_af_file = $base_dir + "NHLBI" + "ESP6500SI-V2-SSA137.GRCh38-liftover.snps_indels.afs.sorted.txt.gz"
  vcfs = Pathname.glob(pancan_vcf_dir + "{snp,indel}/chrs/*/*.tg.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    esp_ta_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".esp_ta.vcf.gz")
    esp_ea_vcf = esp_ta_vcf.dirname + esp_ta_vcf.basename(".gz").sub_ext(".esp_ea.vcf.gz")
    esp_aa_vcf = esp_ta_vcf.dirname + esp_ta_vcf.basename(".gz").sub_ext(".esp_aa.vcf.gz")
    final_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".esp.vcf.gz")
    lsfout = final_vcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -g /germ/esp \\
      -q short -W 12:0 \\
      -o #{lsfout} "
      zcat #{vcf} | vcf-annotate -a #{esp_af_file} -d key=INFO,ID=TA_AF,Number=A,Type=Float,Description=\\"Total American Allele Frequency from NHLBI ESP6500SI-V2-SSA137\\" -c CHROM,POS,REF,ALT,-,-,INFO/TA_AF | bgzip -c > #{esp_ta_vcf} &&
      tabix -p vcf #{esp_ta_vcf} &&
      zcat #{esp_ta_vcf} | vcf-annotate -a #{esp_af_file} -d key=INFO,ID=EA_AF,Number=A,Type=Float,Description=\\"European American Allele Frequency from NHLBI ESP6500SI-V2-SSA137\\" -c CHROM,POS,REF,ALT,INFO/EA_AF,-,- | bgzip -c > #{esp_ea_vcf} &&
      tabix -p vcf #{esp_ea_vcf} &&
      zcat #{esp_ea_vcf} | vcf-annotate -a #{esp_af_file} -d key=INFO,ID=AA_AF,Number=A,Type=Float,Description=\\"African American Allele Frequency from NHLBI ESP6500SI-V2-SSA137\\" -c CHROM,POS,REF,ALT,-,INFO/AA_AF,- | bgzip -c > #{esp_aa_vcf} &&
      rm -rf #{esp_ta_vcf} #{esp_ea_vcf} &&
      rm -rf #{esp_ta_vcf}.tbi #{esp_ea_vcf}.tbi &&
      mv #{esp_aa_vcf} #{final_vcf} &&
      tabix -p vcf #{final_vcf}
      "
    CMD
    submit cmd
  end
end

def annotate_cancer_snp_chr_vcfs_with_cadd_scores
  cancer_type = "BLCA"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  cadd = $base_dir + "CADD/v1.3/whole_genome_SNVs.tsv.gz"
  vcfs = Pathname.glob(cancer_vcf_dir + "snp/chrs/*/*.esp.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    cadd_raw_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".cadd_raw.vcf.gz")
    cadd_phred_vcf = vcf.dirname + cadd_raw_vcf.basename(".gz").sub_ext(".cadd_phred.vcf.gz")
    final_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".cadd.vcf.gz")
    lsfout = final_vcf.sub_ext(".gz.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
    bsub \\
      -g /germ/cadd \\
      -q short -W 12:0 \\
      -o #{lsfout} \\
      "zcat #{vcf} | vcf-annotate -a #{cadd} -d key=INFO,ID=CADD_RAW,Number=A,Type=Float,Description=\\"CADD Raw Score\\" -c CHROM,POS,REF,ALT,INFO/CADD_RAW,- | bgzip -c > #{cadd_raw_vcf} &&
      tabix -p vcf #{cadd_raw_vcf} &&
      zcat #{cadd_raw_vcf} | vcf-annotate -a #{cadd} -d key=INFO,ID=CADD_PHRED,Number=A,Type=Float,Description=\\"CADD Phred Score\\" -c CHROM,POS,REF,ALT,-,INFO/CADD_PHRED | bgzip -c > #{cadd_phred_vcf}
      rm -rf #{cadd_raw_vcf} #{cadd_raw_vcf}.tbi &&
      mv #{cadd_phred_vcf} #{final_vcf} &&
      tabix -p vcf #{final_vcf}"
    CMD
    submit cmd
  end
end

def annotate_pancan_snp_chr_vcfs_with_cadd_scores
  pancan_vcf_dir = $vcf_dir + "PANCAN"
  cadd = $base_dir + "CADD/v1.3/whole_genome_SNVs.tsv.gz"
  vcfs = Pathname.glob(pancan_vcf_dir + "snp/chrs/*/*.esp.vcf.gz").sort
  Parallel.each_with_index(vcfs, :in_threads => 4) do |vcf, vi|
    cadd_raw_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".cadd_raw.vcf.gz")
    cadd_phred_vcf = vcf.dirname + cadd_raw_vcf.basename(".gz").sub_ext(".cadd_phred.vcf.gz")
    final_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".cadd.vcf.gz")
    lsfout = final_vcf.sub_ext(".gz.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
    bsub \\
      -g /germ/cadd \\
      -q i2b2_7d \\
      -o #{lsfout} \\
      "zcat #{vcf} | vcf-annotate -a #{cadd} -d key=INFO,ID=CADD_RAW,Number=A,Type=Float,Description=\\"CADD Raw Score\\" -c CHROM,POS,REF,ALT,INFO/CADD_RAW,- | bgzip -c > #{cadd_raw_vcf} &&
      tabix -p vcf #{cadd_raw_vcf} &&
      zcat #{cadd_raw_vcf} | vcf-annotate -a #{cadd} -d key=INFO,ID=CADD_PHRED,Number=A,Type=Float,Description=\\"CADD Phred Score\\" -c CHROM,POS,REF,ALT,-,INFO/CADD_PHRED | bgzip -c > #{cadd_phred_vcf}
      rm -rf #{cadd_raw_vcf} #{cadd_raw_vcf}.tbi &&
      mv #{cadd_phred_vcf} #{final_vcf} &&
      tabix -p vcf #{final_vcf}"
    CMD
    submit cmd
  end
end


def annotate_cancer_chr_vcfs_with_vep
  num_cores = 4
  vep_bin = Pathname.new "/home/sl279/vep/variant_effect_predictor.pl"
  vep_dir = Pathname.new "/home/sl279/.vep"
  encode_dcc_dir = $base_dir + "UCSC/encodeDCC"
  ucsc_db_dir = $base_dir + "UCSC/database"
  fantom5_dir = $base_dir + "FANTOM5"
  dbsuper_dir = $base_dir + "dbSUPER"
  fdgenome_dir = $base_dir + "4DGenome"
  insituhic_dir = $base_dir + "HiC/GSE63525/suppl"
  roadmap_dir = $base_dir + "Roadmap"
  mtab2323_dir = $base_dir + "E-MTAB-2323"
  snp_vcfs = Pathname.glob($vcf_dir + "CANCER/BLCA/snp/chrs/*/*.cadd.vcf.gz").sort
  puts snp_vcfs.size
  indel_vcfs = Pathname.glob($vcf_dir + "CANCER/BLCA/indel/chrs/*/*.esp.vcf.gz").sort
  vcfs = [snp_vcfs, indel_vcfs].flatten
  Parallel.each(vcfs, :in_threads => 10) do |vcf|
  #vcfs.each do |vcf|
    vep_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".vep.vcf.gz")
    lsfout = vep_vcf.sub_ext(".gz.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
    bsub \\
      -g /germ/vep \\
      -q mcore -W 700:0 \\
      -n #{num_cores} \\
      -R "rusage[mem=10000] span[hosts=1]" -M 10000000 \\
      -o #{lsfout} \\
      "perl #{vep_bin} \\
        --force_overwrite \\
        --fork #{num_cores} \\
        --offline \\
        --no_stats \\
        --everything \\
        --check_existing \\
        --total_length \\
        --allele_number \\
        --no_escape \\
        --per_gene \\
        --symbol \\
        --gencode_basic \\
        --vcf \\
        --assembly GRCh37 \\
        --dir #{vep_dir} \\
        --fasta #{vep_dir}/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \\
        --custom #{encode_dcc_dir}/distalDhsToPromoterDhs.bed.gz,distalDhsToPromoterDhs,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/dhsToGeneExpression.bed.gz,dhsToGeneExpression,bed,overlap,0 \\
        --custom #{fantom5_dir}/fantom5EnhancerTssAssociations.bed.gz,fantom5EnhancerTssAssociations,bed,overlap,0 \\
        --custom #{dbsuper_dir}/all/allDbSuperEnhancerGeneAssociations.bed.gz,allDbSuperEnhancerGeneAssociations,bed,overlap,0 \\
        --custom #{fdgenome_dir}/4dGenomeHomoSapiens.bed.gz,fourDGenomeHomoSapiens,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.annotated.bed.gz,GSE63525_GM12878_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_HeLa_HiCCUPS_looplist.annotated.bed.gz,GSE63525_HeLa_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_HMEC_HiCCUPS_looplist.annotated.bed.gz,GSE63525_HMEC_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_HUVEC_HiCCUPS_looplist.annotated.bed.gz,GSE63525_HUVEC_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_IMR90_HiCCUPS_looplist.annotated.bed.gz,GSE63525_IMR90_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_K562_HiCCUPS_looplist.annotated.bed.gz,GSE63525_K562_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_KBM7_HiCCUPS_looplist.annotated.bed.gz,GSE63525_KBM7_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_NHEK_HiCCUPS_looplist.annotated.bed.gz,GSE63525_NHEK_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_NHEK_HiCCUPS_looplist.annotated.bed.gz,GSE63525_NHEK_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_GM12878_100kb_inter_MAPQGE30_SQRTVC.bed.gz,GSE63525_GM12878_100kb_inter_MAPQGE30_SQRTVC,bed,overlap,0 \\
        --custom #{mtab2323_dir}/TS5_CD34_promoter-other_significant_interactions.bed.gz,TS5_CD34_promoter-other_significant_interactions,bed,overlap,0 \\
        --custom #{mtab2323_dir}/TS5_CD34_promoter-promoter_significant_interactions.bed.gz,TS5_CD34_promoter-promoter_significant_interactions,bed,overlap,0 \\
        --custom #{mtab2323_dir}/TS5_GM12878_promoter-other_significant_interactions.bed.gz,TS5_GM12878_promoter-other_significant_interactions,bed,overlap,0 \\
        --custom #{mtab2323_dir}/TS5_GM12878_promoter-promoter_significant_interactions.bed.gz,TS5_GM12878_promoter-promoter_significant_interactions,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgDnaseMasterSites.bed.gz,wgEncodeAwgDnaseMasterSites,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeRegDnaseClusteredV3.bed.gz,wgEncodeRegDnaseClusteredV3,bed,overlap,0 \\
        --custom #{fantom5_dir}/fantom5PermissiveEnhancers.bed.gz,fantom5PermissiveEnhancers,bed,overlap,0 \\
        --custom #{dbsuper_dir}/all/allDbSuperEnhancers.bed.gz,allDbSuperEnhancers,bed,overlap,0 \\
        --custom #{roadmap_dir}/dnase/roadmapDnaseDyadic.bed.gz,roadmapDnaseDyadic,bed,overlap,0 \\
        --custom #{roadmap_dir}/dnase/roadmapDnaseEnh.bed.gz,roadmapDnaseEnh,bed,overlap,0 \\
        --custom #{roadmap_dir}/dnase/roadmapDnaseProm.bed.gz,roadmapDnaseProm,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz,wgEncodeRegTfbsClusteredWithCellsV3,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedGm12878.bed.gz,wgEncodeAwgSegmentationCombinedGm12878,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedH1hesc.bed.gz,wgEncodeAwgSegmentationChromhmmH1hesc,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedHelas3.bed.gz,wgEncodeAwgSegmentationCombinedHelas3,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedHepg2.bed.gz,wgEncodeAwgSegmentationCombinedHepg2,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedHuvec.bed.gz,wgEncodeAwgSegmentationCombinedHuvec,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedK562.bed.gz,wgEncodeAwgSegmentationCombinedK562,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E001_15_coreMarks_mnemonics.bed.gz,E001_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E002_15_coreMarks_mnemonics.bed.gz,E002_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E003_15_coreMarks_mnemonics.bed.gz,E003_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E004_15_coreMarks_mnemonics.bed.gz,E004_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E005_15_coreMarks_mnemonics.bed.gz,E005_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E006_15_coreMarks_mnemonics.bed.gz,E006_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E007_15_coreMarks_mnemonics.bed.gz,E007_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E008_15_coreMarks_mnemonics.bed.gz,E008_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E009_15_coreMarks_mnemonics.bed.gz,E009_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E010_15_coreMarks_mnemonics.bed.gz,E010_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E011_15_coreMarks_mnemonics.bed.gz,E011_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E012_15_coreMarks_mnemonics.bed.gz,E012_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E013_15_coreMarks_mnemonics.bed.gz,E013_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E014_15_coreMarks_mnemonics.bed.gz,E014_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E015_15_coreMarks_mnemonics.bed.gz,E015_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E016_15_coreMarks_mnemonics.bed.gz,E016_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E017_15_coreMarks_mnemonics.bed.gz,E017_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E018_15_coreMarks_mnemonics.bed.gz,E018_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E019_15_coreMarks_mnemonics.bed.gz,E019_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E020_15_coreMarks_mnemonics.bed.gz,E020_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E021_15_coreMarks_mnemonics.bed.gz,E021_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E022_15_coreMarks_mnemonics.bed.gz,E022_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E023_15_coreMarks_mnemonics.bed.gz,E023_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E024_15_coreMarks_mnemonics.bed.gz,E024_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E025_15_coreMarks_mnemonics.bed.gz,E025_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E026_15_coreMarks_mnemonics.bed.gz,E026_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E027_15_coreMarks_mnemonics.bed.gz,E027_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E028_15_coreMarks_mnemonics.bed.gz,E028_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E029_15_coreMarks_mnemonics.bed.gz,E029_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E030_15_coreMarks_mnemonics.bed.gz,E030_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E031_15_coreMarks_mnemonics.bed.gz,E031_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E032_15_coreMarks_mnemonics.bed.gz,E032_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E033_15_coreMarks_mnemonics.bed.gz,E033_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E034_15_coreMarks_mnemonics.bed.gz,E034_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E035_15_coreMarks_mnemonics.bed.gz,E035_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E036_15_coreMarks_mnemonics.bed.gz,E036_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E037_15_coreMarks_mnemonics.bed.gz,E037_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E038_15_coreMarks_mnemonics.bed.gz,E038_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E039_15_coreMarks_mnemonics.bed.gz,E039_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E040_15_coreMarks_mnemonics.bed.gz,E040_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E041_15_coreMarks_mnemonics.bed.gz,E041_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E042_15_coreMarks_mnemonics.bed.gz,E042_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E043_15_coreMarks_mnemonics.bed.gz,E043_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E044_15_coreMarks_mnemonics.bed.gz,E044_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E045_15_coreMarks_mnemonics.bed.gz,E045_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E046_15_coreMarks_mnemonics.bed.gz,E046_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E047_15_coreMarks_mnemonics.bed.gz,E047_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E048_15_coreMarks_mnemonics.bed.gz,E048_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E049_15_coreMarks_mnemonics.bed.gz,E049_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E050_15_coreMarks_mnemonics.bed.gz,E050_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E051_15_coreMarks_mnemonics.bed.gz,E051_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E052_15_coreMarks_mnemonics.bed.gz,E052_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E053_15_coreMarks_mnemonics.bed.gz,E053_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E054_15_coreMarks_mnemonics.bed.gz,E054_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E055_15_coreMarks_mnemonics.bed.gz,E055_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E056_15_coreMarks_mnemonics.bed.gz,E056_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E057_15_coreMarks_mnemonics.bed.gz,E057_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E058_15_coreMarks_mnemonics.bed.gz,E058_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E059_15_coreMarks_mnemonics.bed.gz,E059_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E061_15_coreMarks_mnemonics.bed.gz,E061_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E062_15_coreMarks_mnemonics.bed.gz,E062_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E063_15_coreMarks_mnemonics.bed.gz,E063_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E065_15_coreMarks_mnemonics.bed.gz,E065_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E066_15_coreMarks_mnemonics.bed.gz,E066_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E067_15_coreMarks_mnemonics.bed.gz,E067_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E068_15_coreMarks_mnemonics.bed.gz,E068_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E069_15_coreMarks_mnemonics.bed.gz,E069_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E070_15_coreMarks_mnemonics.bed.gz,E070_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E071_15_coreMarks_mnemonics.bed.gz,E071_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E072_15_coreMarks_mnemonics.bed.gz,E072_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E073_15_coreMarks_mnemonics.bed.gz,E073_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E074_15_coreMarks_mnemonics.bed.gz,E074_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E075_15_coreMarks_mnemonics.bed.gz,E075_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E076_15_coreMarks_mnemonics.bed.gz,E076_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E077_15_coreMarks_mnemonics.bed.gz,E077_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E078_15_coreMarks_mnemonics.bed.gz,E078_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E079_15_coreMarks_mnemonics.bed.gz,E079_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E080_15_coreMarks_mnemonics.bed.gz,E080_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E081_15_coreMarks_mnemonics.bed.gz,E081_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E082_15_coreMarks_mnemonics.bed.gz,E082_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E083_15_coreMarks_mnemonics.bed.gz,E083_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E084_15_coreMarks_mnemonics.bed.gz,E084_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E085_15_coreMarks_mnemonics.bed.gz,E085_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E086_15_coreMarks_mnemonics.bed.gz,E086_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E087_15_coreMarks_mnemonics.bed.gz,E087_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E088_15_coreMarks_mnemonics.bed.gz,E088_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E089_15_coreMarks_mnemonics.bed.gz,E089_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E090_15_coreMarks_mnemonics.bed.gz,E090_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E091_15_coreMarks_mnemonics.bed.gz,E091_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E092_15_coreMarks_mnemonics.bed.gz,E092_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E093_15_coreMarks_mnemonics.bed.gz,E093_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E094_15_coreMarks_mnemonics.bed.gz,E094_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E095_15_coreMarks_mnemonics.bed.gz,E095_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E096_15_coreMarks_mnemonics.bed.gz,E096_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E097_15_coreMarks_mnemonics.bed.gz,E097_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E098_15_coreMarks_mnemonics.bed.gz,E098_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E099_15_coreMarks_mnemonics.bed.gz,E099_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E100_15_coreMarks_mnemonics.bed.gz,E100_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E101_15_coreMarks_mnemonics.bed.gz,E101_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E102_15_coreMarks_mnemonics.bed.gz,E102_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E103_15_coreMarks_mnemonics.bed.gz,E103_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E104_15_coreMarks_mnemonics.bed.gz,E104_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E105_15_coreMarks_mnemonics.bed.gz,E105_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E106_15_coreMarks_mnemonics.bed.gz,E106_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E107_15_coreMarks_mnemonics.bed.gz,E107_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E108_15_coreMarks_mnemonics.bed.gz,E108_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E109_15_coreMarks_mnemonics.bed.gz,E109_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E110_15_coreMarks_mnemonics.bed.gz,E110_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E111_15_coreMarks_mnemonics.bed.gz,E111_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E112_15_coreMarks_mnemonics.bed.gz,E112_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E113_15_coreMarks_mnemonics.bed.gz,E113_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E114_15_coreMarks_mnemonics.bed.gz,E114_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E115_15_coreMarks_mnemonics.bed.gz,E115_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E116_15_coreMarks_mnemonics.bed.gz,E116_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E117_15_coreMarks_mnemonics.bed.gz,E117_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E118_15_coreMarks_mnemonics.bed.gz,E118_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E119_15_coreMarks_mnemonics.bed.gz,E119_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E120_15_coreMarks_mnemonics.bed.gz,E120_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E121_15_coreMarks_mnemonics.bed.gz,E121_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E122_15_coreMarks_mnemonics.bed.gz,E122_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E123_15_coreMarks_mnemonics.bed.gz,E123_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E124_15_coreMarks_mnemonics.bed.gz,E124_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E125_15_coreMarks_mnemonics.bed.gz,E125_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E126_15_coreMarks_mnemonics.bed.gz,E126_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E127_15_coreMarks_mnemonics.bed.gz,E127_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E128_15_coreMarks_mnemonics.bed.gz,E128_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E129_15_coreMarks_mnemonics.bed.gz,E129_15_coreMarks,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/gwasCatalog.bed.gz,gwasCatalog,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/targetScanS.bed.gz,targetScanS,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/tfbsConsSites.bed.gz,tfbsConsSites,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/wgRna.bed.gz,wgRna,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/phastConsElements46way.bed.gz,phastConsElements46way,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/phastConsElements100way.bed.gz,phastConsElements100way,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/gap.bed.gz,gap,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/genomicSuperDups.bed.gz,genomicSuperDups,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/rmsk.bed.gz,rmsk,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/nestedRepeats.bed.gz,nestedRepeats,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/simpleRepeat.bed.gz,simpleRepeat,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/microsat.bed.gz,microsat,bed,overlap,0 \\
        --input_file #{vcf} \\
        --output_file STDOUT | bgzip -c > #{vep_vcf} && tabix -p vcf #{vep_vcf}"
    CMD
    submit cmd
  end
end

def annotate_pancan_chr_vcfs_with_vep
  num_cores = 4
  vep_bin = Pathname.new "/home/sl279/vep/variant_effect_predictor.pl"
  vep_dir = Pathname.new "/home/sl279/.vep"
  encode_dcc_dir = $base_dir + "UCSC/encodeDCC"
  ucsc_db_dir = $base_dir + "UCSC/database"
  fantom5_dir = $base_dir + "FANTOM5"
  dbsuper_dir = $base_dir + "dbSUPER"
  fdgenome_dir = $base_dir + "4DGenome"
  insituhic_dir = $base_dir + "HiC/GSE63525/suppl"
  roadmap_dir = $base_dir + "Roadmap"
  mtab2323_dir = $base_dir + "E-MTAB-2323"
  #snp_vcfs = Pathname.glob($vcf_dir + "PANCAN/snp/chrs/*/*.cadd.vcf.gz").sort
  indel_vcfs = Pathname.glob($vcf_dir + "PANCAN/indel/chrs/*/*.esp.vcf.gz").sort
  #vcfs = [snp_vcfs, indel_vcfs].flatten
  vcfs = [indel_vcfs].flatten
  Parallel.each(vcfs, :in_threads => 10) do |vcf|
  #vcfs.each do |vcf|
    vep_vcf = vcf.dirname + vcf.basename(".gz").sub_ext(".vep.vcf.gz")
    lsfout = vep_vcf.sub_ext(".gz.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
    bsub \\
      -g /germ/vep \\
      -q mcore -W 700:0 \\
      -n #{num_cores} \\
      -R "rusage[mem=10000] span[hosts=1]" -M 10000000 \\
      -o #{lsfout} \\
      "perl #{vep_bin} \\
        --fork #{num_cores} \\
        --force_overwrite \\
        --offline \\
        --no_stats \\
        --everything \\
        --check_existing \\
        --total_length \\
        --allele_number \\
        --no_escape \\
        --per_gene \\
        --symbol \\
        --gencode_basic \\
        --vcf \\
        --assembly GRCh37 \\
        --dir #{vep_dir} \\
        --fasta #{vep_dir}/homo_sapiens/81_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \\
        --custom #{encode_dcc_dir}/distalDhsToPromoterDhs.bed.gz,distalDhsToPromoterDhs,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/dhsToGeneExpression.bed.gz,dhsToGeneExpression,bed,overlap,0 \\
        --custom #{fantom5_dir}/fantom5EnhancerTssAssociations.bed.gz,fantom5EnhancerTssAssociations,bed,overlap,0 \\
        --custom #{dbsuper_dir}/all/allDbSuperEnhancerGeneAssociations.bed.gz,allDbSuperEnhancerGeneAssociations,bed,overlap,0 \\
        --custom #{fdgenome_dir}/4dGenomeHomoSapiens.bed.gz,fourDGenomeHomoSapiens,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_GM12878_primary+replicate_HiCCUPS_looplist.annotated.bed.gz,GSE63525_GM12878_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_HeLa_HiCCUPS_looplist.annotated.bed.gz,GSE63525_HeLa_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_HMEC_HiCCUPS_looplist.annotated.bed.gz,GSE63525_HMEC_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_HUVEC_HiCCUPS_looplist.annotated.bed.gz,GSE63525_HUVEC_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_IMR90_HiCCUPS_looplist.annotated.bed.gz,GSE63525_IMR90_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_K562_HiCCUPS_looplist.annotated.bed.gz,GSE63525_K562_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_KBM7_HiCCUPS_looplist.annotated.bed.gz,GSE63525_KBM7_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_NHEK_HiCCUPS_looplist.annotated.bed.gz,GSE63525_NHEK_HiCCUPS_looplist,bed,overlap,0 \\
        --custom #{insituhic_dir}/GSE63525_GM12878_100kb_inter_MAPQGE30_SQRTVC.bed.gz,GSE63525_GM12878_100kb_inter_MAPQGE30_SQRTVC,bed,overlap,0 \\
        --custom #{mtab2323_dir}/TS5_CD34_promoter-other_significant_interactions.bed.gz,TS5_CD34_promoter-other_significant_interactions,bed,overlap,0 \\
        --custom #{mtab2323_dir}/TS5_CD34_promoter-promoter_significant_interactions.bed.gz,TS5_CD34_promoter-promoter_significant_interactions,bed,overlap,0 \\
        --custom #{mtab2323_dir}/TS5_GM12878_promoter-other_significant_interactions.bed.gz,TS5_GM12878_promoter-other_significant_interactions,bed,overlap,0 \\
        --custom #{mtab2323_dir}/TS5_GM12878_promoter-promoter_significant_interactions.bed.gz,TS5_GM12878_promoter-promoter_significant_interactions,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgDnaseMasterSites.bed.gz,wgEncodeAwgDnaseMasterSites,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeRegDnaseClusteredV3.bed.gz,wgEncodeRegDnaseClusteredV3,bed,overlap,0 \\
        --custom #{fantom5_dir}/fantom5PermissiveEnhancers.bed.gz,fantom5PermissiveEnhancers,bed,overlap,0 \\
        --custom #{dbsuper_dir}/all/allDbSuperEnhancers.bed.gz,allDbSuperEnhancers,bed,overlap,0 \\
        --custom #{roadmap_dir}/dnase/roadmapDnaseDyadic.bed.gz,roadmapDnaseDyadic,bed,overlap,0 \\
        --custom #{roadmap_dir}/dnase/roadmapDnaseEnh.bed.gz,roadmapDnaseEnh,bed,overlap,0 \\
        --custom #{roadmap_dir}/dnase/roadmapDnaseProm.bed.gz,roadmapDnaseProm,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeRegTfbsClusteredWithCellsV3.bed.gz,wgEncodeRegTfbsClusteredWithCellsV3,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedGm12878.bed.gz,wgEncodeAwgSegmentationCombinedGm12878,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedH1hesc.bed.gz,wgEncodeAwgSegmentationChromhmmH1hesc,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedHelas3.bed.gz,wgEncodeAwgSegmentationCombinedHelas3,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedHepg2.bed.gz,wgEncodeAwgSegmentationCombinedHepg2,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedHuvec.bed.gz,wgEncodeAwgSegmentationCombinedHuvec,bed,overlap,0 \\
        --custom #{encode_dcc_dir}/wgEncodeAwgSegmentationCombinedK562.bed.gz,wgEncodeAwgSegmentationCombinedK562,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E001_15_coreMarks_mnemonics.bed.gz,E001_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E002_15_coreMarks_mnemonics.bed.gz,E002_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E003_15_coreMarks_mnemonics.bed.gz,E003_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E004_15_coreMarks_mnemonics.bed.gz,E004_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E005_15_coreMarks_mnemonics.bed.gz,E005_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E006_15_coreMarks_mnemonics.bed.gz,E006_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E007_15_coreMarks_mnemonics.bed.gz,E007_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E008_15_coreMarks_mnemonics.bed.gz,E008_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E009_15_coreMarks_mnemonics.bed.gz,E009_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E010_15_coreMarks_mnemonics.bed.gz,E010_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E011_15_coreMarks_mnemonics.bed.gz,E011_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E012_15_coreMarks_mnemonics.bed.gz,E012_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E013_15_coreMarks_mnemonics.bed.gz,E013_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E014_15_coreMarks_mnemonics.bed.gz,E014_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E015_15_coreMarks_mnemonics.bed.gz,E015_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E016_15_coreMarks_mnemonics.bed.gz,E016_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E017_15_coreMarks_mnemonics.bed.gz,E017_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E018_15_coreMarks_mnemonics.bed.gz,E018_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E019_15_coreMarks_mnemonics.bed.gz,E019_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E020_15_coreMarks_mnemonics.bed.gz,E020_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E021_15_coreMarks_mnemonics.bed.gz,E021_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E022_15_coreMarks_mnemonics.bed.gz,E022_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E023_15_coreMarks_mnemonics.bed.gz,E023_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E024_15_coreMarks_mnemonics.bed.gz,E024_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E025_15_coreMarks_mnemonics.bed.gz,E025_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E026_15_coreMarks_mnemonics.bed.gz,E026_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E027_15_coreMarks_mnemonics.bed.gz,E027_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E028_15_coreMarks_mnemonics.bed.gz,E028_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E029_15_coreMarks_mnemonics.bed.gz,E029_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E030_15_coreMarks_mnemonics.bed.gz,E030_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E031_15_coreMarks_mnemonics.bed.gz,E031_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E032_15_coreMarks_mnemonics.bed.gz,E032_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E033_15_coreMarks_mnemonics.bed.gz,E033_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E034_15_coreMarks_mnemonics.bed.gz,E034_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E035_15_coreMarks_mnemonics.bed.gz,E035_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E036_15_coreMarks_mnemonics.bed.gz,E036_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E037_15_coreMarks_mnemonics.bed.gz,E037_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E038_15_coreMarks_mnemonics.bed.gz,E038_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E039_15_coreMarks_mnemonics.bed.gz,E039_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E040_15_coreMarks_mnemonics.bed.gz,E040_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E041_15_coreMarks_mnemonics.bed.gz,E041_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E042_15_coreMarks_mnemonics.bed.gz,E042_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E043_15_coreMarks_mnemonics.bed.gz,E043_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E044_15_coreMarks_mnemonics.bed.gz,E044_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E045_15_coreMarks_mnemonics.bed.gz,E045_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E046_15_coreMarks_mnemonics.bed.gz,E046_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E047_15_coreMarks_mnemonics.bed.gz,E047_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E048_15_coreMarks_mnemonics.bed.gz,E048_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E049_15_coreMarks_mnemonics.bed.gz,E049_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E050_15_coreMarks_mnemonics.bed.gz,E050_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E051_15_coreMarks_mnemonics.bed.gz,E051_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E052_15_coreMarks_mnemonics.bed.gz,E052_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E053_15_coreMarks_mnemonics.bed.gz,E053_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E054_15_coreMarks_mnemonics.bed.gz,E054_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E055_15_coreMarks_mnemonics.bed.gz,E055_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E056_15_coreMarks_mnemonics.bed.gz,E056_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E057_15_coreMarks_mnemonics.bed.gz,E057_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E058_15_coreMarks_mnemonics.bed.gz,E058_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E059_15_coreMarks_mnemonics.bed.gz,E059_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E061_15_coreMarks_mnemonics.bed.gz,E061_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E062_15_coreMarks_mnemonics.bed.gz,E062_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E063_15_coreMarks_mnemonics.bed.gz,E063_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E065_15_coreMarks_mnemonics.bed.gz,E065_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E066_15_coreMarks_mnemonics.bed.gz,E066_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E067_15_coreMarks_mnemonics.bed.gz,E067_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E068_15_coreMarks_mnemonics.bed.gz,E068_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E069_15_coreMarks_mnemonics.bed.gz,E069_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E070_15_coreMarks_mnemonics.bed.gz,E070_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E071_15_coreMarks_mnemonics.bed.gz,E071_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E072_15_coreMarks_mnemonics.bed.gz,E072_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E073_15_coreMarks_mnemonics.bed.gz,E073_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E074_15_coreMarks_mnemonics.bed.gz,E074_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E075_15_coreMarks_mnemonics.bed.gz,E075_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E076_15_coreMarks_mnemonics.bed.gz,E076_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E077_15_coreMarks_mnemonics.bed.gz,E077_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E078_15_coreMarks_mnemonics.bed.gz,E078_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E079_15_coreMarks_mnemonics.bed.gz,E079_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E080_15_coreMarks_mnemonics.bed.gz,E080_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E081_15_coreMarks_mnemonics.bed.gz,E081_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E082_15_coreMarks_mnemonics.bed.gz,E082_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E083_15_coreMarks_mnemonics.bed.gz,E083_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E084_15_coreMarks_mnemonics.bed.gz,E084_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E085_15_coreMarks_mnemonics.bed.gz,E085_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E086_15_coreMarks_mnemonics.bed.gz,E086_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E087_15_coreMarks_mnemonics.bed.gz,E087_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E088_15_coreMarks_mnemonics.bed.gz,E088_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E089_15_coreMarks_mnemonics.bed.gz,E089_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E090_15_coreMarks_mnemonics.bed.gz,E090_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E091_15_coreMarks_mnemonics.bed.gz,E091_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E092_15_coreMarks_mnemonics.bed.gz,E092_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E093_15_coreMarks_mnemonics.bed.gz,E093_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E094_15_coreMarks_mnemonics.bed.gz,E094_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E095_15_coreMarks_mnemonics.bed.gz,E095_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E096_15_coreMarks_mnemonics.bed.gz,E096_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E097_15_coreMarks_mnemonics.bed.gz,E097_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E098_15_coreMarks_mnemonics.bed.gz,E098_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E099_15_coreMarks_mnemonics.bed.gz,E099_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E100_15_coreMarks_mnemonics.bed.gz,E100_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E101_15_coreMarks_mnemonics.bed.gz,E101_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E102_15_coreMarks_mnemonics.bed.gz,E102_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E103_15_coreMarks_mnemonics.bed.gz,E103_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E104_15_coreMarks_mnemonics.bed.gz,E104_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E105_15_coreMarks_mnemonics.bed.gz,E105_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E106_15_coreMarks_mnemonics.bed.gz,E106_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E107_15_coreMarks_mnemonics.bed.gz,E107_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E108_15_coreMarks_mnemonics.bed.gz,E108_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E109_15_coreMarks_mnemonics.bed.gz,E109_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E110_15_coreMarks_mnemonics.bed.gz,E110_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E111_15_coreMarks_mnemonics.bed.gz,E111_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E112_15_coreMarks_mnemonics.bed.gz,E112_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E113_15_coreMarks_mnemonics.bed.gz,E113_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E114_15_coreMarks_mnemonics.bed.gz,E114_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E115_15_coreMarks_mnemonics.bed.gz,E115_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E116_15_coreMarks_mnemonics.bed.gz,E116_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E117_15_coreMarks_mnemonics.bed.gz,E117_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E118_15_coreMarks_mnemonics.bed.gz,E118_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E119_15_coreMarks_mnemonics.bed.gz,E119_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E120_15_coreMarks_mnemonics.bed.gz,E120_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E121_15_coreMarks_mnemonics.bed.gz,E121_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E122_15_coreMarks_mnemonics.bed.gz,E122_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E123_15_coreMarks_mnemonics.bed.gz,E123_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E124_15_coreMarks_mnemonics.bed.gz,E124_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E125_15_coreMarks_mnemonics.bed.gz,E125_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E126_15_coreMarks_mnemonics.bed.gz,E126_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E127_15_coreMarks_mnemonics.bed.gz,E127_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E128_15_coreMarks_mnemonics.bed.gz,E128_15_coreMarks,bed,overlap,0 \\
        --custom #{roadmap_dir}/chromhmm/E129_15_coreMarks_mnemonics.bed.gz,E129_15_coreMarks,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/gwasCatalog.bed.gz,gwasCatalog,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/targetScanS.bed.gz,targetScanS,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/tfbsConsSites.bed.gz,tfbsConsSites,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/wgRna.bed.gz,wgRna,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/phastConsElements46way.bed.gz,phastConsElements46way,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/phastConsElements100way.bed.gz,phastConsElements100way,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/gap.bed.gz,gap,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/genomicSuperDups.bed.gz,genomicSuperDups,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/rmsk.bed.gz,rmsk,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/nestedRepeats.bed.gz,nestedRepeats,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/simpleRepeat.bed.gz,simpleRepeat,bed,overlap,0 \\
        --custom #{ucsc_db_dir}/microsat.bed.gz,microsat,bed,overlap,0 \\
        --input_file #{vcf} \\
        --output_file STDOUT | bgzip -c > #{vep_vcf} && tabix -p vcf #{vep_vcf}"
    CMD
    submit cmd
  end
end

def extract_eur_cancer_chr_vcfs
  canTypes = %w[BLCA CESC CRC HNSC LGG LUAD PRAD SKCM STAD THCA UCEC]
  canTypes.each do |canType|
    eurPopFile = $vcf_dir + canType + "EUR_Samples.txt"
    eurPopFile.open('w') do |file|
      popFile = $base_dir + "Table" + "Continental_Ancestry_Prediction_For_TCGA_Individuals_With_Shared_SNPs-#{canType}.txt"
      popFile.readlines[1..-1].each do |line|
        cols = line.chomp.split("\t")
        pop = cols[-1]
        sampleId = cols[2].gsub(".", "-")
        if pop == "EUR"
          file.puts sampleId
        end
      end
    end
    vcfChrFiles = Pathname.glob("#{$vcf_dir}/#{canType}/{snp,indel}/chrs/*/*.vep.vcf.gz").sort
    vcfChrFiles.each do |vcfChrFile|
      vcfChrEurFile = vcfChrFile.dirname + vcfChrFile.basename(".gz").sub_ext(".eur.vcf.gz")
      lsfout = vcfChrEurFile.sub_ext(".gz.lsfout")
      next if lsfout.exist?
      cmd = <<-CMD
        bsub \\
          -q short -W 12:0 \\
          -g /gcc/germ/eur \\
          -R "rusage[mem=10000]" -M 10000000 \\
          -o #{lsfout} \\
          #{$java7_bin} -Xmx10g -jar #{$gatk3_bin} \\
            -T SelectVariants \\
            -R #{$refseq_broad} \\
            --variant #{vcfChrFile} \\
            --sample_file #{eurPopFile} \\
            -o #{vcfChrEurFile}
      CMD
      submit cmd
    end
  end
end

def extract_afs_from_eur_cancer_chr_vcfs
  vcfChrEurFiles = Pathname.glob("#{$vcf_dir}/*/{snp,indel}/chrs/*/*.vep.eur.vcf.gz").sort
  Parallel.each(vcfChrEurFiles, :in_threads => 10) do |vcfChrEurFile|
    puts vcfChrEurFile
    vcfChrFile = Pathname.new(vcfChrEurFile.to_s.gsub("eur.", ""))
    chrEurAfFile = vcfChrEurFile.dirname + vcfChrEurFile.basename(".gz").sub_ext(".afs.txt")
    vcfChrEurAfFile = vcfChrFile.dirname + vcfChrFile.basename(".gz").sub_ext(".eur_af.vcf.gz")
    lsfout = vcfChrEurFile.sub_ext(".txt.gz.lsfout")
    chrEurAfFile.open('w') do |file|
      file.puts %w[#CHROM POS REF ALT EUR_AF].join("\t")
      `zcat #{vcfChrEurFile}`.split("\n").each do |line|
          next if line.start_with?("#")
          cols = line.chomp.split("\t")
          chr = cols[0]
          pos = cols[1]
          ref = cols[3]
          alts = cols[4].split(",")
          info = cols[7]
          afs = info.match(/;AF=(\S+?);/)[1].split(",")
          alts.each_with_index do |alt, i|
            file.puts [chr, pos, ref, alt, afs[i]].join("\t")
          end
        end
    end
    system "bgzip -f #{chrEurAfFile}; tabix -s 1 -b 2 -e 2 #{chrEurAfFile}.gz"
    system "zcat #{vcfChrFile} | vcf-annotate -a #{chrEurAfFile}.gz -c CHROM,FROM,REF,ALT,INFO/EUR_AF -d key=INFO,ID=EUR_AF,Number=A,Type=Float,Description=\"EUR Allele Frequency, for each ALT allele, in the same order as listed\" | bgzip > #{vcfChrEurAfFile} && tabix -p vcf #{vcfChrEurAfFile}"
  end
end

def concat_chr_vep_vcfs
  require 'naturally'
  canType = "UCEC"
  %w[snp indel].each do |varType|
    vcfFiles = Pathname.glob($vcf_dir + "#{canType}/#{varType}/chrs/*/*.vep.eur_af.vcf.gz")
    vcfFiles = Naturally.sort(vcfFiles)
    outVcf = $vcf_dir + "#{canType}/#{varType}/#{vcfFiles[0].basename.to_s.gsub(/chr\S+?\./, "")}"
    lsfout = outVcf.sub_ext(".gz.lsfout")
    cmd = <<-CMD
      bsub \\
        -g /gcc/germ/concat \\
        -q short -W 12:0 \\
        -R "rusage[mem=10000]" -M 10000000 \\
        -o #{lsfout} \\
        "vcf-concat #{vcfFiles.join(" ")} | bgzip -c > #{outVcf} && tabix -p vcf #{outVcf}"
    CMD
    submit cmd
  end
end


def call_variants_using_haplotypercaller
  #gdSession     = GoogleDrive.login("seminlee@gmail.com", "xgmuriuxalzmrgqe")
  #foundPatient  = {}
  #readySamples  = {}
  #totalCount    = 0
  #$cancer_types.each do |cancerType|
    #gdSsheet = gdSession.spreadsheet_by_key($normalWskey)
    #gdWsheet = gdSsheet.worksheet_by_title(cancerType)
    #readySamples[cancerType] = []
    #if gdWsheet.respond_to?("rows")
      #gdWsheet.rows.each do |row|
        #next unless row[0].start_with?("TCGA")
        #nsid, nrid, segSize, cnv, snp = row
        #if snp.to_i == 1
          #pid = nsid.split('-')[0..2].join('-')
          #if foundPatient.has_key?(pid)
            #$logger.error "#{pid} has beend already found in #{cancerType}!"
            #exit
          #else
            #foundPatient[pid] = true
          #end
          #shortNsid = nsid.split('-')[0..3].join('-')
          #key = [shortNsid, nrid].join('___')
          #readySamples[cancerType] << key
        #end
      #end
    #end
    #totalCount += readySamples[cancerType].size
    #$logger.info "[#{cancerType}] #{readySamples[cancerType].size} analysis-ready bams found"
  #end
  #$logger.info "[PanCan#{$cancer_types.size}] #{totalCount} analysis-ready bams found in total"

  ## Create a list of bam files to be fed into GATK callers
  panBamListFile = $panVcfDir + "BamFileListForGATK-PanCan.list"
  #panBamListFile.open('w') do |pfile|
    #readySamples.keys.each do |cancerType|
      #canNomalSamples = readySamples[cancerType]
      #canNomalSamples.each do |canNormalSample|
        #sid, rid    = canNormalSample.split("___")
        #nrid        = rid.gsub(/\_L\d+$/, "")
        #globPattern = nil
        #bams        = if cancerType == "CRC"
                        ##globPattern = "/files/CBMI/parklab/semin/BiO/Research/GCC/CGM/CRC/bam/bwa/b37/*/#{sid}*#{rid}*.bam"
                        #globPattern = "/hms/scratch1/sl279/BiO/Research/GCC/Germline/bam/CRC/#{sid}*#{rid}*.bam"
                        #Pathname.glob(globPattern)
                      #else
                        ##globPattern = "/groups/kucherlapati/GCC/LevelII/#{cancerType}/*/#{sid}*#{rid}*.bam"
                        #globPattern = "/groups/kucherlapati/GCC/LevelII/#{cancerType}/#{nrid}/#{sid}*.bam"
                        #Pathname.glob(globPattern)
                      #end
        #bam         = if bams.size == 0
                        #$logger.warn "[#{cancerType}] Cannot find #{canNormalSample}: #{globPattern} !!!"
                        #exit
                      #elsif bams.size > 1
                        #if (sid == "TCGA-CG-4466-10A" || sid == "TCGA-CG-4477-10A" || sid == "TCGA-A5-A0GE-10A" || sid == "TCGA-A5-A0GV-10A")
                          #bams[0]
                        #else
                          #$logger.warn "[#{cancerType}] Found #{bams.size} bams: #{globPattern} !!!\n#{bams.map { |b| b.to_s }.join("\n")}"
                          #exit
                        #end
                      #else
                        #bams[0]
                      #end
        #pfile.puts bam
      #end
    #end
  #end

  # chunk genomic regions
  chunkSize = 500000
  chromLen  = {}
  $chromLenb37.each_line do |line|
    chrom, len = line.chomp.split(/\s+/)
    len = len.to_i
    next if len == 0
    chrom = "MT" if chrom == "M"
    chromLen[chrom] = len
  end

  $chrs.each do |chr|
    outVcfDir = $panVcfDir + "HaplotypeCaller" + "chrs" + chr.to_s; outVcfDir.mkpath
    slices = (1..chromLen[chr.to_s]).each_slice(chunkSize)
    slices.each_with_index do |slice, si|
      roi     = "#{chr}:#{slice[0]}-#{slice[-1]}"
      roi_num  = "#{(si+1).to_s.rjust(slices.size.to_s.size, '0')}_of_#{slices.size}"
      roi_vcf  = outVcfDir + "PanCan12-WGS-Harvard_GCC.#{roi_num}.vcf"
      next if roi_vcf.exist?
      lsf_out = roi_vcf.sub_ext(".vcf.lsfout")
            #-q park_1d \\
      cmd = <<-CMD
          bsub \\
            -g /gcc/cgm/pan \\
            -q short \\
            -W 12:0 \\
            -M 15000000 \\
            -R "rusage[mem=15000]" \\
            -o #{lsf_out} \\
            #{$java7_bin} -Xmx10G -jar #{$gatk2_bin} \\
              -T HaplotypeCaller \\
              -I #{panBamListFile} \\
              -o #{roi_vcf} \\
              -R #{$refseq_broad} \\
              --dbsnp #{$dbsnp} \\
              -stand_call_conf 10 \\
              -stand_emit_conf 4 \\
              -l INFO \\
              -rf BadCigar \\
              --allow_potentially_misencoded_quality_scores \\
              -L #{roi}
      CMD
      submit cmd
    end
  end
end


def call_variants_using_unifiedgenotyper
  #gdSession     = GoogleDrive.login("seminlee@gmail.com", "xgmuriuxalzmrgqe")
  #foundPatient  = {}
  #readySamples  = {}
  #totalCount    = 0
  #$cancer_types.each do |cancerType|
    #gdSsheet = gdSession.spreadsheet_by_key($normalWskey)
    #gdWsheet = gdSsheet.worksheet_by_title(cancerType)
    #readySamples[cancerType] = []
    #if gdWsheet.respond_to?("rows")
      #gdWsheet.rows.each do |row|
        #next unless row[0].start_with?("TCGA")
        #nsid, nrid, segSize, cnv, snp = row
        #if snp.to_i == 1
          #pid = nsid.split('-')[0..2].join('-')
          #if foundPatient.has_key?(pid)
            #$logger.error "#{pid} has beend already found in #{cancerType}!"
            #exit
          #else
            #foundPatient[pid] = true
          #end
          #shortNsid = nsid.split('-')[0..3].join('-')
          #key = [shortNsid, nrid].join('___')
          #readySamples[cancerType] << key
        #end
      #end
    #end
    #totalCount += readySamples[cancerType].size
    #$logger.info "[#{cancerType}] #{readySamples[cancerType].size} analysis-ready bams found"
  #end
  #$logger.info "[PanCan#{$cancer_types.size}] #{totalCount} analysis-ready bams found in total"

  ## Create a list of bam files to be fed into GATK callers
  panBamListFile = $panVcfDir + "BamFileListForGATK-PanCan.list"
  #panBamListFile.open('w') do |pfile|
    #readySamples.keys.each do |cancerType|
      #canNomalSamples = readySamples[cancerType]
      #canNomalSamples.each do |canNormalSample|
        #sid, rid    = canNormalSample.split("___")
        #nrid        = rid.gsub(/\_L\d+$/, "")
        #globPattern = nil
        #bams        = if cancerType == "CRC"
                        ##globPattern = "/files/CBMI/parklab/semin/BiO/Research/GCC/CGM/CRC/bam/bwa/b37/*/#{sid}*#{rid}*.bam"
                        ##globPattern = "/hms/scratch1/sl279/BiO/Research/GCC/Germline/bam/CRC/#{sid}*#{rid}*.bam"
                        #globPattern = "/groups/kucherlapati/GCC/LevelII/#{cancerType}/#{sid}*#{nrid}*.bam"
                        #Pathname.glob(globPattern)
                      #else
                        ##globPattern = "/groups/kucherlapati/GCC/LevelII/#{cancerType}/*/#{sid}*#{rid}*.bam"
                        #globPattern = "/groups/kucherlapati/GCC/LevelII/#{cancerType}/#{nrid}/#{sid}*.bam"
                        #Pathname.glob(globPattern)
                      #end
        #bam         = if bams.size == 0
                        #$logger.warn "[#{cancerType}] Cannot find #{canNormalSample}: #{globPattern} !!!"
                        #exit
                      #elsif bams.size > 1
                        #if (sid == "TCGA-CG-4466-10A" || sid == "TCGA-CG-4477-10A" || sid == "TCGA-A5-A0GE-10A" || sid == "TCGA-A5-A0GV-10A")
                          #bams[0]
                        #else
                          #$logger.warn "[#{cancerType}] Found #{bams.size} bams: #{globPattern} !!!\n#{bams.map { |b| b.to_s }.join("\n")}"
                          #exit
                        #end
                      #else
                        #bams[0]
                      #end
        #pfile.puts bam
      #end
    #end
  #end

  #exit

  # chunk genomic regions
  chunkSize = 100000
  chromLen  = {}
  $chromLenb37.each_line do |line|
    chrom, len = line.chomp.split(/\s+/)
    len = len.to_i
    next if len == 0
    chrom = "MT" if chrom == "M"
    chromLen[chrom] = len
  end

  $chrs.each do |chr|
    #next if chr != "MT"
    outVcfDir   = $panVcfDir + "UnifiedGenotyper" + "chrs" + chr.to_s; outVcfDir.mkpath
    slices      = (1..chromLen[chr.to_s]).each_slice(chunkSize)
    num_slices  = slices.to_a.size
    slices.each_with_index do |slice, si|
      roi = if chr == "MT"
              "#{chr}:1-16569"
            else
              "#{chr}:#{slice[0]}-#{slice[-1]}"
            end
      roi_num    = "#{(si+1).to_s.rjust(num_slices.to_s.size, '0')}_of_#{num_slices}"
      roi_vcf    = outVcfDir + "PanCan12-WGS-Harvard_GCC.#{chr}.#{roi_num}.vcf.gz"
      roi_vcfIdx = roi_vcf.sub_ext(".idx")
      lsf_out   = roi_vcf.sub_ext(".lsfout")
      if lsf_out.exist?
        lsf_out_txt = lsf_out.read
        if lsf_out_txt.match(/^Exit/)
          warn "GATK for #{roi_vcf} failed, Resubmitting a job..."
          lsf_out.delete
          roi_vcf.delete if roi_vcf.exist?
          roi_vcfIdx.delete if roi_vcfIdx.exist?
        elsif lsf_out_txt.match(/^Successfully completed/)
          #warn "GATK for #{roi_vcf} done. Skip."
          next
        else
          warn "Unknown GATK status in #{lsf_out}. Skip."
          next
          #warn "Unknown GATK status in #{lsf_out}! Resubmitting a job...."
          #roi_vcf.delete if roi_vcf.exist?
          #roi_vcfIdx.delete if roi_vcfIdx.exist?
        end
      else
        #warn("Cannot find #{lsf_out}! Skip.")
        #next
        warn("Cannot find #{lsf_out}! Resubmitting a job....")
        roi_vcf.delete if roi_vcf.exist?
        roi_vcfIdx.delete if roi_vcfIdx.exist?
      end
            #-q short -W 12:0 \\
      cmd = <<-CMD
          bsub \\
            -g /cgm/pan/#{chr} \\
            -q park_7d \\
            -M 6000000 \\
            -R "rusage[mem=6000]" \\
            -o #{lsf_out} \\
            #{$java7_bin} -Xmx5G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk2_bin} \\
              -T UnifiedGenotyper \\
              -glm BOTH \\
              -I #{panBamListFile} \\
              -o #{roi_vcf} \\
              -R #{$refseq_broad} \\
              --dbsnp #{$dbsnp} \\
              -stand_call_conf 10 \\
              -stand_emit_conf 4 \\
              -l INFO \\
              -rf BadCigar \\
              --allow_potentially_misencoded_quality_scores \\
              -L #{roi}
      CMD
      submit cmd
    end
  end
end

def combine_sliced_vcfs
  $chrs.each do |chr|
    chrVcfDir = $panVcfDir + "UnifiedGenotyper" + "chrs" + "#{chr}"
    outVcf = $panVcfDir + "UnifiedGenotyper" + "chrs" + "PanCan12-WGS-Harvard_GCC.#{chr}.vcf.gz"
    vcfs = Pathname.glob(chrVcfDir + "*.vcf.gz").sort
    cmd = <<-CMD
    bsub \\
      -q i2b2_12h \\
      "vcf-concat #{chrVcfDir + "*.vcf.gz"} | bgzip -c > #{outVcf}"
    CMD
    submit cmd
  end
end

def combine_chr_vcfs_vcftools
  vcf = $panVcfDir + "UnifiedGenotyper" + "PanCan12-WGS-Harvard_GCC.vcf.gz"
  combinedVcfDir = $panVcfDir + "UnifiedGenotyper" + "chrs" + "combined"
  cmd = <<-CMD
    bsub \\
      -q park_7d \\
      "vcf-concat #{ $chrs.map { |c| (combinedVcfDir + "PanCan12-WGS-Harvard_GCC.#{c}.vcf.gz").to_s }.join(" ") } | bgzip -c > #{vcf}"
  CMD
  puts cmd
end

def combine_chr_vcfs_gatk
  vcf = $panVcfDir + "UnifiedGenotyper" + "PanCan12-WGS-Harvard_GCC.gatk.vcf.gz"
  combinedVcfDir = $panVcfDir + "UnifiedGenotyper" + "chrs" + "combined"
  cmd = <<-CMD
    bsub \\
      -q park_7d \\
      java -Xmx10g -jar /groups/park/semin/BiO/Install/GenomeAnalysisTK-2.7-2-g6bda569/GenomeAnalysisTK.jar \\
        -R /groups/park/semin/BiO/Install/GATK-bundle/human_g1k_v37_decoy.fasta \\
        -T CombineVariants \\
        #{ $chrs.map { |c| "--variant " + (combinedVcfDir + "PanCan12-WGS-Harvard_GCC.#{c}.vcf.gz").to_s }.join(" ") } \\
        -o #{vcf} \\
        --assumeIdenticalSamples
  CMD
  puts cmd
end

def convert_vcfs_to_bcfs
  vcfs = Pathname.glob($panVcfDir + "UnifiedGenotyper" + "chrs" + "combined" + "*.vcf.gz").sort
  vcfs.each do |vcf|
    bcf = vcf.dirname + vcf.basename(".gz").sub_ext(".bcf")
    cmd = <<-CMD
    bsub \\
      -g /gcc/germ/bcf \\
      -q park_7d \\
      "bcftools view -bS #{vcf} -D #{$ref_dict} > #{bcf}"
    CMD
    puts cmd
  end
end

def recalibrate_bgzipped_snp_vcfs
  #vcfs = Pathname.glob($panVcfDir + "UnifiedGenotyper" + "chrs" + "combined" + "*.vcf.gz").sort
  vcf = $panVcfDir + "UnifiedGenotyper" + "PanCan12-WGS-Harvard_GCC.vcf.gz"
  snpVcf = vcf.dirname + vcf.basename(".gz").sub_ext(".snp.vcf.gz")
  snpVqsrVcf = vcf.dirname + vcf.basename(".gz").sub_ext(".snp.vqsr.vcf.gz")
  snpRecalFile = snpVcf.sub_ext(".recal")
  snpTrancheFile = snpVcf.sub_ext(".traches")
  snpRscriptFile = snpVcf.sub_ext("_plots.R")
  snpLsfOut = snpVcf.sub_ext(".recal.lsfout")
    #-q long -W 700:0 \\
    ##{$java7_bin} -Xmx20G -jar #{$gatk2_bin} \\
        #-T SelectVariants \\
        #-R #{$refseq_broad} \\
        #--variant #{vcf} \\
        #-o #{snpVcf} \\
        #-selectType SNP \\
        #-selectType MNP && \\
    #tabix -p vcf #{snpVcf} &&
  cmd = <<-CMD
  bsub \\
    -g /gcc/germ/vqsr \\
    -q i2b2_unlimited \\
    -R "rusage[mem=25000]" -M 25000000 \\
    -o #{snpLsfOut} \\
    "#{$java7_bin} -Xmx20G -jar #{$gatk2_bin} \\
        -T VariantRecalibrator \\
        -R #{$refseq_broad} \\
        -input #{snpVcf} \\
        -resource:hapmap,known=false,training=true,truth=true,prior=15.0 #{$hapmap} \\
        -resource:omni,known=false,training=true,truth=true,prior=12.0 #{$g1kOmni} \\
        -resource:1000G,known=false,training=true,truth=false,prior=10.0 #{$g1kSnp} \\
        -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 #{$dbsnp} \\
        -an DP -an QD -an FS -an MQRankSum -an ReadPosRankSum \\
        -mode SNP \\
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
        -recalFile #{snpRecalFile} \\
        -tranchesFile #{snpTrancheFile} \\
        -rscriptFile #{snpRscriptFile} && \\
    #{$java7_bin} -jar #{$gatk2_bin} \\
        -T ApplyRecalibration \\
        -R #{$refseq_broad} \\
        -input #{snpVcf} \\
        -mode SNP \\
        --ts_filter_level 99.0 \\
        -recalFile #{snpRecalFile} \\
        -tranchesFile #{snpTrancheFile} \\
        -o #{snpVqsrVcf} && \\
      tabix -p vcf #{snpVqsrVcf}"
  CMD
  submit cmd
end

def recalibrate_bgzipped_indel_vcfs
  vcf = $panVcfDir + "UnifiedGenotyper" + "PanCan12-WGS-Harvard_GCC.vcf.gz"
  indelVcf = vcf.basename(".gz").sub_ext(".indel.vcf.gz")
  indelVqsrVcf = vcf.basename(".gz").sub_ext(".indel.vqsr.vcf.gz")
  indelRecalFile = indelVcf.sub_ext(".recal")
  indelTrancheFile = indelVcf.sub_ext(".traches")
  indelRscriptFile = indelVcf.sub_ext("_plots.R")
  indelLsfOut = indelVcf.sub_ext(".recal.lsfout")
    #-q long -W 700:0 \\
  cmd = <<-CMD
  bsub \\
    -g /gcc/germ/vqsr \\
    -q park_unlimited \\
    -R "rusage[mem=25000]" -M 25000000 \\
    -o #{indelLsfOut} "\\
    #{$java7_bin} -Xmx20G -jar #{$gatk2_bin} \\
        -T SelectVariants \\
        -R #{$refseq_broad} \\
        --variant #{vcf} \\
        -o #{indelVcf} \\
        -selectType INDEL && \\
    tabix -p vcf #{indelVcf} && \\
    #{$java7_bin} -Xmx20G -jar #{$gatk2_bin} \\
        -T VariantRecalibrator \\
        -R #{$refseq_broad} \\
        -input #{indelVcf} \\
        -resource:mills,known=true,training=true,truth=true,prior=12.0 #{$millsIndel} \\
        -an DP \\
        -an FS \\ 
        -an MQRankSum \\
        -an ReadPosRankSum \\
        -mode INDEL \\
        -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
        --maxGaussians 4 \\
        -recalFile #{indelRecalFile} \\
        -tranchesFile #{indelTrancheFile} \\
        -rscriptFile #{indelRscriptFile} && \\
    #{$java7_bin} -jar #{$gatk2_bin} \\
        -T ApplyRecalibration \\
        -R #{$refseq_broad} \\
        -input #{indelVcf} \\
        -mode INDEL \\
        --ts_filter_level 99.0 \\
        -recalFile #{indelRecalFile} \\
        -tranchesFile #{indelTrancheFile} \\
        -o #{indelVqsrVcf} &&
      tabix -p vcf #{indelVqsrVcf}"
  CMD
  puts cmd
end

def extract_chr22_snp_vcfs
  vcfFiles = Pathname.glob($vcf_dir + "CANCER" + "*" + "Harvard_GCC_WGS-*-Normal.snp.vqsr.vcf.gz")
  vcfFiles.each do |vcfFile|
    vcf22File = vcfFile.dirname + vcfFile.basename(".gz").sub_ext(".22.vcf.gz")
    lsfout = vcf22File.sub_ext(".gz.lsfout")
    cmd = <<-CMD
    bsub \\
      -q short -W 12:00 \\ 
      -o #{lsfout} \\
      java -Xmx5G -jar #{$gatk3_bin} \\
        -R /groups/kucherlapati/GCC/LevelII/hg19Ref/Homo_sapiens_assembly19.fasta \\
        -T SelectVariants \\
        --variant #{vcfFile} \\
        -o #{vcf22File} \\
        -L 22
    CMD
    puts cmd
  end
end

def split_1k_chr_vcfs_into_even_sized_chunks
  chunkSize = 10000
  vqsrChrVcfs = Pathname.glob($base_dir + "1000Genomes" + "ALL*.genotypes.vcf").sort
  vqsrChrVcfs.each do |vqsrChrVcf|
    chunkDir = vqsrChrVcf.dirname + "chunks"; chunkDir.mkpath
    chunkStem = chunkDir + vqsrChrVcf.basename.sub_ext(".chunk")
    lsfout = chunkStem.sub_ext(".chunks.lsfout")
    cmd = <<-CMD
      bsub \\
        -q short -W 12:0 \\
        -g /gcc/germ/split \\
        -o #{lsfout} \\
        'grep -vP "^#" #{vqsrChrVcf} | split -l #{chunkSize} -a 5 -d - #{chunkStem}'
    CMD
    submit cmd
  end
end

def add_header_to_even_sized_1k_chr_vcf_chunks
  vqsrChrVcfs = Pathname.glob($base_dir + "1000Genomes" + "ALL*20130502.genotypes.vcf").sort
  vqsrChrVcfs.each do |vqsrChrVcf|
    puts vqsrChrVcf
    vcfHeader = vqsrChrVcf.sub_ext(".vcf.header")
    vcfHeader.open('w') do |file|
      vqsrChrVcf.each_line do |line|
        if line.start_with?("#")
          file.print(line)
        else
          break
        end
      end
    end
    vcfChunks = Pathname.glob($base_dir + "1000Genomes" + "chunks" + "#{vqsrChrVcf.basename('.vcf')}*chunk[0-0][0-9][0-9][0-9][0-9]").sort
    vcfChunks.each do |vcfChunk|
      headeredVcf = Pathname.new("#{vcfChunk}.headered")
      lsfout = "#{headeredVcf}.gz.lsfout"
      cmd = <<-CMD
        bsub \\
          -g /gcc/germ/header \\
          -q short -W 12:0 ' \\
          cat #{vcfHeader} #{vcfChunk} > #{headeredVcf}
          bgzip #{headeredVcf}
          tabix -p vcf #{headeredVcf}.gz'
      CMD
      submit cmd
    end
  end
end

def run_vep_to_overlap_targetscans
  vepVcfFiles = Pathname.glob("*/#{canType}/snp/chrs/*/*.maf05.vcf.gz").sort
  vepVcfFiles.each do |vepVcfFile|
    tsVcfFile = vepVcfFile.dirname + vepVcfFile.basename(".gz").sub_ext(".targetScanS.vcf.gz")
    lsfout = tsVcfFile.sub_ext(".gz.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
      bsub \\
        -q short -W 12:0 \\
        -g /gcc/germ/eur \\
        -R "rusage[mem=10000]" -M 10000000 \\
        -o #{lsfout} \\
          #{$java7_bin} -Xmx10g -jar #{$gatk3_bin} \\
            -T SelectVariants \\
            -R #{$refseq_broad} \\
            --variant #{vcfChrFile} \\
            --sample_file #{eurPopFile} \\
            -o #{vcfChrEurFile}
    CMD
    puts cmd
  end
end


def annotate_targetscan
  vep_bin = Pathname.new "/home/sl279/vep/variant_effect_predictor.pl"
  vep_dir = Pathname.new "/home/sl279/.vep"
  ucsc_db_dir = $base_dir + "UCSC/database"
  vepVcfFiles = Pathname.glob($vcf_dir + "*/snp/chrs/*/*.maf05.vcf.gz").sort
  vepVcfFiles.each do |vepVcfFile|
    vep_vcf = vepVcfFile.dirname + vepVcfFile.basename(".gz").sub_ext(".vep.vcf.gz")
    lsfout = vep_vcf.sub_ext(".gz.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
    bsub \\
      -g /gcc/germ/vep \\
      -q short -W 12:0 \\
      -R "rusage[mem=20000] span[hosts=1]" -M 20000000 \\
      -o #{lsfout} \\
      "perl #{vep_bin} \\
        --force_overwrite \\
        --offline \\
        --no_stats \\
        --vcf \\
        --assembly GRCh37 \\
        --dir #{vep_dir} \\
        --fasta /home/sl279/.vep/homo_sapiens/78_GRCh37/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa \\
        --custom #{ucsc_db_dir}/targetScanS.bed.gz,targetScanS,bed,overlap,0 \\
        --input_file #{vepVcfFile} \\
        --output_file STDOUT | bgzip -c > #{vep_vcf} && tabix -p vcf #{vep_vcf}"
    CMD
    submit cmd
  end
end

def extract_targetscan
  vepVcfFiles = Pathname.glob($vcf_dir + "*/snp/chrs/*/*.maf05.vep.vcf.gz").sort
  vepVcfFiles.each do |vepVcfFile|
    ts_vcf = vepVcfFile.dirname + vepVcfFile.basename(".gz").sub_ext(".targetScanS.vcf.gz")
    lsfout = ts_vcf.sub_ext(".gz.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
    bsub \\
      -g /gcc/germ/vep \\
      -q short -W 12:0 \\
      -o #{lsfout} 'zgrep -P "(^#|targetScanS)" #{vepVcfFile} | bgzip -c > #{ts_vcf} && tabix -p vcf #{ts_vcf}'
    CMD
    system cmd
  end
end

def run_hc_for_each_sample
  cancer_bam_dirs = $bam_dir.children.select { |d| d.directory? }.sort
  cancer_bam_dirs.each do |cancer_bam_dir|
    cancer_type = cancer_bam_dir.basename.to_s
    cancer_vcf_dir = $gvcf_dir + cancer_type;cancer_vcf_dir.mkpath
    cancer_normal_bams = Pathname.glob(cancer_bam_dir + "*.rc.bam")
    Parallel.each(cancer_normal_bams, :in_threads => 10) do |cancer_normal_bam|
    #cancer_normal_bams.each do |cancer_normal_bam|
      out_vcf_file = cancer_vcf_dir + cancer_normal_bam.basename.sub_ext(".bam.g.vcf.gz")
      lsfout = out_vcf_file.sub_ext(".gz.lsfout")
      next if lsfout.exist?
          #-q long -W 700:0 \\
      cmd = <<-CMD
        bsub \\
          -g /germ/hc \\
          -q i2b2_7d \\
          -R "rusage[mem=5500]" -M 5500000 \\
          -o #{lsfout} \\
          #{$java7_bin} -XX:+UseSerialGC -Xmx5G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk3_bin} \\
            -T HaplotypeCaller \\
            -R #{$refseq} \\
            -I #{cancer_normal_bam} \\
            --emitRefConfidence GVCF \\
            --dbsnp #{$dbsnp} \\
            -o #{out_vcf_file}
      CMD
      submit cmd
    end
  end
end

def check_redundant_samples
  cancer_type = "HNSC"
  cancer_gvcf_dir = $gvcf_dir + cancer_type
  cancer_gvcf_files = Pathname.glob(cancer_gvcf_dir + "*.g.vcf.gz").sort
  patient_ids = cancer_gvcf_files.map { |f| f.basename.to_s.match(/(TCGA-\S{2}-\S{4})/)[1] }.uniq
  patient_ids.each do |patient_id|
    patient_gvcf_files = Pathname.glob(cancer_gvcf_dir + "#{patient_id}*.g.vcf.gz").sort
    if patient_gvcf_files.size > 1
      puts "#{patient_gvcf_files.size} files found for #{patient_id}:"
      puts patient_gvcf_files.join("\n")
      #puts "#{patient_gvcf_files[0]} selected."
      #(1...patient_gvcf_files.size).each do |i|
        #patient_gvcf_files[i].rename(patient_gvcf_files[i].to_s.gsub(".gz", ".gz.discarded"))
      #end
    end
  end
end

def combine_gvcfs_for_each_cancer_chr
  cancer_type = "BLCA"
  cancer_gvcf_dir = $gvcf_dir + cancer_type
  cancer_combined_chr_gvcf_dir = cancer_gvcf_dir + "combined/chrs"; cancer_combined_chr_gvcf_dir.mkpath
  cancer_gvcf_files = Pathname.glob(cancer_gvcf_dir + "*.g.vcf.gz").sort
  $chrs.each do |chr|
    cancer_combined_gvcf = cancer_combined_chr_gvcf_dir + "Harvard_GCC-WGS-#{cancer_type}.#{chr}.g.vcf.gz"
    lsfout = cancer_combined_gvcf.sub_ext(".gz.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
      bsub \\
        -g /germ/gvcf \\
        -q i2b2_7d \\
        -R "rusage[mem=5500]" -M 5500000 \\
        -o #{lsfout} \\
        #{$java7_bin} -XX:+UseSerialGC -Xmx5G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk3_bin} \\
          -T CombineGVCFs \\
          -R #{$refseq} \\
          #{ cancer_gvcf_files.map { |v| "--variant #{v}" }.join(' ') } \\
          -L #{chr} \\
          -o #{cancer_combined_gvcf}
    CMD
    submit cmd
  end
end

def combine_gvcfs_for_each_cancer_roi
  #chunkSize = 500000
  chunkSize = 2000000
  chromLen  = {}
  $chromLenb37.each_line do |line|
    chrom, len = line.chomp.split(/\s+/)
    len = len.to_i
    next if len == 0
    chrom = "MT" if chrom == "M"
    chromLen[chrom] = len
  end
  cancer_type = "UCEC"
  cancer_gvcf_dir = $gvcf_dir + cancer_type
  cancer_combined_chr_gvcf_dir = cancer_gvcf_dir + "combined/chrs"; cancer_combined_chr_gvcf_dir.mkpath
  cancer_gvcf_files = Pathname.glob(cancer_gvcf_dir + "*.g.vcf.gz").sort
  $chrs.each do |chr|
    cancer_chr_gvcf_dir = cancer_combined_chr_gvcf_dir + chr.to_s; cancer_chr_gvcf_dir.mkpath
    slices = (1..chromLen[chr.to_s]).each_slice(chunkSize)
    Parallel.each_with_index(slices, :in_threads => 10) do |slice, si|
      roi = "#{chr}:#{slice[0]}-#{slice[-1]}"
      roi_num = "#{(si+1).to_s.rjust(slices.size.to_s.size, '0')}_of_#{slices.size}"
      roi_vcf = cancer_chr_gvcf_dir + "Harvard_GCC-WGS-#{cancer_type}.#{chr}.#{roi_num}.g.vcf.gz"
      lsfout = roi_vcf.sub_ext(".gz.lsfout")
      next if lsfout.exist?
      cmd = <<-CMD
          bsub \\
            -g /germ/cg \\
            -q short -W 12:0 \\
            -R "rusage[mem=5500]" -M 5500000 \\
            -o #{lsfout} \\
              #{$java7_bin} -XX:+UseSerialGC -Xmx5G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk3_bin} \\
              -T CombineGVCFs \\
              -R #{$refseq} \\
              #{ cancer_gvcf_files.map { |v| "--variant #{v}" }.join(' ') } \\
              -L #{roi} \\
              -o #{roi_vcf}
      CMD
      submit cmd
    end
  end
end

def combine_roi_gvcfs_for_each_cancer
  ncore = 4
  cancer_type = "UCEC"
  cancer_gvcf_dir = $gvcf_dir + cancer_type
  cancer_chr_gvcf_dir = cancer_gvcf_dir + "combined/chrs"
  $chrs.each do |chr|
    cancer_ind_chr_gvcf_dir = cancer_chr_gvcf_dir + chr.to_s
    roi_gvcfs = Pathname.glob(cancer_ind_chr_gvcf_dir + "*.g.vcf.gz").sort
    cancer_ind_chr_gvcf = cancer_chr_gvcf_dir + "Harvard_GCC-WGS-#{cancer_type}.#{chr}.g.vcf.gz"
    lsfout = cancer_ind_chr_gvcf.sub_ext(".gz.lsfout")
    next if lsfout.exist?
        #-R "rusage[mem=5500] span[hosts=1]" -M 5500000 \\
    cmd = <<-CMD
      bsub \\
        -g /germ/cv \\
        -q mcore -W 700:0 \\
        -R "rusage[mem=1500] span[hosts=1]" -M 1500000 \\
        -n #{ncore} \\
        -o #{lsfout} \\
        #{$java7_bin} -XX:+UseSerialGC -Xmx5G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk3_bin} \\
          -T CombineVariants \\
          -R #{$refseq} \\
          --assumeIdenticalSamples \\
          -nt #{ncore} \\
          #{ roi_gvcfs.map { |v| "--variant #{v}" }.join(' ') } \\
          -o #{cancer_ind_chr_gvcf}
    CMD
    submit cmd
  end
end

def combine_chr_gvcfs_for_each_cancer
  num_cores = 4
  cancer_type = "UCEC"
  cancer_gvcf_dir = $gvcf_dir + cancer_type
  cancer_combined_gvcf_dir = cancer_gvcf_dir + "combined"
  cancer_combined_chr_gvcf_dir = cancer_combined_gvcf_dir + "chrs"
  cancer_combined_gvcf = cancer_combined_gvcf_dir + "Harvard_GCC-WGS-#{cancer_type}.g.vcf.gz"
  lsfout = cancer_combined_gvcf.sub_ext(".gz.lsfout")
  exit if lsfout.exist?
      #-q mcore -W 700:0 \\
  cmd = <<-CMD
    bsub \\
      -g /germ/cv \\
      -q i2b2_7d \\
      -R "rusage[mem=5000] span[hosts=1]" -M 5000000 \\
      -n #{num_cores} \\
      -o #{lsfout} \\
      #{$java7_bin} -XX:+UseSerialGC -Xmx15G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk3_bin} \\
        -T CombineVariants \\
        -R #{$refseq} \\
        --assumeIdenticalSamples \\
        -nt #{num_cores} \\
        #{ $chrs.map { |c| "--variant #{cancer_combined_chr_gvcf_dir}/Harvard_GCC-WGS-#{cancer_type}.#{c}.g.vcf.gz" }.join(' ') } \\
        -o #{cancer_combined_gvcf}
  CMD
  submit cmd
end

def genotype_gvcfs_for_each_cancer_chr
  num_cores = 4
  cancer_type = "ESCA"
  cancer_gvcf_dir = $gvcf_dir + cancer_type
  cancer_gvcf_files = Pathname.glob(cancer_gvcf_dir + "*.g.vcf.gz").sort
  cancer_chr_vcf_dir = $vcf_dir + "CANCER" + cancer_type + "chrs"; cancer_chr_vcf_dir.mkpath
  $chrs.each do |chr|
    cancer_chr_vcf = cancer_chr_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.#{chr}.vcf.gz"
    lsfout = cancer_chr_vcf.sub_ext(".gz.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
      bsub \\
        -g /germ/vcf \\
        -q mcore -W 700:0 \\
        -R "rusage[mem=5100] span[hosts=1]" -M 5100000 \\
        -n #{num_cores} \\
        -o #{lsfout} \\
        #{$java7_bin} -Xmx20G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk3_bin} \\
          -T GenotypeGVCFs \\
          -R #{$refseq} \\
          --dbsnp #{$dbsnp} \\
          #{ cancer_gvcf_files.map { |v| "--variant #{v}" }.join(' ') } \\
          -L #{chr} \\
          -nt #{num_cores} \\
          -o #{cancer_chr_vcf}
    CMD
    submit cmd
  end
end

def genotype_gvcfs_for_each_cancer_roi
  #chunkSize = 500000
  chunkSize = 2000000
  chromLen  = {}
  $chromLenb37.each_line do |line|
    chrom, len = line.chomp.split(/\s+/)
    len = len.to_i
    next if len == 0
    chrom = "MT" if chrom == "M"
    chromLen[chrom] = len
  end
  cancer_type = "UVM"
  cancer_gvcf_dir = $gvcf_dir + cancer_type
  cancer_gvcf_files = Pathname.glob(cancer_gvcf_dir + "*.g.vcf.gz").sort
  $chrs.each do |chr|
    cancer_chr_vcf_dir = $vcf_dir + "CANCER" + cancer_type + "chrs" + chr.to_s; cancer_chr_vcf_dir.mkpath
    slices = (1..chromLen[chr.to_s]).each_slice(chunkSize)
    Parallel.each_with_index(slices, :in_threads => 10) do |slice, si|
    #slices.each_with_index do |slice, si|
      roi = "#{chr}:#{slice[0]}-#{slice[-1]}"
      roi_num = "#{(si+1).to_s.rjust(slices.size.to_s.size, '0')}_of_#{slices.size}"
      roi_vcf = cancer_chr_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.#{chr}.#{roi_num}.vcf.gz"
      lsfout = roi_vcf.sub_ext(".gz.lsfout")
      next if lsfout.exist?
      cmd = <<-CMD
        bsub \\
          -g /germ/gg \\
          -q short -W 12:0 \\
          -R "rusage[mem=5500]" -M 5500000 \\
          -o #{lsfout} \\
          #{$java7_bin} -XX:+UseSerialGC -Xmx5G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk3_bin} \\
            -T GenotypeGVCFs \\
            -R #{$refseq} \\
            --dbsnp #{$dbsnp} \\
            #{ cancer_gvcf_files.map { |v| "--variant #{v}" }.join(' ') } \\
            -L #{roi} \\
            -o #{roi_vcf}
      CMD
      submit cmd
    end
  end
end

def genotype_gvcfs_for_pancan_roi
  chunkSize = 500000
  chromLen = {}
  $chromLenb37.each_line do |line|
    chrom, len = line.chomp.split(/\s+/)
    len = len.to_i
    next if len == 0
    chrom = "MT" if chrom == "M"
    chromLen[chrom] = len
  end
  cancer_gvcf_files = Pathname.glob($gvcf_dir + "*/combined/*.g.vcf.gz").sort
  $chrs.each do |chr|
    pancan_chr_vcf_dir = $vcf_dir + "PANCAN" + "chrs" + chr.to_s; pancan_chr_vcf_dir.mkpath
    slices = (1..chromLen[chr.to_s]).each_slice(chunkSize)
    Parallel.each_with_index(slices, :in_threads => 10) do |slice, si|
    #slices.each_with_index do |slice, si|
      roi = "#{chr}:#{slice[0]}-#{slice[-1]}"
      roi_num = "#{(si+1).to_s.rjust(slices.size.to_s.size, '0')}_of_#{slices.size}"
      roi_vcf = pancan_chr_vcf_dir + "Harvard_GCC-WGS-PANCAN.#{chr}.#{roi_num}.vcf.gz"
      lsfout = roi_vcf.sub_ext(".gz.lsfout")
      next if lsfout.exist?
      cmd = <<-CMD
        bsub \\
          -g /germ/gg \\
          -q short -W 12:0 \\
          -R "rusage[mem=5500]" -M 5500000 \\
          -o #{lsfout} \\
          #{$java7_bin} -XX:+UseSerialGC -Xmx5G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk3_bin} \\
            -T GenotypeGVCFs \\
            -R #{$refseq} \\
            --dbsnp #{$dbsnp} \\
            #{ cancer_gvcf_files.map { |v| "--variant #{v}" }.join(' ') } \\
            -L #{roi} \\
            -o #{roi_vcf}
      CMD
      submit cmd
    end
  end
end

def combine_roi_vcfs_for_each_cancer
  cancer_type = "UVM"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  cancer_chr_vcf_dir = cancer_vcf_dir + "chrs"
  $chrs.each do |chr|
    cancer_ind_chr_vcf_dir = cancer_chr_vcf_dir + chr.to_s
    roi_vcfs = Pathname.glob(cancer_ind_chr_vcf_dir + "*.vcf.gz").sort
    cancer_ind_chr_vcf = cancer_chr_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.#{chr}.vcf.gz"
    lsfout = cancer_ind_chr_vcf.sub_ext(".gz.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
      bsub \\
        -g /germ/cv \\
        -q short -W 12:0 \\
        -R "rusage[mem=5500]" -M 5500000 \\
        -o #{lsfout} \\
        #{$java7_bin} -XX:+UseSerialGC -Xmx5G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk3_bin} \\
          -T CombineVariants \\
          -R #{$refseq} \\
          --assumeIdenticalSamples \\
          #{ roi_vcfs.map { |v| "--variant #{v}" }.join(' ') } \\
          -o #{cancer_ind_chr_vcf}
    CMD
    submit cmd
  end
end

def combine_roi_vcfs_for_pancan
  num_cores = 4
  pancan_vcf_dir = $vcf_dir + "PANCAN"
  pancan_chr_vcf_dir = pancan_vcf_dir + "chrs"
  $chrs.each do |chr|
    pancan_ind_chr_vcf_dir = pancan_chr_vcf_dir + chr.to_s
    roi_vcfs = Pathname.glob(pancan_ind_chr_vcf_dir + "*.vcf.gz").sort
    pancan_ind_chr_vcf = pancan_chr_vcf_dir + "Harvard_GCC-WGS-PANCAN.#{chr}.vcf.gz"
    lsfout = pancan_ind_chr_vcf.sub_ext(".gz.lsfout")
    next if lsfout.exist?
    cmd = <<-CMD
      bsub \\
        -g /germ/cv \\
        -q i2b2_1d \\
        -R "rusage[mem=5500]" -M 5500000 \\
        -n #{num_cores} \\
        -o #{lsfout} \\
        #{$java7_bin} -XX:+UseSerialGC -Xmx20G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk3_bin} \\
          -T CombineVariants \\
          -R #{$refseq} \\
          -nt #{num_cores} \\
          --assumeIdenticalSamples \\
          #{ roi_vcfs.map { |v| "--variant #{v}" }.join(' ') } \\
          -o #{pancan_ind_chr_vcf}
    CMD
    submit cmd
  end
end

def combine_chr_vcfs_for_each_cancer
  cancer_type = "UVM"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  cancer_chr_vcf_dir = cancer_vcf_dir + "chrs"
  cancer_vcf = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.vcf.gz"
  lsfout = cancer_vcf.sub_ext(".gz.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
    bsub \\
      -g /germ/cv \\
      -q i2b2_12h -W 12:0 \\
      -R "rusage[mem=5500]" -M 5500000 \\
      -o #{lsfout} \\
      #{$java7_bin} -XX:+UseSerialGC -Xmx5G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk3_bin} \\
        -T CombineVariants \\
        -R #{$refseq} \\
        --assumeIdenticalSamples \\
        #{ $chrs.map { |c| "--variant #{cancer_chr_vcf_dir}/Harvard_GCC-WGS-#{cancer_type}.#{c}.vcf.gz" }.join(' ') } \\
        -o #{cancer_vcf}
  CMD
  submit cmd
end

def combine_chr_vcfs_for_pancan
  num_cores = 8
  pancan_vcf_dir = $vcf_dir + "PANCAN"
  pancan_chr_vcf_dir = pancan_vcf_dir + "chrs"
  pancan_vcf = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.vcf.gz"
  lsfout = pancan_vcf.sub_ext(".gz.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
    bsub \\
      -g /germ/cv \\
      -q i2b2_1d \\
      -R "rusage[mem=5500] span[hosts=1]" -M 5500000 \\
      -n #{num_cores} \\
      -o #{lsfout} \\
      #{$java7_bin} -XX:+UseSerialGC -Xmx20G -Djava.io.tmpdir=#{$tmp_dir} -jar #{$gatk3_bin} \\
        -T CombineVariants \\
        -R #{$refseq} \\
        --assumeIdenticalSamples \\
        -nt #{num_cores} \\
        #{ $chrs.map { |c| "--variant #{pancan_chr_vcf_dir}/Harvard_GCC-WGS-PANCAN.#{c}.vcf.gz" }.join(' ') } \\
        -o #{pancan_vcf}
  CMD
  puts cmd
end


def extract_snps_from_cancer_vcf
  cancer_type = "CESC"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  cancer_vcf = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.vcf.gz"
  cancer_snp_vcf = cancer_vcf.dirname + cancer_vcf.basename(".gz").sub_ext(".snp.vcf.gz")
  lsfout = cancer_snp_vcf.sub_ext(".gz.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
    bsub \\
      -g /germ/sv \\
      -q i2b2_12h \\
      -R "rusage[mem=15000]" -M 15000000 \\
      -o #{lsfout} \\
      #{$java7_bin} -XX:+UseSerialGC -Xmx13G -jar #{$gatk3_bin} \\
          -T SelectVariants \\
          -R #{$refseq} \\
          --variant #{cancer_vcf} \\
          -o #{cancer_snp_vcf} \\
          -selectType SNP \\
          -selectType MNP
  CMD
  submit cmd
end

def extract_snps_from_pancan_vcf
  num_cores = 8
  pancan_vcf_dir = $vcf_dir + "PANCAN"
  pancan_vcf = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.vcf.gz"
  pancan_snp_vcf = pancan_vcf.dirname + pancan_vcf.basename(".gz").sub_ext(".snp.vcf.gz")
  lsfout = pancan_snp_vcf.sub_ext(".gz.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
    bsub \\
      -g /germ/sv \\
      -q i2b2_1d \\
      -n #{num_cores} \\
      -R "rusage[mem=5500] span[hosts=1]" -M 5500000 \\
      -o #{lsfout} \\
      #{$java7_bin} -XX:+UseSerialGC -Xmx40G -jar #{$gatk3_bin} \\
        -T SelectVariants \\
        -R #{$refseq} \\
        -nt #{num_cores} \\
        --variant #{pancan_vcf} \\
        -o #{pancan_snp_vcf} \\
        -selectType SNP \\
        -selectType MNP
  CMD
  puts cmd
end

def extract_indels_from_cancer_vcf
  cancer_type = "CESC"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  cancer_vcf = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.vcf.gz"
  cancer_indel_vcf = cancer_vcf.dirname + cancer_vcf.basename(".gz").sub_ext(".indel.vcf.gz")
  lsfout = cancer_indel_vcf.sub_ext(".gz.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
    bsub \\
      -g /germ/sv \\
      -q i2b2_12h \\
      -R "rusage[mem=15000]" -M 15000000 \\
      -o #{lsfout} \\
      #{$java7_bin} -XX:+UseSerialGC -Xmx13G -jar #{$gatk3_bin} \\
          -T SelectVariants \\
          -R #{$refseq} \\
          --variant #{cancer_vcf} \\
          -o #{cancer_indel_vcf} \\
          -selectType INDEL
  CMD
  submit cmd
end

def extract_indels_from_pancan_vcf
  num_cores = 8
  pancan_vcf_dir = $vcf_dir + "PANCAN"
  pancan_vcf = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.vcf.gz"
  pancan_indel_vcf = pancan_vcf.dirname + pancan_vcf.basename(".gz").sub_ext(".indel.vcf.gz")
  lsfout = pancan_indel_vcf.sub_ext(".gz.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
    bsub \\
      -g /germ/sv \\
      -q i2b2_1d \\
      -R "rusage[mem=5500] span[hosts=1]" -M 5500000 \\
      -n #{num_cores} \\
      -o #{lsfout} \\
      #{$java7_bin} -XX:+UseSerialGC -Xmx40G -jar #{$gatk3_bin} \\
          -T SelectVariants \\
          -R #{$refseq} \\
          -nt #{num_cores} \\
          --variant #{pancan_vcf} \\
          -o #{pancan_indel_vcf} \\
          -selectType INDEL
  CMD
  puts cmd
end

def recalibrate_snp_cancer_vcf
  cancer_type = "BLCA"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  cancer_snp_vcf = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.snp.vcf.gz"
  cancer_snp_recal = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.snp.vcf.gz.recal"
  cancer_snp_tranches = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.snp.vcf.gz.tranches"
  cancer_snp_rscript = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.snp.vcf.gz.plots.R"
  lsfout = cancer_snp_recal.sub_ext(".recal.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
  bsub \\
    -g /germ/vqsr \\
    -q i2b2_7d \\
    -R "rusage[mem=15000]" -M 15000000 \\
    -o #{lsfout} \\
    #{$java7_bin} -Xmx10G -jar #{$gatk3_bin} \\
      -T VariantRecalibrator \\
      -R #{$refseq} \\
      -input #{cancer_snp_vcf} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15.0 #{$hapmap} \\
      -resource:omni,known=false,training=true,truth=true,prior=12.0 #{$g1kOmni} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10.0 #{$g1kSnp} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 #{$dbsnp} \\
      -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff \\
      -mode SNP \\
      -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
      -recalFile #{cancer_snp_recal} \\
      -tranchesFile #{cancer_snp_tranches} \\
      -rscriptFile #{cancer_snp_rscript}
  CMD
  puts cmd
end

def recalibrate_snp_pancan_vcf
  num_cores = 12
  pancan_vcf_dir = $vcf_dir + "PANCAN"
  pancan_snp_vcf = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.snp.vcf.gz"
  pancan_snp_recal = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.snp.vcf.gz.recal"
  pancan_snp_tranches = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.snp.vcf.gz.tranches"
  pancan_snp_rscript = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.snp.vcf.gz.plots.R"
  lsfout = pancan_snp_recal.sub_ext(".recal.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
  bsub \\
    -g /germ/vqsr \\
    -q i2b2_7d \\
    -n #{num_cores} \\
    -R "rusage[mem=5500] span[hosts=1]" -M 5500000 \\
    -o #{lsfout} \\
    #{$java7_bin} -XX:+UseSerialGC -Xmx60G -jar #{$gatk3_bin} \\
      -T VariantRecalibrator \\
      -R #{$refseq} \\
      -nt #{num_cores} \\
      -input #{pancan_snp_vcf} \\
      -resource:hapmap,known=false,training=true,truth=true,prior=15.0 #{$hapmap} \\
      -resource:omni,known=false,training=true,truth=true,prior=12.0 #{$g1kOmni} \\
      -resource:1000G,known=false,training=true,truth=false,prior=10.0 #{$g1kSnp} \\
      -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 #{$dbsnp} \\
      -an DP -an QD -an FS -an SOR -an MQ -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff \\
      -mode SNP \\
      -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
      -recalFile #{pancan_snp_recal} \\
      -tranchesFile #{pancan_snp_tranches} \\
      -rscriptFile #{pancan_snp_rscript}
  CMD
  puts cmd
end

def recalibrate_indel_cancer_vcf
  cancer_type = "BLCA"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  cancer_indel_vcf = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.indel.vcf.gz"
  cancer_indel_recal = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.indel.vcf.gz.recal"
  cancer_indel_tranches = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.indel.vcf.gz.tranches"
  cancer_indel_rscript = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.indel.vcf.gz.plots.R"
  lsfout = cancer_indel_recal.sub_ext(".recal.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
  bsub \\
    -g /germ/vqsr \\
    -q i2b2_7d \\
    -R "rusage[mem=15000]" -M 15000000 \\
    -o #{lsfout} \\
    #{$java7_bin} -XX:+UseSerialGC -Xmx10G -jar #{$gatk3_bin} \\
      -T VariantRecalibrator \\
      -R #{$refseq} \\
      -input #{cancer_indel_vcf} \\
      -resource:mills,known=true,training=true,truth=true,prior=12.0 #{$millsIndel} \\
      -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff \\
      -mode INDEL \\
      -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
      --maxGaussians 4 \\
      -recalFile #{cancer_indel_recal} \\
      -tranchesFile #{cancer_indel_tranches} \\
      -rscriptFile #{cancer_indel_rscript}
  CMD
  puts cmd
end

def recalibrate_indel_pancan_vcf
  num_cores = 8
  pancan_vcf_dir = $vcf_dir + "PANCAN"
  pancan_indel_vcf = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.indel.vcf.gz"
  pancan_indel_recal = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.indel.vcf.gz.recal"
  pancan_indel_tranches = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.indel.vcf.gz.tranches"
  pancan_indel_rscript = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.indel.vcf.gz.plots.R"
  lsfout = pancan_indel_recal.sub_ext(".recal.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
  bsub \\
    -g /germ/vqsr \\
    -q i2b2_7d \\
    -n #{num_cores} \\
    -R "rusage[mem=5500] span[hosts=1]" -M 5500000 \\
    -o #{lsfout} \\
    #{$java7_bin} -XX:+UseSerialGC -Xmx40G -jar #{$gatk3_bin} \\
      -T VariantRecalibrator \\
      -R #{$refseq} \\
      -nt #{num_cores} \\
      -input #{pancan_indel_vcf} \\
      -resource:mills,known=true,training=true,truth=true,prior=12.0 #{$millsIndel} \\
      -an QD -an DP -an FS -an SOR -an MQRankSum -an ReadPosRankSum -an InbreedingCoeff \\
      -mode INDEL \\
      -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 \\
      --maxGaussians 4 \\
      -recalFile #{pancan_indel_recal} \\
      -tranchesFile #{pancan_indel_tranches} \\
      -rscriptFile #{pancan_indel_rscript}
  CMD
  puts cmd
end

def apply_recalibration_snp_cancer_vcf
  cancer_type = "BLCA"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  cancer_snp_vcf = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.snp.vcf.gz"
  cancer_snp_recal = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.snp.vcf.gz.recal"
  cancer_snp_tranches = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.snp.vcf.gz.tranches"
  cancer_snp_vqsr_vcf = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.snp.vqsr.vcf.gz"
  lsfout = cancer_snp_vqsr_vcf.sub_ext(".gz.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
  bsub \\
    -g /germ/vqsr \\
    -q i2b2_1d \\
    -R "rusage[mem=15000]" -M 15000000 \\
    -o #{lsfout} \\
    #{$java7_bin} -XX:+UseSerialGC -Xmx10G -jar #{$gatk3_bin} \\
        -T ApplyRecalibration \\
        -R #{$refseq} \\
        -input #{cancer_snp_vcf} \\
        -mode SNP \\
        --ts_filter_level 99.0 \\
        -recalFile #{cancer_snp_recal} \\
        -tranchesFile #{cancer_snp_tranches} \\
        -o #{cancer_snp_vqsr_vcf}
  CMD
  puts cmd
end

def apply_recalibration_snp_pancan_vcf
  num_cores = 12
  pancan_vcf_dir = $vcf_dir + "PANCAN"
  pancan_snp_vcf = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.snp.vcf.gz"
  pancan_snp_recal = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.snp.vcf.gz.recal"
  pancan_snp_tranches = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.snp.vcf.gz.tranches"
  pancan_snp_vqsr_vcf = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.snp.vqsr.vcf.gz"
  lsfout = pancan_snp_vqsr_vcf.sub_ext(".gz.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
  bsub \\
    -g /germ/vqsr \\
    -q i2b2_7d \\
    -R "rusage[mem=5500] span[hosts=1]" -M 5500000 \\
    -o #{lsfout} \\
    #{$java7_bin} -XX:+UseSerialGC -Xmx60G -jar #{$gatk3_bin} \\
      -T ApplyRecalibration \\
      -R #{$refseq} \\
      -nt #{num_cores} \\
      -input #{pancan_snp_vcf} \\
      -mode SNP \\
      --ts_filter_level 99.0 \\
      -recalFile #{pancan_snp_recal} \\
      -tranchesFile #{pancan_snp_tranches} \\
      -o #{pancan_snp_vqsr_vcf}
  CMD
  puts cmd
end

def apply_recalibration_indel_cancer_vcf
  cancer_type = "BLCA"
  cancer_vcf_dir = $vcf_dir + "CANCER" + cancer_type
  cancer_indel_vcf = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.indel.vcf.gz"
  cancer_indel_recal = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.indel.vcf.gz.recal"
  cancer_indel_tranches = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.indel.vcf.gz.tranches"
  cancer_indel_vqsr_vcf = cancer_vcf_dir + "Harvard_GCC-WGS-#{cancer_type}.indel.vqsr.vcf.gz"
  lsfout = cancer_indel_vqsr_vcf.sub_ext(".gz.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
  bsub \\
    -g /germ/vqsr \\
    -q i2b2_1d \\
    -R "rusage[mem=15000]" -M 15000000 \\
    -o #{lsfout} \\
    #{$java7_bin} -XX:+UseSerialGC -Xmx10G -jar #{$gatk3_bin} \\
        -T ApplyRecalibration \\
        -R #{$refseq} \\
        -input #{cancer_indel_vcf} \\
        -mode INDEL \\
        --ts_filter_level 99.0 \\
        -recalFile #{cancer_indel_recal} \\
        -tranchesFile #{cancer_indel_tranches} \\
        -o #{cancer_indel_vqsr_vcf}
  CMD
  puts cmd
end

def apply_recalibration_indel_pancan_vcf
  num_cores = 8
  pancan_vcf_dir = $vcf_dir + "PANCAN"
  pancan_indel_vcf = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.indel.vcf.gz"
  pancan_indel_recal = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.indel.vcf.gz.recal"
  pancan_indel_tranches = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.indel.vcf.gz.tranches"
  pancan_indel_vqsr_vcf = pancan_vcf_dir + "Harvard_GCC-WGS-PANCAN.indel.vqsr.vcf.gz"
  lsfout = pancan_indel_vqsr_vcf.sub_ext(".gz.lsfout")
  exit if lsfout.exist?
  cmd = <<-CMD
  bsub \\
    -g /germ/vqsr \\
    -q i2b2_1d \\
    -R "rusage[mem=5500] span[hosts=1]" -M 5500000 \\
    -n #{num_cores} \\
    -o #{lsfout} \\
    #{$java7_bin} -XX:+UseSerialGC -Xmx40G -jar #{$gatk3_bin} \\
      -T ApplyRecalibration \\
      -R #{$refseq} \\
      -nt #{num_cores} \\
      -input #{pancan_indel_vcf} \\
      -mode INDEL \\
      --ts_filter_level 99.0 \\
      -recalFile #{pancan_indel_recal} \\
      -tranchesFile #{pancan_indel_tranches} \\
      -o #{pancan_indel_vqsr_vcf}
  CMD
  puts cmd
end


if __FILE__ == $0
  ## Call SNPs/indels from GCC data

  ## UnifiedGenotyper pipeline
  #call_snps_using_queue_using_old_bams_normal
  #recalibrate_vcf_using_queue

  ## HaplotypeCaller pipeline
  #run_hc_for_each_sample
  #check_redundant_samples

  # Individual CANCER
  #
  #combine_gvcfs_for_each_cancer_chr
  #combine_gvcfs_for_each_cancer_roi
  #combine_roi_gvcfs_for_each_cancer
  #combine_chr_gvcfs_for_each_cancer

  #genotype_gvcfs_for_each_cancer_chr
  #genotype_gvcfs_for_each_cancer_roi
  #combine_roi_vcfs_for_each_cancer
  #combine_chr_vcfs_for_each_cancer

  #extract_snps_from_cancer_vcf
  #extract_indels_from_cancer_vcf
  #recalibrate_snp_cancer_vcf
  #recalibrate_indel_cancer_vcf
  #apply_recalibration_snp_cancer_vcf
  #apply_recalibration_indel_cancer_vcf

  #split_cancer_vcf_by_chr
  #filter_cancer_chr_vcfs

  ##split_vqsr_cancer_chr_vcfs_into_even_sized_chunks
  ##add_cancer_vcf_header_to_even_sized_chunks
  ##bgzip_headered_even_sized_cancer_vcf_chunks

  #annotate_cancer_chr_vcfs_with_dbsnp_rsids
  #annotate_cancer_chr_vcfs_with_1000genomes_afs
  #annotate_cancer_chr_vcfs_with_esp6500si_afs
  #annotate_cancer_snp_chr_vcfs_with_cadd_scores
  #annotate_cancer_chr_vcfs_with_vep
  #extract_eur_cancer_chr_vcfs
  #extract_afs_from_eur_cancer_chr_vcfs

  # PANCAN
  #
  #genotype_gvcfs_for_pancan_roi
  #combine_roi_vcfs_for_pancan
  #combine_chr_vcfs_for_pancan
  #extract_snps_from_pancan_vcf
  #extract_indels_from_pancan_vcf
  #recalibrate_snp_pancan_vcf
  #recalibrate_indel_pancan_vcf
  #apply_recalibration_snp_pancan_vcf
  #apply_recalibration_indel_pancan_vcf

  #split_pancan_vcf_by_chr
  #filter_pancan_chr_vcfs

  #annotate_pancan_chr_vcfs_with_dbsnp_rsids
  #annotate_pancan_chr_vcfs_with_1000genomes_afs
  #annotate_pancan_chr_vcfs_with_esp6500si_afs
  annotate_pancan_snp_chr_vcfs_with_cadd_scores
  #annotate_pancan_chr_vcfs_with_vep
  #extract_eur_pancan_chr_vcfs
  #extract_afs_from_eur_pancan_chr_vcfs


  # For Nik
  #
  #annotate_targetscan
  #extract_targetscan

  ## ?

  #call_variants_using_haplotypercaller
  #call_variants_using_unifiedgenotyper
  #combine_sliced_vcfs
  #combine_chr_vcfs_vcftools
  #combine_chr_vcfs_gatk
  #recalibrate_bgzipped_snp_vcfs
  #recalibrate_bgzipped_indel_vcfs

  ## Handle 1000 Genomes VCFs

  #split_1k_chr_vcfs_into_even_sized_chunks
  #add_header_to_even_sized_0k_chr_vcf_chunks
  #remove_raw_even_sized_1k_chr_vcf_chunks
end
