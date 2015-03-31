#!/usr/bin/env ruby

require "logger"
require "pathname"
require "fileutils"
require "google_drive"

$logger           = Logger.new(STDOUT)
$logger.level     = Logger::INFO
#$normalWsheetKey  = "0AqyF8I3ZSraYdFl3dmZReEhEOU1iOXdFVWs2TFFmcGc"
$normalWsheetKey  = "0AqyF8I3ZSraYdEN4LTcxdVZtS21PQUtsOUFaeFJDR2c"
$tumorWsheetKey   = "0AqyF8I3ZSraYdGNMRlU3SWpmd0g0LVptSzBGSGxSbGc"
$cancerTypes      = %w[BLCA BRCA CRC ESCA HNSC LGG LUAD PRAD SKCM STAD THCA UCEC]

def download_GDAC_data
  date = "20141017"
  dateUnder = "2014_10_17"
  gdacDir = Pathname.new "/groups/kucherlapati/GCC/Germline/GDAC/CNV_Analysis"; gdacDir.mkpath
  $cancerTypes.each do |cancerType|
    cancerType = "COADREAD" if cancerType == "CRC"
    gdacHttpStem = "http://gdac.broadinstitute.org/runs/stddata__#{dateUnder}/data/#{cancerType}/#{date}/"
    mRNAseqLvl3 = "gdac.broadinstitute.org_#{cancerType}.Merge_rnaseq__illuminahiseq_rnaseq__bcgsc_ca__Level_3__gene_expression__data.Level_3.#{date}00.0.0.tar.gz"
    mRNAseqV2Lvl3 = "gdac.broadinstitute.org_#{cancerType}.Merge_rnaseqv2__illuminahiseq_rnaseqv2__unc_edu__Level_3__RSEM_genes_normalized__data.Level_3.#{date}00.0.0.tar.gz"
    clinicalLvl4 = "gdac.broadinstitute.org_#{cancerType}.Clinical_Pick_Tier1.Level_4.#{date}00.0.0.tar.gz"
    [mRNAseqLvl3, mRNAseqV2Lvl3, clinicalLvl4].each do |file|
    #[clinicalLvl4].each do |file|
      source = gdacHttpStem + file
      target = gdacDir + file
      cmd = "wget -O #{target} #{source}"
      system cmd
      FileUtils.rm(target) if target.zero?
      cmd = "tar xzf #{target} -C #{gdacDir}"
      system cmd unless target.zero?
    end
  end
end

def convert_bicseq_to_cbs(input, output, cf = false)
  lines = input.readlines; lines.shift
  nrow  = lines.size

  # check # of segments
  #if (nrow > segcut || nrow < 1)
  if (nrow < 1)
    abort "Skipped #{bic} (# of segments: #{nrow})"
  end

  # check if there is Inf
  txt = input.read
  if (txt =~ /Inf/m)
    abort "Skipped #{bic} (contains Inf value)"
  end

  output.open('w') do |out|
    hdr = %w[Chromosome Start End Tumor_Read Normal_Read Segment_Mean].join("\t")
    out.puts hdr
    lines.each do |line|
      # control free
      #chrom   start   end     binNum  observed        expected        log2.copyRatio  pvalue
      cols = line.chomp.split(/\s+/)
      chrom     = cols[0]
      pos_start = Integer(cols[1]) + 1
      pos_end   = Integer(cols[2]) + 1
      rc_tumor  = cols[4]
      rc_normal = cf ? cols[5] : cols[6]
      log2ratio = cf ? cols[6] : cols[8]
      elms = [chrom, pos_start, pos_end, rc_tumor, rc_normal, log2ratio]
      out.puts elms.join("\t")
    end
  end
end


def collect_somatic_segments
  pwd           = FileUtils.pwd
  nbicseqDir    = Pathname.new "/groups/kucherlapati/cnv/NBICseq"
  pancanDir     = Pathname.new "/groups/kucherlapati/cnv/NBICseq/PanCan"
  pancanSegDir  = pancanDir + "somatic_segments_hg19"; pancanSegDir.mkpath
  gdSession     = GoogleDrive.login("seminlee@gmail.com", "xgmuriuxalzmrgqe")
  gdSsheet      = gdSession.spreadsheet_by_key($normalWsheetKey)
  $cancerTypes.each do |cancerType|
    $logger.info "Fetching QC-passed #{cancerType} somatic CNVS..."
    gdWsheet      = gdSsheet.worksheet_by_title(cancerType)
    cancerOutDir  = pancanSegDir + cancerType; cancerOutDir.mkpath
    cancerSegDir  = cancerOutDir + "seg"; cancerSegDir.mkpath
    cancerCnvDir  = nbicseqDir + cancerType + "cnv-hg19"
    # convert cnv to seg
    if gdWsheet.respond_to?("rows")
      gdWsheet.rows.each do |row|
        next unless row[0].start_with?("TCGA")
        tsid, trid, nsid, nrid, linkToPlot, ready = row
        cnvProfile = Pathname.glob(cancerCnvDir + "#{tsid}___#{trid}---#{nsid}___#{nrid}*-L3_B1000_NS" + "*.out")[0]
        if ready.to_i == 1
          segProfile  = cancerSegDir + cnvProfile.basename.sub_ext(".seg")
          convert_bicseq_to_cbs(cnvProfile, segProfile)
        end
      end
    end
    # make GISTIC input for each cancer type
    FileUtils.cd cancerOutDir
    mergedSegFile = cancerOutDir + "GCC_#{cancerType}_SOMATIC_CNVS-L3_B1000_NS.hg19.seg"
    mergeBin      = "~/BiO-scratch/Install/merge_dnaseq_segment_files.pl"
    mergeCmd      = "perl #{mergeBin} #{cancerSegDir}/*.seg > #{mergedSegFile}"
    system mergeCmd

    gisticInputBin = "~/BiO-scratch/Install/generate_gistic_input_files_from_dnaseq_segment_data.pl"
    gisticInputCmd = "perl #{gisticInputBin} #{mergedSegFile}"
    system gisticInputCmd
    FileUtils.cd pwd
  end
  # make GISTIC input for pan cancer type
  pancanOutDir  = pancanSegDir + "PANCAN"; pancanOutDir.mkpath
  mergedSegFile = pancanOutDir + "GCC_PANCAN_SOMATIC_CNVS-L3_B1000_NS.hg19.seg"
  mergeBin      = "~/BiO-scratch/Install/merge_dnaseq_segment_files.pl"
  mergeCmd      = "perl #{mergeBin} #{pancanSegDir}/*/seg/*.seg > #{mergedSegFile}"
  FileUtils.cd pancanOutDir
  system mergeCmd

  gisticInputBin = "~/BiO-scratch/Install/generate_gistic_input_files_from_dnaseq_segment_data.pl"
  gisticInputCmd = "perl #{gisticInputBin} #{mergedSegFile}"
  system gisticInputCmd
  FileUtils.cd pwd
end


def collect_normal_segments
  pwd           = FileUtils.pwd
  nbicseqDir    = Pathname.new "/groups/kucherlapati/cnv/NBICseq"
  #pancanDir     = Pathname.new "/groups/kucherlapati/cnv/NBICseq/PanCan"
  pancanDir     = Pathname.new "/groups/park/semin/BiO/Research/GCC/Germline"
  pancanSegDir  = pancanDir + "CNV"; pancanSegDir.mkpath
  gdSession     = GoogleDrive.login("seminlee@gmail.com", "xgmuriuxalzmrgqe")
  gdSsheet      = gdSession.spreadsheet_by_key($normalWsheetKey)
  $cancerTypes.each do |cancerType|
    $logger.info "Fetching QC-passed #{cancerType} normal CNVS..."
    gdWsheet      = gdSsheet.worksheet_by_title(cancerType)
    cancerOutDir  = pancanSegDir + cancerType; cancerOutDir.mkpath
    cancerSegDir  = cancerOutDir + "seg"; cancerSegDir.mkpath
    cancerCnvDir  = nbicseqDir + cancerType + "cnv-cf-hg19"
    # convert cnv to seg
    if gdWsheet.respond_to?("rows")
      gdWsheet.rows.each do |row|
        next unless row[0].start_with?("TCGA")
        nsid, nrid, linkToPlot, ready = row
        #cnvProfile = Pathname.glob(cancerCnvDir + "#{nsid}___#{nrid}*-L3_B1000_NS" + "*.out")[0]
        cnvProfile = Pathname.glob(cancerCnvDir + "#{nsid}___#{nrid}*-L3_B1000" + "*.out")[0]
        if ready.to_i == 1
          segProfile  = cancerSegDir + cnvProfile.basename.sub_ext(".seg")
          convert_bicseq_to_cbs(cnvProfile, segProfile, true)
        end
      end
    end
    # make GISTIC input for each cancer type
    FileUtils.cd cancerOutDir
    #mergedSegFile = cancerOutDir + "GCC-#{cancerType}-NORMAL_CNVS-L3_B1000_NS.hg19.seg"
    mergedSegFile = cancerOutDir + "GCC-#{cancerType}-NORMAL_CNVS-L3_B1000.hg19.seg"
    mergeBin      = "~/BiO/Install/merge_dnaseq_segment_files.pl"
    mergeCmd      = "perl #{mergeBin} #{cancerSegDir}/*.seg > #{mergedSegFile}"
    system mergeCmd

    gisticInputBin = "~/BiO/Install/generate_gistic_input_files_from_dnaseq_segment_data.pl"
    gisticInputCmd = "perl #{gisticInputBin} #{mergedSegFile}"
    system gisticInputCmd
    FileUtils.cd pwd
  end
  # make GISTIC input for pan cancer type
  pancanOutDir  = pancanSegDir + "PANCAN"; pancanOutDir.mkpath
  #mergedSegFile = pancanOutDir + "GCC-PANCAN-NORMAL_CNVS-L3_B1000_NS.hg19.seg"
  mergedSegFile = pancanOutDir + "GCC-PANCAN-NORMAL_CNVS-L3_B1000.hg19.seg"
  mergeBin      = "~/BiO/Install/merge_dnaseq_segment_files.pl"
  mergeCmd      = "perl #{mergeBin} #{pancanSegDir}/*/seg/*.seg > #{mergedSegFile}"
  FileUtils.cd pancanOutDir
  system mergeCmd

  gisticInputBin = "~/BiO/Install/generate_gistic_input_files_from_dnaseq_segment_data.pl"
  gisticInputCmd = "perl #{gisticInputBin} #{mergedSegFile}"
  system gisticInputCmd
  FileUtils.cd pwd
end


def collect_tumor_segments
  pwd           = FileUtils.pwd
  nbicseqDir    = Pathname.new "/groups/kucherlapati/cnv/NBICseq"
  #pancanDir     = Pathname.new "/groups/kucherlapati/cnv/NBICseq/PanCan"
  pancanDir     = Pathname.new "/hms/scratch1/sl279/BiO/Research/GCC/Germline"
  pancanSegDir  = pancanDir + "CNV" + "Tumor"; pancanSegDir.mkpath
  gdSession     = GoogleDrive.login("seminlee@gmail.com", "xgmuriuxalzmrgqe")
  gdSsheet      = gdSession.spreadsheet_by_key($tumorWsheetKey)
  $cancerTypes.each do |cancerType|
    $logger.info "Fetching QC-passed #{cancerType} tumor CNVS..."
    gdWsheet      = gdSsheet.worksheet_by_title(cancerType)
    cancerOutDir  = pancanSegDir + cancerType; cancerOutDir.mkpath
    cancerSegDir  = cancerOutDir + "seg"; cancerSegDir.mkpath
    cancerCnvDir  = nbicseqDir + cancerType + "cnv-cf-hg19"
    # convert cnv to seg
    if gdWsheet.respond_to?("rows")
      gdWsheet.rows.each do |row|
        next unless row[0].start_with?("TCGA")
        nsid, nrid, linkToPlot, ready = row
        cnvProfile = Pathname.glob(cancerCnvDir + "#{nsid}___#{nrid}*-L3_B1000_NS" + "*.out")[0]
        if ready.to_i == 1
          segProfile  = cancerSegDir + cnvProfile.basename.sub_ext(".seg")
          convert_bicseq_to_cbs(cnvProfile, segProfile, true)
        end
      end
    end
    # make GISTIC input for each cancer type
    FileUtils.cd cancerOutDir
    mergedSegFile = cancerOutDir + "GCC-#{cancerType}-TUMOR_CNVS-L3_B1000_NS.hg19.seg"
    mergeBin      = "~/BiO-scratch/Install/merge_dnaseq_segment_files.pl"
    mergeCmd      = "perl #{mergeBin} #{cancerSegDir}/*.seg > #{mergedSegFile}"
    system mergeCmd

    gisticInputBin = "~/BiO-scratch/Install/generate_gistic_input_files_from_dnaseq_segment_data.pl"
    gisticInputCmd = "perl #{gisticInputBin} #{mergedSegFile}"
    system gisticInputCmd
    FileUtils.cd pwd
  end
  # make GISTIC input for pan cancer type
  pancanOutDir  = pancanSegDir + "PANCAN"; pancanOutDir.mkpath
  mergedSegFile = pancanOutDir + "GCC-PANCAN-TUMOR_CNVS-L3_B1000_NS.hg19.seg"
  mergeBin      = "~/BiO-scratch/Install/merge_dnaseq_segment_files.pl"
  mergeCmd      = "perl #{mergeBin} #{pancanSegDir}/*/seg/*.seg > #{mergedSegFile}"
  FileUtils.cd pancanOutDir
  system mergeCmd

  gisticInputBin = "~/BiO-scratch/Install/generate_gistic_input_files_from_dnaseq_segment_data.pl"
  gisticInputCmd = "perl #{gisticInputBin} #{mergedSegFile}"
  system gisticInputCmd
  FileUtils.cd pwd
end

#collect_somatic_segments
#collect_normal_segments
#collect_tumor_segments
#download_GDAC_data
