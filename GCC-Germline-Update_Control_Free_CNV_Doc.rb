#!/usr/bin/env ruby

require "logger"
require "pathname"
require "fileutils"
require "google_drive_v0"


$logger         = Logger.new(STDOUT)
$logger.level   = Logger::INFO
$bicseq_out_dir = Pathname.new("/groups/kucherlapati/cnv/NBICseq")
#$cancer_types   = %w[BLCA BRCA CRC ESCA HNSC LGG LUAD PRAD SKCM STAD THCA UCEC]
$cancer_types   = %w[CESC]
#$normal_wskey   = "0AqyF8I3ZSraYdFl3dmZReEhEOU1iOXdFVWs2TFFmcGc"
#$normal_wskey   = "0AqyF8I3ZSraYdEN4LTcxdVZtS21PQUtsOUFaeFJDR2c"
$normal_wskey   = "1ypZKmoAeE7OVbrp_Wfxasn1yEVpDPQBQcZD5z_zfOJs"
#$tumor_wskey    = "0AqyF8I3ZSraYdGNMRlU3SWpmd0g0LVptSzBGSGxSbGc"
$tumor_wskey    = "0AqyF8I3ZSraYdEZScGxETWI4Ty1tTUxfNEU4c0J4eGc"


def create_doc_files
  lambdas   = [3]
  bin_size  = 1000
  gdSession = GoogleDriveV0.login("seminlee@gmail.com", "xgmuriuxalzmrgqe")
  #%w[normal tumor].each do |sample_type|
  %w[normal].each do |sample_type|
    $cancer_types.each do |cancer_type|
      $logger.debug "Updating #{cancer_type}"
      gdSsheet      = sample_type == "normal" ? gdSession.spreadsheet_by_key($normal_wskey) : gdSession.spreadsheet_by_key($tumor_wskey)
      gdWsheet      = gdSsheet.worksheet_by_title(cancer_type)
      base_dir      = $bicseq_out_dir + cancer_type
      cnv_dir       = base_dir + "cnv-cf-hg19"
      done_samples  = {}
      if gdWsheet.respond_to?("rows")
        gdWsheet.rows.each do |row|
          next unless row[0].start_with?("TCGA")
          nsid, nrid = row
          short_nsid = nsid.split('-')[0..3].join('-')
          key = [short_nsid, nrid].join('___')
          done_samples[key] = true
        end
      end
      doc_dir   = base_dir + "doc-cf-hg19"; doc_dir.mkpath
      time_tag  = Time.now.strftime("%Y-%m-%d")
      #doc_file1 = doc_dir + "GCC-#{cancer_type}-#{sample_type.upcase}_CFCNVs-#{time_tag}.tsv"
      doc_file1 = doc_dir + "GCC-#{cancer_type}-#{sample_type.upcase}_CFCNVs_L3_B1000-#{time_tag}.tsv"
      #headers1  = %w[Sample_ID Run_Folder #_of_Segments_L3_B1000_NS]
      headers1  = %w[Sample_ID Run_Folder #_of_Segments_L3_B1000]
      cnv_cnt_added = 0
      file1     = doc_file1.open('w')
      file1.puts headers1.join("\t")
      samples = Pathname.glob(cnv_dir + "*" + "*.out").map { |o| o.basename.to_s.split('---')[0] }.uniq
      samples = if (sample_type == "normal")
                  samples.select { |s| s.match(/TCGA\-\S{2}\-\S{4}-(\d{2})\S{1}/)[1].to_i >= 10 }
                else
                  samples.select { |s| s.match(/TCGA\-\S{2}\-\S{4}-(\d{2})\S{1}/)[1].to_i < 10 }
                end
      samples.each do |sample|
        #TCGA-DU-7290-10A-01D-2022___130910_SN208_0495_BD2FGNACXX
        next if sample.match(/POOL/)
        sid, rid  = sample.split("___")
        #rid.gsub!(/_L\d+$/, "") # for new STAD renamed sample runs
        key       = [sid.split('-')[0..3].join('-'), rid].join('___')
        if done_samples.has_key?(key)
          $logger.debug "[#{cancer_type}] Skip #{key}: already included in the google doc"
          next
        end
        #l3_done_cnv_tag_file    = cnv_dir + ".#{sample}---L3_B1000_NS.done"
        #l3_failed_cnv_tag_file  = cnv_dir + ".#{sample}---L3_B1000_NS.failed"
        #l3_running_cnv_tag_file = cnv_dir + ".#{sample}---L3_B1000_NS.running"
        l3_done_cnv_tag_file    = cnv_dir + ".#{sample}---L3_B1000.done"
        l3_failed_cnv_tag_file  = cnv_dir + ".#{sample}---L3_B1000.failed"
        l3_running_cnv_tag_file = cnv_dir + ".#{sample}---L3_B1000.running"
        if l3_failed_cnv_tag_file.exist?
          $logger.warn "[#{cancer_type}] Skip #{sample}: cnv file generation was failed"
          next
        end
        if l3_running_cnv_tag_file.exist?
          $logger.warn "[#{cancer_type}] Skip #{sample}: cnv files are being generated"
          next
        end
        unless l3_done_cnv_tag_file.exist?
          $logger.warn "[#{cancer_type}] Skip #{sample}: cannot verify CNV file generation status"
          next
        end
        cnv_cnt_added += 1
        columns1      = []
        columns1      << sid << rid
        #l3_cnv_out    = cnv_dir + "#{sample}---L3_B1000_NS" + "#{sample}---L3_B1000_NS.NBICseq.out"
        l3_cnv_out    = cnv_dir + "#{sample}---L3_B1000" + "#{sample}---L3_B1000.NBICseq.out"
        l3_cnv_plot   = l3_cnv_out.sub_ext(".l2r.png").basename
        l3_nseg       = l3_cnv_out.readlines.size - 1
        l3_plot_link  = "=HYPERLINK(\"http://compbio.med.harvard.edu/smlee/GCC/NBICseq/#{cancer_type}/#{l3_cnv_plot}\", #{l3_nseg})"
        columns1      << l3_plot_link
        file1.puts columns1.join("\t")
      end
      file1.close
      if cnv_cnt_added == 0
        $logger.debug "[#{cancer_type}] Remove #{doc_file1}: no newly added samples"
        FileUtils.rm doc_file1
      end
      $logger.info "[#{cancer_type} #{sample_type}] Total number of CNVs newly added to the google doc (#{time_tag}): #{cnv_cnt_added}"
      $logger.info "[#{cancer_type} #{sample_type}] #{doc_file1} has been created" if cnv_cnt_added > 0
    end
  end
end


if __FILE__ == $0
  create_doc_files
end
