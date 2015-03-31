#!/usr/bin/env ruby

require "logger"
require "pathname"
require "fileutils"


$logger         = Logger.new(STDOUT)
$logger.level   = Logger::INFO
$bicseq_out_dir = Pathname.new("/groups/kucherlapati/cnv/NBICseq")
$scale          = false
$cancer_types   = %w[BLCA BRCA CRC ESCA HNSC LGG LUAD PRAD SKCM STAD THCA UCEC]
$swapped_sample_file = $bicseq_out_dir + "STAD" + "out.mat.csv"
$swapped_samples = {}
$swapped_sample_file.readlines.each_with_index do |line, li|
  next if li == 0
  cols = line.strip.gsub("\"", "").split(",")
  bam_in = Pathname.new cols[7]
  bam_out = Pathname.new cols[9]
  sample_in, run_in = $1, $2 if bam_in.basename.to_s =~ /^(\S+?)_(\S+)_s_\d_rg/
  sample_out, run_out = $1, $2 if bam_out.basename.to_s =~ /^(\S+?)_(\S+)_s_\d_rg/
  key_in = "#{sample_in}___#{run_in}"
  key_out = "#{sample_out}___#{run_out}"
  $swapped_samples[key_in] = key_out
  #puts "#{key_in} --> #{key_out}"
end


def submit(cmd, stdout = false)
  #puts cmd; return
  out = `#{cmd}`.chomp
  if out =~ /Job\s{1}<(\d+)>/
    $logger.info "Job #{$1}: #{cmd} submitted." if stdout
    return Integer($1)
  else
    abort out
  end
end


def create_cnv_files
  lambdas         = [3]
  bin_sizes       = [1000]
  read_len        = 50
  max_dist        = 100000
  nbicseq_script  = Pathname.new("/home/sl279/BiO/Install/NBICseq/BICseq.pl")
  plot_script     = Pathname.new("/home/sl279/BiO/Install/NBICseq/plot.l2r.R")
  gene_script     = Pathname.new("/home/sl279/BiO/Install/NBICseq/annotGene.pl")
  $cancer_types.each do |cancer_type|
    build           = "hg19"
    base_dir        = $bicseq_out_dir + cancer_type
    map_dir         = base_dir + "map"
    bin_dir         = if cancer_type == "CRC"
                        base_dir + "bin-hg19"
                      else
                        base_dir + "bin"
                      end
    cnv_dir         = base_dir + "cnv-hg19"; cnv_dir.mkpath
    tmp_dir         = base_dir + "tmp"; tmp_dir.mkpath
    chrs            = (1..22).to_a << "X" << "Y"
    chrchrs         = chrs.map { |c| "chr#{c}" }
    all_bin_dirs    = bin_dir.children.select { |c| c.directory? }
    tmr_bin_dirs    = all_bin_dirs.select { |d| d.basename.to_s.match(/TCGA\-\S{2}\-\S{4}-(\d{2})\S{1}/)[1].to_i < 10 }
    nrm_bin_dirs    = all_bin_dirs.select { |d| d.basename.to_s.match(/TCGA\-\S{2}\-\S{4}-(\d{2})\S{1}/)[1].to_i >= 10 }
    bin_sizes.each do |bin_size|
      cnt_done, cnt_failed, cnt_running, cnt_submitted, cnt_pairs = 0, 0, 0, 0, 0
      nr_pairs, ind_cnt_pairs = {}, {}
      tmr_bin_dirs.each_with_index do |tmr_bin_dir, ti|
        tmr_bname                 = tmr_bin_dir.basename.to_s
        tmr_pid                   = tmr_bname.match(/(TCGA\-\S{2}\-\S{4})/)[1]
        tmr_prm                   = tmr_bin_dir.sub_ext("_B#{bin_size}.prm")
        tmr_bid, tmr_rid          = tmr_bname.split("___")
        tmr_key                   = "#{tmr_bid}___#{tmr_rid}"
        if cancer_type == "STAD" && $swapped_samples.keys.include?(tmr_key)
          #$logger.warn "#{tmr_key} has been renamed! Skip."
          next
        end
        tmr_done_bin_tag_file     = bin_dir + ".#{tmr_bname}.done"
        tmr_failed_bin_tag_file   = bin_dir + ".#{tmr_bname}.failed"
        tmr_running_bin_tag_file  = bin_dir + ".#{tmr_bname}.running"
        if tmr_failed_bin_tag_file.exist?
          $logger.warn "[#{cancer_type}] Skip #{tmr_bin_dir}: bin file generation failed"
          next
        end
        if tmr_running_bin_tag_file.exist?
          $logger.warn "[#{cancer_type}] Skip #{tmr_bin_dir}: bin files are being generated"
          next
        end
        unless tmr_done_bin_tag_file.exist?
          $logger.warn "[#{cancer_type}] Skip #{tmr_bin_dir}: bin file generation is not done"
          next
        end
        nrm_bin_dirs.each_with_index do |nrm_bin_dir, ni|
          nrm_bname                 = nrm_bin_dir.basename.to_s
          nrm_pid                   = nrm_bname.match(/(TCGA\-\S{2}\-\S{4})/)[1]
          nrm_prm                   = nrm_bin_dir.sub_ext("_B#{bin_size}.prm")
          nrm_bid, nrm_rid          = nrm_bname.split("___")
          nrm_key                   = "#{nrm_bid}___#{nrm_rid}"
          if cancer_type == "STAD" && $swapped_samples.keys.include?(nrm_key)
            #$logger.warn "#{nrm_key} has been renamed! Skip."
            next
          end
          #if cancer_type == "PRAD"
            #if (tmr_pid != nrm_pid && nrm_pid != "TCGA-PC-POOL")
              #if nrm_pid != "TCGA-PC-POOL"
                #ind_cnt_pairs[tmr_bname] = 0 unless ind_cnt_pairs.has_key?(tmr_bname)
                #ind_cnt_pairs[nrm_bname] = 0 unless ind_cnt_pairs.has_key?(nrm_bname)
              #end
              #next
            #end
          #end
          if tmr_pid != nrm_pid
            ind_cnt_pairs[tmr_bname] = 0 unless ind_cnt_pairs.has_key?(tmr_bname)
            ind_cnt_pairs[nrm_bname] = 0 unless ind_cnt_pairs.has_key?(nrm_bname)
            next
          end
          nrm_done_bin_tag_file     = bin_dir + ".#{nrm_bname}.done"
          nrm_failed_bin_tag_file   = bin_dir + ".#{nrm_bname}.failed"
          nrm_running_bin_tag_file  = bin_dir + ".#{nrm_bname}.running"
          if nrm_failed_bin_tag_file.exist?
            $logger.warn "[#{cancer_type}] Skip #{nrm_bin_dir}: bin file generation failed"
            next
          end
          if nrm_running_bin_tag_file.exist?
            $logger.warn "[#{cancer_type}] Skip #{nrm_bin_dir}: bin files are being generated"
            next
          end
          unless nrm_done_bin_tag_file.exist?
            $logger.warn "[#{cancer_type}] Skip #{nrm_bin_dir}: bin file generation is not done"
            next
          end
          cnt_pairs         += 1
          nr_pairs[tmr_pid] = 1
          ind_cnt_pairs.has_key?(tmr_bname) ? ind_cnt_pairs[tmr_bname] += 1 : ind_cnt_pairs[tmr_bname] = 1
          ind_cnt_pairs.has_key?(nrm_bname) ? ind_cnt_pairs[nrm_bname] += 1 : ind_cnt_pairs[nrm_bname] = 1
          pair_id = "#{tmr_bid}___#{tmr_rid}---#{nrm_bid}___#{nrm_rid}"
          lambdas.each do |lam|
            fig_title             = "#{pair_id}---L#{lam}_B#{bin_size}#{ $scale == true ? '' : '_NS' }"
            cnv_pair_dir          = cnv_dir + fig_title
            done_cnv_tag_file     = cnv_dir + ".#{fig_title}.done"
            failed_cnv_tag_file   = cnv_dir + ".#{fig_title}.failed"
            running_cnv_tag_file  = cnv_dir + ".#{fig_title}.running"
            if done_cnv_tag_file.exist?
              $logger.debug "[#{cancer_type}] Skip #{cnv_pair_dir}: cnv files were already generated"
              cnt_done += 1
              next
            end
            if failed_cnv_tag_file.exist?
              $logger.warn "[#{cancer_type}] Remove #{cnv_pair_dir}: cnv file generation was failed"
              FileUtils.rm_rf(cnv_pair_dir)
              FileUtils.rm(failed_cnv_tag_file)
              cnt_failed += 1
            end
            if running_cnv_tag_file.exist?
              $logger.warn "[#{cancer_type}] Skip #{cnv_pair_dir}: cnv files are being generated"
              cnt_running += 1
              next
            end
            $logger.info "[#{cancer_type}] Generate CNV segments and plots for #{cnv_pair_dir}"
            FileUtils.touch running_cnv_tag_file
            cnv_pair_dir.mkpath
            cnt_submitted   += 1
            fig_file        = cnv_pair_dir + "#{fig_title}.NBICseq.seg.png"
            nbicseq_config  = cnv_pair_dir + "#{fig_title}.NBICseq.cfg"
            nbicseq_config.open('w') do |file|
              file.puts %w[chrom tumor normal].join("\t")
              chrchrs.each do |chr|
                cols = %W[
                  #{chr}
                  #{tmr_bin_dir.join("#{chr}.b#{bin_size}.bin")}
                  #{nrm_bin_dir.join("#{chr}.b#{bin_size}.bin")}]
                  file.puts cols.join("\t")
              end
            end
            nbicseq_output      = cnv_pair_dir + "#{fig_title}.NBICseq.out"
            genes_output        = nbicseq_output.sub_ext(".genes")
            nbicseq_lsf_output  = nbicseq_output.sub_ext(".lsf")
            nbicseq_cmd = <<-CMD
              bsub \\
                -g /gcc/#{cancer_type.downcase}/seg \\
                -q short -W 12:0 \\
                -R "rusage[mem=2000]" \\
                -M 2000000 \\
                -o #{nbicseq_lsf_output} "
                xvfb-run -a perl #{nbicseq_script} --bootstrap #{ $scale == true ? "--strict" : "--noscale" } --lambda=#{lam} --tmp=#{tmp_dir} --fig=#{fig_file} --title=#{fig_title} #{nbicseq_config} #{nbicseq_output}
                xvfb-run -a Rscript #{plot_script} #{tmr_bin_dir} #{nrm_bin_dir} #{nbicseq_output} #{bin_size} #{lam} 0 #{fig_title}
                perl #{gene_script} -v #{build} -o #{genes_output} #{nbicseq_output}"
            CMD
            nbicseq_cmd_id = submit(nbicseq_cmd)
            check_cmd = <<-CMD
              bsub \\
                -w 'ended(#{nbicseq_cmd_id})' \\
                -g /gcc/#{cancer_type.downcase}/seg \\
                -q mini \\
                -o /dev/null \\
                ruby -e '
                  require "fileutils";
                  require "pathname";
                  d=Pathname.new("#{cnv_pair_dir}");
                  if (d.children.size == 6 && d.children.all? { |c| c.size > 0 });
                    FileUtils.touch("#{done_cnv_tag_file}");
                    FileUtils.rm("#{failed_cnv_tag_file}") if File.exist?("#{failed_cnv_tag_file}");
                  else;
                    FileUtils.touch("#{failed_cnv_tag_file}");
                  end;
                  FileUtils.rm("#{running_cnv_tag_file}") if File.exist?("#{running_cnv_tag_file}")'
            CMD
            submit(check_cmd, false)
          end
        end
      end
      cnt_singltons = 0
      cnt_tumors    = 0
      cnt_normals   = 0
      if ind_cnt_pairs.values.size > 0
        $logger.warn "[#{cancer_type}] Singleton bin sets:"
        ind_cnt_pairs.each_pair do |k, v|
          if v == 0
            cnt_singltons += 1
            cnt_tumors    += 1 if k.match(/TCGA\-\S{2}\-\S{4}-(\d{2})\S{1}/)[1].to_i < 10
            cnt_normals   += 1 if k.match(/TCGA\-\S{2}\-\S{4}-(\d{2})\S{1}/)[1].to_i >= 10
            $logger.warn "[#{cancer_type}] #{k}"
          end
        end
        $logger.info "[#{cancer_type}] No singletons" if cnt_singltons == 0
        $logger.info "[#{cancer_type}] Total number of sample bin sets: #{ind_cnt_pairs.size}"
        $logger.info "[#{cancer_type}] Total number of tumor sample bin sets: #{tmr_bin_dirs.size}"
        $logger.info "[#{cancer_type}] Total number of normal sample bin sets: #{nrm_bin_dirs.size}"
        $logger.info "[#{cancer_type}] Total number of singletons bin sets: #{cnt_singltons}"
        $logger.info "[#{cancer_type}] Total number of tumor singletons bin sets: #{cnt_tumors}"
        $logger.info "[#{cancer_type}] Total number of normal singletons bin sets: #{cnt_normals}"
        $logger.info "[#{cancer_type}] Total number of bin set pairs: #{cnt_pairs}"
        $logger.info "[#{cancer_type}] Total number of non-redundant bin set pairs: #{nr_pairs.keys.size}"
      end
    end
  end
end


if __FILE__ == $0
  create_cnv_files
end
