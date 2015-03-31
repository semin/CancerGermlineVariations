#!/usr/bin/env ruby

require "logger"
require "pathname"
require "fileutils"


$logger         = Logger.new(STDOUT)
$logger.level   = Logger::INFO
$bicseq_out_dir = Pathname.new("/groups/kucherlapati/cnv/NBICseq")
#$scale          = false
$scale          = true
#$cancer_types   = %w[BLCA BRCA CRC ESCA HNSC LGG LUAD PRAD SKCM STAD THCA UCEC]
$cancer_types   = %w[CESC]


def submit(cmd, stdout = false)
  #puts cmd;return
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
  plot_script     = Pathname.new("/home/sl279/BiO/Install/NBICseq/plot.l2r.cf.R")
  gene_script     = Pathname.new("/home/sl279/BiO/Install/NBICseq/annotGeneSingle.pl")
  $cancer_types.each do |cancer_type|
    build           = "hg19"
    base_dir        = $bicseq_out_dir + cancer_type
    map_dir         = base_dir + "map"
    bin_dir         = if cancer_type == "CRC"
                        base_dir + "bin-hg19"
                      else
                        base_dir + "bin"
                      end
    cnv_dir         = base_dir + "cnv-cf-hg19"; cnv_dir.mkpath
    tmp_dir         = base_dir + "tmp"; tmp_dir.mkpath
    chrs            = (1..22).to_a << "X" << "Y"
    chrchrs         = chrs.map { |c| "chr#{c}" }
    bin_sizes.each do |bin_size|
      all_bin_dirs  = bin_dir.children.select { |c| c.directory?  && c.to_s.end_with?("B#{bin_size}") }
      cnt_done, cnt_failed, cnt_running, cnt_submitted, cnt_indivs = 0, 0, 0, 0, 0
      all_bin_dirs.each_with_index do |ind_bin_dir, ni|
        ind_bname                 = ind_bin_dir.basename.to_s
        ind_pid                   = ind_bname.match(/(TCGA\-\S{2}\-\S{4})/)[1]
        ind_prm                   = ind_bin_dir.sub_ext("_B#{bin_size}.prm")
        ind_bid, ind_rid          = ind_bname.split("___")
        ind_done_bin_tag_file     = bin_dir + ".#{ind_bname}.done"
        ind_failed_bin_tag_file   = bin_dir + ".#{ind_bname}.failed"
        ind_running_bin_tag_file  = bin_dir + ".#{ind_bname}.running"
        if ind_failed_bin_tag_file.exist?
          $logger.warn "[#{cancer_type}] Skip #{ind_bin_dir}: bin file generation failed"
          next
        end
        if ind_running_bin_tag_file.exist?
          $logger.warn "[#{cancer_type}] Skip #{ind_bin_dir}: bin files are being generated"
          next
        end
        unless ind_done_bin_tag_file.exist?
          $logger.warn "[#{cancer_type}] Skip #{ind_bin_dir}: bin file generation is not done"
          next
        end
        cnt_indivs  += 1
        sample_id   = "#{ind_bid}___#{ind_rid}"
        lambdas.each do |lam|
          fig_title             = "#{sample_id}---L#{lam}_B#{bin_size}#{ $scale == true ? '' : '_NS' }"
          cnv_sample_dir        = cnv_dir + fig_title
          done_cnv_tag_file     = cnv_dir + ".#{fig_title}.done"
          failed_cnv_tag_file   = cnv_dir + ".#{fig_title}.failed"
          running_cnv_tag_file  = cnv_dir + ".#{fig_title}.running"
          if done_cnv_tag_file.exist?
            $logger.debug "[#{cancer_type}] Skip #{cnv_sample_dir}: cnv files were already generated"
            cnt_done += 1
            next
          end
          if failed_cnv_tag_file.exist?
            $logger.warn "[#{cancer_type}] Remove #{cnv_sample_dir}: cnv file generation was failed"
            FileUtils.rm_rf(cnv_sample_dir)
            FileUtils.rm(failed_cnv_tag_file)
            cnt_failed += 1
          end
          if running_cnv_tag_file.exist?
            $logger.warn "[#{cancer_type}] Skip #{cnv_sample_dir}: cnv files are being generated"
            cnt_running += 1
            next
          end
          $logger.info "[#{cancer_type}] Generate control-free CNV segments and plots for #{cnv_sample_dir}"
          FileUtils.touch running_cnv_tag_file
          cnv_sample_dir.mkpath
          cnt_submitted   += 1
          fig_file        = cnv_sample_dir + "#{fig_title}.NBICseq.seg.png"
          nbicseq_config  = cnv_sample_dir + "#{fig_title}.NBICseq.cfg"
          nbicseq_config.open('w') do |file|
            file.puts %w[chrom bin].join("\t")
            chrchrs.each do |chr|
              cols = %W[
                #{chr}
                #{ind_bin_dir.join("#{chr}.b#{bin_size}.bin")}]
                file.puts cols.join("\t")
            end
          end
          nbicseq_output      = cnv_sample_dir + "#{fig_title}.NBICseq.out"
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
                xvfb-run -a /opt/R-3.1.2/bin/Rscript #{plot_script} #{ind_bin_dir} #{nbicseq_output} #{bin_size} #{lam} 0 #{fig_title}
                perl #{gene_script} -v #{build} -o #{genes_output} #{nbicseq_output}"
          CMD
          #puts nbicseq_cmd;exit
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
                  d=Pathname.new("#{cnv_sample_dir}");
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
      $logger.info "[#{cancer_type}] Total number of bin sets: #{cnt_indivs}"
    end
  end
end


if __FILE__ == $0
  create_cnv_files
end
