#!/usr/bin/env ruby

require 'pathname'
require 'parallel'

rscript = Pathname.new("/groups/kucherlapati/GCC/Germline/Scripts/GCC-Germline-Normalize_HiC_Contact_Matrices.R")
rawObsFiles = Pathname.glob("/groups/kucherlapati/GCC/Germline/HiC/GSE63525/suppl/*/*intra*/*/*/*.RAWobserved")
Parallel.each(rawObsFiles, :in_threads => 8) do |rawObsFile|
#rawObsFiles.each do |rawObsFile|
  #"chr1_10kb.RAWobserved"
  binPrefix, binPostfix = $1, $2 if rawObsFile.basename.to_s =~ /.*\_(\d+)(\w+)b\.RAWobserved/
  binPostfixNum = case binPostfix
                  when "k"
                    1000
                  when "m"
                    1000000
                  end
  binSize = binPrefix.to_i * binPostfixNum
  normTypes = %w[KR SQRTVC VC]
  normTypes.each do |normType|
    normFile = rawObsFile.sub_ext(".#{normType}norm")
    if (normFile.size == 0)
      warn "#{normFile}'s size is 0!"
      next
    end
    expFile = rawObsFile.sub_ext(".#{normType}expected")
    expFile = "NA" unless expFile.exist?
    outFile = rawObsFile.sub_ext(".#{normType}normalized.RData")
    lsfout = outFile.sub_ext(".RData.lsfout")
    next if lsfout.exist?

  #-q short -W 12:0 \\
    cmd = <<-CMD
bsub \\
  -g /gcc/germ/hic \\
  -q mini \\
  -o #{lsfout} \\
  /opt/R-3.1.2/bin/Rscript #{rscript} #{binSize} #{rawObsFile} #{normFile} #{expFile} #{outFile}
    CMD
    system cmd
  end
end
