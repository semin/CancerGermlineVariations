#!/usr/bin/env ruby

require 'pathname'

$baseDir = Pathname.new("/groups/kucherlapati/GCC/Germline")
$scriptDir = $baseDir + "Scripts"

def filter_vcfs
  #cancerTypes = %w[BLCA CESC CRC HNSC LGG LUAD PRAD SKCM STAD THCA UCEC]
  cancerTypes = %w[LUAD]
  rScript = $scriptDir + "GCC-Germline-Filter_VCFs.R"
  cancerTypes.each do |cancerType|
    lsfout = rScript.sub_ext(".R.#{cancerType}.lsfout")
    cmd =<<-CMD
      bsub \\
        -g /gcc/germ/filter \\
        -q i2b2_12h -W 12:0 \\
        -R "rusage[mem=70000] span[hosts=1]" \\
        -M 70000000 \\
        -o #{lsfout} \\
        xvfb-run -a /opt/R-3.1.2/bin/Rscript #{rScript} #{cancerType}
    CMD
    system cmd
  end
end

filter_vcfs
