#!/usr/bin/env ruby

require 'pathname'

$baseDir = Pathname.new("/groups/kucherlapati/GCC/Germline")
$scriptDir = $baseDir + "Scripts"

def integer_tcga_eur_vcfs
  cancerTypes = %w[BLCA CESC CRC HNSC LGG LUAD PRAD SKCM STAD THCA UCEC]
  rScript = $scriptDir + "GCC-Germline-Integer_TCGA_EUR_VCFs.R"
  cancerTypes.each do |cancerType|
    lsfout = rScript.sub_ext(".R.#{cancerType}.lsfout")
    cmd =<<-CMD
      bsub \\
        -g /gcc/germ/integer \\
        -q short -W 12:0 \\
        -R "rusage[mem=40000] span[hosts=1]" \\
        -M 40000000 \\
        -o #{lsfout} \\
        xvfb-run -a /opt/R-3.1.2/bin/Rscript #{rScript} #{cancerType}
    CMD
    system cmd
  end
end

integer_tcga_eur_vcfs

