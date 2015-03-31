#!/usr/bin/env ruby

require 'pathname'

$baseDir = Pathname.new("/groups/kucherlapati/GCC/Germline")
$scriptDir = $baseDir + "Scripts"

def postprocess_broad_gdac_gdac_methylation_data
  cancerTypes = %w[BLCA CESC CRC HNSC LGG LUAD PRAD SKCM STAD THCA UCEC]
  rScript = $scriptDir + "GCC-Germline-Postprocess_Broad_GDAC_Methylation450K_Data.R"
  cancerTypes.each do |cancerType|
    lsfout = rScript.sub_ext(".R.#{cancerType}.lsfout")
    cmd =<<-CMD
      bsub \\
        -g /gcc/germ/methyl450 \\
        -q short -W 12:0 \\
        -R "rusage[mem=80000] span[hosts=1]" \\
        -M 80000000 \\
        -o #{lsfout} \\
        xvfb-run -a /opt/R-3.1.2/bin/Rscript #{rScript} #{cancerType}
    CMD
    system cmd
  end
end

postprocess_broad_gdac_gdac_methylation_data

