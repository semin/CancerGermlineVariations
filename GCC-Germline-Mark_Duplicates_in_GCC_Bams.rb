#!/usr/bin/env ruby

require "pathname"
require "fileutils"

def run_mark_duplicates
  outBaseDir = Pathname.new("/home/sl279/BiO/Research/GCC/LevelII")
  lvl2CancerDir = Pathname.new("/groups/kucherlapati/GCC/LevelII/UVM")
  lvl2CancerOutDir = outBaseDir + lvl2CancerDir.basename; lvl2CancerOutDir.mkpath
  lvl2CancerRunDirs = lvl2CancerDir.children.select { |c| c.directory? && c.basename.to_s.start_with?("1") }.sort
  lvl2CancerRunDirs.each do |lvl2CancerRunDir|
    runId = lvl2CancerRunDir.basename.to_s
    lvl2CancerRunOutDir = lvl2CancerOutDir + runId; lvl2CancerRunOutDir.mkpath
    bams = lvl2CancerRunDir.children.select { |f| f.basename.to_s.match(/TCGA/) && f.to_s.end_with?("bam") }.sort
    bams.each do |bam|
      nbam = lvl2CancerRunOutDir + bam.basename.sub_ext(".dedup.bam")
      metrics = lvl2CancerRunOutDir + bam.basename.sub_ext(".dedup.metrics")
      lsfout = lvl2CancerRunOutDir + bam.basename.sub_ext(".dedup.lsfout")
      cmd = <<-CMD
        bsub \\
          -g /gcc/germ/md \\
          -q short -W 12:0 \\
          -R "rusage[mem=20000]" \\
          -M 20000000 \\
          -o #{lsfout} \\
          /opt/java/jdk7/bin/java -Xmx18G -jar /home/sl279/BiO/Install/picard-tools-1.129/picard.jar \\
            MarkDuplicates \\
            INPUT=#{bam} \\
            OUTPUT=#{nbam} \\
            METRICS_FILE=#{metrics} \\
            CREATE_INDEX=true \\
            TMP_DIR=/home/sl279/BiO/Temp \\
            MAX_RECORDS_IN_RAM=250000 \\
            VALIDATION_STRINGENCY=LENIENT
      CMD
      system cmd
    end
  end
end

run_mark_duplicates
