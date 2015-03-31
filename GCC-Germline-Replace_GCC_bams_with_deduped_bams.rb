#!/usr/bin/env ruby

require "pathname"
require "fileutils"

def replace_with_new_bams
  outBaseDir = Pathname.new("/home/sl279/BiO/Research/GCC/LevelII")
  lvl2CancerDir = Pathname.new("/groups/kucherlapati/GCC/LevelII/UVM")
  lvl2CancerOutDir = outBaseDir + lvl2CancerDir.basename; lvl2CancerOutDir
  lvl2CancerRunDirs = lvl2CancerDir.children.select { |c| c.directory? && c.basename.to_s.start_with?("1") }.sort
  lvl2CancerRunDirs.each do |lvl2CancerRunDir|
    runId = lvl2CancerRunDir.basename.to_s
    lvl2CancerRunOutDir = lvl2CancerOutDir + runId
    bams = lvl2CancerRunDir.children.select { |f| f.basename.to_s.match(/TCGA/) && f.to_s.end_with?("bam") }.sort
    bams.each do |bam|
      bai = bam.sub_ext(".bam.bai")
      nbam = lvl2CancerRunOutDir + bam.basename.sub_ext(".dedup.bam")
      nbai = nbam.sub_ext(".bai")
      lsfout = bam.sub_ext(".bam.lsfout")
      if !nbam.exist?
        die "#{nbam} cannot be found!"
      end
      cmd =<<-CMD
        bsub \\
          -g /gcc/germ/cp \\
          -q short \\
          -W 12:0 \\
          -o #{lsfout} \\
          "cp #{nbam} #{bam}; cp #{nbai} #{bai}"
      CMD
      system cmd
    end
  end
end

replace_with_new_bams
