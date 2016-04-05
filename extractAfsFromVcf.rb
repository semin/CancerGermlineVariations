#!/usr/bin/env ruby

puts %w[#CHROM POS REF ALT TCGA_EUR_DEDUP_AF].join("\t")
ARGF.each_line do |line|
  next if line.start_with?("#")
  cols = line.chomp.split("\t")
  chr = cols[0]
  pos = cols[1]
  ref = cols[3]
  alts = cols[4].split(",")
  info = cols[7]
  afs = info.match(/;AF=(\S+?);/)[1].split(",")
  alts.each_with_index do |alt, i|
    puts [chr, pos, ref, alt, afs[i]].join("\t")
  end
end
