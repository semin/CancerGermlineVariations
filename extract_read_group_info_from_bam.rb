#!/usr/bin/env ruby

readGroups = {}

while line = gets
  #FCD1LAJACXX:5:1112:8803:86730#CGCGGTGA
  #FCD1J3AACXX:5:1308:3175:24699#
  readName, barcode = line.split("\t").first.split("#")
  flowcellId, laneNo, tileNo, xPos, yPos = readName.split(":")[0..4]
  #puts [flowcellId, laneNo, tileNo, xPos, yPos, barcode].join("\t")
  readGroup = "#{flowcellId}:#{laneNo}:#{barcode}"
  unless readGroups.has_key?(readGroup)
    readGroups[readGroup] = true
    puts [flowcellId, laneNo, barcode].join("\t")
  end
end
