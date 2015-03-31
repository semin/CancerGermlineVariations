#!/usr/bin/env ruby

readGroups = {}

#@FCD1J3AACXX:5:1101:1408:1987#/1
#@FCD1J27ACXX:6:1101:1198:2090#CGCGGTGA/1
while line = gets
  if line.start_with?("@")
    readName = line.gsub("@", "").split("\/").first
    runName, laneNo = readName.split(":")[0..1]
    barcode = readName.split("#")[1]
    readGroup = "#{runName}:#{laneNo}:#{barcode}"
    unless readGroups.has_key?(readGroup)
      readGroups[readGroup] = true
      puts [runName, laneNo, barcode].join("\t")
    end
  end
end

