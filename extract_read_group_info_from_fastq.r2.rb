#!/usr/bin/env ruby

readGroups = {}

#@FCD1J3AACXX:5:1101:1408:1987#/1
#@FCD1J27ACXX:6:1101:1198:2090#CGCGGTGA/1
sleep(600)

while line = gets
  if line.start_with?("@")
    readName = line.gsub("@", "").split("\/").first
    runName, laneNo = readName.split(":")[0..1]
    barcode = readName.split("#")[1]
    readGroup = "#{runName}:#{laneNo}:#{barcode}"
    if readGroups.has_key?(readGroup)
      readGroups[readGroup] += 1
      #puts [runName, laneNo, barcode].join("\t")
    else
      readGroups[readGroup] = 1
    end
  end
end

readGroups.each_pair do |key, value|
  puts [key.split(":"), value].flatten.join("\t")
end

