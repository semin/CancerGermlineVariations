#!/usr/bin/env ruby

formatIdx = nil

ARGF.each_line do |line|
  next if line.start_with?('##')
  line.chomp!
  if line.start_with?('#')
    headers = line.gsub("#", "").chomp.split("\t")
    formatIdx = headers.index("FORMAT")
    #headers = headers.map { |h| h.gsub("-", "_").gsub("1465_", "").gsub("_WGSb", "") }
    #headers[(formatIdx+1)..-1] = headers[(formatIdx+1)..-1].map { |h| %w[GT AD1 AD2 DP GQ PL1 PL2 PL3].map { |f| "#{h}_#{f}" } }.flatten
    puts headers.join("\t")
  else
    next if line =~ /\.\/\./
    fields = line.split("\t")
    #0/1:173,141:282:99:255,0,255
    #fields[(formatIdx+1)..-1] = fields[(formatIdx+1)..-1].map { |f| f.start_with?("./.") ? ["./.", nil, nil, nil, nil, nil, nil, nil] : f.split(':').map { |s| s.split(",") } }.flatten
    fields[(formatIdx+1)..-1] = fields[(formatIdx+1)..-1].map { |f| f.split(':')[0] }
    puts fields.join("\t")
  end
end
