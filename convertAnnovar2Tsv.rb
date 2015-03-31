#!/usr/bin/env ruby

formatIdx = nil
prev_fields_size = nil

ARGF.each_line do |line|
  line.chomp!
  if line.start_with?('#')
    headers = line[1..-1].split("\t")
    formatIdx = headers.index("FORMAT")
    headers = headers.map { |h| h.gsub("-", "_").gsub("1465_", "").gsub("_WGSb", "") }
    #headers[(formatIdx+1)..-1] = headers[(formatIdx+1)..-1].map { |h| %w[GT AD1 AD2 DP GQ PL1 PL2 PL3].map { |f| "#{h}_#{f}" } }.flatten
    headers[(formatIdx+1)..-1] = headers[(formatIdx+1)..-1].map { |h| %w[GT AD1 AD2 AD3 DP GQ PL1 PL2 PL3 PL4 PL5 PL6].map { |f| "#{h}_#{f}" } }.flatten
    puts headers.join("\t")
    #warn headers.size
  else
    fields = line.split("\t")
    if prev_fields_size.nil?
      prev_fields_size = fields.size
    elsif prev_fields_size != fields.size
      abort "Cannot recognize this line: #{line}"
    end
    #0/1:173,141:282:99:255,0,255
    fields[(formatIdx+1)..-1] = fields[(formatIdx+1)..-1].map { |f|
      f.start_with?("./.") ? ["./.", nil, nil, nil, nil, nil, nil, nil, nil, nil, nil, nil] : f.split(':').each_with_index.map { |s, si|
        elems = s.split(",")
        if si == 1
          elems[2] = elems[2]
          elems[0..2]
        elsif si == 4
          elems[5] = elems[5]
          elems[0..5]
        else
          elems
        end
      }
    }.flatten
    #warn fields.size
    puts fields.join("\t")
  end
end
