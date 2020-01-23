#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# RM2bed
RM2BEDVER = "0.2.2"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

if ARGV[0].nil?
	puts "\033[1mRM2bed " + RM2BEDVER + "\033[0m"
	puts "\nUsage: ruby RM2bed.rb <in_RM.txt> > <out.bed>"
else
	start = false
	File.open(ARGV[0]) do |f1|
		while line = f1.gets
			if start
				line_arr = line.split(" ")
				puts line_arr[4] + "\t" + (line_arr[5].to_i - 1).to_s + "\t" + line_arr[6]
			elsif line == "\n"
				start = true
			end
		end
	end
end