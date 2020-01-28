#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# FilterGM
FILTERGMVER = "0.1.0"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

if ARGV[0].nil?
	puts "\033[1mFilterGM " + FILTERGMVER + "\033[0m"
	puts "\nUsage: ruby FilterGM.rb <in_GenMap.bed> <cutoff> > <out.bed>"
else
	File.open(ARGV[0]) do |f1|
		while line = f1.gets
			line_arr = line.split
			puts line if line_arr[4].to_f >= ARGV[1].to_f
		end
	end
end