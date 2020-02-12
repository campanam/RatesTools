#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# filterGM
FILTERGMVER = "0.3.0"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

require_relative 'denovolib'

if ARGV[0].nil?
	format_splash('filterGM', FILTERGMVER, '<in_GenMap.bed[.gz]> <cutoff> [exclude] > <out.bed>')
else
	gz_file_open(ARGV[0]).open(ARGV[0]) do |f1|
		while line = f1.gets
			line_arr = line.split
			if ARGV[2] == "exclude"
				puts line if line_arr[4].to_f < ARGV[1].to_f
			else
				puts line if line_arr[4].to_f >= ARGV[1].to_f
			end
		end
	end
end