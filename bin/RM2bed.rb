#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# RM2bed
RM2BEDVER = "0.4.1"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

# This script converts a RepeatMasker out file to a BED for exclusion of repeat regions using VCFtools.

require_relative 'denovolib'

if ARGV[0].nil?
	# If no input RepeatMasker output, output usage using denovolib: format_splash
	format_splash('RM2bed', RM2BEDVER, '<in_RM.txt[.gz]> > <out.bed>')
else
	# If provided a RepeatMasker out file, convert to a BED file for exclusion
	start = false
	gz_file_open(ARGV[0]).open(ARGV[0]) do |f1|
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
