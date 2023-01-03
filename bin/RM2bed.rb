#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# RM2bed
RM2BEDVER = "0.4.3"
# Michael G. Campana and Ellie E. Armstrong, 2020-2023
# Smithsonian Institution and Stanford University

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to RatesTools;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.
 
# We politely request that this work be cited as:
# Armstrong, E.E. & M.G. Campana. 2022. RatesTools: a Nextflow pipeline for detecting
# de novo germline mutations in pedigree sequence data. Bioinformatics. btac784.
# 10.1093/bioinformatics/btac784.
#----------------------------------------------------------------------------------------

# This script converts a RepeatMasker out file to a BED for exclusion of repeat regions using VCFtools.

require_relative 'denovolib'

if ARGV[0].nil?
	# If no input RepeatMasker output, output usage using denovolib: format_splash
	format_splash('RM2bed', RM2BEDVER, '<in_RM.txt[.gz]> > <out.bed>')
else
	# If provided a RepeatMasker out file, convert to a BED file for exclusion
	start = false
	gz_file_open(ARGV[0]) do |f1|
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
