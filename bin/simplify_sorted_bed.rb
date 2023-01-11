#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# simplify_sorted_bed
SIMSORTBEDVER = "0.1.2"
# Michael G. Campana and Ellie E. Armstrong, 2020-2023
# Smithsonian Institution and Stanford University

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to RatesTools;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.
 
# We politely request that this work be cited as:
# Armstrong, E.E. & M.G. Campana. 2023. RatesTools: a Nextflow pipeline for detecting
# de novo germline mutations in pedigree sequence data. Bioinformatics. 39: btac784.
# 10.1093/bioinformatics/btac784.
#----------------------------------------------------------------------------------------

# Script to detect and combine overlapping entries in a coordinate-sorted BED file
# Fed to simplify_bed for more complex overlaps

require_relative 'denovolib'

if ARGV[0].nil?
	# If no parameters passed, print help screen
	format_splash('simplify_sorted_bed', SIMSORTBEDVER, '<in.bed[.gz]> > <out.bed>')
else
	# Detect and combine overlapping BED entries
	@current_contig = nil # Current contig for co-ordinates
	@startsite = nil # Set start site
	@endsite = nil  # Set end site
	gz_file_open(ARGV[0]) do |f1|
		while line = f1.gets
			line_arr = line.split
			if @current_contig != line_arr[0] # Write output if contig changes
				unless @current_contig.nil?
					puts @current_contig + "\t" + @startsite.to_s + "\t" + @endsite.to_s # Puts last entry
				end
				@current_contig = line_arr[0]
				@startsite = line_arr[1].to_i
				@endsite = line_arr[2].to_i
			elsif @endsite != line_arr[1].to_i # Only extend region if EXACTLY one bp greater (contiguous). simplify_bed will handle more complex overlaps
				puts @current_contig + "\t" + @startsite.to_s + "\t" + @endsite.to_s # Puts last entry
				@startsite = line_arr[1].to_i
				@endsite = line_arr[2].to_i
			else
				@endsite = line_arr[2].to_i # If everything ok, extend region
			end
		end
		puts @current_contig + "\t" + @startsite.to_s + "\t" + @endsite.to_s # Puts last entry
	end
end
