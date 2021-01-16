#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# simplify_bed
SIMBEDVER = "0.2.1"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

# Script to sort, detect and combine overlapping entries in an unsorted BED file

require_relative 'denovolib'

if ARGV[0].nil?
	# If no parameters passed, print help screen
	format_splash('simplify_bed', SIMBEDVER, '<in.bed[.gz]> > <out.bed>')
else
	# Detect and combine overlapping BED entries
	$stderr.puts 'Reading BED'
	$contig_hash = {} # Hash of contigs and annotated BED regions
	gz_file_open(ARGV[0]) do |f1|
		while line = f1.gets
			line_arr = line.split
			if $contig_hash.keys.include?(line_arr[0])
				$contig_hash[line_arr[0]].push([line_arr[1].to_i, line_arr[2].to_i - 1])
			else
				$contig_hash[line_arr[0]] = [[line_arr[1].to_i, line_arr[2].to_i - 1]]
			end
		end
	end
	$stderr.puts 'Simplifying BED'
	for hash in $contig_hash.keys
		$stderr.puts "Simplifying #{hash}"
		# Sort the contig by coordinates
		$contig_hash[hash].sort_by! { |val| val[1] } # Sort by end coordinate
		$contig_hash[hash].sort_by! { |val| val[0] } # Sort end-sorted results by first coordinate
		bed_sites = $contig_hash[hash][0] # Initialize array of output bed values
		$contig_hash[hash].shift
		if $contig_hash[hash].empty? # Write out results if only one entry for contig
			puts hash + "\t" + bed_sites[0].to_s + "\t" + (bed_sites[1]+1).to_s # Adjust end-site for BED semi-open format
		else
			for val in $contig_hash[hash]
				if val[0] > bed_sites[1] + 1 # Write out results if not overlapping or adjacent (hence +1) ; Reset current region
					puts hash + "\t" + bed_sites[0].to_s + "\t" + (bed_sites[1]+1).to_s # Adjust end-site for BED semi-open format
					bed_sites = val
				elsif val[1] > bed_sites[1] # Extend overlapping region
					bed_sites[1] = val[1]
				end
			end
			puts hash + "\t" + bed_sites[0].to_s + "\t" + (bed_sites[1]+1).to_s # Puts last entry for contig
		end
	end
end
