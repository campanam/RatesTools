#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# simplify_bed
SIMBEDVER = "0.1.2"
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
	gz_file_open(ARGV[0]).open(ARGV[0]) do |f1|
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
		bed_sites = [] # Array of output bed values
		for val in $contig_hash[hash]
			for site in val[0] .. val[1]
				bed_sites.push(site) unless bed_sites.include?(site)
			end
		end
		bed_sites = bed_sites.sort!.uniq
		start_site = nil
		end_site = nil
		unless bed_sites.empty?
			for i in 0 ... bed_sites.size
				start_site ||= bed_sites[i]
				end_site ||= start_site + 1
				if bed_sites[i] > end_site
					puts hash + "\t" + start_site.to_s + "\t" + (end_site).to_s
					start_site = bed_sites[i]
					end_site = start_site + 1
				elsif bed_sites[i] == end_site # Accounting for bizarre semi-open BED format
					end_site += 1
				end
			end
			puts hash + "\t" + start_site.to_s + "\t" + (end_site).to_s # Puts last entry
		end
	end
end
