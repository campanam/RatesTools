#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# indels2bed
INDELS2BEDVER = "0.1.2"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

if ARGV[0].nil?
	puts "\033[1mindels2bed " + INDELS2BEDVER + "\033[0m"
	puts "\nUsage: ruby indels2bed.rb <indels.vcf> <bp_to_exclude_upstream/downstream> > <out.bed>"
else
	$contigs = {} # Hash of contig lengths
	start = false
	File.open(ARGV[0]) do |f1|
		while line = f1.gets
			if line[0..5] == "#CHROM"
				start = true
			elsif start
				line_arr = line.split("\t")
				chrom = line_arr[0]
				chromStart = line_arr[1].to_i - ARGV[1].to_i - 1 # -1 due to BED open 0-indexed
				chromStart = 0 if chromStart < 0
				alleles = [line_arr[3],line_arr[4].split(",")].flatten
				allele_lengths = alleles.map { |x| x.length}
				max_allele_length = allele_lengths.max
				chromEnd = line_arr[1].to_i + ARGV[1].to_i + max_allele_length - 1 # 1 bp of indel is reference base
				chromEnd = $contigs[chrom] if chromEnd > $contigs[chrom]
				puts chrom + "\t" + chromStart.to_s + "\t" + chromEnd.to_s
			elsif line[0..12] == "##contig=<ID="
				key = line.split(",")[0][13..-1]
				value = line.split("=")[-1][0..-3].to_i #Remove final line break and closing >
				$contigs[key] = value
			end
		end
	end
end 