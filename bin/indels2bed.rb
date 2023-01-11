#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# indels2bed
INDELS2BEDVER = "0.3.3"
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

# Script to identify indel regions from a VCF, exclude a set number of bps around indels, and output BED for VCFtools exclusion

require_relative 'denovolib'

if ARGV[0].nil?
	# If no parameters passed, print basic help screen
	format_splash('indels2bed', INDELS2BEDVER, '<indels.vcf[.gz]> <bp_to_exclude_upstream/downstream> > <out.bed>')
else
	# If given VCF input, open file and look for indels (including those not annotated as indels) and generate exclusion BED
	$contigs = {} # Hash of contig lengths
	start = false
	gz_file_open(ARGV[0]) do |f1|
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
				if max_allele_length > 1
					chromEnd = line_arr[1].to_i + ARGV[1].to_i + max_allele_length - 1 # 1 bp of indel is reference base
					chromEnd = $contigs[chrom] if chromEnd > $contigs[chrom]
					puts chrom + "\t" + chromStart.to_s + "\t" + chromEnd.to_s
				elsif alleles.include?("*")
					chromEnd = line_arr[1].to_i + ARGV[1].to_i
					chromEnd = $contigs[chrom] if chromEnd > $contigs[chrom]
					puts chrom + "\t" + chromStart.to_s + "\t" + chromEnd.to_s
				end
			elsif line[0..12] == "##contig=<ID="
				key = line.split(",")[0][13..-1]
				value = line.split("=")[-1][0..-3].to_i #Remove final line break and closing >
				$contigs[key] = value
			end
		end
	end
end 
