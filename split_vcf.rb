#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# split_vcf
SPLITVCFVER = "0.1.0"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

require_relative 'denovolib'

def split_vcf
	@vcfs = [] # Array of VCF names. Exclude filepath from names
	start = false
	header = ""
	outfile = ""
	outlines = "" # Lines to store until writing to disk
	writecycles = 0 # Number of write cycles passed
	gz_file_open($options.infile).open($options.infile) do |f1|
		while line = f1.gets
			if start
				contig = line.split("\t")[0]
				if contig != outfile
					if outfile != "" # Write remaining lines stored in memory
						File.open($options.outdir + "/chr" + outfile + ".vcf", 'a') do |write|
							write << outlines
						end
					end
					writecycles = 0 # Reset write cycles count
					outlines = "" # Reset outlines
					outfile = contig
					@vcfs.push("chr" + outfile) # Need to add chr to avoid qsub issue with contigs IDed by only a number
					File.open($options.outdir + "/chr" + outfile + ".vcf", 'w') do |write|
						write.puts header
						write.puts line
					end
				else
					writecycles += 1
					if writecycles == $options.writecycles
						File.open($options.outdir + "/chr" + outfile + ".vcf", 'a') do |write|
							write << outlines + line
						end
						writecycles = 0 # Reset for next round
						outlines = "" # Reset for next round
					else
						outlines << line
					end
				end
			elsif line[0..5] == "#CHROM"
				start = true
				header = line
			end
		end
	end
	File.open($options.outdir + "/chr" + outfile + ".vcf", 'a') do |write| # Output final line
		write << outlines
	end
	return @vcfs
end