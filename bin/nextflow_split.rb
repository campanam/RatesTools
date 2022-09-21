#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# nextflow_split
SPLITNFVCFVER = "0.1.3"
# Michael G. Campana and Ellie E. Armstrong, 2020
# Smithsonian Institution and Stanford University

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to RatesTools;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.
 
# We politely request that this work be cited as:
# Armstrong, E.E. & M.G. Campana. 2022. RatesTools: a Nextflow pipeline for detecting
# de novo germline mutations in pedigree sequence data. *bioRxiv*.
# doi: 10.1101/2022.07.18.500472.
#----------------------------------------------------------------------------------------

# Script to split VCFs for parallel nextflow, bypassing calculations in other Ruby scripts

require_relative 'denovolib'

#----------------------------------------------------------------------------------------
class Nextflow_Parser # Reduced parser including only input/output and bypassing other Ruby scripts
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		args.infile = "" # Input file
		args.outdir = "./" # Output directory
		args.writecycles = 1000000 # number of cycles before writing to disk
		opt_parser = OptionParser.new do |opts|
			opts.on("-i","--input [FILE]", String, "Input VCF") do |infile|
				args.infile = File.expand_path(infile) if infile != nil
			end
			opts.on("-o","--output [DIRECTORY]", String, "Output Directory (Default is current directory)") do |outdir|
				args.outdir = File.expand_path(outdir) if outdir != nil
			end
			opts.on("-W", "--writecycles [VALUE]", Integer, "Number of variants to read before writing to disk (Default = 1000000)") do |wrt|
				args.writecycles = wrt if wrt != nil
				args.writecycles = 1 if args.writecycles < 1
			end
			opts.on_tail("-h","--help", "Show help") do
				puts opts
				exit
			end
		end
		opt_parser.parse!(options)
		return args
	end
end
#----------------------------------------------------------------------------------------
ARGV[0] ||= "-h" # Print help if no parameters passed
$options = Nextflow_Parser.parse(ARGV)
Dir.mkdir($options.outdir) if !FileTest.directory?($options.outdir)
split_vcf(true)
