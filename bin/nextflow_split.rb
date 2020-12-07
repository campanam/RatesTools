#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# split_nextflow
SPLITNFVCFVER = "0.1.1"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

require_relative 'denovolib'

#----------------------------------------------------------------------------------------
class Nextflow_Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		args.infile = "" # Input file
		args.outdir = "./" # Output directory
		opt_parser = OptionParser.new do |opts|
			opts.on("-i","--input [FILE]", String, "Input VCF") do |infile|
				args.infile = File.expand_path(infile) if infile != nil
			end
			opts.on("-o","--output [DIRECTORY]", String, "Output Directory (Default is current directory)") do |outdir|
				args.outdir = File.expand_path(outdir) if outdir != nil
			end
		end
		opt_parser.parse!(options)
		return args
	end
end
#----------------------------------------------------------------------------------------
$options = Nextflow_Parser.parse(ARGV)
Dir.mkdir($options.outdir) if !FileTest.directory?($options.outdir)
split_vcf(true)
