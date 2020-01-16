#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# denovolib
DENOVOLIBVER = "0.1.0"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

require 'ostruct'
require 'optparse'
require 'zlib'

#-----------------------------------------------------------------------------------------------
# From BaitsTools 1.6.6: Campana et al. 2018
def gz_file_open(file)
	if file[-3..-1] == ".gz"
		return Zlib::GzipReader
	else
		return File
	end
end
#-----------------------------------------------------------------------------------------
class Parser
	def self.parse(options, parallel = false)
		# Set defaults
		args = OpenStruct.new
		args.infile = "" # Input file
		args.outdir = "./" # Output directory
		args.sire = "" # Sire name 
		args.dam = "" # Dam name
		args.window = 1000000 # Window length for bootstrapping
		args.minbslen = 1000000 # Minimum window length for bootstrapping
		args.minwindows = 10 # Minimum number of bootstrap windows
		args.step = 1000000 # Window step for bootrapping
		args.bootstrap = 0 # Number of bootstrap replicates
		args.gvcf = false # Whether input VCF is a gVCF
		args.rng = srand # Random number seed
		if parallel
			args.writecycles = 1000000 # Number of variants to read before writing to disk
			args.memory = "1G" # Memory reserved
			args.himem = false # High memory flag
			args.lopri = false # Low priority falg
			args.queue = "sThC.q" # Qsub queue
			args.email = "" # Email address to notify
			args.restart = false # Restart if partitioned VCFs already exit
		end
		opt_parser = OptionParser.new do |opts|
			if parallel
				opts.banner = "\033[1mparallel_denovo " + PARDENOVOVER + "\033[0m"
			else
				opts.banner = "\033[1mcalc_denovo_mutation_rate " + CALCDENOVOVER + "\033[0m"
			end
			opts.separator ""
			if parallel
				opts.separator "Command-line usage: ruby parallel_denovo.rb [options]"
			else
				opts.separator "Command-line usage: ruby calc_denovo_mutation_rate.rb [options]"
			end
			opts.separator ""
			opts.separator "calc_denovo_mutation_rate options:"
			opts.on("-i","--input [FILE]", String, "Input VCF") do |infile|
				args.infile = File.expand_path(infile) if infile != nil
			end
			opts.on("-s","--sire [NAME]", String, "Sire's name in VCF") do |sire|
				args.sire = sire if sire != nil
			end
			opts.on("-d", "--dam [NAME]", String, "Dam's name in VCF") do |dam|
				args.dam = dam if dam != nil
			end
			opts.on("-w", "--window [VALUE]", Integer, "Sequence window length (bp) for bootstrapping (Default = 1000000)") do |window|
				args.window = window if window != nil
			end
			opts.on("-S", "--step [VALUE]", Integer, "Window step (bp) for bootstrapping (Default = 1000000)") do |step|
				args.step = step if step != nil
			end
			opts.on("-b", "--bootstrap [VALUE]", Integer, "Number of bootstrap replicates (Default = 0)") do |bs|
				args.bootstrap = bs if bs != nil
			end
			opts.on("-l", "--minbootstraplength [VALUE]", Integer, "Minimum bootstrap window length (bp) to retain (Default = 1000000)") do |minbslen|
				args.minbslen = minbslen if minbslen != nil
			end
			opts.on("-M", "--minwindows [VALUE]", Integer, "Minimum number of bootstrap windows to retain contig (Default = 10)") do |minwindows|
				args.minwindows = minwindows if minwindows != nil
				args.minwindows = 1 if args.minwindows < 1
			end
			opts.on("-g", "--gvcf", "Input is a gVCF (Default = false)") do |gvcf|
				args.gvcf = true
			end
			opts.on("--rng [VALUE]", Integer, "Random number seed") do |rng|
				args.rng = rng if rng != nil
			end
			if parallel
				opts.separator ""
				opts.separator "parallel_denovo Options:"
				opts.on("-o","--output [DIRECTORY]", String, "Output Directory (Default is current directory)") do |outdir|
					args.outdir = File.expand_path(outdir) if outdir != nil
				end
				opts.on("-W", "--writecycles [VALUE]", Integer, "Number of variants to read before writing to disk (Default = 1000000)") do |wrt|
					args.writecycles = wrt if wrt != nil
					args.writecycles = 1 if args.writecycles < 1
				end
				opts.on("-r", "--restart", "Restart from previously split VCFs (Default = false)") do |restart|
					args.restart = true
				end
				opts.separator ""
				opts.separator "SI/HPC Options:"
				opts.on("-q", "--queue [VALUE]", String, "Qsub queue to use (Default = sThC.q") do |queue|
					args.queue = queue if queue != nil
				end
				opts.on("-m", "--memory [VALUE]", String, "Reserved memory (Default = 1G)") do |memory|
					args.memory = memory if memory != nil
				end
				opts.on("-H", "--himem", "Use high-memory queue (Default is false)") do |himem|
					args.himem = true
				end
				opts.on("-L", "--lopri", "Use low priority queue (Default is false)") do |lopri|
					args.lopri = true
				end
				opts.on("-e", "--email [VALUE]", String, "E-mail address to notify") do |email|
					args.email = email if email != nil
				end
			end
			opts.separator ""
			opts.separator "Program information:"
			opts.on("-v", "--version", "Show program version") do
				if parallel
					puts "parallel_denovo version " + PARDENOVOVER
				else 
					puts "calc_denovo_mutation rate version " + CALCDENOVOVER
				end
				exit
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