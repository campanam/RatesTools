#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# denovolib
DENOVOLIBVER = "0.7.0"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

require 'ostruct'
require 'optparse'
require 'zlib'

def split_vcf(gzip = false)
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
						`gzip #{$options.outdir + "/chr" + outfile + ".vcf"}` if gzip # If outputting compressed vcfs
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
	`gzip #{$options.outdir + "/chr" + outfile + ".vcf"}` if gzip # If outputting compressed vcfs
	return @vcfs
end
#-----------------------------------------------------------------------------------------------
# From BaitsTools 1.6.6: Campana 2018
def gz_file_open(file)
	if file[-3..-1] == ".gz"
		return Zlib::GzipReader
	else
		return File
	end
end
#-----------------------------------------------------------------------------------------
def format_splash(cmd, version, cmdline)
	puts "\033[1m#{cmd} #{version}\033[0m"
	puts "\nUsage: ruby #{cmd}.rb #{cmdline}"
end
#-----------------------------------------------------------------------------------------
def print_results # Method to print basic results
	puts "\nTotal Sample Results:"
	puts "Total number of retained sites: " + $total_sites.to_s
	puts "\nTotal numbers of observed de novo mutations:"
	puts "Offspring\tSingle-Forward\tDouble-Forward\tBackward"
	for offspr in $total_denovo.keys
		puts offspr + "\t" + $total_denovo[offspr].join("\t")
	end
	puts "\nInferred mutation rates:"
	puts "Offspring\tAllsites\tSingle-ForwardOnly"
	for offspr in $total_denovo.keys
		sitedenom = 2 * $total_sites.to_f
		puts offspr + "\t" + ($total_denovo[offspr].sum.to_f/sitedenom).to_s + "\t" + ($total_denovo[offspr][0].to_f/sitedenom).to_s
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
		args.parhom = false # Flag to require homozygous for sire/dam DNMs
		args.minAD1 = false # Flag for exclusion of DNMs that parents have any alleles for
		args.minAF = nil # Value for minimum frequency to include allele. Nil = off
		if parallel
			args.writecycles = 1000000 # Number of variants to read before writing to disk
			args.memory = "1G" # Memory reserved
			args.himem = false # High memory flag
			args.lopri = false # Low priority falg
			args.queue = "sThC.q" # Qsub queue
			args.email = "" # Email address to notify
			args.restart = false # Restart if partitioned VCFs already exit
			args.submit = false # Submit previously generated jobs
			args.nosubmit = false # Generate split VCFs and jobs, but do not submit
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
			opts.on("--parhom", "Require parents to be homozoygous at DNM sites") do |parhom|
				args.parhom = true
			end
			opts.on("--minAD1", "Discard DNMs if parents have DNM alleles even if not called") do |minAD1|
				args.minAD1 = true
			end
			opts.on("--minAF [VALUE]", Float, "Filter alleles by minimum frequency") do |minAF|
				args.minAF = minAF
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
				opts.on("--nosubmit", "Generate split VCFs and jobs, but do not submit them") do |nosubmit|
					args.nosubmit = true
				end
				opts.on("-r", "--restart", "Restart from previously split VCFs (Default = false)") do |restart|
					args.restart = true
				end
				opts.on("--submit", "Submit previously generated jobs and split VCFs (Implies -r)") do |submit|
					args.submit = true
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
