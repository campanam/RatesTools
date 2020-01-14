#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# parallel_denovo
PARDENOVOVER = "0.2.0"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

require 'ostruct'
require 'optparse'
require 'zlib'

def split_vcf
	@vcfs = [] # Array of VCF names. Exclude filepath from names
	start = false
	header = ""
	outfile = ""
	eval(gz_file_open($options.infile)).open($options.infile) do |f1|
		while line = f1.gets
			if start
				contig = line.split("\t")[0]
				if contig != outfile
					outfile = contig
					@vcfs.push("chr" + outfile) # Need to add chr to avoid qsub issue with contigs IDed by only a number
					File.open($options.outdir + "/chr" + outfile + ".vcf", 'w') do |write|
						write.puts header
						write.puts line
					end
				else
					File.open($options.outdir + "/chr" + outfile + ".vcf", 'a') do |write|
						write << line
					end
				end
			elsif line[0..5] == "#CHROM"
				start = true
				header = line
			end
		end
	end
	return @vcfs
end
#-----------------------------------------------------------------------------------------------
# From BaitsTools 1.6.5: Campana et al. 2018 -- Some redundant code
def gz_file_open(file)
	if file[-3..-1] == ".gz"
		return "Zlib::GzipReader"
	else
		return "File"
	end
end
#-----------------------------------------------------------------------------------------
def execute_qsub(file)
	system("qsub #{$options.outdir}/#{file}.job")
end
#-----------------------------------------------------------------------------------------
def write_qsub(vcf)
	header = "#!/bin/sh\n#$ -S /bin/sh\n#$ -q #{$options.queue}\n#$ -l mres=#{$options.memory},h_data=#{$options.memory},h_vmem=#{$options.memory}"
	header << ",himem" if $options.himem
	header << ",lopri" if $options.lopri
	header << "\n#$ -j y\n#$ -cwd\n#$ -N #{vcf}\n"
	header << "#$ -m bea\n#$ -M #{$options.email}\n" if $options.email != ""
	header << "module load bioinformatics/ruby/2.6.3\n"
	header << "ruby calc_denovo_mutation_rate.rb -i #{$options.outdir}/#{vcf}.vcf -s #{$options.sire} -d #{$options.dam} -w #{$options.window} -S #{$options.step} -b #{$options.bootstrap} --rng #{$options.rng}"
	header << " -g" if $options.gvcf
	header << " > #{$options.outdir}/#{vcf}.log"
	File.open("#{$options.outdir}/#{vcf}.job", 'w') do |qsub|
		qsub.puts header
	end
end
#-----------------------------------------------------------------------------------------
class ParallelParser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		args.infile = "" # Input file
		args.outdir = "./" # Output directory
		args.sire = "" # Sire name 
		args.dam = "" # Dam name
		args.window = 1000000 # Window length for bootstrapping
		args.minbslen = 1000000 # Minimum window length for bootstrapping
		args.step = 1000000 # Window step for bootrapping
		args.bootstrap = 0 # Number of bootstrap replicates
		args.gvcf = false # Whether input VCF is a gVCF
		args.rng = srand # Random number seed
		args.memory = "1G" # Memory reserved
		args.himem = false # High memory flag
		args.lopri = false # Low priority falg
		args.queue = "sThC.q" # Qsub queue
		args.email = "" # Email address to notify
		opt_parser = OptionParser.new do |opts|
			opts.banner = "\033[1mparallel_denovo " + PARDENOVOVER + "\033[0m"
			opts.separator ""
			opts.separator "Command-line usage: ruby parallel_denovo.rb [options]"
			opts.separator ""
			opts.separator "calc_denovo_mutation_rate options:"
			opts.on("-i","--input [FILE]", String, "Input VCF") do |infile|
				args.infile = File.expand_path(infile) if infile != nil
			end
			opts.on("-o","--output [DIRECTORY]", String, "Output Directory (Default is current directory)") do |outdir|
				args.outdir = File.expand_path(outdir) if outdir != nil
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
			opts.on("-g", "--gvcf", "Input is a gVCF (Default = false)") do |gvcf|
				args.gvcf = true
			end
			opts.on("--rng [VALUE]", Integer, "Random number seed") do |rng|
				args.rng = rng if rng != nil
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
			opts.separator ""
			opts.separator "Program information:"
			opts.on("-v", "--version", "Show program version") do
				puts "parallel_denovo version " + PARDENOVOVER
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
#-----------------------------------------------------------------------------------------
ARGV[0] ||= "-h"
$options = ParallelParser.parse(ARGV)
Dir.mkdir($options.outdir) if !FileTest.directory?($options.outdir)
vcfs = split_vcf
for vcf in vcfs
	write_qsub(vcf)
	execute_qsub(vcf)
end