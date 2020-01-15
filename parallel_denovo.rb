#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# parallel_denovo
PARDENOVOVER = "0.5.0"
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
	outlines = "" # Lines to store until writing to disk
	writecycles = 0 # Number of write cycles passed
	eval(gz_file_open($options.infile)).open($options.infile) do |f1|
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
#-----------------------------------------------------------------------------------------------
def get_vcfs
	@vcfs = []
	Dir.foreach($options.outdir + "/") do |f1|
		@vcfs.push(f1[0..-5]) if f1[-3..-1] == "vcf"
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
	system("qsub -N #{file} -v SCAFFOLD=#{file} #{$options.outdir}/calc_denovo_mutation_rate.job")
end
#-----------------------------------------------------------------------------------------
def write_qsub
	header = "#!/bin/sh\n#$ -S /bin/sh\n#$ -q #{$options.queue}\n#$ -l mres=#{$options.memory},h_data=#{$options.memory},h_vmem=#{$options.memory}"
	header << ",himem" if $options.himem
	header << ",lopri" if $options.lopri
	header << "\n#$ -j y\n#$ -cwd\n"
	header << "#$ -m bea\n#$ -M #{$options.email}\n" if $options.email != ""
	header << "module load bioinformatics/ruby/2.6.3\n"
	header << "ruby calc_denovo_mutation_rate.rb -i #{$options.outdir}/${SCAFFOLD}.vcf -s #{$options.sire} -d #{$options.dam} -w #{$options.window} -S #{$options.step} -b #{$options.bootstrap} --rng #{$options.rng}"
	header << " -g" if $options.gvcf
	header << " > #{$options.outdir}/${SCAFFOLD}.log"
	File.open("#{$options.outdir}/calc_denovo_mutation_rate.job", 'w') do |qsub|
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
		args.minwindows = 10 # Minimum number of bootstrap windows
		args.step = 1000000 # Window step for bootrapping
		args.bootstrap = 0 # Number of bootstrap replicates
		args.gvcf = false # Whether input VCF is a gVCF
		args.rng = srand # Random number seed
		args.writecycles = 1000000 # Number of variants to read before writing to disk
		args.memory = "1G" # Memory reserved
		args.himem = false # High memory flag
		args.lopri = false # Low priority falg
		args.queue = "sThC.q" # Qsub queue
		args.email = "" # Email address to notify
		args.restart = false # Restart if partitioned VCFs already exit
		opt_parser = OptionParser.new do |opts|
			opts.banner = "\033[1mparallel_denovo " + PARDENOVOVER + "\033[0m"
			opts.separator ""
			opts.separator "Command-line usage: ruby parallel_denovo.rb [options]"
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
$options.restart ? vcfs = get_vcfs : vcfs = split_vcf
write_qsub
for vcf in vcfs
	execute_qsub(vcf)
end