#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# denovolib
DENOVOLIBVER = "0.8.0"
# Michael G. Campana and Ellie E. Armstrong, 2021
# Smithsonian Institution and Stanford University

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to RatesTools;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.
 
# We politely request that this work be cited as:
# Campana, M.G. & E.E. Armstrong. 2020. RatesTools: Pipeline to calculate de novo
# mutation rates from parent-offspring trios. Smithsonian Institution and Stanford
# University. <https://github.com/campanam/RatesTools>.
#----------------------------------------------------------------------------------------

require 'ostruct'
require 'optparse'
require 'zlib'

# Library of methods accessed by other RatesTools ruby scripts

def split_vcf(gzip = false) # Split an input VCF into chromosome pieces for parallelization using paralle_denovo.rb or nextflow_split.rb
	@vcfs = [] # Array of VCF names. Exclude filepath from names
	start = false
	header = ""
	outfile = ""
	outlines = "" # Lines to store until writing to disk
	writecycles = 0 # Number of write cycles passed
	gz_file_open($options.infile) do |f1|
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
# From BaitsTools 1.6.8: Campana 2018
def gz_file_open(file) # Determine whether input file is gzipped or not and set method to open it
	if file[-3..-1] == ".gz"
		yield Zlib::GzipReader.open(file)
	else
		yield File.open(file)
	end
end
#-----------------------------------------------------------------------------------------
def filter_exit(message, snp_record) # Method to exit program if tag missing or data unsuited to specified filter
	$stderr.puts message
	$stderr.puts snp_record
	exit
end
#-----------------------------------------------------------------------------------------
def format_splash(cmd, version, cmdline) # Format output for basic script help screens
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
class Parser # Parse input for calc_denovo_mutation_rate.rb and post-hoc kochDNP.rb
	def self.parse(options, kochDNp = false)
		# Set defaults
		args = OpenStruct.new
		args.infile = "" # Input file
		args.outdir = "./" # Output directory
		args.sire = "" # Sire name 
		args.dam = "" # Dam name
		args.mu = 1e-6
		args.theta = 0.008
		args.cutoff = 0.3
		unless kochDNp
			args.window = 1000000 # Window length for bootstrapping
			args.minbslen = 1000000 # Minimum window length for bootstrapping
			args.minwindows = 10 # Minimum number of bootstrap windows
			args.step = 1000000 # Window step for bootrapping
			args.bootstrap = 0 # Number of bootstrap replicates
			args.gvcf = false # Whether input VCF is a gVCF
			args.rng = srand # Random number seed
			args.parhom = false # Flag to require homozygous for sire/dam DNMs
			args.minAD1 = false # Flag for exclusion of DNMs that parents have any alleles for candidate DNMs
			args.minAF = nil # Value for minimum frequency to include allele. Nil = off
			args.kochDNp = false # Flag to use Koch DNp filter
		end
		opt_parser = OptionParser.new do |opts|
			if kochDNp
				opts.banner = opts.banner = "\033[1mkochDNp " + KOCHDNPVER + "\033[0m"
				opts.separator ""
				opts.separator "Command-line usage: ruby kochDNp.rb [options] > <outfile>"
				opts.separator ""
				opts.separator "kochDNp options:"
			else
				opts.banner = "\033[1mcalc_denovo_mutation_rate " + CALCDENOVOVER + "\033[0m"
				opts.separator ""
				opts.separator "Command-line usage: ruby calc_denovo_mutation_rate.rb [options] > <outfile>"
				opts.separator ""
				opts.separator "calc_denovo_mutation_rate options:"
			end			
			opts.on("-i","--input [FILE]", String, "Input VCF") do |infile|
				args.infile = File.expand_path(infile) if infile != nil
			end
			opts.on("-s","--sire [NAME]", String, "Sire's name in VCF") do |sire|
				args.sire = sire if sire != nil
			end
			opts.on("-d", "--dam [NAME]", String, "Dam's name in VCF") do |dam|
				args.dam = dam if dam != nil
			end
			unless kochDNp
				opts.on("-k", "--kochDNp", "Use Koch DNp statistic to filter DNM sites") do |kochDNP|
					args.kochDNp = true
				end
			end
			opts.on("-m", "--mu [VALUE]", Float, "Mutation rate (mu) for Koch DNp (Default = 1e-6)") do |mu|
				args.mu = mu if mu != nil
			end
			opts.on("-t", "--theta [VALUE]", Float, "Heterozygosity (theta) for Koch DNp (Default = 0.008)") do |theta|
				args.theta = theta if theta != nil
			end
			opts.on("-c", "--cutoff [VALUE]", Float, "Koch DNp cutoff (Default = 0.3)") do |cutoff|
				args.cutoff = cutoff if cutoff != nil
			end
			unless kochDNp
				opts.on("--parhom", "Require parents to be homozygous at DNM sites") do |parhom|
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
			end
			opts.separator ""
			opts.separator "Program information:"
			opts.on("-v", "--version", "Show program version") do
				if kochDNp
					puts "kochDNp version " + KOCHDNPVER
				else 
					puts "calc_denovo_mutation_rate version " + CALCDENOVOVER
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
