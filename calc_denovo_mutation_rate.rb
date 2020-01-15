#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# calc_denovo_mutation_rate
CALCDENOVOVER = "0.6.0"
# Michael G. Campana, 2019-2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

require 'ostruct'
require 'optparse'
require 'zlib'

$previous_snp = nil # Global variable storing previous SNP information if phasing is known. nil when new contig/scaffold or new phasing block.
$total_sites = 0 # Total number of retained sites
$total_denovo = {} # Hash of denovo variants keyed by offspring
$sample_index = {} # Hash of offspring names keyed by index
$windows = [] # Array of completed bootstrap windows
$current_windows = []  # Array of active uncompleted bootstrap windows
$bootstraps = {} # Hash of bootstrap replicate summary values keyed by offspring
$mutations = "\nIdentified de novo mutations:\n" # Reduced VCF of ided mutations
#-----------------------------------------------------------------------------------------
class Snp
	attr_accessor :sire, :dam, :offspring, :alleles, :denovo
	def initialize(sire = "", dam = "", offspring = {}, alleles = []) # Offspring is an array because multiple offspring
		@alleles = alleles # Array of alleles sorted by allele ID (so allele 0 is in array slot 0, allele 1 in slot 1, etc)
		@sire = sire # Sire genotype, Maintenance of phasing information as string
		@dam = dam # Dam genotype, Maintenance of phasing information as string
		@offspring = offspring # Hash of offspring genotypes, Maintenance of phasing information as string
		@denovo = {} # Hash of denovo mutation presence statuses and types sorted by offspring
	end
	#-------------------------------------------------------------------------------------
	def split_alleles(alleles)
		alleles.include?("/") ? split_chr = "/" : split_chr = "|"
		return alleles.split(split_chr).sort, alleles.include?("|") # Sorting here for simplicity of later comparisons. Returning the splitting for phasing information
	end
	#-------------------------------------------------------------------------------------
	def classify_denovo(sire_alleles, dam_alleles, offspr_alleles, offspr) # Determine whether denovo is forward or reverse mutation(s)
		sire_alleles.uniq.size == 1 ? sire_hom = true : sire_hom = false # Determine sire homozygosity
		dam_alleles.uniq.size == 1 ? dam_hom = true : dam_hom = false # Determine dam homozygosity
		allele_pool = (sire_alleles + dam_alleles).flatten.uniq # Get all parental alleles
		offspr_alleles.uniq.size == 1 ? offspr_hom = true : offspr_hom = false # Determine offspring homozygosity
		mutation_class = "dummy" # Dummy variable to prevent compile error while working out logic
		if !allele_pool.include?(offspr_alleles[0]) and !allele_pool.include?(offspr_alleles[1]) # Test for a double-forward mutation
			mutation_class = "double-forward"
		elsif allele_pool.include?(offspr_alleles[0]) and !allele_pool.include?(offspr_alleles[1]) # Test for a de-novo single forward
			mutation_class = "single-forward"
		elsif !allele_pool.include?(offspr_alleles[0]) and allele_pool.include?(offspr_alleles[1]) # Inverse of previous test
			mutation_class = "single-forward"
		elsif offspr_alleles == sire_alleles or offspr_alleles == dam_alleles # Test for a reverse mutation to parental type
			mutation_class = "backward"
		elsif offspr_hom and sire_alleles.include?(offspr_alleles[0]) and !dam_alleles.include?(offspr_alleles[0]) # Test for a reverse mutation to parental type
			mutation_class = "backward"
		elsif offspr_hom and !sire_alleles.include?(offspr_alleles[0]) and dam_alleles.include?(offspr_alleles[0]) # Inverse of previous test
			mutation_class = "backward"
		end
		case mutation_class
		when "single-forward"
			update_denovo_counts(offspr, 0)
		when "double-forward"
			update_denovo_counts(offspr, 1)
		when "backward"
			update_denovo_counts(offspr, 2)
		end
		return mutation_class
	end
	#-------------------------------------------------------------------------------------
	def update_denovo_counts(offspr, type)
		$total_denovo[offspr][type] += 1
		for window in $current_windows
			window.denovo[offspr][type] += 1
		end
	end
	#-------------------------------------------------------------------------------------
	def check_phasing(poss_genotypes)
		# NEED TO FIGURE OUT HOW TO USE GENOTYPING PHASING HERE
		# Compare offspring genotypes against parents' where known
		# Possible going forward, but also in reverse direction?
		unless $previous_snp.nil?
		
		end
		return poss_genotypes
	end
	#-------------------------------------------------------------------------------------
	def identify_denovo
		# Calculate possible parental alleles without mutation
		poss_genotypes_unphased = [] # Array holding possible unphased genotypes without mutation
		sire_alleles, sire_phased = split_alleles(@sire)
		dam_alleles, dam_phased = split_alleles(@dam)
		poss_genotypes_unphased.push([sire_alleles[0],dam_alleles[0]].sort) # Resorting because could change during recombination
		poss_genotypes_unphased.push([sire_alleles[0],dam_alleles[1]].sort)
		poss_genotypes_unphased.push([sire_alleles[1],dam_alleles[0]].sort)
		poss_genotypes_unphased.push([sire_alleles[1],dam_alleles[1]].sort)
		poss_genotypes_unphased = poss_genotypes_unphased.uniq # If only one possible, uniq! will return nil
		for offspr in offspring.keys
			offspr_alleles, offspr_phased = split_alleles(offspring[offspr])
			offspr_phased && dam_phased && sire_phased ? genopool = check_phasing(poss_genotypes_unphased) : genopool = poss_genotypes_unphased # Determine whether to use phased or unphased genotypes	
			denovo_status = [false,""] # Default status is no de novo mutation
			unless genopool.include?(offspr_alleles)
				mutation_class = classify_denovo(sire_alleles, dam_alleles, offspr_alleles, offspr)
				denovo_status = [true,mutation_class]
			end
			@denovo[offspr] = denovo_status
		end
	end
end
#-----------------------------------------------------------------------------------------
class Bootstrap_Window
	attr_accessor :startbp, :endbp, :denovo # Window start, end coordinates. Hash of observed mutations
	def initialize(startbp)
		@startbp = startbp
		@denovo = {}
		for key in $total_denovo.keys
			@denovo[key] = [0,0,0]
		end
	end
	def length # Get length of bootstrap window
		return @endbp - @startbp + 1
	end
end
#-----------------------------------------------------------------------------------------------
# From BaitsTools 1.6.5: Campana et al. 2018
def gz_file_open(file)
	if file[-3..-1] == ".gz"
		return "Zlib::GzipReader"
	else
		return "File"
	end
end
#-----------------------------------------------------------------------------------------
def read_vcf # Method to read vcf
	collect_data = false # Flag to start collecting data
	next_window_site = 1 # Site to start new window
	next_window_close_site = next_window_site + $options.window - 1 # Site to close window and calculate rate
	eval(gz_file_open($options.infile)).open($options.infile) do |f1|
		while line = f1.gets
			if line[0..5] == "#CHROM"
				$mutations << line
				collect_data = true
				header_arr = line[0..-2].split("\t")
				for sample in header_arr[9..-1]
					if sample == $options.sire
						sire_index = header_arr.index(sample)
					elsif sample == $options.dam
						dam_index = header_arr.index(sample)
					else
						$total_denovo[sample] = [0,0,0] # Total single-forward, double-forward, backward
					end
					$sample_index[header_arr.index(sample)] = sample # Push the name into the index for identification
				end
			elsif collect_data
				snp_array = line[0..-2].split("\t")
				$total_sites += 1
				$stderr.puts $total_sites.to_s + " processed" if $total_sites % 1000000 == 0
				$total_sites += snp_array[7].split("=")[1].to_i - snp_array[1].to_i if $options.gvcf # Adjust for gVCF blocks
				while $total_sites >= next_window_site # Allow multiple windows to pass if gVCF block sufficiently long
					$current_windows.push(Bootstrap_Window.new(next_window_site))
					next_window_site += $options.step
				end
				while $total_sites >= next_window_close_site
					$current_windows[0].endbp = next_window_close_site
					$windows.push($current_windows.shift)
					next_window_close_site += $options.step
				end
				alleles = ([snp_array[3]] + snp_array[4].split(",")).flatten.uniq # Get alleles
				alleles.delete("<NON_REF>") if alleles.include?("<NON_REF>")
				if alleles.size > 1 # Only process sites with actual SNPs in them
					snp = Snp.new
					snp.alleles = alleles
					for i in 9...snp_array.size
						genotype = snp_array[i].split(":")[0] # Assuming GT in first slot for now
						if i == sire_index
							snp.sire = genotype
						elsif i == dam_index
							snp.dam = genotype
						else
							snp.offspring[$sample_index[i]] = genotype
						end
					end
					snp.identify_denovo
					printline = false # check whether to print site
					for offspr in snp.offspring.keys
						if snp.denovo[offspr][0]
							printline = true 
							break # maybe save a bit of processing
						end
					end
					$mutations << line if printline
				end
			end
		end
	end
	# Finish remaining uncompleted windows
	for window in $current_windows
		window.endbp = $total_sites
		$windows.push(window) if window.length >= $options.minbslen # Discard undersized windows
	end
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
		puts offspr + "\t" + ($total_denovo[offspr].sum.to_f/$total_sites.to_f).to_s + "\t" + ($total_denovo[offspr][0].to_f/$total_sites.to_f).to_s
	end
	puts $mutations
end
#-----------------------------------------------------------------------------------------
def mean(values)
	return values.sum.to_f/values.size.to_f
end
#-----------------------------------------------------------------------------------------
def conf95(values)
	meanval = mean(values)
	sdnum = values.map { |x| (x - meanval) * (x + meanval) }
	stdev = Math.sqrt(sdnum.sum/(values.size.to_f - 1.0))
	critval = 1.96 * stdev/Math.sqrt(values.size.to_f)
	finalstring = meanval.to_s + "\t" + (meanval - critval).to_s + "-" + (meanval + critval).to_s
	return finalstring
end
#-----------------------------------------------------------------------------------------
def bootstrap_results
	for offspr in $total_denovo.keys
		$bootstraps[offspr] = [[],[]]
	end
	for i in 1..$options.bootstrap
		puts "\nBootstrap Replicate " + i.to_s + ":"
		current_bootstrap_length = 0
		bootstrap_denovo = {}
		for offspr in $total_denovo.keys
			bootstrap_denovo[offspr] = [0,0,0]
		end
		while current_bootstrap_length < $total_sites
			selected_window = $windows[rand($windows.size)]
			current_bootstrap_length += selected_window.length
			for offspr in $total_denovo.keys
				for val in 0..2
					bootstrap_denovo[offspr][val] += selected_window.denovo[offspr][val]
				end
			end
		end
		puts "Total number of retained sites: " + current_bootstrap_length.to_s
		puts "\nTotal numbers of observed de novo mutations:"
		puts "Offspring\tSingle-Forward\tDouble-Forward\tBackward"
		for offspr in bootstrap_denovo.keys
			puts offspr + "\t" + bootstrap_denovo[offspr].join("\t")
		end
		puts "\nInferred mutation rates:"
		puts "Offspring\tAllsites\tSingle-ForwardOnly"
		for offspr in bootstrap_denovo.keys
			puts offspr + "\t" + (bootstrap_denovo[offspr].sum.to_f/current_bootstrap_length.to_f).to_s + "\t" + (bootstrap_denovo[offspr][0].to_f/current_bootstrap_length.to_f).to_s
			$bootstraps[offspr][0].push(bootstrap_denovo[offspr].sum.to_f/current_bootstrap_length.to_f)
			$bootstraps[offspr][1].push(bootstrap_denovo[offspr][0].to_f/current_bootstrap_length.to_f)
		end
	end
	puts "\nBootstrapped Estimates:"
	puts "Offspring\tMean_AllSites\t95%C.I._Allsites\tMean_Single-ForwardOnly\t95%C.I._Single-ForwardOnly:"
	for offspr in $bootstraps.keys
		puts offspr + "\t" + conf95($bootstraps[offspr][0]) + "\t" + conf95($bootstraps[offspr][1])
	end
end
#-----------------------------------------------------------------------------------------
def print_options
	puts "calc_denovo_mutation_rate " + CALCDENOVOVER + " started with parameters:"
	cmdline = "-i " + $options.infile + " -s " + $options.sire + " -d " + $options.dam + " -w " + $options.window.to_s + " -S " + $options.step.to_s + " -l " + $options.minbslen.to_s + " -M " + $options.minwindows.to_s
	cmdline << " -g" if $options.gvcf
	cmdline << " --rng " + $options.rng.to_s
	puts cmdline
end
#-----------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		args.infile = "" # Input file
		args.sire = "" # Sire name 
		args.dam = "" # Dam name
		args.window = 1000000 # Window length for bootstrapping
		args.minbslen = 1000000 # Minimum window length for bootstrapping
		args.step = 1000000 # Window step for bootrapping
		args.minwindows = 10 # Minimum number of bootstrap windows
		args.bootstrap = 0 # Number of bootstrap replicates
		args.gvcf = false # Whether input VCF is a gVCF
		args.rng = srand # Random number seed
		opt_parser = OptionParser.new do |opts|
			opts.banner = "\033[1mcalc_denovo_mutation_rate " + CALCDENOVOVER + "\033[0m"
			opts.separator ""
			opts.separator "Command-line usage: ruby calc_denovo_mutation_rate.rb [options]"
			opts.on("-i","--input [FILE]", String, "Input file") do |infile|
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
			opts.on("-v", "--version", "Show program version") do
				puts "calc_denovo_mutation_rate version " + CALCDENOVOVER
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
$options = Parser.parse(ARGV)
$options.minbslen = $options.window if $options.minbslen > $options.window # Prevent nonsensical results
srand($options.rng)
print_options
read_vcf
print_results
bootstrap_results if ($options.bootstrap > 0 && $windows.size >= $options.minwindows)
