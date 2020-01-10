#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# calc_denovo_mutation_rate
CALCDENOVOVER = "0.1.0"
# Michael G. Campana, 2019-2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

require 'ostruct'
require 'optparse'

$previous_snp = nil # Global variable storing previous SNP information if phasing is known. nil when new contig/scaffold or new phasing block.
$total_sites = 0 # Total number of retained sites
$total_denovo = {} # Hash of denovo variants keyed by offspring
$sample_index = {} # Hash of offspring names keyed by index
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
			$total_denovo[offspr][0] += 1
		when "double-forward"
			$total_denovo[offspr][1] += 1
		when "backward"
			$total_denovo[offspr][2] += 1
		end
		return mutation_class
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
def read_vcf # Method to read vcf
	collect_data = false # Flag to start collecting data
	File.open($options.infile) do |f1|
		while line = f1.gets
			if line[0..5] == "#CHROM"
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
				$total_sites += snp_array[7].split("=")[1].to_i - snp_array[1].to_i if $options.gvcf # Adjust for gVCF blocks
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
				end
			end
		end
	end
end
#-----------------------------------------------------------------------------------------
def print_results # Method to print basic results
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
end
#-----------------------------------------------------------------------------------------
class Parser
	def self.parse(options)
		# Set defaults
		args = OpenStruct.new
		args.infile = "" # Input file
		args.sire = "" # Sire name 
		args.dam = "" # Dam name
		args.gvcf = false # Whether input VCF is a gVCF
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
			opts.on("-g", "--gvcf", "Input is a gVCF (Default = false)") do |gvcf|
				args.gvcf = true
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
read_vcf
print_results
