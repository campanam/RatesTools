#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# calc_denovo_mutation_rate
CALCDENOVOVER = "0.11.0"
# Michael G. Campana and Ellie E. Armstrong, 2019-2021
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

# This script calculates the de novo mutation rate from trios from a provided filtered VCF

require_relative 'denovolib'

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
	# Class to define SNPs and classify them
	attr_accessor :sire, :dam, :offspring, :alleles, :denovo, :sire_pl, :dam_pl, :offspring_pl
	def initialize(sire = "", dam = "", offspring = {}, alleles = [], sire_pl = [], dam_pl = [], offspring_pl = {}) # Offspring is an array because multiple offspring
		@alleles = alleles # Array of alleles sorted by allele ID (so allele 0 is in array slot 0, allele 1 in slot 1, etc)
		@sire = sire # Sire genotype, Maintenance of phasing information as string
		@dam = dam # Dam genotype, Maintenance of phasing information as string
		@offspring = offspring # Hash of offspring genotypes, Maintenance of phasing information as string
		@denovo = {} # Hash of denovo mutation presence statuses and types sorted by offspring
		@sire_pl = sire_pl # Array of PL values for sire. Used for Koch DNp.
		@dam_pl = dam_pl # Array of PL values for dam. Used for Koch DNp.
		@offspring_pl = offspring_pl # Hash of arrays of PL values for offspring. Used for Koch DNp.
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
		type == 1 ? count = 2 : count = 1 # Add two mutations if a double-forward, one mutation is otherwise
		$total_denovo[offspr][type] += count
		for window in $current_windows
			window.denovo[offspr][type] += count
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
		parhom_filter = false # Do not classify as denovo if parents heterozyogous and options parhom
		if $options.parhom
			parhom_filter = true if sire_alleles.uniq.size != 1
			parhom_filter = true if dam_alleles.uniq.size != 1
		end
		poss_genotypes_unphased.push([sire_alleles[0],dam_alleles[0]].sort) # Resorting because could change during recombination
		poss_genotypes_unphased.push([sire_alleles[0],dam_alleles[1]].sort)
		poss_genotypes_unphased.push([sire_alleles[1],dam_alleles[0]].sort)
		poss_genotypes_unphased.push([sire_alleles[1],dam_alleles[1]].sort)
		poss_genotypes_unphased = poss_genotypes_unphased.uniq # If only one possible, uniq! will return nil
		for offspr in offspring.keys
			offspr_alleles, offspr_phased = split_alleles(offspring[offspr])
			offspr_phased && dam_phased && sire_phased ? genopool = check_phasing(poss_genotypes_unphased) : genopool = poss_genotypes_unphased # Determine whether to use phased or unphased genotypes	
			denovo_status = [false,""] # Default status is no de novo mutation
			unless parhom_filter # Do not classify as denovo if parents heterozyogous and options parhom
				unless genopool.include?(offspr_alleles)
					retain_snp = true # Remove SNPs that do not pass Koch DNp filter if option is used
					retain_snp = filter_candidate(offspring_pl, dam_pl, sire_pl, cutoff) if $options.kochDNp
					if retain_snp # Classify SNPs that passed Koch DNp filter (or all SNPs if Koch DNp filter not used)
						mutation_class = classify_denovo(sire_alleles, dam_alleles, offspr_alleles, offspr)
						denovo_status = [true,mutation_class]
					end
				end
			end
			@denovo[offspr] = denovo_status
		end
	end
end
#-----------------------------------------------------------------------------------------
class Bootstrap_Window
	# Class to define windows for bootstrapping
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
#-----------------------------------------------------------------------------------------
def depth_to_alleles(adepth, comparator) # Method to convert allele coverages/frequencies to alleles
	genotype = []
	for all in 0 ... adepth.size
		genotype.push(all.to_s) if adepth[all].to_f >= comparator
	end
	genotype.push(genotype[0]) if genotype.size == 1 # Make homozygotes
	return genotype.join("/")
end
#-----------------------------------------------------------------------------------------
def read_vcf # Method to read vcf
	collect_data = false # Flag to start collecting data
	next_window_site = 1 # Site to start new window
	next_window_close_site = next_window_site + $options.window - 1 # Site to close window and calculate rate
	gz_file_open($options.infile) do |f1|
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
				alleles.delete(".") if alleles.include?(".") # Ignore for non-polymorphic sites
				if alleles.size > 1 # Only process sites with actual SNPs in them
					snp = Snp.new
					snp.alleles = alleles
					tags = snp_array[8].split(":") # Get indexes of AD and GT tags
					ad = gt = pl = nil # Reset from previous cycle
					ad = tags.index("AD") if tags.include?("AD")
					pl = tags.index("PL") if tags.include?("PL")
					gt = tags.index("GT")
					for i in 9...snp_array.size
						genotype = snp_array[i].split(":")[gt]
						if $options.minAD1
							if ad.nil?
								filter_exit("AD tag required for minAD1 filter. Exiting.\nError found here:", line)
							else
								adepth = snp_array[i].split(":")[ad].split(",") # Convert allele coverages to alleles
								genotype = depth_to_alleles(adepth, 1.0)
							end
						elsif !$options.minAF.nil?
							if ad.nil?
								filter_exit("AD tag required for minAF filter. Exiting.\nError found here:", line)
							else
								adepth = snp_array[i].split(":")[ad].split(",") # Convert allele coverages to alleles
								adepth.map! { |x| x.to_f }
								total_depth = adepth.sum
								adepth.map! { |x| x/total_depth } # Convert allele coverages to frequencies
								genotype = depth_to_alleles(adepth, $options.minAF)
							end
						end
						pl_array = []
						if $options.kochDNp
							if pl.nil?
								filter_exit("PL tag required for Koch_DNp filter. Exiting.\nError found here:", line)
							else
								pl_array = snp_array[i].split(":")[pl].split(',')
								filter_exit("Biallelic SNPs required for Koch_DNp filter. Exiting.\nError found here:", snp_record) if pl_array.size != 3
							end
						end
						if i == sire_index
							snp.sire = genotype
							snp.sire_pl = pl_array
						elsif i == dam_index
							snp.dam = genotype
							snp.dam_pl = pl_array
						else
							snp.offspring[$sample_index[i]] = genotype
							snp.offspring_pl[$sample_index[i]] = pl_array
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
def mean(values) # Calculate mean of an array of values
	return values.sum.to_f/values.size.to_f
end
#-----------------------------------------------------------------------------------------
def conf95(values) # Calculate 95% confidence interval for an array of values
	meanval = mean(values)
	sdnum = values.map { |x| (x - meanval) * (x + meanval) }
	stdev = Math.sqrt(sdnum.sum/(values.size.to_f - 1.0))
	critval = 1.96 * stdev/Math.sqrt(values.size.to_f)
	finalstring = meanval.to_s + "\t" + (meanval - critval).to_s + ".." + (meanval + critval).to_s
	return finalstring
end
#-----------------------------------------------------------------------------------------
def bootstrap_results # Calculate bootstrapped results using predefined bootstrap windows and print results
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
			bootdenom = 2 * current_bootstrap_length.to_f
			# Get sum of de novo mutations. Double-forward mutations (index 1) count as two mutations.
			boot_dnm_sum = (bootstrap_denovo[offspr][0] + 2 * bootstrap_denovo[offspr][1] + bootstrap_denovo[offspr][2]).to_f
			puts offspr + "\t" + (boot_dnm_sum/bootdenom).to_s + "\t" + (bootstrap_denovo[offspr][0].to_f/bootdenom).to_s
			$bootstraps[offspr][0].push(bootstrap_denovo[offspr].sum.to_f/bootdenom)
			$bootstraps[offspr][1].push(bootstrap_denovo[offspr][0].to_f/bootdenom)
		end
	end
	puts "\nBootstrapped Estimates:"
	puts "Offspring\tMean_AllSites\t95%C.I._Allsites\tMean_Single-ForwardOnly\t95%C.I._Single-ForwardOnly:"
	for offspr in $bootstraps.keys
		puts offspr + "\t" + conf95($bootstraps[offspr][0]) + "\t" + conf95($bootstraps[offspr][1])
	end
end
#-----------------------------------------------------------------------------------------
def print_options # Print options used at startup to output
	puts "calc_denovo_mutation_rate " + CALCDENOVOVER + " started with parameters:"
	cmdline = "-i " + $options.infile + " -s " + $options.sire + " -d " + $options.dam + " -w " + $options.window.to_s + " -S " + $options.step.to_s + " -l " + $options.minbslen.to_s + " -M " + $options.minwindows.to_s
	cmdline << " -k -m " + $options.mu.to_s + " -t " + $options.theta.to_s + " -c " + $options.cutoff.to_s if $options.kochDNp
	cmdline << " --parhom" if $options.parhom
	cmdline << " --minAD1" if $options.minAD1
	cmdline << " --minAF " + $options.minAF.to_s unless $options.minAF.nil?
	cmdline << " -g" if $options.gvcf
	cmdline << " --rng " + $options.rng.to_s
	puts cmdline
end
#-----------------------------------------------------------------------------------------
ARGV[0] ||= "-h" # If no parameters passed, print help screen
# Parse options and run script below
$options = Parser.parse(ARGV)
$options.minbslen = $options.window if $options.minbslen > $options.window # Prevent nonsensical results
srand($options.rng)
print_options
setup_kochDNp($options.mu, $options.theta) if $options.kochDNp # Precalculate probabilities if using DNp filter
read_vcf
print_results
setup_kochDNp($options.mu, $options.theta) if $options.kochDNp # Precalculate probabilities if using DNp filter
bootstrap_results if ($options.bootstrap > 0 && $windows.size >= $options.minwindows) # Bootstrap results if possible
puts $mutations # Print list of mutations extracted from VCF
