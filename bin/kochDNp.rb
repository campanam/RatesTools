#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# kochDNp
KOCHDNPVER = "0.1.1"
# Michael G. Campana and Ellie E. Armstrong, 2022
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

require 'bigdecimal'
require_relative 'denovolib'

def precalc_parental_heterozygosity(theta) # Build parental heterozygosity array and store for single calculation
	# Calculate probability of parental heterozygosity from theta using 1 = p^2 + 2pq + q^2 where q = 1 - p
	prob_both_het = BigDecimal((theta ** 2).to_s) # Probability both parents are heterozygous
	prob_one_het = BigDecimal(((1-theta) * theta).to_s) # Probability one of the parents is heterozygous (split between dam/sire)
	prob_neither_het = BigDecimal((1 - prob_both_het - 2 * prob_one_het).to_s) # Probability neither is heterozygous
	het_array = [prob_neither_het, prob_both_het, prob_one_het].map { |x| Math.log10(x) } # Return log-transformed 
	return het_array 
end
#-----------------------------------------------------------------------------------------
def precalc_transmission(mu) # Precalculate transmission probabilities based on Koch model

	# Attempting to work-out transition probability matrix used by Koch et al.
	# Modelling as the probability of a mutation occurred in one of 4 parental alleles:
	#     (1-mu)^4 = (1 - 2mu + mu^2)(1 - 2mu + mu^2) = 1 + -2mu + mu^2 + -2mu + 4mu^2 + -2mu^3 + mu^2 + - 2mu^3 + mu^4
    #               	= 1 + - 4mu + 6mu^2 + -4mu^3 + mu^4      	
    # Probability of no mutation occurring in 4 alleles roughly equals 1 - 4mu if mu^2 and higher exponents rounded to 0
    # Probability of a mutation occurring in 4 alleles therefore roughly equals 4mu  
    # Even if mutation occurred, probability of unmutated alleles being passed down is 0.5
    # Probability of non-mutated alleles being passed down is 2mu
    # Total probability of no mutation is therefore 1 - 4mu + 2mu = 1 - 2mu (possibility of no mutations + possibility of mutations not passed)
    # If a single mutation occured, the parent has 3 alleles, so there is 1/3 probability of passing the mutated allele
    # Probability when single mutation must occur is therefore: (1 - 4mu) * 0 + (4mu)*(0.5 * 0 + 0.33 * 0.5) = (2/3) * mu
	
	# Return all probabilities in log scale
	$p_nomutation = Math.log10(BigDecimal((1.0 - 2.0 * mu).to_s)) # Probability of no mutation at a site
	$p_double = -307 # Probability of required double mutation. Use lowest exponent for Ruby to prevent underflow which is default -307
	$p_singleforward_parents_hom = Math.log10(BigDecimal((2.0/3.0 * mu).to_s)) # Probability of a single forward mutation at a 2 parental homozygous site
	$p_singlebackward_doublebackward = Math.log10(BigDecimal((1.0/6.0 * mu).to_s)) # Probability of a site that could either be a single or double backward mutation
	$p_singlebackward_het_nomutation = Math.log10(BigDecimal((0.5 - (2.0/3.0) * mu).to_s)) # Probability of a site that could either be a single backward mutation to a heterozygous state or no mutation
	$p_singlebackward_hom_nomutation = Math.log10(BigDecimal((0.5 - (5.0/6.0) * mu).to_s)) # Probability of a site that could either be a single backward mutation to a homozygous state or no mutation
	$p_singlebackward = Math.log10(BigDecimal((1.0/3.0 * mu).to_s)) # Probability of a single backward mutation
	$p_singleforward_nomutation = Math.log10(BigDecimal((0.25 - (1.0/3.0) * mu).to_s)) # Probability of a single forward mutation or no mutation
end
#-----------------------------------------------------------------------------------------
def calculate_pp_array(offspring_pl, dam_pl, sire_pl)
	pp_array = []
	for i in 0 .. 2 # Offspring pl states
		for j in 0 .. 2 # Dam pl states
			for k in 0 .. 2 # Sire pl states
				pl_state = [offspring_pl[i],dam_pl[j],sire_pl[k]].map { |x| BigDecimal((x.to_f / -10.0).to_s) } # Convert PL to error probabilitiies
				trans_prob = nil
				firstterm = pl_state.reduce(:+) # Calculate genotype likelihoods of Koch DNp statistic in log scale
				# Get precalculated heterozygosity probabilities
				if (j == 0 && k == 0) || (j == 2 && k == 2) || (j == 0 && k == 2) || (j == 2 && k == 0) # Parents both homozygous REF or ALT
					pr_het = $het_array[0]
					if j == k # If parents both homozygous and same
						# If parents both Hom and same, probabilities of hom same offspring is 1 assuming no mutation
						# With mutation, probability is 1 - 2p(1-p) - p^2 = 1 - 2p + p^2
						if i == j # If offspring same as parents
							trans_prob = $p_nomutation		
						elsif i == 1 # If offspring heterozygous, case of homozygous parents and single-forward mutation
								trans_prob = $p_singleforward_parents_hom
						else # If offspring homozygous for double mutation
							trans_prob = $p_double
						end		
					else # If parents homozygous but different
						# If parents both Hom and different, probability of het offspring with mutation is 1 - 2p(1-p) - p^2 = 1 - 2p + p^2
						if i == 1 # if offspring heterozygous
							trans_prob = $p_nomutation		
						else # Probability of a single backward mutation
							trans_prob = $p_singlebackward
						end	
					end
				elsif j == 1 && k == 1 # Parents both heterozygous
					pr_het = $het_array[1]
					if i == 1 # if offspring heterozygous, probability of genotype without mutation is 0.5
						# Could either be single mutation to a heterozygous state, double mutation (p = 0 with approximation) or single mutation
						trans_prob = $p_singlebackward_het_nomutation
					else # Could be either single-forward or no mutaiton
						trans_prob = $p_singleforward_nomutation
					end
				else # Remaining cases are one parent heterozygous
					pr_het = $het_array[2]
					if (i != 1 && j != i && k == 1) || (i != 1 && k != i && j == 1) # Case can be either single- or double-backward
						trans_prob = $p_singlebackward_doublebackward 
					elsif i == 1 # Case can be either single backward to a heterozygous state or no mutation
						trans_prob = $p_singlebackward_het_nomutation
					else # Can be single backward to a homozygous state or no mutation
						trans_prob = $p_singlebackward_hom_nomutation
					end
				end
				pl_state = firstterm + pr_het + trans_prob
				pp_array.push(pl_state)
			end
		end
	end
	return pp_array
end
#-----------------------------------------------------------------------------------------
def calculate_DNp(pp_array)
	# Translate log-sum-exp trick into Ruby
	scalefactor = pp_array.max
	pp_array_scaled = pp_array.map { |x| 10 ** (x - scalefactor) }
	#print pp_array_scaled.sum
    logsum = scalefactor + Math.log10(pp_array_scaled.sum)
    dnp_array = pp_array.map { |x| 10 ** (x - logsum) }
	koch_DNp = dnp_array[9] + dnp_array[17]
	return koch_DNp
end
#-----------------------------------------------------------------------------------------
def get_pl_arrays(snp_record, offspring_index, dam_index, sire_index) # Method to get PL values directly from VCF records
	# Permits use of code in both post-hoc analysis and as a filter in calc_denovo_mutation_rate
	snp_array = snp_record.split
	tags = snp_array[8].split(":") # Get indexes of PL scores
	pl_index = nil
	pl_index = tags.index("PL")
	# Error handling for absence of PL tag.
	filter_exit("PL tag required for Koch_DNp filter. Exiting.\nError found here:", snp_record) if pl_index.nil? # Exit out of code if PL tag does not exist
	# Get individual PL scores and transform into arrays
	offspring_pl = snp_array[offspring_index].split(":")[pl_index].split(',')
	dam_pl = snp_array[dam_index].split(":")[pl_index].split(',')
	sire_pl = snp_array[sire_index].split(":")[pl_index].split(',')
	# Error handling for non-biallelic sites
	filter_exit("Biallelic SNPs required for Koch_DNp filter. Exiting.\nError found here:", snp_record) if offspring_pl.size != 3 || dam_pl.size != 3 or sire_pl.size != 3
	return offspring_pl, dam_pl, sire_pl
end
#-----------------------------------------------------------------------------------------
def setup_kochDNp(mu,theta) # Control method to execute once to calculate DNp global variables
	BigDecimal.limit(16)
	precalc_transmission(mu)
	$het_array = precalc_parental_heterozygosity(theta)
end
#-----------------------------------------------------------------------------------------
def filter_candidate(offspring_pl, dam_pl, sire_pl, cutoff)
	# Method to calculate and filter by DNp during calc_denovo_mutation_rate using pre-extracted values
	pp_array = calculate_pp_array(offspring_pl, dam_pl, sire_pl)
	koch_DNp = calculate_DNp(pp_array)
	koch_DNp >= cutoff ? koch_filter = true : koch_filter = false
	return koch_filter
end
#-----------------------------------------------------------------------------------------
def posthoc_filter # Post-hoc filter previously generated de novo mutation candidates
	# Precalculate the various probabilities
	setup_kochDNp($options.mu, $options.theta)
	# Input is a RatesTools final DNM log
	collect_data = false # Don't do anything until header line found
	gz_file_open($options.infile) do |f1|
		while line = f1.gets
			if line[0..5] == "#CHROM"
				print line
				collect_data = true
				line_arr = line.split
				# Get indexes of sire, dam and offspring
				sire_index = line_arr.index($options.sire)
				dam_index = line_arr.index($options.dam)
				offspring_indexes = [] # Could be multiple offspring
				for i in 9...line_arr.size
					offspring_indexes.push(i) if i != sire_index and i != dam_index
				end
			elsif collect_data
				line_arr = line.split
				dnm_state = false # Variable controlling whether site is considered a likely DNM in at least one offspring
				for offspring_index in offspring_indexes
					offspring_pl, dam_pl, sire_pl = get_pl_arrays(line, offspring_index, dam_index, sire_index)
					offspring_dnm_state = filter_candidate(offspring_pl, dam_pl, sire_pl, $options.cutoff) # Is this a DNM for this offspring
					dnm_state = true if offspring_dnm_state # If at least one offspring has DNM, dnm_state overall is true
				end
				print line if dnm_state
			end
		end
	end
end
#-----------------------------------------------------------------------------------------
unless ARGV.include?("-k") or ARGV.include?("--kochDNp") # Run post-hoc code only when specified, not when incorporated into calc_denovo_mutation_rate
	$options = Parser.parse(ARGV, true)
	posthoc_filter
end
