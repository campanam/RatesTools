#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# dnm_summary_stats
DNMSUMSTATSVER = "1.2.0"
# Michael G. Campana and Ellie E. Armstrong, 2022-2024
# Smithsonian Institution and Stanford University

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to RatesTools;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.
 
# We politely request that this work be cited as:
# Armstrong, E.E. & M.G. Campana. 2023. RatesTools: a Nextflow pipeline for detecting
# de novo germline mutations in pedigree sequence data. Bioinformatics. 39: btac784.
# 10.1093/bioinformatics/btac784.
#----------------------------------------------------------------------------------------

# Script to summarize output from ratestools pipeline

require_relative 'denovolib'

def extract_site_count(count, logpattern, indiv)
	Dir.foreach(ARGV[0] + "/") do |f1|
		val = logpattern.length * -1
		if f1[val..-1] == logpattern && f1[0...indiv.length] == indiv
			File.open(ARGV[0] + "/" + f1) do |f2|
				while line = f2.gets
					line_arr = line.split
					if count[indiv].nil?
						count[indiv] = line_arr[3].to_i
					else
						count[indiv] += line_arr[3].to_i
					end
				end
			end
		end
	end
end
#----------------------------------------------------------------------------------------
def update_counts(class_hash,mutclass,snpsite,outindiv) # Update number and individuals for types of mutations
	case class_hash
	when 'sf'
		if $mutclasses[mutclass].nil? 
			$otherscnt += 1
		else
			$mutclasses[mutclass] += 1
			$candidates[snpsite][0].push([outindiv,mutclass]) # Add to single-forward array
		end
	when 'df'
		if $dfclasses[mutclass].nil? # Dumps all other mutations into other.
			$otherscnt += 2 
		else 
			$dfclasses[mutclass] += 2 
			$candidates[snpsite][1].push([outindiv,mutclass]) # Add to double-forward array
		end
	when 'bk'
		if $backclasses[mutclass].nil? # Dumps all other mutations into other
			$otherscnt += 1
		else
			$backclasses[mutclass] += 1
			$candidates[snpsite][2].push([outindiv,mutclass]) # Add to backward array
		end
	end
end
#----------------------------------------------------------------------------------------
def print_spectrum(outindiv) # Print mutational spectra
	puts "\n" + outindiv + " Mutation Classes"
	puts "Single-Forward,Count"
	for mut in $spectra[outindiv][0].keys
		puts mut + "," + $spectra[outindiv][0][mut].to_s
	end
	puts "\nDouble-Forward,Count"
	for mut in $spectra[outindiv][1].keys
		puts mut + "," + $spectra[outindiv][1][mut].to_s
	end
	puts "\nBackward,Count"
	for mut in $spectra[outindiv][2].keys
		puts mut + "," + $spectra[outindiv][2][mut].to_s
	end
end
#----------------------------------------------------------------------------------------
def classify_sites(outindiv)
	# Single-Forward Mutations
	$mutclasses = { "A->T" => 0, "A->C" => 0, "A->G" => 0, "T->A" => 0, "T->C" => 0,
					"T->G" => 0, "C->T" => 0, "C->G" => 0, "C->A" => 0, "G->T" => 0,
					"G->A" => 0, "G->C" => 0 }
	# Double-forward mutations
	$dfclasses = { "A->T" => 0, "A->C" => 0, "A->G" => 0, "T->A" => 0, "T->C" => 0,
					"T->G" => 0, "C->T" => 0, "C->G" => 0, "C->A" => 0, "G->T" => 0,
					"G->A" => 0, "G->C" => 0 }
	# Backward mutations
	$backclasses = { "A->T" => 0, "A->C" => 0, "A->G" => 0, "T->A" => 0, "T->C" => 0,
					"T->G" => 0, "C->T" => 0, "C->G" => 0, "C->A" => 0, "G->T" => 0,
					"G->A" => 0, "G->C" => 0 }
	$otherscnt = 0
	start = false # Flag to collect data
	getmutationcnts = false # Flag to collect mutation count data
	getbootcnts = false # Flag to collect bootstrap count data
	$bootstrapped_data = false # Flag that data is bootstrapped
	File.open(ARGV[0] + '/' + ARGV[1] + "_offspring" + outindiv + "_summary.log") do |f3|
		while line = f3.gets
			if start
				snp_array = line[0..-2].split("\t")
				snpsite = snp_array[0] + ':' + snp_array[1]
				if $candidates[snpsite].nil? # Structure is array of [[individuals for single-forward mutation],[double-forward individuals],[back mutations individuals], number of individuals sharing this site]
					$candidates[snpsite] = [[],[],[],1]
				else
					$candidates[snpsite][3]  += 1
				end
				alleles = ([snp_array[3]] + snp_array[4].split(",")).flatten.uniq # Get alleles
				alleles.delete("<non_ref>") if alleles.include?("<non_ref>")
				alleles.delete(".") if alleles.include?(".") # Ignore for non-polymorphic sites
				tags = snp_array[8].split(":") # Get indexes of GT tags
				gt = tags.index("DNM_GT")
				par_genotypes = [] # Array of parental genotypes
				for i in 9...snp_array.size # There should always be 3 samples if using RatesTools pipeline. Will break if using calc_denovo_mutation_rate within more individuals.
					genotype = snp_array[i].split(":")[gt].gsub("|", "/").split("/").uniq.map { |x| x.to_i }
					if i == @off_index
						off_genotype = genotype
					else
						par_genotypes.push(genotype)
					end
				end
				# Handling for SNP class mutational spectra. Dumps other types of mutations into 'other'
				if par_genotypes[0] == [0] && par_genotypes[1] == [0] && off_genotype == [0,1]
					mutclass = "#{alleles[0]}->#{alleles[1]}" # Dumps all other mutations into other
					update_counts('sf',mutclass,snpsite,outindiv)
				elsif par_genotypes[0] == [1] && par_genotypes[1] == [1] && off_genotype == [0,1]
					mutclass = "#{alleles[1]}->#{alleles[0]}"
					update_counts('sf',mutclass,snpsite,outindiv)
				elsif par_genotypes[0] == [0] && par_genotypes[1] == [0] && off_genotype == [1]
					mutclass = "#{alleles[0]}->#{alleles[1]}"
					update_counts('df',mutclass,snpsite,outindiv)
				elsif par_genotypes[0] == [1] && par_genotypes[1] == [1] && off_genotype == [0]
					mutclass = "#{alleles[1]}->#{alleles[0]}"
					update_counts('df',mutclass,snpsite,outindiv)
				elsif par_genotypes[0] == [0,1] && par_genotypes[1] && off_genotype == [0]
					mutclass = "#{alleles[1]}->#{alleles[0]}"
					update_counts('bk',mutclass,snpsite,outindiv)
				elsif par_genotypes[0] == [0] && par_genotypes[0,1] && off_genotype == [1]
					mutclass = "#{alleles[0]}->#{alleles[1]}"
					update_counts('bk',mutclass,snpsite,outindiv)
				elsif par_genotypes[0] == [1] && par_genotypes[1] == [0] && off_genotype == [1]
					mutclass = "#{alleles[0]}->#{alleles[1]}"
					update_counts('bk',mutclass,snpsite,outindiv)
				elsif par_genotypes[0] == [1] && par_genotypes[1] == [0] && off_genotype == [0]
					mutclass = "#{alleles[1]}->#{alleles[0]}"
					update_counts('bk',mutclass,snpsite,outindiv)
				elsif par_genotypes[0] == [0] && par_genotypes[1] == [1] && off_genotype == [1]
					mutclass = "#{alleles[0]}->#{alleles[1]}"
					update_counts('bk',mutclass,snpsite,outindiv)
				elsif par_genotypes[0] == [0] && par_genotypes[1] == [1] && off_genotype == [0]
					mutclass = "#{alleles[1]}->#{alleles[0]}"
					update_counts('bk',mutclass,snpsite,outindiv)
				end
			elsif line[0..30] == "Total number of retained sites:"
				$totalbases[outindiv] = [line.split(":")[1].to_i,0,0] # Total callable bases, total mutations, single-forward mutations
			elsif line == "Offspring\tSingle-Forward\tDouble-Forward\tBackward\n"
				getmutationcnts = true
				$bootstrapped_data = true
			elsif getmutationcnts
				cnts = line.strip.split.map { |x| x.to_i }
				getmutationcnts = false
				$totalbases[outindiv][1] = cnts[1..3].sum
				$totalbases[outindiv][2] = cnts[1]
			elsif line == "Offspring\tMean_AllSites\t95%C.I._Allsites\tMean_Single-ForwardOnly\t95%C.I._Single-ForwardOnly:\n"
				getbootcnts = true
			elsif getbootcnts
				getbootcnts = false
				meanallsites = line.split[1].to_f
				stderrallsites = (line.split[2].split('..')[1].to_f - meanallsites)/1.96
				meansingleforward = line.split[3].to_f
				stderrsingleforward = (line.split[4].split('..')[1].to_f - meansingleforward)/1.96
				$bootstrapmeans[outindiv] = [meanallsites,meansingleforward,stderrallsites,stderrsingleforward]
			elsif line[0..5] == "#CHROM"
				header_arr = line[0..-2].split("\t")
				@off_index = header_arr.index(outindiv)
				start = true
			end
		end	
	end
	$spectra[outindiv] = [$mutclasses.dup, $dfclasses.dup, $backclasses.dup]
	print_spectrum(outindiv)
	puts "\nIndels/Other,Count"
	puts "Total: " + $otherscnt.to_s
end
#----------------------------------------------------------------------------------------
if ARGV[0].nil?
	# If no parameters passed, print help screen
	format_splash('dnm_summary_stats', DNMSUMSTATSVER, '<logs_directory> <output_prefix> <DNM_clump_range> > <out.csv> 2> <retained sites.tsv>')
else
	$individuals = [] # Array of individual names
	$allsites = nil # Total number of sites before filtration
	$chrsites = nil # Number of sites on chromosomes
	$triosites = {} # Hash of site counts in each trio
	$vcf_filtsites = {} # Hash of site counts after VCFtools site filtering
	$gatk_filtsites = {} # Hash of site counts after GATK site filtering
	$regionsites = {} # Hash of site counts after region filtering
	$candidates = {} # Hash of mutation site counts to identify sibling overlaps
	$totalbases = {} # Hash indexing total site counts to individuals
	$bootstrapmeans = {} # Hash indexing boostrap mean values by individuals
	$dnmclump = ARGV[2].to_i # Number of bases to search for clumped DNM candidates. Default of 0.
	$spectra = {} # Hash of mutational spectra indexed by individuals
	
	# Get individual names
	Dir.foreach(ARGV[0] + "/") do |f1|
		if f1[-12..-1] == "_summary.log"
			$individuals.push(f1[0..-13])
		end
	end
	outindivs = []
	for indiv in $individuals
		outindiv = indiv.gsub(ARGV[1] + "_offspring", "")
		outindivs.push(outindiv)
	end
	$total_removed = {} # Hash of total overlapping individuals and clumped sites per individual
	$sfindv = {} # Hash of single-forward counts per individual
	print "\nRaw Mutational Spectra:"
	for outindiv in outindivs
		classify_sites(outindiv)
		$total_removed[outindiv] = 0 # Initialize total removed site count to zero
		$sfindv[outindiv] = 0 # Initialize removed single-forward site count to zero
	end		
	if $dnmclump > 0 # Don't bother if no clumping
		puts "\nClumped Candidate DNM Sites"
		sorted_candidates = {} # Hash to put in sites keyed by chr
		for key in $candidates.keys
			key_chr = key.split(":")[0]
			key_snp = key.split(":")[1].to_i
			if sorted_candidates[key_chr].nil?
				sorted_candidates[key_chr] = [key_snp]
			else
				sorted_candidates[key_chr].push(key_snp)
			end
		end
		for key in sorted_candidates.keys
			sorted_sites = sorted_candidates[key].sort
			prev_site = sorted_sites[0]
			if sorted_sites.size > 1 # Ignore case when clumped sites impossible
				for i in 1 ... sorted_sites.size
					if sorted_sites[i] <= prev_site + $dnmclump
						removed_site = key + ":" + prev_site.to_s
						puts removed_site
						for sfindv in $candidates[removed_site][0] # Code to count number of single-forward mutations
							$sfindv[sfindv[0]] +=1
							$total_removed[sfindv[0]] += 1
							$spectra[sfindv[0]][0][sfindv[1]] -= 1 # Update mutational spectra
						end
						for tindv in $candidates[removed_site][1] # Code to count total number of double-forward removed sites
							$total_removed[tindv[0]] += 2
							$spectra[tindv[0]][1][tindv[1]] -= 2 # Update mutational spectra
						end
						for tindv in $candidates[removed_site][2] # Code to count total number of backward removed sites
							$total_removed[tindv[0]] += 1
							$spectra[tindv[0]][2][tindv[1]] -= 1 # Update mutational spectra
						end
						$candidates.delete(removed_site)
						if i == sorted_sites.size - 1 # Add last removed site if goes to end
							removed_site = key + ":" + sorted_sites[i].to_s
							puts removed_site
							for sfindv in $candidates[removed_site][0]
								$sfindv[sfindv[0]] += 1
								$total_removed[sfindv[0]] += 1
								$spectra[sfindv[0]][0][sfindv[1]] -= 1 # Update mutational spectra
							end
							for tindv in $candidates[removed_site][1]
								$total_removed[tindv[0]] += 2
								$spectra[tindv[0]][1][tindv[1]] -= 2 # Update mutational spectra
							end
							for tindv in $candidates[removed_site][2]
								$total_removed[tindv[0]] += 1
								$spectra[tindv[0]][2][tindv[1]] -= 1 # Update mutational spectra
							end
							$candidates.delete(removed_site)
						end
					elsif (sorted_sites[i-1] - sorted_sites[i-2]).abs <= $dnmclump # Handling for numerical gap
						removed_site = key + ":" + prev_site.to_s
						puts removed_site
						for sfindv in $candidates[removed_site][0] # Code to count number of single-forward mutations
							$sfindv[sfindv[0]] +=1
							$total_removed[sfindv[0]] += 1
							$spectra[sfindv[0]][0][sfindv[1]] -= 1 # Update mutational spectra
						end
						for tindv in $candidates[removed_site][1] # Code to count total number of double-forward removed sites
							$total_removed[tindv[0]] += 2
							$spectra[tindv[0]][1][tindv[1]] -= 2 # Update mutational spectra
						end
						for tindv in $candidates[removed_site][2] # Code to count total number of backward removed sites
							$total_removed[tindv[0]] += 1
							$spectra[tindv[0]][2][tindv[1]] -= 1 # Update mutational spectra
						end
						$candidates.delete(removed_site)
					end
					prev_site = sorted_sites[i]
				end
			end
		end
	end
	puts "\nRemaining Candidate Sites Overlapping Between Offspring"
	for key in $candidates.keys
		if $candidates[key][3] > 1
			puts key
			for sfindv in $candidates[key][0] # Code to count number of single-forward mutations
				$sfindv[sfindv[0]] +=1
				$total_removed[sfindv[0]] += 1
				$spectra[sfindv[0]][0][sfindv[1]] -= 1 # Update mutational spectra
			end
			for tindv in $candidates[key][1] # Code to count number of double-forward mutations
				$total_removed[tindv[0]] +=2
				$spectra[tindv[0]][1][tindv[1]] # Update mutational spectra
			end
			for tindv in $candidates[key][2] # Code to count number of backward mutations
				$total_removed[tindv[0]] +=1
				$spectra[tindv[0]][2][tindv[1]] -= 1 # Update mutational spectra
			end
			$candidates.delete(key)
		end
	end
	puts "\nSurviving Candidate DNM Sites (Across Individuals)"
	for key in $candidates.keys
		puts key
		$stderr.puts key.gsub(':',"\t") # Creates a list of retained SNPs for VCFtools site filtering
	end
	print "\nFiltered Mutational Spectra:"
	for key in $spectra.keys
		print_spectrum(key)
	end
	puts "\nOffspring,Single-ForwardRemovedSites,TotalRemovedSites,RemainingSingle-ForwardSites,RemainingTotalSites,RecalcSingle-ForwardRate,Single-Forward95%BinomialConfidence,RecalcAllsitesRate,AllSites95%BinomialConfidence"
	for key in $sfindv.keys
		scount = $totalbases[key][2]-$sfindv[key] # Number of single-forward sites
		srate = scount.to_f/$totalbases[key][0].to_f/2.to_f # recalculate single-forward rate
		sconf = `Rscript -e 'library(Hmisc, quietly = TRUE);binconf(x = #{scount}, n = #{$totalbases[key][0] * 2}, alpha = 0.05)'` # Get single-forward binomial confidence interval
		sconf2 = sconf.split("\n")[-1].split[1..2].join('...')
		acount = $totalbases[key][1]-$total_removed[key] # Number of all-sites mutations
		arate = acount.to_f/$totalbases[key][0].to_f/2.to_f # recalculate all mutation rate rate
		aconf = `Rscript -e 'library(Hmisc, quietly = TRUE);binconf(x = #{acount}, n = #{$totalbases[key][0] * 2}, alpha = 0.05)'` # Get all-sites binomial confidence interval
		aconf2 = aconf.split("\n")[-1].split[1..2].join('...')
		puts key + "," + $sfindv[key].to_s + "," + $total_removed[key].to_s + "," + ($totalbases[key][2]-$sfindv[key]).to_s + "," + ($totalbases[key][1]-$total_removed[key]).to_s + "," + srate.to_s + "," + sconf2 + "," + arate.to_s + "," + aconf2
	end
end

if $bootstrapped_data # Do not bother if no bootstrapped data
	puts "\nOffspring,Single-ForwardCorrectedMean,Single-ForwardCorrectedSE,Single-ForwardCorrected95%C.I.,AllSitesCorrectedMean,AllSitesCorrectedSE,AllSitesCorrected95%C.I."
	for key in $bootstrapmeans.keys
		sf_correction_factor = ($totalbases[key][2]-$sfindv[key]).to_f/$totalbases[key][2].to_f
		sf_corrected_mean = sf_correction_factor * $bootstrapmeans[key][1]
		sf_corrected_se = $bootstrapmeans[key][3] * sf_correction_factor
		sf_corrected_crit = 1.96 * sf_corrected_se
		sf_corrected_ci = (sf_corrected_mean - sf_corrected_crit).to_s + '...' + (sf_corrected_mean + sf_corrected_crit).to_s
		all_correction_factor = ($totalbases[key][1]-$total_removed[key]).to_f/$totalbases[key][1].to_f
		all_corrected_mean = all_correction_factor * $bootstrapmeans[key][0]
		all_corrected_se = $bootstrapmeans[key][2] * all_correction_factor
		all_corrected_crit = 1.96 * all_corrected_se
		all_corrected_ci = (all_corrected_mean - all_corrected_crit).to_s + '...' + (all_corrected_mean + all_corrected_crit).to_s
		puts key + "," + sf_corrected_mean.to_s + "," + sf_corrected_se.to_s + "," + sf_corrected_ci + "," + all_corrected_mean.to_s + ',' + all_corrected_se.to_s + ',' + all_corrected_ci
	end
end