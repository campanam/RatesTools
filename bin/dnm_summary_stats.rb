#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# dnm_summary_stats
DNMSUMSTATSVER = "1.1.0"
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
def classify_sites(outindiv)
	# Single-Forward Mutations
	mutclasses = { "A->T" => 0, "A->C" => 0, "A->G" => 0, "T->A" => 0, "T->C" => 0,
					"T->G" => 0, "C->T" => 0, "C->G" => 0, "C->A" => 0, "G->T" => 0,
					"G->A" => 0, "G->C" => 0 }
	# Double-forward mutations
	dfclasses = { "A->T" => 0, "A->C" => 0, "A->G" => 0, "T->A" => 0, "T->C" => 0,
					"T->G" => 0, "C->T" => 0, "C->G" => 0, "C->A" => 0, "G->T" => 0,
					"G->A" => 0, "G->C" => 0 }
	# Backward mutations
	backclasses = { "A->T" => 0, "A->C" => 0, "A->G" => 0, "T->A" => 0, "T->C" => 0,
					"T->G" => 0, "C->T" => 0, "C->G" => 0, "C->A" => 0, "G->T" => 0,
					"G->A" => 0, "G->C" => 0 }
	otherscnt = 0
	start = false # Flag to collect data
	getmutationcnts = false # Flag to collect mutation count data
	File.open(ARGV[1] + "_offspring" + outindiv + "_summary.log") do |f3|
		while line = f3.gets
			if start
				snp_array = line[0..-2].split("\t")
				snpsite = snp_array[0] + ':' + snp_array[1]
				$candidates[snpsite].nil? ? $candidates[snpsite] = [1,[]] :  $candidates[snpsite][0] += 1 # Structure is array of [total_number, [individuals for single-forward mutation]]
				alleles = ([snp_array[3]] + snp_array[4].split(",")).flatten.uniq # Get alleles
				alleles.delete("<non_ref>") if alleles.include?("<non_ref>")
				alleles.delete(".") if alleles.include?(".") # Ignore for non-polymorphic sites
				tags = snp_array[8].split(":") # Get indexes of GT tags
				gt = tags.index("GT")
				par_genotypes = [] # Array of parental genotypes
				for i in 9...snp_array.size # There should always be 3 samples if using RatesTools pipeline. Will break if using calc_denovo_mutation_rate within more individuals.
					genotype = snp_array[i].split(":")[gt].gsub("|", "/").split("/").uniq.map { |x| x.to_i }
					if i == @off_index
						off_genotype = genotype
					else
						par_genotypes.push(genotype)
					end
				end
				# Basic handling assuming biallelic SNPs, single-forward, parental homozygous. Dumps all other types into "other".
				if par_genotypes[0] == [0] && par_genotypes[1] == [0] && off_genotype == [0,1]
					mutclass = "#{alleles[0]}->#{alleles[1]}" # Dumps all other mutations into other
					if mutclasses[mutclass].nil? 
						otherscnt += 1
					else
						mutclasses[mutclass] += 1
						$candidates[snpsite][1].push(outindiv) # Add to single-forward array
					end
				elsif par_genotypes[0] == [1] && par_genotypes[1] == [1] && off_genotype == [0,1]
					mutclass = "#{alleles[1]}->#{alleles[0]}"
					if mutclasses[mutclass].nil?  # Dumps all other mutations into other
						otherscnt += 1
					else
						mutclasses[mutclass] += 1
						$candidates[snpsite][1].push(outindiv) # Add to single-forward array
					end
				elsif par_genotypes[0] == [0] && par_genotypes[1] == [0] && off_genotype == [1]
					mutclass = "#{alleles[0]}->#{alleles[1]}"
					dfclasses[mutclass].nil? ? otherscnt += 1 : dfclasses[mutclass] += 1 # Dumps all other mutations into other
				elsif par_genotypes[0] == [1] && par_genotypes[1] == [1] && off_genotype == [0]
					mutclass = "#{alleles[1]}->#{alleles[0]}"
					dfclasses[mutclass].nil? ? otherscnt += 1 : dfclasses[mutclass] += 1 # Dumps all other mutations into other
				elsif par_genotypes[0] == [0,1] && par_genotypes[1] && off_genotype == [0]
					mutclass = "#{alleles[1]}->#{alleles[0]}"
					backclasses[mutclass].nil? ? otherscnt += 1 : backclasses[mutclass] += 1 # Dumps all other mutations into other
				elsif par_genotypes[0] == [0] && par_genotypes[0,1] && off_genotype == [1]
					mutclass = "#{alleles[0]}->#{alleles[1]}"
					backclasses[mutclass].nil? ? otherscnt += 1 : backclasses[mutclass] += 1 # Dumps all other mutations into other
				elsif par_genotypes[0] == [1] && par_genotypes[1] == [0] && off_genotype == [1]
					mutclass = "#{alleles[0]}->#{alleles[1]}"
					backclasses[mutclass].nil? ? otherscnt += 1 : backclasses[mutclass] += 1 # Dumps all other mutations into other
				elsif par_genotypes[0] == [1] && par_genotypes[1] == [0] && off_genotype == [0]
					mutclass = "#{alleles[1]}->#{alleles[0]}"
					backclasses[mutclass].nil? ? otherscnt += 1 : backclasses[mutclass] += 1 # Dumps all other mutations into other
				elsif par_genotypes[0] == [0] && par_genotypes[1] == [1] && off_genotype == [1]
					mutclass = "#{alleles[0]}->#{alleles[1]}"
					backclasses[mutclass].nil? ? otherscnt += 1 : backclasses[mutclass] += 1 # Dumps all other mutations into other
				elsif par_genotypes[0] == [0] && par_genotypes[1] == [1] && off_genotype == [0]
					mutclass = "#{alleles[1]}->#{alleles[0]}"
					backclasses[mutclass].nil? ? otherscnt += 1 : backclasses[mutclass] += 1 # Dumps all other mutations into other
				end
			elsif line[0..30] == "Total number of retained sites:"
				$totalbases[outindiv] = [line.split(":")[1].to_i,0,0] # Total callable bases, total mutations, single-forward mutations
			elsif line == "Offspring\tSingle-Forward\tDouble-Forward\tBackward\n"
				getmutationcnts = true
			elsif getmutationcnts
				cnts = line.strip.split.map { |x| x.to_i }
				getmutationcnts = false
				$totalbases[outindiv][1] = cnts[1..3].sum
				$totalbases[outindiv][2] = cnts[1]
			elsif line[0..5] == "#CHROM"
				header_arr = line[0..-2].split("\t")
				@off_index = header_arr.index(outindiv)
				start = true
			end
		end	
	end
	puts "\n" + outindiv + " Mutation Classes\nSingle-Forward,Count"
	for mut in mutclasses.keys
		puts mut + "," + mutclasses[mut].to_s
	end
	puts "\nDouble-Forward,Count"
	for mut in dfclasses.keys
		puts mut + "," + dfclasses[mut].to_s
	end
	puts "\nBackward,Count"
	for mut in backclasses.keys
		puts mut + "," + backclasses[mut].to_s
	end
	puts "\nIndels/Other,Count"
	puts "Total: " + otherscnt.to_s
end
#----------------------------------------------------------------------------------------
if ARGV[0].nil?
	# If no parameters passed, print help screen
	format_splash('dnm_summary_stats', DNMSUMSTATSVER, '<logs_directory> <output_prefix> <DNM_clump_range> > <out.csv>')
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
	$dnmclump = ARGV[2] # Number of bases to search for clumped DNM candidates. Default of 0.
	
	# Get individual names
	Dir.foreach(ARGV[0] + "/") do |f1|
		if f1[-12..-1] == "_summary.log"
			$individuals.push(f1[0..-13])
		end
	end
	# Get allsites and chrsites if chrfilt.log exists.
	if File.exist?(ARGV[0] + "/chrfilt.log")
		File.open(ARGV[0] + "/chrfilt.log") do |f1|
			while line = f1.gets
				line_arr = line.split
				$chrsites = line_arr[3]
				$allsites = line_arr[-2]
			end
		end
	end
	# Get trio sites and fill in allsites/chrsites if chrfilt.log does not exist
	Dir.foreach(ARGV[0] + "/") do |f2|
		if f2[-9..-1] == '_trio.log'
			@individual = f2[0..-10]
			File.open(ARGV[0] + "/" + f2) do |f3|
				while line = f3.gets
					line_arr = line.split
					if $allsites.nil?
						$allsites = line_arr[-2]
						$chrsites = line_arr[-2]
					end
					$triosites[@individual] = line_arr[3]
				end
			end
		end
	end
	puts "Individual,AllSites,ChrSites,TrioSites,VCFtoolsFiltSites,GATKFiltSites,RegionFiltSites"
	outindivs = []
	for indiv in $individuals
		extract_site_count($vcf_filtsites, '.sitefilt.log', indiv)
		extract_site_count($gatk_filtsites, '.gatksitefilt.log', indiv)
		if Dir.glob(ARGV[0]+ "/*regionfilt.log").any? # Handling for when region filters are turned off
			extract_site_count($regionsites, '.regionfilt.log', indiv)
		else
			$regionsites = $gatk_filtsites
		end
		outindiv = indiv.gsub(ARGV[1] + "_offspring", "")
		outindivs.push(outindiv)
		puts outindiv + "," + $allsites + "," + $chrsites + "," + $triosites[indiv].to_s + "," + $vcf_filtsites[indiv].to_s + "," + $gatk_filtsites[indiv].to_s + "," + $regionsites[indiv].to_s
		
	end
	for outindiv in outindivs
		classify_sites(outindiv)
	end
	total_overlap = 0 # Both overlapping individuals and overlapping sites
	$sfindv = {} # Hash of single-forward counts per individual
	if $dnmclump > 0 # Don't bother if no clumping
		puts "\nClumped Candidate DNM Sites"
		sorted_candidates = {} # Hash to put in sites keyed by chr
		for key in $candidates.keys
			key_chr = key.split(":")[0]
			key_snp = key.split(":")[1].to_i
			if sorted_candidates[key].nil?
				sorted_candidates[key] = [key_snp]
			else
				sorted_candidates[key].push(key_snp)
			end
		end
		for key in sorted_candidates.key
			sorted_sites = sorted_candidates[key].sort
			if sorted_sites.size < 2 # Ignore case when clumped sites impossible
				prev_site = sorted_sites[0]
				for i in 1 ... sorted_sites.size
					if sorted_sites[i] <= prev_site
						removed_site = key + ":" + prev_site.to_s
						total_overlap += 1
						prev_site = sorted_sites[i]
						puts removed_site
						for sfindv in $candidates[removed_site][1] # Code to count number of single-forward mutations
							$sfindv[sfindv].nil? ? $sfindv[sfindv] = 1 : $sfindv[sfindv] +=1
						end
						$candidates.delete(removed_site)
					end
				end
				if sorted_sites[-1] == prev_site # Check the last site for being in a clump
					removed_site = key + ":" + prev_site.to_s
					total_overlap += 1
					prev_site = sorted_sites[i]
					puts removed_site
					for sfindv in $candidates[removed_site][1] # Code to count number of single-forward mutations
						$sfindv[sfindv].nil? ? $sfindv[sfindv] = 1 : $sfindv[sfindv] +=1
					end
					$candidates.delete(removed_site)
				end
			end
		end
	end
	puts "\nRemaining Candidate Sites Overlapping Between Offspring"
	for key in $candidates.keys
		if $candidates[key][0] > 1
			puts key
			total_overlap +=1
			for sfindv in $candidates[key][1] # Code to count number of single-forward mutations
				$sfindv[sfindv].nil? ? $sfindv[sfindv] = 1 : $sfindv[sfindv] +=1
			end
		end
	end
	puts "\nOffspring,Single-ForwardOverlappingSites,TotalOverlappingSites,RecalcSingle-ForwardRate,RecalcAllsitesRate"
	for key in $sfindv.keys
		srate = ($totalbases[key][2]-$sfindv[key]).to_f/$totalbases[key][0].to_f/2.to_f # recalculate single-forward rate
		arate = ($totalbases[key][1]-total_overlap).to_f/$totalbases[key][0].to_f/2.to_f # recalculate all mutation rate rate
		puts key + "," + $sfindv[key].to_s + "," + total_overlap.to_s + "," + srate.to_s + "," + arate.to_s
	end
end
