#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# dnm_summary_stats
DNMSUMSTATSVER = "0.1.0"
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

if ARGV[0].nil?
	# If no parameters passed, print help screen
	format_splash('filterGM', DNMSUMSTATSVER, '<logs_directory> > <out.csv>')
else
	$individuals = [] # Array of individual names
	$allsites = nil # Total number of sites before filtration
	$chrsites = nil # Number of sites on chromosomes
	$triosites = {} # Hash of site counts in each trio
	$vcf_filtsites = {} # Hash of site counts after VCFtools site filtering
	$gatk_filtsites = {} # Hash of site counts after GATK site filtering
	$regionsites = {} # Hash of site counts after region filtering
	
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
	for indiv in $individuals
		extract_site_count($vcf_filtsites, '_sitefilt.log', indiv)
		extract_site_count($gatk_filtsites, '_gatksitefilt.log', indiv)
		extract_site_count($regionsites, '_regionfilt.log', indiv)
		puts indiv + "," + $allsites + "," + $chrsites + "," + $triosites[indiv].to_s + "," + $vcf_filtsites[indiv].to_s + "," + $gatk_filtsites[indiv].to_s + "," + $regionsites[indiv].to_s
	end
end
