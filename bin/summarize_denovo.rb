#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# summarize_denovo
SUMMARIZEDENOVOVER = "0.5.6"
# Michael G. Campana and Ellie E. Armstrong, 2020-2024
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

# Script to calculate genome-wide denovo mutation rates from a directory of previously split DNM logs

require_relative 'denovolib'

def mean_ci(mean, crit) # Calculate mean and confidence interval for bootstrapped results and return formatted output line
	bootdenom = $bootstrap_bp.to_f # bootdenom does NOT need to be multiplied by 2 (already done during previous bootstrapped calculation)
	return (mean/bootdenom).to_s + "\t" + ((mean-crit)/bootdenom).to_s + ".." + ((mean+crit)/bootdenom).to_s
end
#----------------------------------------------------------------------------------------
def print_summary # Print results from unbootstrapped and then bootstrapped results
	print_results
	puts "\nBootstrapped Estimates:"
	puts "Bootstrapped retained sites: " + $bootstrap_bp.to_s
	puts "Offspring\tMean_AllSites\t95%C.I._Allsites\tMean_Single-ForwardOnly\t95%C.I._Single-ForwardOnly:"
	for offspr in $total_boot_rate.keys
		puts offspr + "\t" + mean_ci($total_boot_rate[offspr][0],$total_boot_rate[offspr][1]) + "\t" + mean_ci($total_boot_rate[offspr][2],$total_boot_rate[offspr][3])
	end
	puts $mutations
end
#----------------------------------------------------------------------------------------
if ARGV[0].nil?
	# If no parameters passed, print help screen
	format_splash('summarize_denovo', SUMMARIZEDENOVOVER, '<directory> > <out.txt>')
else
	# If input logs found, calculate the overall DNM rates
	$total_sites = 0 # Total length of observed bp across all chromosomes
	$bootstrap_bp = 0 # Total length of bootstrapped bp across all chromosomes
	$total_denovo = {} # Total mutations across all chromosomes keyed by offspring
	$total_boot_rate = {} # Total bootstrapped mutation rates across all chromosomes keyed by offspring
	$mutations = "\nIdentified de novo mutations:\n" # Reduced VCF of ided mutations
	vcf_header_collected = false # Whether reduced VCF header has been collected
	Dir.foreach(ARGV[0] + "/") do |f1|
		if f1[-3..-1] == "log"
			File.open(ARGV[0] + "/" + f1) do |f2|
				collection_stage = 0 # 0: no collection, 1: collect original data, 2: collect boostrapped data
				collect_offspring = false # Get offspring data
				collect_vcf = false # Get VCF data
				while line = f2.gets
					if line == "Total Sample Results:\n" # If never triggered, prevents inclusion of bad data
						collection_stage = 1
					elsif collection_stage == 1
						if line == "Offspring\tSingle-Forward\tDouble-Forward\tBackward\n"
							collect_offspring = true
						elsif collect_offspring
							if line == "\n"
								collect_offspring = false
							else
								line_arr = line.split("\t")
								unless $total_denovo.keys.include?(line_arr[0])
									$total_denovo[line_arr[0]] = line_arr[1..3].map { |x| x.to_i }
									$total_boot_rate[line_arr[0]] = [0,0,0,0]
								else
									for i in 1..3
										$total_denovo[line_arr[0]][i-1] += line_arr[i].to_i
									end
								end 
							end
						elsif line[0..5] == "#CHROM"
							collect_vcf = true
							unless vcf_header_collected
								$mutations << line
								vcf_header_collected = true
							end
						elsif collect_vcf
							if line == "\n"
								collect_vcf = false
							else
								$mutations << line
							end
						elsif line[0..8] == "Bootstrap"
							collection_stage = 2
							$bootstrap_bp += @current_contig_bp # Only add this value if the contig is actually bootstrapped
						elsif line[0..31] == "Total number of retained sites: "
							@current_contig_bp = line.split(" ")[-1].to_i
							$total_sites += @current_contig_bp
						end
					elsif collection_stage == 2
						if line == "Offspring\tMean_AllSites\t95%C.I._Allsites\tMean_Single-ForwardOnly\t95%C.I._Single-ForwardOnly:\n"
							collect_offspring = true
						elsif collect_offspring
							if line == "\n"
								collect_offspring = false
								collection_stage = 1 # need to reset stage to collect unbootstrapped mutations
							else
								line_arr = line.split("\t")
								for i in 1..4
									case i
									when 1,3
										@current_mean = line_arr[i].to_f
										$total_boot_rate[line_arr[0]][i-1] += @current_mean * @current_contig_bp.to_f
									when 2,4
										crit = @current_mean - line_arr[i].split("..")[0].to_f
										$total_boot_rate[line_arr[0]][i-1] += crit * @current_contig_bp.to_f
									end
								end
							end
						end
					end
				end
			end
		end
	end
	print_summary
end 
