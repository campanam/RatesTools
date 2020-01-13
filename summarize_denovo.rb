#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# summarize_denovo
SUMMARIZEDENOVOVER = "0.1.0"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

#----------------------------------------------------------------------------------------
def mean_ci(mean, crit)
	return (mean/$total_bp.to_f).to_s + "\t" + ((mean-crit)/$total_bp.to_f).to_s + "-" + ((mean+crit)/$total_bp.to_f).to_s
end
#----------------------------------------------------------------------------------------
def print_summary
	puts "\nTotal Sample Results:"
	puts "Total number of retained sites: " + $total_bp.to_s
	puts "\nTotal numbers of observed de novo mutations:"
	puts "Offspring\tSingle-Forward\tDouble-Forward\tBackward"
	for offspr in $total_mutations.keys
		puts offspr + "\t" + $total_mutations[offspr].join("\t")
	end
	puts "\nInferred mutation rates:"
	puts "Offspring\tAllsites\tSingle-ForwardOnly"
	for offspr in $total_mutations.keys
		puts offspr + "\t" + ($total_mutations[offspr].sum.to_f/$total_bp.to_f).to_s + "\t" + ($total_mutations[offspr][0].to_f/$total_bp.to_f).to_s
	end
	puts "\nBootstrapped Estimates:"
	puts "Offspring\tMean_AllSites\t95%C.I._Allsites\tMean_Single-ForwardOnly\t95%C.I._Single-ForwardOnly:"
	for offspr in $total_boot_rate.keys
		puts offspr + "\t" + mean_ci($total_boot_rate[offspr][0],$total_boot_rate[offspr][1]) + "\t" + mean_ci($total_boot_rate[offspr][2],$total_boot_rate[offspr][3])
	end
end
#----------------------------------------------------------------------------------------
if ARGV[0].nil?
	puts "\033[1msummarize_denovo " + SUMMARIZEDENOVOVER + "\033[0m"
	puts "\nUsage: ruby summarize_denovo.rb <directory> > <out.txt>"
else
	$total_bp = 0 # Total length of observed bp across all chromosomes
	$total_mutations = {} # Total mutations across all chromosomes keyed by offspring
	$total_boot_rate = {} # Total bootstrapped mutation rates across all chromosomes keyed by offspring
	Dir.foreach(ARGV[0] + "/") do |f1|
		if f1[-3..-1] == "log"
			File.open(ARGV[0] + "/" + f1) do |f2|
				collection_stage = 0 # 0: no collection, 1: collect original data, 2: collect boostrapped data
				collect_offspring = false # Get offspring data
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
								unless $total_mutations.keys.include?(line_arr[0])
									$total_mutations[line_arr[0]] = line_arr[1..3].map { |x| x.to_i }
									$total_boot_rate[line_arr[0]] = [0,0,0,0]
								else
									for i in 1..3
										$total_mutations[line_arr[0]][i-1] += line_arr[i].to_i
									end
								end 
							end
						elsif line[0..8] == "Bootstrap"
							collection_stage = 2
						elsif line[0..31] == "Total number of retained sites: "
							@current_contig_bp = line.split(" ")[-1].to_i
							$total_bp += @current_contig_bp
						end
					elsif collection_stage == 2
						if line == "Offspring\tMean_AllSites\t95%C.I._Allsites\tMean_Single-ForwardOnly\t95%C.I._Single-ForwardOnly:\n"
							collect_offspring = true
						elsif collect_offspring
							line_arr = line.split("\t")
							for i in 1..4
								case i
								when 1,3
									@current_mean = line_arr[i].to_f
									$total_boot_rate[line_arr[0]][i-1] += @current_mean * @current_contig_bp.to_f
								when 2,4
									crit = @current_mean - line_arr[i].split("-")[0].to_f
									$total_boot_rate[line_arr[0]][i-1] += crit * @current_contig_bp.to_f
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
