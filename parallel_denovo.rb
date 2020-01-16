#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# parallel_denovo
PARDENOVOVER = "0.6.0"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

require_relative 'denovolib'

def split_vcf
	@vcfs = [] # Array of VCF names. Exclude filepath from names
	start = false
	header = ""
	outfile = ""
	outlines = "" # Lines to store until writing to disk
	writecycles = 0 # Number of write cycles passed
	gz_file_open($options.infile).open($options.infile) do |f1|
		while line = f1.gets
			if start
				contig = line.split("\t")[0]
				if contig != outfile
					if outfile != "" # Write remaining lines stored in memory
						File.open($options.outdir + "/chr" + outfile + ".vcf", 'a') do |write|
							write << outlines
						end
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
	return @vcfs
end
#-----------------------------------------------------------------------------------------------
def get_vcfs
	@vcfs = []
	Dir.foreach($options.outdir + "/") do |f1|
		@vcfs.push(f1[0..-5]) if f1[-3..-1] == "vcf"
	end
	return @vcfs
end
#-----------------------------------------------------------------------------------------
def execute_qsub(file)
	system("qsub -N #{file} -v SCAFFOLD=#{file} #{$options.outdir}/calc_denovo_mutation_rate.job")
end
#-----------------------------------------------------------------------------------------
def summarize_vcfs(vcfs)
	# Clean up job output
	#for vcf in vcfs
	#	system("rm #{vcf}.o*")
	#end
	system("qsub -N summarize_mutation_rate -hold_jid #{vcfs.join(",")} #{$options.outdir}/summarize_mutation_rate.job")
end
#-----------------------------------------------------------------------------------------
def write_qsub
	header = "#!/bin/sh\n#$ -S /bin/sh\n#$ -q #{$options.queue}\n#$ -l mres=#{$options.memory},h_data=#{$options.memory},h_vmem=#{$options.memory}"
	header << ",himem" if $options.himem
	header << ",lopri" if $options.lopri
	header << "\n#$ -j y\n#$ -cwd\n"
	header << "#$ -m bea\n#$ -M #{$options.email}\n" if $options.email != ""
	header << "module load bioinformatics/ruby/2.6.3\n"
	File.open($options.outdir + "/summarize_mutation_rate.job", 'w') do |f1|
		f1.puts header
		f1.puts "ruby summarize_denovo.rb #{$options.outdir} > #{$options.outdir}/summary.log"
	end
	header << "ruby calc_denovo_mutation_rate.rb -i #{$options.outdir}/${SCAFFOLD}.vcf -s #{$options.sire} -d #{$options.dam} -w #{$options.window} -S #{$options.step} -b #{$options.bootstrap} --rng #{$options.rng}"
	header << " -g" if $options.gvcf
	header << " > #{$options.outdir}/${SCAFFOLD}.log"
	File.open("#{$options.outdir}/calc_denovo_mutation_rate.job", 'w') do |qsub|
		qsub.puts header
	end
end
#-----------------------------------------------------------------------------------------
ARGV[0] ||= "-h"
$options = Parser.parse(ARGV, true)
Dir.mkdir($options.outdir) if !FileTest.directory?($options.outdir)
$options.restart ? vcfs = get_vcfs : vcfs = split_vcf
write_qsub
for vcf in vcfs
	execute_qsub(vcf)
end
summarize_vcfs(vcfs)
