#!/usr/bin/env ruby

#----------------------------------------------------------------------------------------
# parallel_denovo
PARDENOVOVER = "0.10.1"
# Michael G. Campana, 2020
# Smithsonian Conservation Biology Institute
#----------------------------------------------------------------------------------------

require_relative 'denovolib'

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
	header << " --parhom" if $options.parhom
	header << " --minAD1" if $options.minAD1
	header << " --minAF " + $options.minAF.to_s unless $options.minAF.nil?
	header << " > #{$options.outdir}/${SCAFFOLD}.log"
	File.open("#{$options.outdir}/calc_denovo_mutation_rate.job", 'w') do |qsub|
		qsub.puts header
		qsub.puts "gzip #{$options.outdir}/${SCAFFOLD}.vcf"
	end
end
#-----------------------------------------------------------------------------------------
ARGV[0] ||= "-h"
$options = Parser.parse(ARGV, true)
Dir.mkdir($options.outdir) if !FileTest.directory?($options.outdir)
$options.restart ? vcfs = get_vcfs : vcfs = split_vcf
write_qsub unless $options.submit
unless $options.nosubmit
	for vcf in vcfs
		execute_qsub(vcf)
	end
	summarize_vcfs(vcfs)
end
