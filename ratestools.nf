#!/usr/bin/env nextflow

/* RatesTools version 0.5.11
Michael G. Campana and Ellie E. Armstrong, 2020-2023
Smithsonian Institution and Stanford University

CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
University have waived all copyright and related or neighboring rights to RatesTools;
this work is published from the United States. You should have received a copy of the
CC0 legal code along with this work. If not, see 
<http://creativecommons.org/publicdomain/zero/1.0/>.
 
We politely request that this work be cited as:
Armstrong, E.E. & M.G. Campana. 2023. RatesTools: a Nextflow pipeline for detecting
de novo germline mutations in pedigree sequence data. Bioinformatics. 39: btac784.
doi: 10.1093/bioinformatics/btac784. */

nextflow.enable.dsl=1

readpairs_ch = Channel.fromFilePairs(params.reads) // Group read pairs
// algorithm for BWA index. Next lines make into command-line option for prepareRef
if ( params.bwa_alg == "" ) {
	bwa_alg = ""
} else {
	bwa_alg = "-a " + params.bwa_alg + " "
}	
chr_file = file(params.chr_file)

// Generate Picard and GATK executable commands
if ( params.picard_conda ) {
	picard = "picard " + params.picard_java
} else {
	picard = "java " + params.picard_java + " -jar " + params.picard
}
if ( params.gatk_conda ) {
	if ( params.gatk_build == 3 ) {
		gatk = "gatk3 " + params.gatk_java
	} else if (params.gatk_build == 4 ) {
		gatk = 'gatk --java-options "' + params.gatk_java + '" '
	}
} else {
	gatk = "java " + params.gatk_java + " -jar " + params.gatk
}

process prepareRef {

	// Prepare reference sequence for downstream processing
	
	label 'bwa'
	label 'samtools'
	errorStrategy 'finish'
	
	input:
	path refseq from params.refseq
	
	output:
	path "${refseq.baseName}*.{amb,ann,bwt,pac,sa,fai,dict}" into bwa_index_ch // Channel to find these files again
	
	"""
	bwa index ${bwa_alg}${refseq}
	samtools faidx ${refseq}
	samtools dict ${refseq} > ${refseq.baseName}.dict
	"""

}

process alignSeqs {

	// Align fastqs against reference sequence
	
	label 'bwa'
	label 'samtools'
	errorStrategy 'finish'
	
	input:
	path refseq from params.refseq
	tuple val(pair_id), path(reads) from readpairs_ch
	path "*" from bwa_index_ch
	
	output:
	path "${pair_id}_${refseq.simpleName}.bam" into raw_bam_ch
	val pair_id into filt_sample_ch, stats_sample_ch // sample name for downstream use
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	bwa mem -t ${task.cpus} ${refseq} ${reads} | samtools view -@ ${samtools_extra_threads} -bS -o ${pair_id}_${refseq.simpleName}.bam -  
	"""
	
}

process sortBAM {

	// Sort aligned bams

	label 'samtools'
	errorStrategy 'finish'

	input:
	path raw_bam from raw_bam_ch
	
	output:
	path "${raw_bam.simpleName}.srt.bam" into sorted_bam_ch
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	samtools sort -@ ${samtools_extra_threads} ${raw_bam} > ${raw_bam.simpleName}.srt.bam
	"""
	
}

process markDuplicates {

	// Mark duplicates using sambamba or picard
	
	label 'sambamba'
	label 'picard'
	errorStrategy 'finish'
	
	input:
	path sorted_bam from sorted_bam_ch
	val markDuplicates from params.markDuplicates
	
	output:
	path "${sorted_bam.simpleName}.markdup.bam" into markdup_bam_ch
	
	script:
	if ( markDuplicates == "sambamba" )
		"""
		sambamba markdup ${sorted_bam} ${sorted_bam.simpleName}.markdup.bam
		"""
	else
		"""
		$picard MarkDuplicates I=${sorted_bam} O=${sorted_bam.simpleName}.markdup.bam M=${sorted_bam.simpleName}.markdup.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
		"""
		
}

process fixReadGroups {

	// Fix read groups using picard
	
	label 'picard'
	errorStrategy 'finish'
	
	input:
	path markdup_bam from markdup_bam_ch
	path refseq from params.refseq
	
	output:
	path "${markdup_bam.simpleName}.rg.bam" into rg_bam_ch
	
	script:
	pair_id = "${markdup_bam.simpleName}".minus("_${refseq.simpleName}")
	"""
	$picard AddOrReplaceReadGroups I=${markdup_bam} O=${markdup_bam.simpleName}.rg.bam RGID=$pair_id RGLB=$pair_id RGPL=illumina RGPU=$pair_id RGSM=$pair_id
	"""

}

process realignIndels {

	// GATK RealignerTargetCreator. Index bam with picard.
	
	label 'gatk'
	label 'picard'
	errorStrategy 'finish'
	
	input:
	path rg_bam from rg_bam_ch
	path refseq from params.refseq
	path "*" from bwa_index_ch
	
	output:
	path "${rg_bam.simpleName}.realn.bam" into realn_bam_ch
	path "${rg_bam.simpleName}.realn.bai" into realn_bai_ch
	path "${rg_bam.baseName}.bai"
	
	script:
	if (params.gatk_build == 3)
		"""
		$picard BuildBamIndex I=${rg_bam}
		$gatk -T RealignerTargetCreator -R ${refseq} -I ${rg_bam} -o ${rg_bam.baseName}.intervals
		$gatk -T IndelRealigner -R ${refseq} --filter_bases_not_stored -I ${rg_bam} -targetIntervals ${rg_bam.baseName}.intervals -o ${rg_bam.simpleName}.realn.bam
		"""
	else if (params.gatk_build == 4)
		"""
		$picard BuildBamIndex I=${rg_bam}
		$gatk LeftAlignIndels -R ${refseq} -I $rg_bam -O ${rg_bam.simpleName}.realn.bam
		"""

}

process filterBAMs {

	// Filter BAMs using GATK
	
	label 'gatk'
	errorStrategy 'finish'
	
	input:
	path refseq from params.refseq
	path realn_bam from realn_bam_ch
	path realn_bai from realn_bai_ch
	path "*" from bwa_index_ch
	
	output:
	path "${realn_bam.simpleName}.filt.bam" into filt_bam_ch
	path "${realn_bam.simpleName}.filt.bai" into filt_bai_ch
	
	script:
	if (params.gatk_build == 3)
		"""
		$gatk -R ${refseq} -T PrintReads -I ${realn_bam} -o ${realn_bam.simpleName}.filt.bam -nct ${task.cpus} --read_filter BadCigar --read_filter DuplicateRead --read_filter FailsVendorQualityCheck --read_filter HCMappingQuality --read_filter MappingQualityUnavailable --read_filter NotPrimaryAlignment --read_filter UnmappedRead --filter_bases_not_stored --filter_mismatching_base_and_quals
		"""
	else if (params.gatk_build == 4)
		"""
		$gatk PrintReads -I ${realn_bam} -O ${realn_bam.simpleName}.filt.bam --read-filter GoodCigarReadFilter --read-filter NotDuplicateReadFilter --read-filter PassesVendorQualityCheckReadFilter --read-filter MappingQualityReadFilter --read-filter MappingQualityAvailableReadFilter --read-filter PrimaryLineReadFilter --read-filter MappedReadFilter --read-filter NotOpticalDuplicateReadFilter --read-filter ProperlyPairedReadFilter
		"""
}

process fixMate {

	// Fix mate information using Picard
	
	label 'picard'
	publishDir "$params.outdir/01_FinalBAMs", mode: 'copy'
	errorStrategy 'finish'
		
	input:
	path filt_bam from filt_bam_ch
	path filt_bai from filt_bai_ch
	
	output:
	path "${filt_bam.simpleName}.fix.bam" into fix_bam_ch
	path "${filt_bam.simpleName}.fix.bai" into fix_bai_ch
	
	"""
	$picard FixMateInformation I=${filt_bam} O=${filt_bam.simpleName}.fix.bam ADD_MATE_CIGAR=true
	$picard BuildBamIndex I=${filt_bam.simpleName}.fix.bam
	"""
	
}

process callVariants {

	// Index final bam and call variants using GATK
	
	label 'gatk'
	publishDir "$params.outdir/02_gVCFs", mode: 'copy'
	errorStrategy 'finish'
		
	input:
	path refseq from params.refseq
	path fix_bam from fix_bam_ch
	path fix_bai from fix_bai_ch
	path "*" from bwa_index_ch
	
	output:
	path "${fix_bam.simpleName}.g.vcf.*" into var_vcf_ch
	
	script:
	if (params.gatk_build == 3)
		"""
		$gatk -T HaplotypeCaller -nct ${task.cpus} -R ${refseq} -A DepthPerSampleHC -A Coverage -A HaplotypeScore -A StrandAlleleCountsBySample -I ${fix_bam} -o ${fix_bam.simpleName}.g.vcf.gz -ERC GVCF -out_mode EMIT_ALL_SITES
		"""
	else if (params.gatk_build == 4)
		"""
		$gatk HaplotypeCaller -R $refseq -I $fix_bam -O ${fix_bam.simpleName}.g.vcf.gz -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation
		"""

}

process genotypegVCFs {

	// Joint genotype gVCFs using GATK
	// Uncompressed combined VCF because GATK 3.8-1 inflate/deflate is glitched for this tool
	// Also uncompressed combined VCF because GATK4 cannot generate index with gzipped output
	
	label 'gatk'
	label 'gzip'
	publishDir "$params.outdir/03_CombinedVCF", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	path refseq from params.refseq
	path "*" from bwa_index_ch
	path "*" from var_vcf_ch.collect()
	val prefix from params.prefix
	
	output:
	path "${prefix}_combined.vcf*"
	path "${prefix}_combined.vcf.gz" into combined_vcf_ch, combined_indels_ch
	
	script:
	if (params.gatk_build == 3)
		"""
		VARPATH=""
		for file in *.vcf.gz; do VARPATH+=" --variant \$file"; done
		$gatk -T GenotypeGVCFs -R ${refseq} --includeNonVariantSites\$VARPATH -o ${prefix}_combined.vcf
		gzip ${prefix}_combined.vcf
		"""
	else if (params.gatk_build == 4)
		"""
		VARPATH=""
		for file in *.vcf.gz; do VARPATH+=" --variant \$file"; done
		$gatk CombineGVCFs -R $refseq -O tmp.g.vcf.gz --convert-to-base-pair-resolution\$VARPATH
		$gatk GenotypeGVCFs -R $refseq --include-non-variant-sites -V tmp.g.vcf.gz -O ${prefix}_combined.vcf
		gzip ${prefix}_combined.vcf
		"""
}

process genMapIndex {

	// Generate GenMap index
	
	label 'genmap'
	errorStrategy 'finish'
	
	input:
	path refseq from params.refseq
	val gm_tmpdir from params.gm_tmpdir
	
	output:
	path "${refseq.simpleName}_index" into genmap_index_ch
	path "${refseq.simpleName}_index/*" into genmap_index_files_ch
	
	"""
	export TMPDIR=${gm_tmpdir}
	if [ ! -d ${gm_tmpdir} ]; then mkdir ${gm_tmpdir}; fi
	genmap index -F ${refseq} -I ${refseq.simpleName}_index
	"""

}

process genMapMap {

	// Calculate mappability using GenMap and filter using filterGM
	
	label 'genmap'
	label 'ruby'
	errorStrategy 'finish'
	
	input:
	path refseq from params.refseq
	path genmap_index from genmap_index_ch
	path '*' from genmap_index_files_ch
	
	output:
	path "${refseq.simpleName}_genmap.1.0.bed" into genmap_ch
	
	"""
	genmap map -K 30 -E 2 -T ${task.cpus} -I ${refseq.simpleName}_index/ -O ${refseq.simpleName}_genmap -b
	filterGM.rb ${refseq.simpleName}_genmap.bed 1.0 exclude > ${refseq.simpleName}_genmap.1.0.bed
	"""
}

 
process repeatMask {

	// Mask repeats using RepeatMasker and default RM libraries
	
	label 'repeatmasker'
	errorStrategy 'finish'
	publishDir "$params.outdir/04_RepeatMasking", mode: 'copy'
	
	input:
	path refseq from params.refseq
	val rm_species from params.rm_species
	
	output:
	path "${refseq}.masked" into rm_ref_ch
	path "${refseq}.out" into rm_out_ch
	
	"""
	RepeatMasker -pa ${task.cpus} -gccalc -nolow -species ${rm_species} ${refseq}
	if [ ! -f ${refseq}.masked ]; then # Handling for no repeats detected
		ln -s ${refseq} ${refseq}.masked
	fi
	"""

}

process repeatModeler {

	// RepeatModeler on soft-masked reference
	
	label 'repeatmodeler'
	errorStrategy 'finish'
	publishDir "$params.outdir/04_RepeatMasking", mode: 'copy'
	
	input:
	path refseq from params.refseq
	path refseq_masked from rm_ref_ch
	
	output:
	path "**consensi.fa.classified" into rm_lib_ch
	path refseq_masked into rm_ref_ch2
	
	script:
	// RepeatModeler adds an extra thread for each core for rmblastn
	rm_pa = task.cpus / 2
	"""
	BuildDatabase -name ${refseq.baseName}-soft ${refseq_masked}
	RepeatModeler -pa ${rm_pa} -database ${refseq.baseName}-soft
	if [ ! -f */consensi.fa.classified ]; then 
		mkdir dummy # For fake library
		touch dummy/consensi.fa.classified
	fi
	"""

}

process repeatMaskRM {

	// RepeatMask using custom RepeatModeler repeat library
	
	label 'repeatmasker'
	label 'ruby'
	errorStrategy 'finish'
	publishDir "$params.outdir/04_RepeatMasking", mode: 'copy'
	
	input:
	path refseq from params.refseq
	path rm_lib from rm_lib_ch
	path rm_out from rm_out_ch
	path refseq_masked from rm_ref_ch2
	
	output:
	path "${refseq}.masked.*" optional true
	path "${refseq}.RM.bed" into rm_bed_ch
	
	"""
	if [ "\$(wc -l < consensi.fa.classified)" -eq 0 ]; then
	# If no output from RepeatModeler, use original RepeatMasker results
		RM2bed.rb ${refseq}.out > ${refseq}.RM.bed
	else
		RepeatMasker -pa ${task.cpus} -gccalc -nolow -lib consensi.fa.classified ${refseq_masked}
		# Convert out file into BED for downstream
		RM2bed.rb ${refseq_masked}.out > ${refseq}.RM.bed
	fi
	"""

}

process maskIndels {

	// Mask indel-affected regions using indels2bed
	
	label 'ruby'
	errorStrategy 'finish'
	
	input:
	path combo_vcf from combined_indels_ch
	val indelpad from params.indelpad
	
	output:
	path "${combo_vcf.simpleName}_indels.bed" into indels_ch
	
	"""
	indels2bed.rb ${combo_vcf} $indelpad > ${combo_vcf.simpleName}_indels.bed
	"""

}

process simplifyBed {

	// Reduce number of unique BED entries using simplify_bed
	
	label 'ruby'
	errorStrategy 'finish'
	publishDir "$params.outdir/05_ExcludedRegions", mode: 'copy'
	
	input:
	path indel_bed from indels_ch
	path rm_bed from rm_bed_ch
	path gm_bed from genmap_ch
	val prefix from params.prefix
	
	output:
	path "${prefix}_excluded_reduced.bed" into exclude_bed_ch
	
	"""
	#!/usr/bin/env bash
	if [ ! "\$(wc -l < ${indel_bed})" -eq 0 ]; then
		simplify_sorted_bed.rb ${indel_bed} > ${indel_bed.baseName}_sorted.bed
	else
		ln -s ${indel_bed} ${indel_bed.baseName}_sorted.bed
	fi
	if [ ! "\$(wc -l < ${rm_bed})" -eq 0 ]; then
		simplify_sorted_bed.rb ${rm_bed} > ${rm_bed.baseName}_sorted.bed
	else
		ln -s ${rm_bed} ${rm_bed.baseName}_sorted.bed
	fi
	if [ ! "\$(wc -l < ${gm_bed})" -eq 0 ]; then
		simplify_sorted_bed.rb ${gm_bed} > ${gm_bed.baseName}_sorted.bed
	else
		ln -s ${gm_bed} ${gm_bed.baseName}_sorted.bed
	fi
	cat ${indel_bed.baseName}_sorted.bed ${rm_bed.baseName}_sorted.bed ${gm_bed.baseName}_sorted.bed > ${prefix}_excluded.bed
	if [ ! "\$(wc -l < ${prefix}_excluded.bed)" -eq 0 ]; then
		simplify_bed.rb ${prefix}_excluded.bed > ${prefix}_excluded_reduced.bed
	else
		ln -s ${prefix}_excluded.bed ${prefix}_excluded_reduced.bed
	fi
	"""

}

process filterChr {
	
	// Optionally include only specific chromsomes
	
	label 'vcftools'
	label 'bcftools'
	label  'gzip'
	publishDir "$params.outdir/06_FilterChrVCFs", mode: 'copy', pattern: '*chrfilt.recode.vcf.gz'
	errorStrategy 'finish'
	
	input:
	path comb_vcf from combined_vcf_ch
	val prefix from params.prefix
	path chrs from chr_file
	
	output:
	path "${prefix}.chrfilt.recode.vcf.gz" into chrfilt_vcf_ch, chrfilt_stats_ch
	tuple path('chrfilt.tmp'), path(comb_vcf), path("${prefix}.chrfilt.recode.vcf.gz") optional true into chrfilt_log_ch 
	
	script:
	if (chrs.name == "NULL")
		"""
		ln -s $comb_vcf ${prefix}.chrfilt.recode.vcf.gz
		"""
	else
		"""
		chr_line=`echo '--chr '`; chr_line+=`awk 1 ORS=' --chr ' ${chrs}`; chr_line=`echo \${chr_line% --chr }` # Awkwardly make into a --chr command-list
		vcftools --gzvcf $comb_vcf --recode -c \$chr_line | gzip > ${prefix}.chrfilt.recode.vcf.gz
		cp .command.log chrfilt.tmp
		"""
}

filt_sample_ch2 = filt_sample_ch.filter { it != params.sire && it != params.dam } // Need new channel after filtering this one to remove dam and sire from offspring lists

process splitTrios {
	
	// Split samples into trios for analysis
	
	label 'vcftools'
	label 'gzip'
	label 'bcftools'
	publishDir "$params.outdir/07_SplitTrioVCFs", mode: 'copy', pattern: "${params.prefix}_offspring*.chrfilt.recode.vcf.gz"
	errorStrategy 'finish'
	
	input:
	path chr_vcf from chrfilt_vcf_ch
	val dam from params.dam
	val sire from params.sire
	val prefix from params.prefix
	val pair_id from filt_sample_ch2
	path chrs from chr_file
	
	output:
	path "${prefix}_offspring${pair_id}.chrfilt.recode.vcf.gz" into triosplit_vcf_ch, candidate_dnms_header_ch
	tuple path("${prefix}_offspring${pair_id}_trio.tmp"), path(chr_vcf), path("${prefix}_offspring${pair_id}.chrfilt.recode.vcf.gz") into triosplit_log_ch
	
	"""
	vcftools --gzvcf $chr_vcf --recode -c --indv ${dam} --indv ${sire} --indv ${pair_id} | gzip > ${prefix}_offspring${pair_id}.chrfilt.recode.vcf.gz
	cp .command.log ${prefix}_offspring${pair_id}_trio.tmp
	"""

}

process pullDPGQ {

	// Extract DP/GQ values from autosome scaffolds to look at the distributions of DP and GQ for variant call sites

	label 'bcftools'
	publishDir "$params.outdir/08_gVCFs_DP_GQ", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	path chrfilt from chrfilt_stats_ch
	val pair_id from stats_sample_ch
	
	output:
	path "${chrfilt.simpleName}_ind${pair_id}.variants.txt" into dp_gq_ch
	
	"""
	bcftools view -v snps ${chrfilt} -s ${pair_id} | bcftools query -f \"%CHROM %POS [ %DP] [ %GQ]\\n\" -o ${chrfilt.simpleName}_ind${pair_id}.variants.txt
	"""

}

process plotDPGQ {

	// Plot DP and GQ distributions
	
	label 'R'
	publishDir "$params.outdir/08_gVCFs_DP_GQ", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	path "*.txt" from dp_gq_ch.collect()
	
	output:
	path "${params.prefix}_*.png"
	path "${params.prefix}_depth_ratestools.csv"
	path "${params.prefix}_qual_ratestools.csv"
	
	"""
	plotDPGQ.R $params.prefix
	"""
	
}

process splitVCFs {

	// Split VCFs by contig/chromosome/scaffold etc
	
	label 'ruby'
	label 'bgzip'
	publishDir "$params.outdir/09_SplitVCFs", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	path chr_vcf from triosplit_vcf_ch
	
	output:
	path "${chr_vcf.simpleName}_split/*vcf.gz" into split_vcfs_ch mode flatten
	
	"""
	nextflow_split.rb -i ${chr_vcf} -o ${chr_vcf.simpleName}_split
	cd ${chr_vcf.simpleName}_split
	for file in *vcf.gz; do mv \$file ${chr_vcf.simpleName}_\${file}; done
	cd ..
	"""
	
}

process vcftoolsFilterSites {

	// Filter sites using VCFtools
	
	label 'vcftools'
	label 'bcftools'
	label 'bgzip'
	publishDir "$params.outdir/10_VCFtoolsSiteFilteredVCFs", mode: 'copy', pattern: '*sitefilt.recode.vcf.gz'
	errorStrategy 'finish'
	
	input:
	path split_vcf from split_vcfs_ch
	val site_filters from params.vcftools_site_filters
	
	output:
	tuple path("${split_vcf.simpleName}_sitefilt.tmp"), path(split_vcf), path("${split_vcf.simpleName}.sitefilt.recode.vcf.gz") into sitefilt_log_ch
	
	script:
	if (site_filters == "NULL")
		"""
		cp -P $split_vcf ${split_vcf.simpleName}.sitefilt.recode.vcf.gz
		vcftools --gzvcf $split_vcf
		cp .command.log ${split_vcf.simpleName}_sitefilt.tmp
		"""
	else
		"""
		vcftools --gzvcf ${split_vcf} --recode -c ${site_filters} | bgzip > ${split_vcf.simpleName}.sitefilt.recode.vcf.gz
		cp .command.log ${split_vcf.simpleName}_sitefilt.tmp
		"""

}

process sanityCheckLogsVcftools {

	// Sanity check logs for VCFtools site filtering and remove too short contigs

	label 'gzip'
	errorStrategy 'finish'

	input:
	tuple path(logfile), path(allvcflog), path(filtvcflog) from sitefilt_log_ch
	val min_contig_length from params.min_contig_length
	val min_filt_contig_length from params.min_filt_contig_length
	
	output:
	path "${logfile.simpleName}.log" into sitefilt_log_sanity_ch
	path "${filtvcflog.simpleName}.sitefilt.recode.OK.vcf.gz" optional true into sitefilt_vcf_ch
	
	"""
	logstats.sh $logfile $allvcflog $filtvcflog $min_contig_length $min_filt_contig_length > ${logfile.simpleName}.log
	"""
	
}

process gatkFilterSites {

	// Apply GATK-only site filters
	
	label 'gatk'
	label 'tabix'
	label 'vcftools'
	label 'bcftools'
	publishDir "$params.outdir/11_GATKSiteFilteredVCFs", mode: 'copy', pattern: '*gatksitefilt.vcf.gz'
	errorStrategy 'finish'
	
	input:
	path refseq from params.refseq
	path "*" from bwa_index_ch
	path site_vcf from sitefilt_vcf_ch
	val site_filters from params.gatk_site_filters
	
	output:
	tuple path("${site_vcf.simpleName}_gatksitefilt.tmp"), path(site_vcf), path("${site_vcf.simpleName}.gatksitefilt.vcf.gz") into gatk_sitefilt_log_ch
	
	script:
	if (site_filters == "NULL")
		"""
		ln -s $site_vcf ${site_vcf.simpleName}.gatksitefilt.vcf.gz
		vcftools --gzvcf $site_vcf
		cp .command.log  ${site_vcf.simpleName}_gatksitefilt.tmp
		"""
	else if (params.gatk_build == 3)
		"""
		tabix $site_vcf
		$gatk -T VariantFiltration -V $site_vcf -o tmp.vcf -R $refseq $site_filters
		$gatk -T SelectVariants -V tmp.vcf -o ${site_vcf.simpleName}.gatksitefilt.vcf.gz -R $refseq --excludeFiltered
		vcftools --gzvcf ${site_vcf.simpleName}.gatksitefilt.vcf.gz
		cp .command.log  ${site_vcf.simpleName}_gatksitefilt.tmp
		"""
	else if (params.gatk_build == 4)
		"""
		tabix $site_vcf
		$gatk VariantFiltration -R $refseq -V $site_vcf -O tmp.vcf.gz $site_filters
		$gatk SelectVariants -R $refseq -V $site_vcf -O ${site_vcf.simpleName}.gatksitefilt.vcf.gz --exclude-filtered
		vcftools --gzvcf ${site_vcf.simpleName}.gatksitefilt.vcf.gz
		cp .command.log  ${site_vcf.simpleName}_gatksitefilt.tmp
		"""

}

process sanityCheckLogsGatk {
	
	// Sanity check logs for GATK site filtering and remove too short contigs
	// Dummy value of 1 for min_contig_length since already evalutated and no longer accurate

	label 'gzip'
	errorStrategy 'finish'

	input:
	tuple path(logfile), path(allvcflog), path(filtvcflog) from gatk_sitefilt_log_ch
	val min_filt_contig_length from params.min_filt_contig_length
	
	output:
	path "${logfile.simpleName}.log" into gatk_sitefilt_log_sanity_ch
	path "${filtvcflog.simpleName}.gatksitefilt.OK.vcf.gz" optional true into gatk_sitefilt_vcf_ch
	
	"""
	logstats.sh $logfile $allvcflog $filtvcflog 1 $min_filt_contig_length > ${logfile.simpleName}.log
	"""
	
}

process filterRegions {

	// Filter regions using BEDTools. If failure, tries to fix using zcat.
	// Attempts BCFtools filtration if BEDTools fails
	// Defaults to VCFtools in case of unrecoverable error 
	
	label 'bedtools'
	label 'bcftools'
	label 'vcftools'
	label 'tabix'
	label 'gzip'
	publishDir "$params.outdir/12_RegionFilteredVCFs", mode: 'copy', pattern: '*regionfilt.vcf.gz'
	errorStrategy 'retry'
	maxRetries 3
	
	input:
	path site_vcf from gatk_sitefilt_vcf_ch
	path exclude_bed from exclude_bed_ch
	
	output:
	tuple path("${site_vcf.simpleName}_regionfilt.tmp"), path(site_vcf), path("${site_vcf.simpleName}.regionfilt.vcf.gz") into regionfilt_log_ch
	
	script:
	chr = site_vcf.simpleName.split('_chr')[1]
	if (task.attempt == 1)
		"""
		#!/usr/bin/env bash
		grep ${chr} ${exclude_bed} > tmp.bed
		if [ ! "\$(wc -l < tmp.bed)" -eq 0 ]; then
			bedtools subtract -a ${site_vcf} -b tmp.bed -header | gzip > ${site_vcf.simpleName}.regionfilt.vcf.gz
			if [[ `grep -n 'Error: Invalid record' .command.log | cut -d ':' -f 1` -eq 0 ]]; then
				vcftools --gzvcf ${site_vcf.simpleName}.regionfilt.vcf.gz
				cp .command.log ${site_vcf.simpleName}_regionfilt.tmp
			else
				rm ${site_vcf.simpleName}.regionfilt.vcf.gz
			fi
		else
			vcftools --gzvcf ${site_vcf} --recode -c | gzip > ${site_vcf.simpleName}.regionfilt.vcf.gz
			cp .command.log ${site_vcf.simpleName}_regionfilt.tmp
		fi	
		"""
	else if (task.attempt == 2)
		"""
		#!/usr/bin/env bash
		grep ${chr} ${exclude_bed} > tmp.bed
		if [ ! "\$(wc -l < tmp.bed)" -eq 0 ]; then
			zcat ${site_vcf} | bedtools subtract -a stdin -b tmp.bed -header | gzip > ${site_vcf.simpleName}.regionfilt.vcf.gz
			if [[ `grep -n 'Error: Invalid record' .command.log | cut -d ':' -f 1` -eq 0 ]]; then
				vcftools --gzvcf ${site_vcf.simpleName}.regionfilt.vcf.gz
				cp .command.log ${site_vcf.simpleName}_regionfilt.tmp
			else
				rm ${site_vcf.simpleName}.regionfilt.vcf.gz
			fi
		else
			vcftools --gzvcf ${site_vcf} --recode -c | gzip > ${site_vcf.simpleName}.regionfilt.vcf.gz
			cp .command.log ${site_vcf.simpleName}_regionfilt.tmp
		fi	
		"""
	else if (task.attempt == 3)
		"""
		#!/usr/bin/env bash
		grep ${chr} ${exclude_bed} > tmp.bed
		if [ ! "\$(wc -l < tmp.bed)" -eq 0 ]; then
			tabix ${site_vcf}
			bcftools view -R tmp.bed -Ob -o tmp.bcf ${site_vcf}
			tabix tmp.bcf
			# bcftools isec gives the target sites not included in the bed
			bcftools isec -C -O v ${site_vcf} tmp.bcf | cut -f1,2 > ${site_vcf.simpleName}.targets
			# Use the isec output to get the output. Needs to stream (-T) rather than index jump (-R) for efficiency.
			bcftools view -T ${site_vcf.simpleName}.targets -Ov ${site_vcf} | gzip > ${site_vcf.simpleName}.regionfilt.vcf.gz
			vcftools --gzvcf ${site_vcf.simpleName}.regionfilt.vcf.gz
			cp .command.log ${site_vcf.simpleName}_regionfilt.tmp
		else
			vcftools --gzvcf ${site_vcf} --recode -c | gzip > ${site_vcf.simpleName}.regionfilt.vcf.gz
			cp .command.log ${site_vcf.simpleName}_regionfilt.tmp
		fi
		"""
	else
		"""
		#!/usr/bin/env bash
		grep ${chr} ${exclude_bed} > tmp.bed
		if [ ! "\$(wc -l < tmp.bed)" -eq 0 ]; then
			vcftools --gzvcf ${site_vcf} --recode -c --exclude-bed tmp.bed | gzip  > ${site_vcf.simpleName}.regionfilt.vcf.gz
			cp .command.log ${site_vcf.simpleName}_regionfilt.tmp
		else
			vcftools --gzvcf ${site_vcf} --recode -c | gzip > ${site_vcf.simpleName}.regionfilt.vcf.gz
			cp .command.log ${site_vcf.simpleName}_regionfilt.tmp
		fi
		"""

}

process sanityCheckLogsRegions {

	// Sanity check logs for region filtering and remove too short contigs
	// Dummy value of 1 for min_contig_length since already evalutated and no longer accurate

	label 'gzip'
	errorStrategy 'finish'

	input:
	tuple path(logfile), path(allvcflog), path(filtvcflog) from regionfilt_log_ch
	val min_filt_contig_length from params.min_filt_contig_length
	
	output:
	path "${logfile.simpleName}.log" into regionfilt_log_sanity_ch
	path "${filtvcflog.simpleName}.regionfilt.OK.vcf.gz" optional true into regionfilt_vcf_ch
	
	"""
	logstats.sh $logfile $allvcflog $filtvcflog 1 $min_filt_contig_length > ${logfile.simpleName}.log
	"""
	
}

process calcDNMRate {

	// Calculate de novo mutations using calc_denovo_mutation_rate
	
	label 'ruby'
	publishDir "$params.outdir/13_SplitCalcDNMLogs", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	path splitvcf from regionfilt_vcf_ch
	val sire from params.sire
	val dam from params.dam
	val dnm_opts from params.dnm_opts
	
	output:
	path "${splitvcf.simpleName}.log" into split_logs_ch
	
	"""
	calc_denovo_mutation_rate.rb -i ${splitvcf} -s ${sire} -d ${dam} ${dnm_opts} > ${splitvcf.simpleName}.log
	"""

}

process summarizeDNM {

	// Calculate genome-wide DNM rate using summarize_denovo
	
	label 'ruby'
	label 'bcftools'
	label 'gzip'
	publishDir "$params.outdir/14_SummarizeDNMLogs", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	path "*" from split_logs_ch.collect()
	path "*" from candidate_dnms_header_ch.collect()
	
	output:
	path "${params.prefix}*_summary.log" into summary_log_ch
	path "${params.prefix}*_candidates.vcf.gz" into candidates_vcf_ch
	
	"""
	#!/usr/bin/env bash
	for file in ${params.prefix}*.log; do 
		if [ ! -d \${file%_chr*.log} ]; then
			mkdir \${file%_chr*.log}
		fi
		mv \$file \${file%_chr*.log}/
	done
	for outdir in ${params.prefix}*; do 
		if [ -d \$outdir ]; then
			summarize_denovo.rb \$outdir > \${outdir}_summary.log
		fi
	done
	for sumlog in *summary.log; do
		bcftools view -h \${sumlog/_summary.log/.chrfilt.recode.vcf.gz} > header.txt
		val=`grep -n \'#CHROM\' \$sumlog | cut -d \':\' -f 1`
		total=`wc -l \$sumlog | cut -d \' \' -f 1`
		let lncount=\$total-\$val
		if [ \$lncount -gt 0 ]; then
			tail -n \$lncount \$sumlog > tmp.txt
			cat header.txt tmp.txt | gzip > \${sumlog/_summary.log/_candidates.vcf.gz}
		else
			mv header.txt \${sumlog/_summary.log/_candidates.vcf}
			gzip \${sumlog/_summary.log/_candidates.vcf}
		fi
	done
	"""

}

all_logs_ch = triosplit_log_ch.mix(chrfilt_log_ch)

process sanityCheckLogs {

	// Sanity check logs for Chromosome filtering and trio-splitting, but perform no filtering

	label 'gzip'
	errorStrategy 'finish'

	input:
	tuple path(logfile), path(allvcflog), path(filtvcflog) from all_logs_ch
	
	output:
	path "${logfile.simpleName}.log" into logs_sanity_ch
	
	"""
	logstats.sh $logfile $allvcflog $filtvcflog 0 0 > ${logfile.simpleName}.log
	"""
	
}

all_logs_sanity_ch = logs_sanity_ch.mix(regionfilt_log_sanity_ch, gatk_sitefilt_log_sanity_ch, sitefilt_log_sanity_ch, summary_log_ch)

process generateSummaryStats {

	label 'ruby'
	publishDir "$params.outdir/15_SummaryStats", mode: "copy"
	errorStrategy 'finish'
	
	input:
	path "*" from all_logs_sanity_ch.collect()
	
	output:
	path "summary_stats.csv"
	
	"""
	dnm_summary_stats.rb . ${params.prefix} > summary_stats.csv
	"""

}

workflow.onError {
	println "RatesTools pipeline encountered an error. Error message: $workflow.errorMessage."
	if (params.email != "NULL") {
		sendMail(to: params.email, subject: "RatesTools error.", body: "RatesTools pipeline encountered an error. Error message: $workflow.errorMessage.")
	}
}
workflow.onComplete {
	// workflow.onComplete is also run when an error is encountered, but is triggered when all processes finish.
	// Running error messages using onError so that the user gets them at first possible opportunity
	if (workflow.success) {
		println "RatesTools pipeline completed successfully at $workflow.complete!"
		if (params.email != "NULL") {
			sendMail(to: params.email, subject: 'RatesTools successful completion', body: "RatesTools pipeline completed successfully at $workflow.complete!")
		}
	} else {
		println "RatesTools pipeline terminated with errors at $workflow.complete.\nError message: $workflow.errorMessage"
		if (params.email != "NULL") {
			sendMail(to: params.email, subject: 'RatesTools terminated with errors', body: "RatesTools pipeline terminated with errors at $workflow.complete.\nError message: $workflow.errorMessage")
		}
	}
}
