#!/usr/bin/env nextflow

/* RatesTools version 0.5
Michael G. Campana and Ellie E. Armstrong, 2020-2022
Smithsonian Institution and Stanford University

CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
University have waived all copyright and related or neighboring rights to RatesTools;
this work is published from the United States. You should have received a copy of the
CC0 legal code along with this work. If not, see 
<http://creativecommons.org/publicdomain/zero/1.0/>.
 
We politely request that this work be cited as:
Campana, M.G. & E.E. Armstrong. 2020. RatesTools: Pipeline to calculate de novo
mutation rates from parent-offspring trios. Smithsonian Institution and Stanford
University. <https://github.com/campanam/RatesTools>. */

readpairs_ch = Channel.fromFilePairs(params.reads) // Group read pairs
// algorithm for BWA index. Next lines make into command-line option for prepareRef
if ( params.bwa_alg == "" ) {
	bwa_alg = ""
} else {
	bwa_alg = "-a " + params.bwa_alg + " "
}	
chr_file = file(params.chr_file)

process prepareRef {

	// Prepare reference sequence for downstream processing
	
	label 'bwa'
	label 'samtools'
	errorStrategy 'finish'
	
	input:
	path refseq from params.refseq
	
	output:
	file "${refseq.baseName}*.{amb,ann,bwt,pac,sa,fai,dict}" into bwa_index_ch // Channel to find these files again
	
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
	val bwa_threads from params.bwa_threads
	file "*" from bwa_index_ch
	
	output:
	file "${pair_id}_${refseq.simpleName}.bam" into raw_bam_ch
	val pair_id into filt_sample_ch, stats_sample_ch // sample name for downstream use
	val samtools_extra_threads into samtools_threads_ch
	
	
	script:
	samtools_extra_threads = bwa_threads - 1
	"""
	bwa mem -t ${bwa_threads} ${refseq} ${reads} | samtools view -@ ${samtools_extra_threads} -bS -o ${pair_id}_${refseq.simpleName}.bam -  
	"""
	
}

process sortBAM {

	// Sort aligned bams

	label 'samtools'
	errorStrategy 'finish'

	input:
	file raw_bam from raw_bam_ch
	val samtools_extra_threads from samtools_threads_ch
	
	output:
	file "${raw_bam.simpleName}.srt.bam" into sorted_bam_ch
	
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
	file sorted_bam from sorted_bam_ch
	val markDuplicates from params.markDuplicates
	path picard from params.picard
	val picard_java from params.picard_java
	
	output:
	file "${sorted_bam.simpleName}.markdup.bam" into markdup_bam_ch
	
	script:
	if ( markDuplicates == "sambamba" )
		"""
		sambamba markdup ${sorted_bam} ${sorted_bam.simpleName}.markdup.bam
		"""
	else
		"""
		java ${picard_java} -jar ${picard} MarkDuplicates I=${sorted_bam} O=${sorted_bam.simpleName}.markdup.bam M=${sorted_bam.simpleName}.markdup.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
		"""
		
}

process fixReadGroups {

	// Fix read groups using picard
	
	label 'picard'
	errorStrategy 'finish'
	
	input:
	path markdup_bam from markdup_bam_ch
	path picard from params.picard
	val picard_java from params.picard_java
	path refseq from params.refseq
	
	output:
	file "${markdup_bam.simpleName}.rg.bam" into rg_bam_ch
	
	script:
	pair_id = "${markdup_bam.simpleName}".minus("_${refseq.simpleName}")
	"""
	java ${picard_java} -jar ${picard} AddOrReplaceReadGroups I=${markdup_bam} O=${markdup_bam.simpleName}.rg.bam RGID=$pair_id RGLB=$pair_id RGPL=illumina RGPU=$pair_id RGSM=$pair_id
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
	path picard from params.picard
	val picard_java from params.picard_java
	val gatk from params.gatk
	val gatk_java from params.gatk_java
	file "*" from bwa_index_ch
	
	output:
	file "${rg_bam.simpleName}.realn.bam" into realn_bam_ch
	file "${rg_bam.simpleName}.realn.bai" into realn_bai_ch
	file "${rg_bam.baseName}.bai"
	
	script:
	if (params.gatk_build == 3)
		"""
		java ${picard_java} -jar ${picard} BuildBamIndex I=${rg_bam}
		java ${gatk_java} -jar ${gatk} -T RealignerTargetCreator -R ${refseq} -I ${rg_bam} -o ${rg_bam.baseName}.intervals
		java ${gatk_java} -jar ${gatk} -T IndelRealigner -R ${refseq} --filter_bases_not_stored -I ${rg_bam} -targetIntervals ${rg_bam.baseName}.intervals -o ${rg_bam.simpleName}.realn.bam
		"""
	else if (params.gatk_build == 4)
		"""
		java ${picard_java} -jar ${picard} BuildBamIndex I=${rg_bam}
		java ${gatk_java} -jar ${gatk} LeftAlignIndels -R ${refseq} -I $rg_bam -O ${rg_bam.simpleName}.realn.bam
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
	val gatk from params.gatk
	val gatk_java from params.gatk_java
	val gatk_nct from params.gatk_nct
	file "*" from bwa_index_ch
	
	output:
	file "${realn_bam.simpleName}.filt.bam" into filt_bam_ch
	file "${realn_bam.simpleName}.filt.bai" into filt_bai_ch
	
	script:
	if (params.gatk_build == 3)
		"""
		java ${gatk_java} -jar ${gatk} -R ${refseq} -T PrintReads -I ${realn_bam} -o ${realn_bam.simpleName}.filt.bam -nct ${gatk_nct} --read_filter BadCigar --read_filter DuplicateRead --read_filter FailsVendorQualityCheck --read_filter HCMappingQuality --read_filter MappingQualityUnavailable --read_filter NotPrimaryAlignment --read_filter UnmappedRead --filter_bases_not_stored --filter_mismatching_base_and_quals
		"""
	else if (params.gatk_build == 4)
		"""
		java ${gatk_java} -jar ${gatk} PrintReads -I ${realn_bam} -O ${realn_bam.simpleName}.filt.bam --read-filter GoodCigarReadFilter --read-filter NotDuplicateReadFilter --read-filter PassesVendorQualityCheckReadFilter --read-filter MappingQualityReadFilter --read-filter MappingQualityAvailableReadFilter --read-filter PrimaryLineReadFilter --read-filter MappedReadFilter --read-filter NotOpticalDuplicateReadFilter --read-filter NotProperlyPairedReadFilter
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
	path picard from params.picard
	val picard_java from params.picard_java
	
	output:
	file "${filt_bam.simpleName}.fix.bam" into fix_bam_ch
	file "${filt_bam.simpleName}.fix.bai" into fix_bai_ch
	
	"""
	java ${picard_java} -jar ${picard} FixMateInformation I=${filt_bam} O=${filt_bam.simpleName}.fix.bam ADD_MATE_CIGAR=true
	java ${picard_java} -jar ${picard} BuildBamIndex I=${filt_bam.simpleName}.fix.bam
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
	path picard from params.picard
	val picard_java from params.picard_java
	val gatk from params.gatk
	val gatk_java from params.gatk_java
	val gatk_nct from params.gatk_nct
	file "*" from bwa_index_ch
	
	output:
	file "${fix_bam.simpleName}.g.vcf.*" into var_vcf_ch
	
	script:
	if (params.gatk_build == 3)
		"""
		java ${gatk_java} -jar ${gatk} -T HaplotypeCaller -nct ${gatk_nct} -R ${refseq} -A DepthPerSampleHC -A Coverage -A HaplotypeScore -A StrandAlleleCountsBySample -I ${fix_bam} -o ${fix_bam.simpleName}.g.vcf.gz -ERC GVCF -out_mode EMIT_ALL_SITES
		"""
	else if (params.gatk_build == 4)
		"""
		java ${gatk_java} -jar ${gatk} HaplotypeCaller -R $refseq -I $fix_bam -O ${fix_bam.simpleName}.g.vcf.gz -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation
		"""

}

process genotypegVCFs {

	// Joint genotype gVCFs using GATK
	// Uncompressed combined VCF because GATK 3.8-1 inflate/deflate is glitched for this tool
	
	label 'gatk'
	publishDir "$params.outdir/03_CombinedVCF", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	path refseq from params.refseq
	file "*" from bwa_index_ch
	file "*" from var_vcf_ch.collect()
	val gatk from params.gatk
	val gatk_java from params.gatk_java
	val prefix from params.prefix
	
	output:
	file "${prefix}_combined.vcf*"
	file "${prefix}_combined.vcf.gz" into combined_vcf_ch, combined_indels_ch
	
	script:
	if (params.gatk_build == 3)
		"""
		VARPATH=""
		for file in *.vcf.gz; do VARPATH+=" --variant \$file"; done
		java ${gatk_java} -jar ${gatk} -T GenotypeGVCFs -o ${prefix}_combined.vcf -R ${refseq} --includeNonVariantSites\$VARPATH
		gzip ${prefix}_combined.vcf # Clean up giant VCF before going on with filtering
		"""
	else if (params.gatk_build == 4)
		"""
		VARPATH=""
		for file in *.vcf.gz; do VARPATH+=" --variant \$file"; done
		java ${gatk_java} -jar ${gatk} CombineGVCFs -R $refseq -O tmp.g.vcf.gz --convert-to-base-pair-resolution\$VARPATH
		java ${gatk_java} -jar ${gatk} GenotypeGVCFs -R $refseq --include-non-variant-sites -I tmp.g.vcf.gz -O ${prefix}_combined.vcf.gz
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
	file "${refseq.simpleName}_index/*" into genmap_index_files_ch
	
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
	val gm_threads from params.gm_threads
	path genmap_index from genmap_index_ch
	file '*' from genmap_index_files_ch
	
	output:
	file "${refseq.simpleName}_genmap.1.0.bed" into genmap_ch
	
	"""
	genmap map -K 30 -E 2 -T ${gm_threads} -I ${refseq.simpleName}_index/ -O ${refseq.simpleName}_genmap -b
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
	val rm_pa from params.rm_pa
	
	output:
	file "${refseq}.masked" into rm_ref_ch
	file "${refseq}.out" into rm_out_ch
	
	"""
	RepeatMasker -pa ${rm_pa} -gccalc -nolow -species ${rm_species} ${refseq}
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
	file refseq_masked from rm_ref_ch
	val rm_pa from params.rm_pa
	
	output:
	file "**consensi.fa.classified" into rm_lib_ch
	file refseq_masked into rm_ref_ch2
	
	"""
	BuildDatabase -name ${refseq.baseName}-soft ${refseq_masked}
	RepeatModeler -pa ${rm_pa} -database ${refseq.baseName}-soft
	if [ ! -f */consensi.fa.classified ]; then 
		mkdir dummy # For fake library
		mkfile -n 0 dummy/consensi.fa.classified
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
	file rm_lib from rm_lib_ch
	file rm_out from rm_out_ch
	file refseq_masked from rm_ref_ch2
	val rm_pa from params.rm_pa
	
	output:
	file "${refseq}.masked.*" optional true
	file "${refseq}.RM.bed" into rm_bed_ch
	
	"""
	if [ "\$(wc -l < consensi.fa.classified)" -eq 0 ]; then
	# If no output from RepeatModeler, use original RepeatMasker results
		RM2bed.rb ${refseq}.out > ${refseq}.RM.bed
	else
		RepeatMasker -pa ${rm_pa} -gccalc -nolow -lib consensi.fa.classified ${refseq_masked}
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
	file combo_vcf from combined_indels_ch
	
	output:
	file "${combo_vcf.simpleName}_indels.bed" into indels_ch
	
	"""
	indels2bed.rb ${combo_vcf} 5 > ${combo_vcf.simpleName}_indels.bed
	"""

}

process simplifyBed {

	// Reduce number of unique BED entries using simplify_bed
	
	label 'ruby'
	errorStrategy 'finish'
	publishDir "$params.outdir/05_ExcludedRegions", mode: 'copy'
	
	input:
	file indel_bed from indels_ch
	file rm_bed from rm_bed_ch
	file gm_bed from genmap_ch
	val prefix from params.prefix
	
	output:
	file "${prefix}_excluded_reduced.bed" into exclude_bed_ch
	
	"""
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
	publishDir "$params.outdir/06_FilterChrVCFs", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	file comb_vcf from combined_vcf_ch
	val prefix from params.prefix
	file chrs from chr_file
	
	output:
	file "${prefix}.chrfilt.recode.vcf.gz" into chrfilt_vcf_ch, chrfilt_stats_ch
	
	script:
	if (chrs.name == "NULL")
		"""
		cp $comb_vcf ${prefix}.chrfilt.recode.vcf.gz
		"""
	else
		"""
		chr_line=`echo '--chr '`; chr_line+=`awk 1 ORS=' --chr ' ${chrs}`; chr_line=`echo \${chr_line% --chr }` # Awkwardly make into a --chr command-list
		vcftools --gzvcf $comb_vcf --recode --out ${prefix}.chrfilt \$chr_line
		gzip ${prefix}.chrfilt.recode.vcf
		"""
}

filt_sample_ch2 = filt_sample_ch.filter { it != params.sire && it != params.dam } // Need new channel after filtering this one to remove dam and sire from offspring lists

process splitTrios {
	
	// Split samples into trios for analysis
	
	label 'vcftools'
	publishDir "$params.outdir/07_SplitTrioVCFs", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	file chr_vcf from chrfilt_vcf_ch
	val dam from params.dam
	val sire from params.sire
	val prefix from params.prefix
	val pair_id from filt_sample_ch2
	file chrs from chr_file
	
	output:
	file "${prefix}_offspring${pair_id}.chrfilt.recode.vcf.gz" into triosplit_vcf_ch
	
	"""
	vcftools --gzvcf $chr_vcf --recode --out ${prefix}_offspring${pair_id}.chrfilt --indv ${dam} --indv ${sire} --indv ${pair_id}
	gzip ${prefix}_offspring${pair_id}.chrfilt.recode.vcf
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
	file "${params.prefix}_*.png"
	file "${params.prefix}_depth_ratestools.csv"
	file "${params.prefix}_qual_ratestools.csv"
	
	"""
	plotDPGQ.R $params.prefix
	"""
	
}

process splitVCFs {

	// Split VCFs by contig/chromosome/scaffold etc
	
	label 'ruby'
	publishDir "$params.outdir/09_SplitVCFs", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	file chr_vcf from triosplit_vcf_ch
	
	output:
	file "${chr_vcf.simpleName}_split/*vcf.gz" into split_vcfs_ch mode flatten
	
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
	label 'bgzip'
	publishDir "$params.outdir/10_VCFtoolsSiteFilteredVCFs", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	file split_vcf from split_vcfs_ch
	val site_filters from params.vcftools_site_filters
	
	output:
	file "${split_vcf.simpleName}.sitefilt.recode.vcf.gz" into sitefilt_vcf_ch
	
	script:
	if (site_filters == "NULL")
		"""
		bgzip -c $split_vcf > ${split_vcf.simpleName}.sitefilt.recode.vcf.gz
		"""
	else
		"""
		vcftools --gzvcf ${split_vcf} --recode --out ${split_vcf.simpleName}.sitefilt ${site_filters}
		bgzip ${split_vcf.simpleName}.sitefilt.recode.vcf
		"""

}

process gatkFilterSites {

	// Apply GATK-only site filters
	
	label 'gatk'
	label 'bgzip'
	label 'tabix'
	publishDir "$params.outdir/11_GATKSiteFilteredVCFs", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	path refseq from params.refseq
	file "*" from bwa_index_ch
	file site_vcf from sitefilt_vcf_ch
	val site_filters from params.gatk_site_filters
	val gatk from params.gatk
	val gatk_java from params.gatk_java
	
	output:
	file "${site_vcf.simpleName}.gatksitefilt.vcf.gz" into gatk_sitefilt_vcf_ch
	
	script:
	if (site_filters == "NULL")
		"""
		ln -s $site_vcf ${site_vcf.simpleName}.gatksitefilt.vcf.gz
		"""
	else if (gatk_build == 3)
		"""
		tabix $site_vcf
		java ${gatk_java} -jar ${gatk} -T VariantFiltration -V $site_vcf -o tmp.vcf -R $refseq $site_filters
		java ${gatk_java} -jar ${gatk} -T SelectVariants -V tmp.vcf -o ${site_vcf.simpleName}.gatksitefilt.vcf -R $refseq --excludeFiltered
		bgzip ${site_vcf.simpleName}.gatksitefilt.vcf 
		"""
	else if (gatk_build == 4)
		"""
		tabix $site_vcf
		java ${gatk_java} -jar ${gatk} VariantFiltration -R $refseq -V $site_vcf -O tmp.vcf.gz $site_filters
		java ${gatk_java} -jar ${gatk} SelectVariants -R $refseq -V $site_vcf -O ${site_vcf.simpleName}.gatksitefilt.vcf.gz --exclude-filtered
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
	publishDir "$params.outdir/12_RegionFilteredVCFs", mode: 'copy'
	errorStrategy 'retry'
	maxRetries 3
	
	input:
	path site_vcf from gatk_sitefilt_vcf_ch
	file exclude_bed from exclude_bed_ch
	
	output:
	file "${site_vcf.simpleName}.regionfilt.vcf.gz" into regionfilt_vcf_ch
	
	script:
	chr = site_vcf.simpleName.split('_chr')[1]
	if (task.attempt == 1)
		"""
		bedtools intersect -a ${site_vcf} -b ${exclude_bed} -v -header | gzip > ${site_vcf.simpleName}.regionfilt.vcf.gz
		"""
	else if (task.attempt == 2)
		"""
		zcat ${site_vcf} | bedtools intersect -a stdin -b ${exclude_bed} -v -header | gzip > ${site_vcf.simpleName}.regionfilt.vcf.gz
		"""
	else if (task.attempt == 3)
		"""
		grep ${chr} ${exclude_bed} > tmp.bed 
		tabix ${site_vcf}
		bcftools view -R tmp.bed -Ob -o tmp.bcf ${site_vcf}
		tabix tmp.bcf
		# bcftools isec gives the target sites not included in the bed
		bcftools isec -C -O v -o ${site_vcf.simpleName}.targets ${site_vcf} tmp.bcf
		# Use the isec output to get the output. Needs to stream (-T) rather than index jump (-R) for efficiency.
		bcftools view -T ${site_vcf.simpleName}.targets -Ov ${site_vcf} | gzip > ${site_vcf.simpleName}.regionfilt.vcf.gz
		"""
	else
		"""
		grep ${chr} ${exclude_bed} > tmp.bed 
		vcftools --gzvcf ${site_vcf} --recode --out ${site_vcf.simpleName}.regionfilt --exclude-bed tmp.bed
		gzip -c ${site_vcf.simpleName}.regionfilt.recode.vcf > ${site_vcf.simpleName}.regionfilt.vcf.gz
		"""

}

process calcDNMRate {

	// Calculate de novo mutations using calc_denovo_mutation_rate
	
	label 'ruby'
	publishDir "$params.outdir/13_SplitCalcDNMLogs", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	file splitvcf from regionfilt_vcf_ch
	val sire from params.sire
	val dam from params.dam
	val dnm_opts from params.dnm_opts
	
	output:
	file "${splitvcf.simpleName}.log" into split_logs_ch
	
	"""
	calc_denovo_mutation_rate.rb -i ${splitvcf} -s ${sire} -d ${dam} ${dnm_opts} > ${splitvcf.simpleName}.log
	"""

}

process summarizeDNM {

	// Calculate genome-wide DNM rate using summarize_denovo
	
	label 'ruby'
	publishDir "$params.outdir/14_SummarizeDNMLogs", mode: 'copy'
	errorStrategy 'finish'
	
	input:
	file "*" from split_logs_ch.collect()
	
	output:
	file "${params.prefix}*_summary.log"
	
	"""
	for file in ${params.prefix}*.log; do 
		if [ ! -d \${file%_chr*.log} ]; then
			mkdir \${file%_chr*.log}
		fi
		mv \$file \${file%_chr*.log}/
	done
	for outdir in ${params.prefix}*; do summarize_denovo.rb \$outdir > \${outdir}_summary.log; done
	"""

}
