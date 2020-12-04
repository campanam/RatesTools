#!/usr/bin/env nextflow

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
	
	input:
	path refseq from params.refseq
	tuple val(pair_id), path(reads) from readpairs_ch
	val bwa_threads from params.bwa_threads
	file "*" from bwa_index_ch
	
	output:
	file "${pair_id}_${refseq.baseName}.srt.bam" into sorted_bam_ch
	val pair_id into (sample_ch,filt_sample_ch) // sample name for downstream use
	
	"""
	bwa mem -t ${bwa_threads} ${refseq} ${reads} | samtools view -bS - | samtools sort > ${pair_id}_${refseq.baseName}.srt.bam
	"""
	
}

process markDuplicates {

	// Mark duplicates using sambamba or picard
	
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
		java ${picard_java} -jar ${picard} MarkDuplicates I=${sorted_bam} O=${sorted_bam.simpleName}.markdup.bam M=${sorted_bam.simpleName}.markdup.txt
		"""
		
}

process fixReadGroups {

	// Fix read groups using picard
	
	input:
	path markdup_bam from markdup_bam_ch
	path picard from params.picard
	val picard_java from params.picard_java
	val pair_id from sample_ch
	
	output:
	file "${markdup_bam.simpleName}.rg.bam" into rg_bam_ch
	
	"""
	java ${picard_java} -jar ${picard} AddOrReplaceReadGroups I=${markdup_bam} O=${markdup_bam.simpleName}.rg.bam RGID=$pair_id RGLB=$pair_id RGPL=illumina RGPU=$pair_id RGSM=$pair_id
	"""

}

process realignIndels {

	// GATK RealignerTargetCreator. Index bam with picard.
	
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

	"""
	java ${picard_java} -jar ${picard} BuildBamIndex I=${rg_bam}
	java ${gatk_java} -jar ${gatk} -T RealignerTargetCreator -R ${refseq} -I ${rg_bam} -o ${rg_bam.baseName}.intervals
	java ${gatk_java} -jar ${gatk} -T IndelRealigner -R ${refseq} --filter_bases_not_stored -I ${rg_bam} -targetIntervals ${rg_bam.baseName}.intervals -o ${rg_bam.simpleName}.realn.bam
	"""

}

process filterBAMs {

	// Filter BAMs using GATK
	
	input:
	path refseq from params.refseq
	path realn_bam from realn_bam_ch
	val gatk from params.gatk
	val gatk_java from params.gatk_java
	val gatk_nct from params.gatk_nct
	file "*" from bwa_index_ch
	
	output:
	file "${realn_bam.simpleName}.filt.bam" into filt_bam_ch
	
	"""
	java ${gatk_java} -jar ${gatk} -R ${refseq} -T PrintReads -I ${realn_bam} -o ${realn_bam.simpleName}.filt.bam -nct ${gatk_nct} --read_filter BadCigar --read_filter DuplicateRead --read_filter FailsVendorQualityCheck --read_filter HCMappingQuality --read_filter MappingQualityUnavailable --read_filter NotPrimaryAlignment --read_filter UnmappedRead --filter_bases_not_stored --filter_mismatching_base_and_quals
	"""
	
}

process fixMate {

	// Fix mate information using Picard
	publishDir "$params.outdir/FinalBAMs"
		
	input:
	path filt_bam from filt_bam_ch
	path picard from params.picard
	val picard_java from params.picard_java
	
	output:
	file "${filt_bam.simpleName}.fix.bam" into fix_bam_ch
	
	"""
	java ${picard_java} -jar ${picard} FixMateInformation I=${filt_bam} O=${filt_bam.simpleName}.fix.bam ADD_MATE_CIGAR=true
	"""
	
}

process callVariants {

	// Index final bam and call variants using GATK
	publishDir "$params.outdir/gVCFs"
		
	input:
	path refseq from params.refseq
	path fix_bam from fix_bam_ch
	path picard from params.picard
	val picard_java from params.picard_java
	val gatk from params.gatk
	val gatk_java from params.gatk_java
	val gatk_nct from params.gatk_nct
	file "*" from bwa_index_ch
	
	output:
	file "${fix_bam.simpleName}.g.vcf.*" into var_vcf_ch 
	
	"""
	java ${picard_java} -jar ${picard} BuildBamIndex I=${fix_bam}
	java ${gatk_java} -jar ${gatk} -T HaplotypeCaller -nct ${gatk_nct} -R ${refseq} -A DepthPerSampleHC -A Coverage -A HaplotypeScore -A StrandAlleleCountsBySample -I ${fix_bam} -o ${fix_bam.simpleName}.g.vcf.gz -ERC GVCF -out_mode EMIT_ALL_SITES
	"""

}

process genotypegVCFs {

	// Joint genotype gVCFs using GATK
	// Uncompressed combined VCF because GATK 3.8-1 inflate/deflate is glitched for this tool
	
	publishDir "$params.outdir/CombinedVCF"
	
	input:
	path refseq from params.refseq
	file "*" from bwa_index_ch
	file "*" from var_vcf_ch.collect()
	val gatk from params.gatk
	val gatk_java from params.gatk_java
	val prefix from params.prefix
	
	output:
	file "${prefix}_combined.vcf*" into combined_vcf_ch 
	file "${prefix}_combined.vcf.gz" into combined_indels_ch
	
	"""
	VARPATH=""
	for file in *.vcf.gz; do VARPATH+=" --variant \$file"; done
	java ${gatk_java} -jar ${gatk} -T GenotypeGVCFs -o ${prefix}_combined.vcf -R ${refseq} --includeNonVariantSites\$VARPATH
	gzip ${prefix}_combined.vcf # Clean up giant VCF before going on with filtering
	"""
}

process genMap {

	// Calculate mappability using genMap and filter using filterGM
	
	input:
	path refseq from params.refseq
	
	output:
	file "${refseq.simpleName}_genmap.1.0.bed" into genmap_ch
	
	"""
	genmap index -F ${refseq} -I ${refseq.simpleName}_index
	genmap map -K 30 -E 2 -I ${refseq.simpleName}_index/ -O ${refseq.simpleName}_genmap -b
	filterGM.rb ${refseq.simpleName}_genmap.bed 1.0 exclude > ${refseq.simpleName}_genmap.1.0.bed
	"""
}

process repeatMask {

	// Mask repeats using RepeatMasker and RepeatModeler
	
	input:
	path refseq from params.refseq
	val rm_species from params.rm_species
	val rm_pa from params.rm_pa
	
	output:
	file "${refseq}.RM.bed" into rm_bed_ch
	
	"""
	repeatmasker -pa ${rm_pa} -gccalc -nolow -species ${rm_species} ${refseq}
	if ( ! test -f ${refseq}.masked); then # Handling for no repeats detected
		cp ${refseq} ${refseq}.masked
	fi
	BuildDatabase -name ${refseq.baseName}-soft ${refseq}.masked
	RepeatModeler -pa ${rm_pa} -database ${refseq.baseName}-soft
	
	# If no output, make dummy values
	if ( ! test -f */consensi.fa.classified); then
		cp ${refseq} ${refseq}.masked.masked
		cp ${refseq}.out ${refseq}.masked.out
	else
		mv */consensi.fa* ./
		repeatmasker -pa ${rm_pa} -gccalc -nolow -lib consensi.fa.classified ${refseq}.masked
	fi
	
	# Convert out file into BED for downstream
	RM2bed.rb ${refseq}.masked.out > ${refseq}.RM.bed
	"""

}

process maskIndels {

	// Mask indel-affected regions using indels2bed
	
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
	
	publishDir "$params.outdir/ExcludedRegions"
	
	input:
	file indel_bed from indels_ch
	file rm_bed from rm_bed_ch
	file gm_bed from genmap_ch
	val prefix from params.prefix
	
	output:
	file "${prefix}_excluded_reduced.bed" into exclude_bed_ch
	
	"""
	cat ${indel_bed} ${rm_bed} ${gm_bed} > ${prefix}_excluded.bed
	simplify_bed.rb ${prefix}_excluded.bed > ${prefix}_excluded_reduced.bed
	"""

}

filt_sample_ch2 = filt_sample_ch.filter { it != params.sire && it != params.dam } // Need new channel after filtering this one to remove dam and sire from offspring lists

process filterSites {

	// Filter sites using VCFtools
	
	publishDir "$params.outdir/SiteFilteredVCFs"
	
	input:
	file "*" from combined_vcf_ch
	val site_filters from params.site_filters
	val dam from params.dam
	val sire from params.sire
	val prefix from params.prefix
	val pair_id from filt_sample_ch2
	
	output:
	file "${prefix}_offspring*.sitefilt.recode.vcf.gz" into sitefilt_vcf_ch
	
	"""
	vcftools --gzvcf *vcf.gz --recode --out ${prefix}_offspring${pair_id}.sitefilt ${site_filters} --indv ${dam} --indv ${sire} --indv ${pair_id}
	gzip ${prefix}_offspring*.sitefilt.recode.vcf
	"""

}

process filterRegions {

	// Filter regions using VCFtools
	
	publishDir "$params.outdir/RegionFilteredVCFs"
	
	input:
	file site_vcf from sitefilt_vcf_ch
	file exclude_bed from exclude_bed_ch
	
	output:
	file "${site_vcf.simpleName}.regionfilt.recode.vcf.gz" into regionfilt_vcf_ch
	
	"""
	vcftools --gzvcf ${site_vcf} --recode --out ${site_vcf.simpleName}.regionfilt --exclude-bed ${exclude_bed}
	gzip ${site_vcf.simpleName}.regionfilt.recode.vcf
	"""

}

process filterChr {
	
	// Optionally include only specific chromsomes
	
	publishDir "$params.outdir/FilterChrVCFs"
	
	input:
	file region_vcf from regionfilt_vcf_ch
	file chrs from chr_file
	
	output:
	file "${region_vcf.simpleName}.chrfilt.recode.vcf.gz" into chrfilt_vcf_ch
	
	script:
	if (chrs.name == "NULL")
		"""
		cp $region_vcf ${region_vcf.simpleName}.chrfilt.recode.vcf.gz
		"""
	else
		"""
		chr_line=`echo '--chr '`; chr_line+=`awk 1 ORS=' --chr ' ${chrs}`; chr_line=`echo \${chr_line% --chr }` # Awkwardly make into a --chr command-list
		vcftools --gzvcf $region_vcf --recode --out ${region_vcf.simpleName}.chrfilt \$chr_line
		gzip ${region_vcf.simpleName}.chrfilt.recode.vcf
		"""
	
}

process splitVCFs {

	// Split VCFs by contig/chromosome/scaffold etc
	
	publishDir "$params.outdir/SplitVCFs"
	
	input:
	file filtvcf from chrfilt_vcf_ch
	
	output:
	file "${filtvcf.simpleName}_split/*vcf.gz" into split_vcfs_ch
	
	"""
	nextflow_split.rb -i ${filtvcf} -o ${filtvcf.simpleName}_split
	cd ${filtvcf.simpleName}_split
	for file in *vcf.gz; do mv \$file ${filtvcf.simpleName}_\${file}; done
	cd ..
	"""
	
}

process calcDNMRate {

	// Calculate de novo mutations using calc_denovo_mutation_rate
	
	publishDir "$params.outdir/SplitCalcDNMLogs"
	
	input:
	file splitvcf from split_vcfs_ch
	val sire from params.sire
	val dam from params.dam
	val dnm_opts from params.dnm_opts
	
	output:
	file "${splitvcf.simpleName}.log" into split_logs_ch
	
	"""
	calc_denovo_mutation_rate.rb -i ${splitvcf} -s ${sire} -d ${dam} > ${splitvcf.simpleName}.log
	"""
}

process summarizeDNM {

	// Calculate genome-wide DNM rate using summarize_denovo
	
	publishDir "$params.outdir/SummarizeDNMLogs"
	
	input:
	file "*" from split_logs_ch.collect()
	
	output:
	file "*_summary.log"
	
	"""
	for file in *.log; do 
		if ( ! \${file%_chr*.log} -d ); then
			mkdir \${file%_chr*.log}
		fi
		mv \$file \${file%_chr*.log}/
	done
	for outdir in *; do summarize_denovo.rb \$outdir > \${outdir}_summary.log; done
	"""

}

/* process seqdictfai {

	cpus 1
	executor 'sge'
	queue 'sThC.q'
	clusterOptions '-l mres=2G,h_data=2G,h_vmem=2G -S /bin/bash'
	module 'bioinformatics/samtools/1.9'
	
	"""
	samtools faidx $refseq
	samtools dict $refseq > ${refseq.getParent()}/$dict
	"""
	#for i in *.g.vcf; do gzip `readlink \$i`; mv \$i \${i}.gz; mv $baseDir/${params.outdir}/gVCFs/\$i $baseDir/${params.outdir}/gVCFs/\${i}.gz; done
} */
