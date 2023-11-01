#!/usr/bin/env nextflow

/* RatesTools version 1.0.0
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

// algorithm for BWA index for prepareRef

if (params.bwa_alg == "") { bwa_alg = "" } else { bwa_alg = "-a " + params.bwa_alg + " " }
// Generate Picard and GATK executable commands
if (params.picard_conda) { picard = "picard " + params.picard_java} else { picard = "java " + params.picard_java + " -jar " + params.picard }
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
		
	input:
	path refseq
	
	output:
	path "${refseq.baseName}*.{amb,ann,bwt,pac,sa,fai,dict}"
	
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
		
	input:
	tuple val(sample), val(pair_id), path(reads1), path(reads2), val(rg)
	path refseq
	path "*"
	
	output:
	tuple path("${pair_id}_${refseq.simpleName}.bam"), val(sample)
	
	script:
	samtools_extra_threads = task.cpus - 1
	"""
	bwa mem -t ${task.cpus} -R '${rg}' ${refseq} ${reads1} ${reads2} | samtools fixmate -@ ${samtools_extra_threads} -m - - | samtools sort -@ ${samtools_extra_threads} -o ${pair_id}_${refseq.simpleName}.bam - 
	"""
	
}

process markDuplicates {

	// Mark duplicates using sambamba, samtools or picard
	
	label 'sambamba'
	label 'samtools'
	label 'picard'
		
	input:
	tuple path(sorted_bam), val(sample)
	
	output:
	tuple path("${sorted_bam.simpleName}.markdup.bam"), val(sample)
	
	script:
	samtools_extra_threads = task.cpus - 1
	if ( params.markDuplicates == "sambamba" )
		"""
		sambamba markdup ${sorted_bam} ${sorted_bam.simpleName}.markdup.bam
		"""
	else if ( params.markDuplicates == "samtools" )
		"""
		samtools markdup -@ ${samtools_extra_threads} ${sorted_bam} ${sorted_bam.simpleName}.markdup.bam
		"""
	else
		"""
		$picard MarkDuplicates I=${sorted_bam} O=${sorted_bam.simpleName}.markdup.bam M=${sorted_bam.simpleName}.markdup.txt MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000
		"""
		
}

process mergeLibraries {

	// Merge libraries by their sample IDs using SAMtools merge
	
	input:
	tuple path(bam), val(sample)
	
	output:
	path "${sample}.merged.bam"
	
	script:
	samtools_extra_threads = task.cpus - 1
	bams = 0
	bamlist = ""
	for (i in bam) {
		bams++
		bamlist = bamlist + " " + i
	}
	if (bams == 1) // Skip merging single libraries
		"""
		ln -s $bamlist ${sample}.merged.bam
		"""
	else
		"""
		samtools merge -@ ${samtools_extra_threads} -o ${sample}.merged.bam $bamlist
		"""
} 

process realignIndels {

	// GATK RealignerTargetCreator. Index bam with picard.
	
	label 'gatk'
	label 'picard'
	
	if (params.filter_bams == false) {
		publishDir  "$params.outdir/01_FinalBAMs", mode: 'copy'
	}
		
	input:
	path rg_bam
	path refseq
	path "*"
	
	output:
	tuple path("${rg_bam.simpleName}.realn.bam"), path("${rg_bam.simpleName}.realn.bai")
	
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
		
	input:
	tuple path(realn_bam), path(realn_bai)
	path refseq
	path "*"
	
	output:
	tuple path("${realn_bam.simpleName}.filt.bam"), path("${realn_bam.simpleName}.filt.bai")
	
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
			
	input:
	tuple path(filt_bam), path(filt_bai)
	
	output:
	tuple path("${filt_bam.simpleName}.fix.bam"), path("${filt_bam.simpleName}.fix.bai")
	
	"""
	$picard FixMateInformation I=${filt_bam} O=${filt_bam.simpleName}.fix.bam ADD_MATE_CIGAR=true
	$picard BuildBamIndex I=${filt_bam.simpleName}.fix.bam
	"""
	
}

process callVariants {

	// Index final bam and call variants using GATK
	
	label 'gatk'
	publishDir "$params.outdir/02_gVCFs", mode: 'copy'
			
	input:
	tuple path(fix_bam), path(fix_bai)
	path refseq
	path "*"
	
	output:
	path "${fix_bam.simpleName}.g.vcf.*"
	
	script:
	if (params.gatk_build == 3)
		"""
		$gatk -T HaplotypeCaller -nct ${task.cpus} -R ${refseq} -A DepthPerSampleHC -A Coverage -A HaplotypeScore -A StrandAlleleCountsBySample -I ${fix_bam} -o ${fix_bam.simpleName}.g.vcf.gz -ERC GVCF -out_mode EMIT_ALL_SITES
		"""
	else if (params.gatk_build == 4)
		"""
		$gatk HaplotypeCaller -R ${refseq} -I $fix_bam -O ${fix_bam.simpleName}.g.vcf.gz -ERC GVCF -G StandardAnnotation -G AS_StandardAnnotation
		"""

}

process genotypegVCFs {

	// Joint genotype gVCFs using GATK
	// Uncompressed combined VCF because GATK 3.8-1 inflate/deflate is glitched for this tool
	// Also uncompressed combined VCF because GATK4 cannot generate index with gzipped output
	
	label 'gatk'
	label 'gzip'
	publishDir "$params.outdir/03_CombinedVCF", mode: 'copy'
		
	input:
	path "*"
	path refseq
	path "*"
	
	output:
	path "${params.prefix}_combined.vcf.gz"
	
	script:
	if (params.gatk_build == 3)
		"""
		VARPATH=""
		for file in *.vcf.gz; do VARPATH+=" --variant \$file"; done
		$gatk -T GenotypeGVCFs -R ${refseq} --includeNonVariantSites\$VARPATH -o ${params.prefix}_combined.vcf
		gzip ${params.prefix}_combined.vcf
		"""
	else if (params.gatk_build == 4)
		"""
		VARPATH=""
		for file in *.vcf.gz; do VARPATH+=" --variant \$file"; done
		$gatk CombineGVCFs -R ${refseq} -O tmp.g.vcf.gz --convert-to-base-pair-resolution\$VARPATH
		$gatk GenotypeGVCFs -R ${refseq} --include-non-variant-sites -V tmp.g.vcf.gz -O ${params.prefix}_combined.vcf
		gzip ${params.prefix}_combined.vcf
		"""
}

process genMapIndex {

	// Generate GenMap index
	
	label 'genmap'
		
	input:
	path refseq
	val gm_tmpdir
	
	output:
	tuple path("$refseq"), path("${refseq.simpleName}_index"), path("${refseq.simpleName}_index/*")
	
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
		
	input:
	tuple path(refseq), path(genmap_index), path("*")
	
	output:
	path "${refseq.simpleName}_genmap.1.0.bed"
	
	"""
	genmap map ${params.gm_opts} -T ${task.cpus} -I ${refseq.simpleName}_index/ -O ${refseq.simpleName}_genmap -b
	filterGM.rb ${refseq.simpleName}_genmap.bed 1.0 exclude > ${refseq.simpleName}_genmap.1.0.bed
	"""
}

 
process repeatMask {

	// Mask repeats using RepeatMasker and default RM libraries
	
	label 'repeatmasker'
		publishDir "$params.outdir/04_RepeatMasking", mode: 'copy'
	
	input:
	path refseq
	val rm_species
	
	output:
	path("${refseq}.masked"), emit: rm1
	path("${refseq}.out"), emit: rm1_out
	path "${refseq}.tbl"
	
	"""
	RepeatMasker -pa ${task.cpus} ${params.rm_mask_opts} -species ${rm_species} ${refseq}
	if [ ! -f ${refseq}.masked ]; then # Handling for no repeats detected
		ln -s ${refseq} ${refseq}.masked
	fi
	"""

}

process repeatModeler {

	// RepeatModeler on soft-masked reference
	// Requires RepeatModeler 2.0.5
	
	label 'repeatmodeler'
		publishDir "$params.outdir/04_RepeatMasking", mode: 'copy'
	
	input:
	path(refseq_masked)
	
	output:
	path "**consensi.fa.classified"
	
	"""
	BuildDatabase -name ${refseq_masked.baseName}-soft ${refseq_masked}
	RepeatModeler -threads ${task.cpus} ${params.rm_model_opts} -database ${refseq_masked.baseName}-soft
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
		publishDir "$params.outdir/04_RepeatMasking", mode: 'copy'
	
	input:
	path refseq_masked
	path rm_out
	path rm_lib
	
	output:
	path "${refseq_masked}.masked.*", optional: true
	path "${refseq_masked.simpleName}.RM.bed", emit: RMbed
	
	"""
	if [ "\$(wc -l < consensi.fa.classified)" -eq 0 ]; then
	# If no output from RepeatModeler, use original RepeatMasker results
		RM2bed.rb ${rm_out} > ${refseq_masked.simpleName}.RM.bed
	else
		RepeatMasker -pa ${task.cpus} ${params.rm_mask_opts} -lib consensi.fa.classified ${refseq_masked}
		# Convert out file into BED for downstream
		RM2bed.rb ${rm_out} > tmp.out
		RM2bed.rb ${refseq_masked}.out > tmp2.out
		cat tmp.out tmp2.out | sort  -k1,1 -k2,2n > ${refseq_masked.simpleName}.RM.bed
	fi
	"""

}

process maskIndels {

	// Mask indel-affected regions using indels2bed
	
	label 'ruby'
		
	input:
	path combo_vcf
	
	output:
	path "${combo_vcf.simpleName}_indels.bed"
	
	"""
	indels2bed.rb ${combo_vcf} ${params.indelpad} > ${combo_vcf.simpleName}_indels.bed
	"""

}

process simplifyBed {

	// Reduce number of unique BED entries using simplify_bed
	
	label 'bedtools'
		publishDir "$params.outdir/05_ExcludedRegions", mode: 'copy'
	
	input:
	path indel_bed
	path rm_bed
	path gm_bed
	
	output:
	path "${prefix}_excluded_reduced.bed", optional: true
	
	"""
	#!/usr/bin/env bash
	cat ${indel_bed} ${rm_bed} ${gm_bed} | sort -k1,1 -k2,2n > tmp.bed
	if [ ! "\$(wc -l < tmp.bed)" -eq 0 ]; then
		bedtools merge -i tmp.bed > ${params.prefix}_excluded_reduced.bed
	else
		touch ${params.prefix}_excluded_reduced.bed
	fi
	"""

}

process filterChr {
	
	// Optionally include only specific chromsomes
	
	label 'vcftools'
	label 'bcftools'
	label  'gzip'
	publishDir "$params.outdir/06_FilterChrVCFs", mode: 'copy', pattern: '*chrfilt.vcf.gz'
		
	input:
	path comb_vcf
	path chrs
	
	output:
	path "${params.prefix}.chrfilt.vcf.gz", emit: chr_vcf
	path 'chrfilt.tmp', emit: chr_tmp

	"""
	chr_line=`echo '--chr '`; chr_line+=`awk 1 ORS=' --chr ' ${chrs}`; chr_line=`echo \${chr_line% --chr }` # Awkwardly make into a --chr command-list
	vcftools --gzvcf $comb_vcf --recode -c \$chr_line | gzip > ${params.prefix}.chrfilt.vcf.gz
	cp .command.log chrfilt.tmp
	"""

}

process splitTrios {
	
	// Split samples into trios for analysis
	
	label 'vcftools'
	label 'gzip'
	label 'bcftools'
	publishDir "$params.outdir/07_SplitTrioVCFs", mode: 'copy', pattern: "${params.prefix}_offspring*.vcf.gz"
		
	input:
	tuple path(chr_vcf), val(sample_id)
	
	output:
	path "${params.prefix}_offspring${sample_id}.vcf.gz", emit: trio_vcf
	path "${params.prefix}_offspring${sample_id}_trio.tmp", emit: trio_tmp
	
	"""
	vcftools --gzvcf $chr_vcf --recode -c --indv ${params.dam} --indv ${params.sire} --indv ${sample_id} | gzip > ${params.prefix}_offspring${sample_id}.vcf.gz
	cp .command.log ${params.prefix}_offspring${sample_id}_trio.tmp
	"""

}

process pullDPGQ {

	// Extract DP/GQ values from autosome scaffolds to look at the distributions of DP and GQ for variant call sites

	label 'bcftools'
	publishDir "$params.outdir/08_gVCFs_DP_GQ", mode: 'copy'
		
	input:
	tuple path(chrfilt), val(sample_id)
	
	output:
	path "${chrfilt.simpleName}_ind${sample_id}.variants.txt"
	
	"""
	bcftools view -v snps ${chrfilt} -s ${sample_id} | bcftools query -f \"%CHROM %POS [ %DP] [ %GQ]\\n\" -o ${chrfilt.simpleName}_ind${sample_id}.variants.txt
	"""

}

process plotDPGQ {

	// Plot DP and GQ distributions
	
	label 'R'
	publishDir "$params.outdir/08_gVCFs_DP_GQ", mode: 'copy'
		
	input:
	path "*.txt"
	
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
		
	input:
	path chr_vcf
	
	output:
	path "${chr_vcf.simpleName}_split/*vcf.gz"
	
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
		
	input:
	path split_vcf
	
	output:
	tuple path("${split_vcf.simpleName}.sitefilt.tmp"), path(split_vcf), path("${split_vcf.simpleName}.sitefilt.vcf.gz")
	
	script:
	if (site_filters == "NULL")
		"""
		cp -P $split_vcf ${split_vcf.simpleName}.sitefilt.recode.vcf.gz
		vcftools --gzvcf $split_vcf
		cp .command.log ${split_vcf.simpleName}.sitefilt.tmp
		"""
	else
		"""
		vcftools --gzvcf ${split_vcf} --recode -c ${params.vcftools_site_filters} | bgzip > ${split_vcf.simpleName}.sitefilt.vcf.gz
		cp .command.log ${split_vcf.simpleName}.sitefilt.tmp
		"""

}

process gatkFilterSites {

	// Apply GATK-only site filters
	
	label 'gatk'
	label 'tabix'
	label 'vcftools'
	label 'bcftools'
	publishDir "$params.outdir/11_GATKSiteFilteredVCFs", mode: 'copy', pattern: '*gatksitefilt.vcf.gz'
		
	input:
	path site_vcf
	path "*"
	
	output:
	tuple path("${site_vcf.simpleName}.gatksitefilt.tmp"), path(site_vcf), path("${site_vcf.simpleName}.gatksitefilt.vcf.gz")
	
	script:
	if (site_filters == "NULL")
		"""
		ln -s $site_vcf ${site_vcf.simpleName}.gatksitefilt.vcf.gz
		vcftools --gzvcf $site_vcf
		cp .command.log  ${site_vcf.simpleName}.gatksitefilt.tmp
		"""
	else if (params.gatk_build == 3)
		"""
		tabix $site_vcf
		$gatk -T VariantFiltration -V $site_vcf -o tmp.vcf -R ${params.refseq} ${params.gatk_site_filters}
		$gatk -T SelectVariants -V tmp.vcf -o ${site_vcf.simpleName}.gatksitefilt.vcf.gz -R ${params.refseq} --excludeFiltered
		vcftools --gzvcf ${site_vcf.simpleName}.gatksitefilt.vcf.gz
		tail .command.log > ${site_vcf.simpleName}.gatksitefilt.tmp
		rm tmp.vcf
		"""
	else if (params.gatk_build == 4)
		"""
		tabix $site_vcf
		$gatk VariantFiltration -R ${params.refseq} -V $site_vcf -O tmp.vcf.gz ${params.gatk_site_filters}
		$gatk SelectVariants -R ${params.refseq} -V tmp.vcf.gz -O ${site_vcf.simpleName}.gatksitefilt.vcf.gz --exclude-filtered
		vcftools --gzvcf ${site_vcf.simpleName}.gatksitefilt.vcf.gz
		tail .command.log > ${site_vcf.simpleName}.gatksitefilt.tmp
		rm tmp.vcf.gz
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
	
	input:
	path site_vcf
	path exclude_bed from exclude_bed_ch
	
	output:
	tuple path("${site_vcf.simpleName}_regionfilt.tmp"), path(site_vcf), path("${site_vcf.simpleName}.regionfilt.vcf.gz")
	
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

process calcDNMRate {

	// Calculate de novo mutations using calc_denovo_mutation_rate
	
	label 'ruby'
	publishDir "$params.outdir/13_SplitCalcDNMLogs", mode: 'copy'
		
	input:
	path splitvcf
	
	output:
	path "${splitvcf.simpleName}.log"
	
	"""
	calc_denovo_mutation_rate.rb -i ${splitvcf} -s ${params.sire} -d ${params.dam} ${param.dnm_opts} > ${splitvcf.simpleName}.log
	"""

}

process summarizeDNM {

	// Calculate genome-wide DNM rate using summarize_denovo
	
	label 'ruby'
	label 'bcftools'
	label 'gzip'
	publishDir "$params.outdir/14_SummarizeDNMLogs", mode: 'copy'
		
	input:
	path "*"
	path "*"
	
	output:
	path "${params.prefix}*_summary.log", emit: log
	path "${params.prefix}*_candidates.vcf.gz", emit: vcf
	
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

process sanityCheckLogs {

	// Sanity check filtering logs and remove too short contigs as needed

	label 'gzip'
	
	input:
	path logfile
	path allvcflog
	path filtvcflog
	val min_contig_length
	val min_filt_contig_length
	
	output:
	path "${logfile.simpleName}.log",  emit: log
	path "${filtvcflog.baseName}.OK.vcf.gz", optional: true, emit: ok_vcf
	
	"""
	logstats.sh $logfile $allvcflog $filtvcflog $min_contig_length $min_filt_contig_length  > ${logfile.simpleName}.log
	"""
	
}

process generateSummaryStats {

	label 'ruby'
	publishDir "$params.outdir/15_SummaryStats", mode: "copy"
		
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

workflow logSanityTrio {
	// Sanity Check logs from trio splitting
	take:
		tmpfile
		rawvcf
		filtvcf
	main:
		sanityCheckLogs(tmpfile, rawvcf, filtvcf, 0, 0)
	emit:
		sanelog = sanityCheckLogs.out.log
		ok_vcf = sanityCheckLogs.out.ok_vcf
}

workflow logVcftoolsSanity {
	// Sanity check logs from VCFtools site filtering
	take:
		data
	main:
		sanityCheckLogs(data[0], data[1], data[2], params.min_contig_length, params.min_filt_contig_length)
	emit:
		sanelog = sanityCheckLogs.out.log
		ok_vcf = sanityCheckLogs.out.ok_vcf
}

workflow logGatkSanity {
	// Sanity check logs for GATK site filtering and remove too short contigs
	// Dummy value of 1 for min_contig_length since already evalutated and no longer accurate
	take:
		data
	main:
		sanityCheckLogs(data[0], data[1], data[2], 1, params.min_filt_contig_length)
	emit:
		sanelog = sanityCheckLogs.out.log
		ok_vcf = sanityCheckLogs.out.ok_vcf
}

workflow logRegionSanity {
	// Sanity check logs for region filtering and remove too short contigs
	// Dummy value of 1 for min_contig_length since already evalutated and no longer accurate
	take:
		data
	main:
		sanityCheckLogs(data[0], data[1], data[2], 1, params.min_filt_contig_length)
	emit:
		sanelog = sanityCheckLogs.out.log
		ok_vcf = sanityCheckLogs.out.ok_vcf
}

workflow {
	main:
		prepareRef(params.refseq)
		read_data = Channel.fromPath(params.libraries).splitCsv(header:true).map { row -> tuple(row.Sample, row.Library, file(params.readDir + row.Read1), file(params.readDir + row.Read2), '@RG\\tID:' + row.Library + '\\tSM:' + row.Sample + '\\tLB:ILLUMINA\\tPL:ILLUMINA') }
		alignSeqs(read_data, params.refseq, prepareRef.out) | markDuplicates
		mergeLibraries(markDuplicates.out.groupTuple(by: 1)) // Need unique samples matched with their file paths
		realignIndels(mergeLibraries.out, params.refseq, prepareRef.out)
		if (params.filter_bams) {
			filterBAMs(realignIndels.out, params.refseq, prepareRef.out) | fixMate
			callVariants(fixMate.out, params.refseq, prepareRef.out)
		} else {
			callVariants(realignIndels.out, params.refseq, prepareRef.out)
		}
		genotypegVCFs(callVariants.out.collect(), params.refseq, prepareRef.out) | maskIndels
		if (params.region_filter) {
			genMapIndex(params.refseq, params.gm_tmpdir) | genMapMap
			repeatMask(params.refseq, params.rm_species)
			repeatModeler(repeatMask.out.rm1)
			repeatMaskRM(repeatMask.out.rm1, repeatMask.out.rm1_out, repeatModeler.out)
			simplifyBed(genMapMap.out, maskIndels.out, repeatMaskRM.out.RMbed)
		}
		all_samples = read_data.map {it -> it[0]}.unique()
		trio_samples = all_samples.filter { it != params.sire && it != params.dam } // Need new channel after filtering this one to remove dam and sire from offspring lists
		if ( params.chr_file != 'NULL') {
			filterChr(genotypegVCFs.out, channel.fromPath(params.chr_file))
			sanityCheckLogs(filterChr.out.chr_tmp, genotypegVCFs.out, filterChr.out.chr_vcf, 0, 0)
			trio = filterChr.out.chr_vcf.combine(trio_samples)
			splitTrios(trio)
			logSanityTrio(splitTrios.out.trio_tmp, filterChr.out.chr_vcf, splitTrios.out.trio_vcf)
			log_trio_sanity = sanityCheckLogs.out.log.mix(logSanityTrio.out.sanelog)
			allDPGQ = filterChr.out.chr_vcf.combine(all_samples)			
		} else {
			trio = genotypegVCFs.out.combine(trio_samples)
			splitTrios(trio)
			logSanityTrio(splitTrios.out.trio_tmp, genotypegVCFs.out, splitTrios.out.trio_vcf)
			log_trio_sanity = logSanityTrio.out.sanelog
			allDPGQ = genotypegVCFs.out.combine(all_samples)
		}
		pullDPGQ(allDPGQ)
		plotDPGQ(pullDPGQ.out.collect())
		splitVCFs(splitTrios.out.trio_vcf) }/* | flatten | vcftoolsFilterSites | logVcftoolsSanity

		gatkFilterSites(logVcftoolsSanity.out.ok_vcf, prepareRef.out) | logGatkSanity
		if (params.region_filter) { 
			filterRegions(logGatkSanity.out.ok_vcf) | logRegionSanity
			calcDNMRate(logRegionSanity.out.ok_vcf)
			summarizeDNM(calcDNMRate.out.collect(),splitTrios.out.trio_vcf.collect())
			all_logs_sanity = log_trio_sanity.mix(logRegionSanity.out.sanelog, logGatkSanity.out.sanelog, logVcftoolsSanity.out.sanelog, summarizeDNM.out.log).collect()
		} else {
			calcDNMRate(logGatkSanity.out.ok_vcf)
			summarizeDNM(calcDNMRate.out.collect(),splitTrios.out.trio_vcf.collect())
			all_logs_sanity = log_trio_sanity.mix(logGatkSanity.out.sanelog, logVcftoolsSanity.out.sanelog, summarizeDNM.out.log).collect()
		}
		generateSummaryStats(all_logs_sanity)
} */