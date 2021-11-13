# RatesTools Pipeline Details  

__Michael G. Campana & Ellie E. Armstrong, 2019-2020__  
Smithsonian Conservation Biology Institute  
Stanford University  

The following document provides detailed descriptions of the steps included in the RatesTools pipeline. Please note that the commands shown below only list non-standard options for clarity. We omit basic input/output required for these commands. See the documentation for these programs for operation details.  

## prepareRef  
The prepareRef process indexes the reference sequence using `bwa index` [1] (and optionally the specified indexing algorithm) and `samtools faidx` [2]. It also generates a sequence dictionary using `samtools dict`.

## alignSeqs  
The alignSeqs process aligns read pairs (in two separate fastq format files) against the indexed reference sequencing using `bwa mem` and converts to bam format using `samtools view` to conserve disk space.  

## sortBAM  
The sortBAM process sorts the bam files from alignSeqs using `samtools sort`. This process was separated from the alignSeqs step to minimize redundant calculations should there be an error necessitating pipeline resumption.  

## markDuplicates  
The markDuplicates process marks PCR duplicates using either `sambamba markdup` [3] or Picard [4] MarkDuplicates (`java -jar picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000`).  

## fixReadGroups  
The fixReadGroups process adds sample read group information to the bam files using Picard AddOrReplaceReadGroups (`java -jar picard.jar AddOrReplaceReadGroups`).  

## realignIndels  
The realignIndels process first builds a bam index (.bai) for the fixReadGroups output bam using Picard BuildBamIndex (`java -jar picard.jar BuildBamIndex`). It then identifies intervals for indel realignment using Genome Analysis Toolkit (GATK) [5] RealignerTargetCreator (`java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator`) and realigns indels using GATK IndelRealigner (`java -jar GenomeAnalysisTk.jar -T IndelRealigner --filter_bases_not_stored`).  

## filterBAMs  
The filterBAM process strictly filters aligned reads using GATK PrintReads following [6] (`java -jar GenomeAnalysisTK.jar -T PrintReads --read_filter BadCigar --read_filter DuplicateRead --read_filter FailsVendorQualityCheck --read_filter HCMappingQuality --read_filter MappingQualityUnavailable --read_filter NotPrimaryAlignment --read_filter UnmappedRead --filter_bases_not_stored --filter_mismatching_base_and_quals`).  

## fixMate  
The fixMate process corrects sequence mate pair tags using Picard FixMateInformation (`java -jar picard.jar FixMateInformation ADD_MATE_CIGAR=true`). It then builds a bam index (.bai) for the output bam file using Picard BuildBamIndex (`java -jar picard.jar BuildBamIndex`).  

## callVariants  
The callVariants process calls individual sequence variants (gVCF format) using GATK HaplotypeCaller (`java -jar GenomeAnalysisTK.jar -T HaplotypeCaller -A DepthPerSampleHC -A Coverage -A HaplotypeScore -A StrandAlleleCountsBySample -ERC GVCF -out_mode EMIT_ALL_SITES`).  

## genotypegVCFs  
The genotypegVCFs process performs joint-genotyping and generates an all-sites VCF for all samples using GATK GenotypeGVCFs (`java -jar GenomeAnalysisTK.jar -T GenotypeGVCFs --includeNonVariantSites`).  It then compresses the resulting multi-sample VCF using `gzip`.  

## genMapIndex  
The genMapIndex process generates the GenMap [7] index for the reference sequence for downstream mappability calculations (`genmap index`).  

## genMapMap  
The genMapMap process calculates the mappability for the reference sequence using GenMap (`genmap map -K 30 -E 2 -b`). It then filters the GenMap results using [`filterGM.rb`](ruby_r_scripts.md#filterGMrb) to exclude any sequences with mappability < 1.0 (`filterGM.rb <raw_genmap.bed> 1.0 exclude > <genmap.1.0.bed>`).  

## repeatMask  
The repeatMask process uses RepeatMasker [8] to soft-mask repeat regions in the reference sequence based on a pre-defined target species (`RepeatMasker -gccalc -nolow -species <specified species>`). The soft-masked reference sequence is passed to the repeatModeler process.  

## repeatModeler  
The repeatModeler process uses RepeatModeler [9] and the RepeatMasker soft-masked reference sequence to generate a species-specific repeat library for the reference sequence. First, a database is built using `BuildDatabase`, then the library is built using `RepeatModeler`. The resulting library (`consensi.fa.classified`) RepeatMasker -pa ${rm_pa} -gccalc -nolow -lib consensi.fa.classified is passed to the repeatMaskRM process.  

## repeatMaskRM  
The repeatMaskRM process uses RepeatMasker and the RepeatModeler custom repeat library to soft mask the reference sequence (`RepeatMasker -gccalc -nolow -lib consensi.fa.classified`). The RepeatMasker out file (.out) is then converted to bed using [`RM2bed.rb`](ruby_r_scripts.md#RM2bedrb).  

## maskIndels  
Using [`indels2bed.rb`](ruby_r_scripts.md#indels2bedrb), the maskIndels process scans the all-sites VCF generated by the genotypeGVCFs process for indels and generates a bed file excluding all sites 5 bp upstream/downstream of the identified indels (`indels2bed.rb <all-sites.vcf> 5 > <indels.bed>`).  

## simplifyBed  
The simplifyBed process identifies and merges overlapping bed entries in the bed results from the repeatMaskRM, genMapMap and maskIndels processes. First, overlapping bed entries are merged for each of the three chromosome-sorted and coordinate-sorted bed files using [`simplify_sorted_bed.rb`](ruby_r_scripts.md#simplify_sorted_bedrb). The three merged bed files are then concatenated using `cat`. Overlapping entries in the resulting, unsorted bed file are merged using [`simplify_bed.rb`](ruby_r_scripts.md#simplify_bedrb).  

## filterChr  
Using VCFtools [10], the filterChr process generates all-sites VCFs for each offspring and its parents (`vcftools --recode --indv <sire> --indv <dam> --indv <offspring>`). The number of resulting VCFs will thus equal the number of offspring in the dataset. If a list of target chromosomes was provided, this process also filters the output VCFs to include only the specified target regions using the `--chr` option. The resulting VCFs are compressed using `gzip`.  

## pullDPGQ
Using BCFtools [11], the pullGQDP process extracts the DP and GQ information from the chromosome-filtered VCFs (`bcftools view -v snps <chromosome_vcf> | bcftools query -f "%CHROM %POS [ %DP] [ %GQ]\n"`).  

## plotDPGQ
The plotDPGQ process plots the DP and GQ information from pullDPGQ using [`plotDPGQ.R`](ruby_r_scripts.md#plotDPGQR).  

## splitVCFs  
The splitVCFs process splits each VCF generated during the filterChr process by chromosome/contig name for parallelization of downstream processes using [`nextflow_split.rb`](ruby_r_scripts.md#nextflow_splitrb). The resulting VCFs are compressed using `gzip`.  

## filterSites  
The filterSites process filters the split VCF files from the splitVCFs process using VCFtools and the site filters provided in the config file (`vcftools --recode <site_filters>`). The resulting VCFs are compressed using `bgzip` [11].  

## filterRegions  
The filterRegions process removes low-reliability regions (repeat regions, indel-affected sites, and regions of non-unique mappability) from the site-filtered VCFs from the filterSites process. The low-reliability regions are specified in the output merged bed from the simplifyBed process. On the first pass, the process uses BEDTools [12] intersect on the site-filtered VCF (`bedtools intersect -a <site-filtered_vcf> -b <out_cat.bed>`). Should this fail (e.g. due to a compression issue), the process repeats using `zcat` to uncompress the VCF and pass it via stdin to BEDTools. If neither BEDTools run is successful, the process attempts to removed the filtered regions using BCFtools view (`bcftools view -R <out_cat.bed> -Ob -o tmp.bcf`), BCFtools isec (`bcftools isec -C <site-filtered_vcf> <tmp.bcf>`) and BCFtools view (`bcftools view -T <isec_out> <site-filtered_vcf>`) [11]. If the BCFtools attempt fails, the process defaults to using VCFtools (`vcftools --recode --exclude-bed <out_cat.bed>`). The resulting VCFs are compressed using `gzip`.  

## calcDNMRate 
Using the region-filtered VCFs output from the filterRegions process, the calcDNMRate process calculates the per-chromosome mutation rate using [`calc_denovo_mutation_rate.rb`](ruby_r_scripts.md#calc_denovo_mutation_raterb) and the options specified in the config file (`calc_denovo_mutation_rate.rb -i <chr.vcf> -s <sire> -d <dam> <denovo_mutation_options> > chr.log)`).  

## summarizeDNM  
Using [`summarize_denovo.rb`](ruby_r_scripts.md#summarize_denovorb), the summarizeDNM process combines the per-chromosome mutation rate results from the calcDNMRate process to obtain the genomic mutation rate.  

## References  
1. Li, H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv*, [1303.3997v2](https://arxiv.org/abs/1303.3997).  
2. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., 1000 Genome Project Data Processing Subgroup (2009) The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25, 2078-2079. DOI: [10.1093/bioinformatics/btp352](https://academic.oup.com/bioinformatics/article/25/16/2078/204688).  
3. Tarasov, A., Vilella, A.J., Cuppen, E., Nijman, I.J., Prins, P. (2015) Sambamba: fast processing of NGS alignment formats. *Bioinformatics*, __31__, 2032–2034. DOI: [10.1093/bioinformatics/btv098](https://academic.oup.com/bioinformatics/article/31/12/2032/214758).  
4. Broad Institute (2020). Picard v. 2.23.8 (https://broadinstitute.github.io/picard/).  
5. McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., DePristo, M.A. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res*, __20__, 1297-1303. DOI: [10.1101/gr.107524.110](https://genome.cshlp.org/content/20/9/1297.abstract).  
6. Besenbacher, S., Hvilsom, C., Marques-Bonet, T., Mailund, T., Schierup, M.H. (2019) Direct estimation of mutations in great apes reconciles phylogenetic dating. *Nat Ecol Evol*, __3__, 286-292. DOI: [10.1038/s41559-018-0778-x](https://www.nature.com/articles/s41559-018-0778-x).  
7. Pockrandt, C., Alzamel, M., Iliopoulos, C.S., Reinert, K. (2020) GenMap: ultra-fast computation of genome mappability. *Bioinformatics*, __36__, 3687–3692, doi: [10.1093/bioinformatics/btaa222](https://academic.oup.com/bioinformatics/article/36/12/3687/5815974?login=true).  
8. Smit, A.F.A., Hubley, R., Green, P. (2013-2015) *RepeatMasker Open-4.0*. (http://www.repeatmasker.org).  
9. Flynn, J.M., Hubley, R., Goubert, C., Rosen, J. Clark,. A.G., Feschotte, C., Smit, A.F. (2020) RepeatModeler2 for automated genomic discovery of transposable element families. *Proc Natl Acad Sci U S A*, __117__, 9451-9457. DOI: [10.1073/pnas.1921046117](https://www.pnas.org/content/117/17/9451.short).  
10. Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., McVean, G., Durbin, R. (2011) The variant call format and VCFtools. *Bioinformatics*, __27__, 2156–2158. DOI: [10.1093/bioinformatics/btr330](https://academic.oup.com/bioinformatics/article/27/15/2156/402296).  
11. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H. (2021) Twelves years of SAMtools and BCFtools. *GigaScience*, __10__, giab008. DOI: [10.1093/gigascience/giab008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722).  
12. Quinlan, A.R., Hall, I.M. (2010) BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, __26__, 841-842, doi: [10.1093/bioinformatics/btq0333](https://academic.oup.com/bioinformatics/article/26/6/841/244688).  
