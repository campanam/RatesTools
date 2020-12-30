# RatesTools Pipeline Details  

__Michael G. Campana & Ellie E. Armstrong, 2019-2020__  
Smithsonian Conservation Biology Institute  
Stanford University  

The following document provides detailed descriptions of the steps included in the RatesTools pipeline. Please note that the commands shown below only list non-standard options for clarity. We omit basic input/output required for these commands. See the documentation for these programs for operation details.  

## prepareRef  
The prepareRef process indexes the reference sequence using `bwa index` (and optionally the specified indexing algorithm) and `samtools faidx`. It also generates a sequence dictionary using `samtools dict`.

## alignSeqs  
The alignSeqs process aligns read pairs (in two separate fastq format files) against the indexed reference sequencing using `bwa mem` and converts to bam format using `samtools view` to conserve disk space.  

## sortBAM  
The sortBAM process sorts the bam files from alignSeqs using `samtools sort`. This process was separated from the alignSeqs step to minimize redundant calculations should there be an error necessitating pipeline resumption.  

## markDuplicates  
The markDuplicates process marks PCR duplicates using either `sambamba markdup` or Picard MarkDuplicates (`java -jar picard.jar MarkDuplicates MAX_FILE_HANDLES_FOR_READ_ENDS_MAP=1000`).  

## fixReadGroups  
The fixReadGroups process adds sample read group information to the bam files using Picard AddOrReplaceReadGroups (`java -jar picard.jar AddOrReplaceReadGroups`).  

## realignIndels  
The realignIndels process first builds a bam index (.bai) for the fixReadGroups output bam using Picard BuildBamIndex (`java -jar picard.jar BuildBamIndex`). It then identifies intervals for indel realignment using GATK RealignerTargetCreator (`java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator`) and realigns indels using GATK IndelRealigner (`java -jar GenomeAnalysisTk.jar -T IndelRealigner --filter_bases_not_stored`).  

## filterBAMs  
The filterBAM process strictly filters BAM alignments using GATK PrintReads following [1] (`java -jar GenomeAnalaysisTK.jar -T PrintReads --read_filter BadCigar --read_filter DuplicateRead --read_filter FailsVendorQualityCheck --read_filter HCMappingQuality --read_filter MappingQualityUnavailable --read_filter NotPrimaryAlignment --read_filter UnmappedRead --filter_bases_not_stored --filter_mismatching_base_and_quals`).  

## fixMate  
The fixMate process corrects sequence mate pair tags using Picard FixMateInformation (`java -jar picard.jar FixMateInformation ADD_MATE_CIGAR=true`). It then builds a bam index (.bai) for the output bam file using Picard BuildBamIndex (`java -jar picard.jar BuildBamIndex`).  

## callVariants  

## genotypegVCFs  

## genMapIndex  

## genMapMap  

## repeatMask  

## repeatModeler  

## repeatMaskRM  

## maskIndels  

## simplifyBed  

## filterChr  

## splitVCFs  

## filterSites  

## filterRegions  

## calcDNMRate 

## summarizeDNM  

## References  
1. Besenbacher, S., Hvilsom, C., Marques-Bonet, T., Mailund, T., Schierup, M.H. (2019) Direct estimation of mutations in great apes reconciles phylogenetic dating. *Nat Ecol Evol*, __3__, 286-292. DOI:[10.1038/s41559-018-0778-x](https://www.nature.com/articles/s41559-018-0778-x).  

