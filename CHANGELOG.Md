# CHANGE LOG  
Michael G. Campana & Ellie E. Armstrong, 2019-2021  
Smithsonian Institution  
Stanford University  

## Table of Contents  
1. [calc_denovo_mutation_rate](#calc_denovo_mutation_rate)  
2. [denovolib](#denovolib)  
3. [filterGM](#filtergm)  
4. [indels2bed](#indels2bed)  
5. [kochDNp](#kochdnp)  
6. [nextflow_split](#nextflow_split)  
7. [nextflow.config](#nextflowconfig)  
8. [ratestools](#ratestools)  
9. [RM2bed](#rm2bed)  
10. [simplify_bed](#simplify_bed)  
11. [simplify_sorted_bed](#simplify_sorted_bed)  
12. [summarize_denovo](#summarize_denovo)  
13. [Deprecated](#deprecated)  


## calc_denovo_mutation_rate  
### Version 0.11.0  
Added handling for Koch DNp statistic filter  
Generalized ad_exit for all filtering errors (method filter_exit) and moved to denovolib  
Fixed double-forward mutation count bug  

### Version 0.10.3  
Improved gz_file_open using yield  

### Version 0.10.2  
Commented code  

### Version 0.10.1  
conf95 uses ".." to separate CI values to permit negative exponents  

### Version 0.10.0  
--minAF filter  

### Version 0.9.2  
--minAD1 and --parhom now shown in interpreted command output  
--minAD1 regenotyping now applies to offspring too  
--parhom flag fixed  

### Version 0.9.1  
Fixed minAD1 option  
Alleles no longer includes "." as valid allele  

### Version 0.9.0  
Added minAD1 and parhom options  

### Version 0.8.1  
Corrected site denominator for diploidy in bootstrap_results  

### Version 0.8.0  
print_results moved to denovolib  

### Version 0.7.0  
Redundant methods (Parser, gz_open_file) moved to denovolib  

### Version 0.6.0  
Prints progress updates to STDERR  

### Version 0.5.0  
Can specify minimum bootstrapping windows to perform analysis  

### Version 0.4.0  
Can read gzipped input VCFs (Imported from BaitsTools 1.6.5)  
Outputs VCF lines of identified de novo mutations  

### Version 0.3.0  
Random number seed setting  
Interpreted options printed at start-up  

### Version 0.2.0  
Bootstrapping  

### Version 0.1.0  
Initial script to calculate point estimate de novo mutation rate  

## denovolib  
### Version 0.8.0  
Added handling for Koch DNp statistic filter  
Generalized ad_exit from calc_denovo_mutation_rate as denovolib method filter_exit  
Removed parsing options for deprecated parallel_denovo  
Fixed double-forward mutation count bug  
split_vcf retains complete VCF header for bedtools intersect  

### Version 0.7.3  
Improved gz_file_open using yield  

### Version 0.7.2  
Commented code  
Fixed missing ')' in help for -q option  

### Version 0.7.1  
Fixed 'homozogyous' typo in --parhom option description  

### Version 0.7.0  
split_vcf moved to denovolib  

### Version 0.6.0  
--minAF filter  

### Version 0.5.1  
Site denominator bug fix  

### Version 0.5.0  
Added minAD1 and parhom options  

### Version 0.4.1  
Corrected site denominator for diploidy in print_results  

### Version 0.4.0  
Option '--nosubmit' splits VCFs and writes jobs, but does not submit them  
Option '--submit' submits previously generated jobs and split VCFs  

### Version 0.3.0  
Method 'format_splash' simplifies basic script method output  

### Version 0.2.0  
print_results moved to denovolib  

### Version 0.1.0  
Redundant methods (Parser, gz_open_file) moved to denovolib  
gz_open_file improved to version in BaitsTools 1.6.6  

## filterGM  
### Version 0.3.2  
Improved gz_file_open using yield  

### Version 0.3.1  
Code commented  

### Version 0.3.0  
Can read gzipped files  
Optional "exclude" parameter produces inverse BED file for exclusion  

### Version 0.2.0  
Uses denovolib 'format_splash' method  
Script capitalization consistency correction  

### Version 0.1.0  
Filters GenMap output by specified mappabilities  

## indels2bed  
### Version 0.3.2  
Improved gz_file_open using yield  

### Version 0.3.1  
Code commented  

### Version 0.3.0  
Can read gzipped input VCFs  
Identifies deletions coded by * character  
Identifies indels from allele lengths without need for external previous filter  

### Version 0.2.0  
Uses denovolib 'format_splash' method  

### Version 0.1.2  
Outputs headerless BED rather than incorrect BED header  

### Version 0.1.1  
Output BED header line starts with # as appropriate  

### Version 0.1.0  
Initial script to convert VCF of indels into BED file for variant site exclusion  

## kochDNp  
### Version 0.1.0  
Calculates and filters by Koch DNp statistic  

## nextflow_split  
### Version 0.1.3  
Fixed bug where $options.writecycles was not defined  
Prints help when no options/-h  

### Version 0.1.2  
Code commented  

### Version 0.1.1  
Requires denovolib rather than split_vcf  

### Version 0.1.0  
Accessor script for split_vcf for nextflow   

## nextflow.config  
### Version 0.3.0  
Fixed RepeatModeler additional thread glitch  
Added module options for bedtools, bcftools, bgzip, tabix  
Updated process defaults  

## plotDPGQ
### Version 0.2.0  
First version integrated into pipeline  

### Version 0.1.0
Initial script (plotDPGQ-clean.R)  

## ratestools  
### Version 0.3.0  
Region filtration now uses BEDTools intersect first, then BCFtools, then VCFtools  
Updated some steps to use bgzip when needed by BCFtools  
Fixed RepeatModeler additional thread glitch  
Added pullDPGQ process  

### Version 0.2.0  
Reordered vcf filtration so that split_vcf occurs first to maximize parallelization  
Fixed major bug in fixReadGroups that assigned read groups to the incorrect BAM file  
Alignments all use refseq simpleName rather than baseName for consistency of file naming  
Removed duplicate BuildBamIndex step in callVariants  

### Version 0.1.0  
Nextflow pipeline for de novo mutation discovery. Previously panthera_dnm.nf during construction   

## RM2bed  
### Version 0.4.2  
Improved gz_file_open using yield  

### Version 0.4.1  
Code commented  

### Version 0.4.0  
Can read gzipped files  

### Version 0.3.0  
Uses denovolib 'format_splash' method  

### Version 0.2.2  
Outputs headerless BED rather than incorrect BED header  

### Version 0.2.1  
Output BED header line starts with # as appropriate  

### Version 0.2.0  
Outputs usage if no parameters specified  

### Version 0.1.0  
Initial script to convert [RepeatMasker](http://www.repeatmasker.org/) output to BED format  

## simplify_bed  
### Version 0.2.1  
Improved gz_file_open using yield  

### Version 0.2.0  
Complete rewrite of sorting/merging algorithm to use intervals (~8x faster)  

### Version 0.1.2  
Status output to stderr  
Redundant sites not added to array to hopefully improve performance  

### Version 0.1.1  
Code commented  

### Version 0.1.0  
Initial script to combine overlapping BED entries for region inclusion/exclusion  

## simplify_sorted_bed  
### Version 0.1.1  
Improved gz_file_open using yield  

### Version 0.1.0  
Script to merge adjacent BED regions in a coordinate-sorted BED file to reduce simplify_bed memory use and compute time  

## summarize_denovo  
### Version 0.5.5  
Code commented  

### Version 0.5.4  
Site denominator returned to previous value as bootstrapped estimates already multiplied by two  
Resolved issue with mean_ci confidence interval splits using - as separator  

### Version 0.5.3  
Corrected site denominator for diploidy in mean_ci  

### Version 0.5.2  
Retained SNP collation bug fixed  

### Version 0.5.1  
Fixed collect_offspring bug  

### Version 0.5.0  
Uses denovolib 'format_splash' method  

### Version 0.4.0  
Redundant print_summary code moved to denovolib print_results  

### Version 0.3.0  
Unbootstrapped contigs are now excluded from calculations of bootstrap sequence lengths  

### Version 0.2.0  
Concatenates headerless VCF results of identified de novo mutations  

### Version 0.1.0  
Initial script to summarize parallelized de novo rate calculations  

## Deprecated 
## parallel_denovo  
Removed in RatesTools pipeline version 0.3  

### Version 0.10.2  
Code commented  

### Version 0.10.1  
split_vcf.rb no longer required  

### Version 0.10.0  
split_vcf moved to split_vcf.rb for separate access by nextflow  

### Version 0.9.0  
--minAF filter  

### Version 0.8.0  
--parhom and --minAD1 options added to write_qsub  

### Version 0.7.1  
'write_qsub' calc_denovo_mutation_rate.job file now gzips analyzed split VCFs  

### Version 0.7.0  
Option '--nosubmit' splits VCFs and writes jobs, but does not submit them  
Option '--submit' submits previously generated jobs and split VCFs  

### Version 0.6.0  
Redundant methods (Parser, gz_open_file) moved to denovolib  
summarize_denovo automatically called upon completion of all split jobs  on cluster  

### Version 0.5.0  
Template qsub written once and arguments passed to it  
Can restart processing from a previously split VCF  

### Version 0.4.0  
Can specify the number of variants to read while splitting before writing to disk  

### Version 0.3.0  
Can specify minimum bootstrapping windows to perform analysis  

### Version 0.2.0  
Can read gzipped input VCFs (Imported from BaitsTools 1.6.5)  

### Version 0.1.0  
Initial script to parallelize de novo mutation rate calcuations on Grid Engine HPC  

## split_vcf  
### Version 0.2.0  
Now option to gzip output for file saving  

### Version 0.1.0  
split_vcf moved to split_vcf.rb for separate access by nextflow  
