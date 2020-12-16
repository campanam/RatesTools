# CHANGE LOG  
Michael G. Campana, 2019-2020  
Smithsonian Conservation Biology Institute  

## calc_denovo_mutation_rate  
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
### Version 0.7.2  
Commented code  

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
### Version 0.3.0  
Can read gzipped files  
Optional "exclude" parameter produces inverse BED file for exclusion  

### Version 0.2.0  
Uses denovolib 'format_splash' method  
Script capitalization consistency correction  

### Version 0.1.0  
Filters GenMap output by specified mappabilities  

## indels2bed  
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

## nextflow_split  
### Version 0.1.1  
Requires denovolib rather than split_vcf  

### Version 0.1.0  
Accessor script for split_vcf for nextflow   

## parallel_denovo  
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

## ratestools.nf  
### Version 0.1.0  
Nextflow pipeline for de novo mutation discovery. Previously panthera_dnm.nf during construction   

## RM2bed  
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
### Version 0.1.0  
Initial script to combine overlapping BED entries for region inclusion/exclusion  

## summarize_denovo  
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

# DEFUNCT  
## split_vcf  
### Version 0.2.0  
Now option to gzip output for file saving  

### Version 0.1.0  
split_vcf moved to split_vcf.rb for separate access by nextflow  

## References
Campana, MG. 2018. BaitsTools: software for hybridization capture bait design. *Molecular Ecology Resources* __18__, 356-361. ( https://doi.org/10.1111/1755-0998.12721)  
Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0. 2013-2015 (http://www.repeatmasker.org).  
