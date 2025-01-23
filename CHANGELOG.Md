# CHANGE LOG  
Michael G. Campana & Ellie E. Armstrong, 2019-2024  
Smithsonian Institution  
Stanford University  

## Table of Contents  
1. [calc_denovo_mutation_rate](#calc_denovo_mutation_rate)  
2. [configure.sh](#configuresh)  
3. [denovolib](#denovolib)  
4. [dnm_summary_stats](#dnm_summary_stats)
5. [filterGM](#filtergm)  
6. [indels2bed](#indels2bed)  
7. [kochDNp](#kochdnp)  
8. [logstats](#logstats)  
9. [nextflow_split](#nextflow_split)  
10. [nextflow.config](#nextflowconfig)  
11. [plotDPGQ](#plotDPGQ)  
12. [ratestools](#ratestools)  
13. [RM2bed](#rm2bed)  
14. [simplify_bed](#simplify_bed)  
15. [simplify_sorted_bed](#simplify_sorted_bed)  
16. [summarize_denovo](#summarize_denovo)  
17. [Deprecated](#deprecated)  


## calc_denovo_mutation_rate  
### Version 1.2.0  
Fixed minor spacing bug in output of command line options  
Add DNM_GT tag for genotypes used after minAF and minAD1 processing  
Fixed bug that doubled double-forward mutations in the rate estimates  

### Version 1.1.2  
Fixed bugs in calculating the bootstrapped 95% confidence interval  

### Version 0.12.3  
Fixed missing bootstrap parameter in command-line output in log  

### Version 0.12.2  
Update file headers  

### Version 0.12.1  
Handling for selfing parent  

### Version 0.12.0  
Applying both --minAD1 and --minAF no longer defaults to --minAD1. --minAD1 is applied to parents, --minAF to offspring.  

### Version 0.11.2  
Fixed GATK4 total AD = 0 bug  
More accurate minAD1 help description  
Fixed snp_record (line 239) glitch  

### Version 0.11.1  
Fixed missing kochDNp require_relative statement  
Fixed cutoff instead of $options.cutoff bug  
Fixed bug that called whole offspring_pl hash rather than individual PL array  

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

## configure.sh  
### Version 1.2.3  
Updated ReadPosRankSum default filter  

### Version 1.2.0  
Removal of redundant dnm_bootstraps parameter  

### Version 1.1.1  
Updates for revised dnm_summary_stats script  

### Version 1.1.0  
trio-phasing configuration  

### Version 1.0.0  
Updated for RatesTools v. 1.0.0  

### Version 0.5.10  
Update file headers  
Contig length filtration  

### Version 0.5.8  
Conda installation handling  

### Version 0.5.7  
$launchDir instead of $baseDir  
Removed thread parameters  

## denovolib  
### Version 1.2.0  
Fixed bug that doubled double-forward mutations in the rate estimates  
Added binomial confidence intervals using Hmisc  

### Version 0.8.2  
Update file headers  

### Version 0.8.1  
split_vcf now uses bgzip rather than gzip for GATK site filtering  
filter_exit returns exit status 1  

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

## dnm_summary_stats  
### Version 1.2.0  
Removal of redundant dnm_bootstraps parameter  
Uses DNM_GT tag to correctly count number of sites after filtering  
Correctly counts double-forward sites as two mutations  
Calculates binomial confidence intervals using R package Hmisc  
Bug fixes in number of all-sites removed  
Outputs a TSV list of sites retained to standard error for VCF processing  
Outputs filtered mutational spectra  

### Version 1.1.2  
Fixes bug that incorrectly calculated adjusted standard errors  

### Version 1.1.1  
Automatically removes clumped DNM candidates  
Recalculates adjusted confidence intervals  

### Version 1.0.0  
Updated for ratestools 1.0.0  
Added handling for more mutation classes  
Added code to identify shared mutations between offspring and recalculate mutation rates  

### Version 0.1.1  
Update file headers  

### Version 0.1.0  
First script to calculate summary statistics  

## filterGM  
### Version 0.3.3  
Update file headers  

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
### Version 0.3.3  
Update file headers  

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
### Version 0.1.2  
Update file headers  

### Version 0.1.1  
Fixed bug that ran command-line parser when run from calc_denovo_mutation_rate

### Version 0.1.0  
Calculates and filters by Koch DNp statistic  

## nextflow_split  
### Version 0.1.4  
Update file headers  

### Version 0.1.3  
Fixed bug where $options.writecycles was not defined  
Prints help when no options/-h  

### Version 0.1.2  
Code commented  

### Version 0.1.1  
Requires denovolib rather than split_vcf  

### Version 0.1.0  
Accessor script for split_vcf for nextflow   

## logstats  
### Version 1.0.0  
Added error handling when the run time line is not recorded  

### Version 0.1.3  
Fixed bug that discarded the unfiltered site total for contigs where all sites were removed  

### Version 0.1.2  
Update file headers  
Contig length filtration  

### Version 0.1.1  
Improved sanity checks  

### Version 0.1.0  
Script to correct vcftools site overflow glitch using bcftools stats  

## nextflow.config  
### Version 1.2.4  
Added BCFtools to conda configuration  

### Version 1.2.3  
Updated ReadPosRankSum default filter  

### Version 1.2.2  
Updated R channel to conda-forge  

### Version 1.2.0  
Removal of redundant dnm_bootstraps parameter  
Update conda dependencies for generateSummaryStatsm, calcDNMRate, summarizeDNM  

### Version 1.1.1  
DNM clump option  

### Version 1.1.0  
phaseTrio process configuration  

### Version 1.0.0  
Updated for RatesTools v. 1.0.0  

### Version 0.5.14  
Set repeatModeler task.cpus = 2 for standard and conda profiles  

### Version 0.5.12  
Added conda.enabled = true  

### Version 0.5.10  
Update file headers  
Contig length filtration  

### Version 0.5.8  
Conda profile and handling  

### Version 0.5.7  
Removed 'hydra' example profile  
Removed storeDir and maxForks from local profile  
Removal of thread parameters  
$launchDir instead of $baseDir  

### Version 0.5.5  
Updated summarizeDNM modules for VCF output  

### Version 0.5.2  
process sanityCheckLogs for logstats.sh  

### Version 0.5.1  
Modules updated for logstats.sh  
Fixed local profile RepeatModeler -> repeatModeler  
Closures rather that objects for Nextflow 22.04.04 compatibility  

### Version 0.5.0  
GATK4 compatibility  
Email errors/completion status  
Variable indel padding  

### Version 0.4.0  
Split site_filters into vcftools_site_filters and gatk_site_filters  
Updated process list in both profiles for new processes  
Fixed some missing processes from v0.3.0 in default profile  

### Version 0.3.0  
Fixed RepeatModeler additional thread glitch  
Added module options for bedtools, bcftools, bgzip, tabix  
Updated process defaults  

## plotDPGQ  
### Version 0.3.1  
Update file headers  

### Version 0.3.0  
Fixed bug that misread sample names  
Corrected output so that all PNGs properly named  
Fixed table coding so that all individual samples calculated in CSVs  

### Version 0.2.0  
First version integrated into pipeline  

### Version 0.1.0
Initial script (plotDPGQ-clean.R)  

## ratestools  
### Version 1.2.1  
Fixed region filter bypass bgzip bug  

### Version 1.2.0  
Removal of redundant dnm_bootstraps parameter  
Fixed header bug on candidates VCF  
Outputs clump- and sibling-filtered candidates VCF  

### Version 1.1.1  
Automatic DNM clump removal  
Automatic recalculation of confidence intervals  
Output summary stats includes file prefix  

### Version 1.1.0  
Haplotype phasing using WhatsHap  
Fixed bug in repeatMaskRM that did not retain new masked files  

### Version 1.0.0  
Upgraded to nextflow DSL2  
Readgroups added and alignment sorted during alignment process  
Input reads specified via a CSV  
Automated sample library merging  
Added samtools markdup as a duplicate option  
Bumped pipeline software versions  
Added GenMap mapping option parameter  
Added RepeatMasker and RepeatModeler options parameters  
BED simplification now uses BEDtools  
simplify_bed and simplify_sorted_bed removed  
Added option to ignore region filtration  
errorStrategy calls moved to configuration file for easier manipulation  
Added option to bypass BAM filtering  
Conda installs now use Mamba by default  
RepeatModeler process updated for RepeatMasker 2.0.5  

### Version 0.5.15  
Fixed bug in GATK4 pipeline that did not apply GATK site filters appropriately  

### Version 0.5.14  
Added handling for RepeatModeler non-even positive thread requests  

### Version 0.5.13  
Added tbl file of first round of RepeatMasker to retained output  
Fixed RM2bed concatenation glitch in repeatMaskRM

### Version 0.5.12  
Reduced hard-disk footprint of gatkFilterSites temporary files  
RepeatMasker soft-masks rather than hard-masks  

### Version 0.5.11  
Fixed glitch in summarizeDNM when no candidate sites found  

### Version 0.5.10  
Update file headers  
Contig length filtration  

### Version 0.5.9  
Handling for selfing parent in calc_denovo_mutation_rate.rb  

### Version 0.5.8  
Conda handling  
Use touch rather than mkfile for RepeatModeler dummy file  

### Version 0.5.7  
DSL1 specified  
Used task.cpus rather than thread parameters  

### Version 0.5.6  
Fixed bug in region filtration when empty region bed  
bedtools subtract rather than interact for second attempt  

### Version 0.5.5  
summarizeDNM outputs VCF of candidate DNMs  

### Version 0.5.4  
Fixed logging glitches that copied VCF output into downstream results  
Fixed glitch in BCFtools region-filtering that removed sites unexpectedly  
bedtools subtract rather than intersect -v  
Fixed glitch where bedtools error was not detected and region-filtered VCFs were truncated  

### Version 0.5.2  
process sanityCheckLogs for logstats.sh  

### Version 0.5.1  
Site count logging uses logstats.sh to address VCFtools overflow  

### Version 0.5.0  
GATK4 compatibility  
Improved use of symbolic links/writing to disk  
Completion messages and emailing of completion status/errors  
Variable indel padding  
Use path rather than file qualifiers  
Site count logging and summary statistics  

### Version 0.4.0  
Updated plotDPGQ output PNGs  
split filterSites in vcftoolsFilterSites and gatkFilterSites  
Added 'NULL' option to site filtrations to bypass unnecessary processing  
Divided filterChr into filterChr and splitTrios  
filterChr (new version) now only does chromosome filtering  
splitTrios does trio splitting after chromosome filtration  
Fixed pullDPGQ input so that all individuals calculated from chr-filtered VCFs  
Numbered output directories  

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
### Version 0.4.3  
Update file headers  

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

## summarize_denovo  
### Version 0.5.6  
Update file headers  

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

## simplify_bed  
Removed in RatesTools version 1.0.0

### Version 0.2.2  
Update file headers  

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
Removed in RatesTools version 1.0.0  

### Version 0.1.2  
Update file headers  

### Version 0.1.1  
Improved gz_file_open using yield  

### Version 0.1.0  
Script to merge adjacent BED regions in a coordinate-sorted BED file to reduce simplify_bed memory use and compute time  

## makefile  
Removed in RatesTools version 0.5.7  

### Version 0.5.3  
Updated makefile for Linux and macOS  

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
