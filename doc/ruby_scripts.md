# RatesTools Ruby Scripts  

__Michael G. Campana & Ellie E. Armstrong, 2019-2020__  
Smithsonian Conservation Biology Institute  
Stanford University  

Here we document the usage and functions of the Ruby scripts included in the RatesTools package.  

## calc_denovo_mutation_rate.rb  
The calc_denovo_mutation_rate.rb script calculates the genomic de novo mutation (DNM) rate from a multi-individual all-sites VCF. The script can optionally perform block bootstrapping to estimate the confidence interval for the estimated DNM rates. All individuals that are not specified as either the 'sire' or 'dam' are assumed to be offspring of the specified individuals. Results are printed to STDOUT.  

Basic usage is: `calc_denovo_mutation_rate.rb [options]`. Help is available using `calc_denovo_mutation_rate.rb -h`.  

The following calc_denovo_mutation_rate.rb options are available:  
`-i, --input [FILE]`: Input VCF (Required). Files with the final extension '.gz' are assumed to be gzip-compressed.  
`-s, --sire [NAME]`: Sire's name in VCF (Required).  
`-d, --dam [NAME]`: Dam's name in VCF (Required).  
`--parhom`: Require parents to be homozygous at candidate DNM sites. Parental heterozygosity forces the candidate site(s) to be discarded.  
`--minAD1`: Discard candidate DNMs if parents have DNM alleles present (even if the parents' alleles are not called). Requires the 'AD' tag to be specified in the VCF.  
`--minAF [VALUE]`: Filter alleles by minimum frequency.  
`-w, --window [VALUE]`: Sequence window length (bp) for bootstrapping (Default = 1000000).  
`-S, --step [VALUE]`: Window step (bp) for bootstrapping (Default = 1000000).  
`-b, --bootstrap [VALUE]`: Number of bootstrap replicates (Default = 0).  
`-l, --minbootstraplength [VALUE]`: Minimum bootstrap window length (bp) to retain (Default = 1000000).  
`-M, --minwindows [VALUE]`: Minimum number of bootstrap windows to retain chromosome/contig (Default = 10).  
 `-g, --gvcf`: Input is a gVCF (Default = false).  
`--rng [VALUE]`: Random number seed.  

Program information options (Do not use other options with the following):  
`-v, --version`: Show program version.  
`-h, --help`: Show help.  

## denovolib.rb  
This script provides a library of methods and classes used by the remaining Ruby scripts in the RatesTools pipeline.  

## filterGM.rb  
The filterGM.rb script filters a GenMap mappability bed file. The user specifies a mappability cut-off above which to retain (default behavior). Optionally, the user can output regions below the cut-off (e.g. for subsequent removal) by appending 'exclude' to the command line. Input files with the final extension '.gz' are assumed to be gzip-compressed.  

Usage is: `filterGM.rb <in_GenMap.bed[.gz]> <cutoff> [exclude] > <out.bed>`.  

## indels2bed.rb  
The indels2bed.rb script identifies indels in the unfiltered all-sites, multi-sample VCF. It then generates a bed file that specifies a set window upstream/downstream around the indel for later exclusion in the pipeline. Input files with the final extension '.gz' are assumed to be gzip-compressed.  

Usage is: `indels2bed.rb <indels.vcf[.gz]> <bp_to_exclude_upstream/downstream> > <out.bed>`.  

## nextflow_split.rb  
The nextflow_split.rb script splits an all-sites VCF by chromosomes/contigs for parallelization of filtering and de novo mutation rate calculations using the RatesTools pipeline.

Basic usage is: `nextflow_split.rb -i <in.VCF> -o <outdir>`. Help is available using `nextflow_split.rb -h`.  

Options available:  
`-i, --input [FILE]`: Input VCF.  
`-o, --output [DIRECTORY]`: Output Directory (Default is current directory).  
`-W, --writecycles [VALUE]`: Number of variants to read before writing to disk (Default = 1000000).  
`-h, --help`: Show help (Do not use with other options).  

## parallel_denovo.rb  
*WARNING: parallel_denovo.rb is deprecated in favor of ratestools.nf and nextflow_split.rb. It is included for reference.*  

parallel_denovo.rb splits a previously filtered, all-sites VCF by chromosome, parallelizes the calculation of the per-chromosome mutation rate using calc_denovo_mutation_rate, and the calculates the genomic mutation rate from the per-chromosome rates. It is written for the Univa Grid Engine (UGE)-based Smithsonian Institution High Performance Computing (SI/HPC) cluster 'Hydra' and will require manual revision for deployment on other systems.  

Basic usage is: `parallel_denovo.rb [options]`. Help is available using `parallel_denovo.rb -h`.  

All calc_denovo_mutation_rate.rb options are available in parallel_denovo.rb. See [calc_denovo_mutation_rate.rb](#calc_denovo_mutation_raterb) for details. Additionally, the following parallel_denovo.rb-specific options are available:  

parallel_denovo.rb Options:  
`-o, --output [DIRECTORY]`: Output Directory (Default is current directory).  
`-W, --writecycles [VALUE]`: Number of variants to read before writing to disk (Default = 1000000).  
 `--nosubmit`: Generate split VCFs and job files, but do not submit them. This allows the splitting to be performed as a job (rather than on the head node), even in systems that do not permit subjobs.  
 `-r, --restart`: Restart from previously split VCFs (Default = false). This permits revision of the calc_denovo_mutation_rate parameters without having to re-split the starting VCF file.  
`--submit`: Submit previously generated jobs and split VCFs (Implies -r).  

SI/HPC Options:  
`-q, --queue [VALUE]`: Qsub queue to use (Default = sThC.q).  
`-m, --memory [VALUE]`: Reserved memory (Default = 1G).
`-H, --himem`: Use high-memory queue (Default is false).  
 `-L, --lopri`: Use low priority queue (Default is false).  
 `-e, --email [VALUE]`: E-mail address to notify. 

## RM2bed.rb  

## simplify_bed.rb  

## simplify_sorted_bed.rb  

## summarize_denovo.rb  
