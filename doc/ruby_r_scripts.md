# RatesTools Ruby and R Scripts  

__Michael G. Campana & Ellie E. Armstrong, 2019-2022__  
Smithsonian's National Zoo and Conservation Biology Institute  
Stanford University  

Here we document the usage and functions of the Ruby scripts included in the RatesTools package.  

## calc_denovo_mutation_rate.rb  
The calc_denovo_mutation_rate.rb script calculates the genomic de novo mutation (DNM) rate from a multi-individual all-sites VCF. The script can optionally perform block bootstrapping to estimate the confidence interval for the estimated DNM rates. All individuals that are not specified as either the 'sire' or 'dam' are assumed to be offspring of the specified individuals. Results are printed to STDOUT.  

Basic usage is: `calc_denovo_mutation_rate.rb [options]`. Help is available using `calc_denovo_mutation_rate.rb -h`.  

The following calc_denovo_mutation_rate.rb options are available:  
`-i, --input [FILE]`: Input VCF (Required). Files with the final extension '.gz' are assumed to be gzip-compressed.  
`-s, --sire [NAME]`: Sire's name in VCF (Required).  
`-d, --dam [NAME]`: Dam's name in VCF (Required).  
`-k, --kochDNp`: Use Koch et al. DNp statistic [1] to filter candidate DNM sites.  
`-m, --mu [VALUE]`: Mutation rate (mu) for Koch DNp (Default = 1e-6)  
`-t, --theta [VALUE]`: Heterozygosity (theta) for Koch DNp (Default = 0.008)  
`-c, --cutoff [VALUE]`: Koch DNp cutoff (Default = 0.3)  
`--parhom`: Require parents to be homozygous at candidate DNM sites. Parental heterozygosity forces the candidate site(s) to be discarded.  
`--minAD1`: Call all alleles at a site with an allelic depth of at least 1 (even if not called by the genotyper). Requires the 'AD' tag to be specified in the VCF. If the total allelic depth is 0 (e.g. if no informative reads are included but the total read depth is greater than the required threshold), the site is discarded as a potential DNM. If both `--minAD1` and `--minAF` options specified, `--minAD1` is applied to the parents and `--minAF` is applied to the offspring for maximally conservative DNM calling.  
`--minAF [VALUE]`: Filter alleles by minimum frequency. Requires the 'AD' tag to be specified in the VCF. If both `--minAD1` and `--minAF` options specified, `--minAD1` is applied to the parents and `--minAF` is applied to the offspring for maximally conservative DNM calling.  
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

## dnm_summary_stats.rb  
The dnm_summary_stats.rb script summarizes retained site counts from the various filtration steps of the RatesTools pipeline. It also classifies single-forward DNM by their substitution class. All other candidate DNMs are counted as "Other".  

Usage is: `dnm_summary_stats.rb <logs_directory> <output_prefix> > <out.csv>`.  

## filterGM.rb  
The filterGM.rb script filters a GenMap [2] mappability bed file. The user specifies a mappability cut-off above which to retain sites (default behavior). Optionally, the user can output regions below the cut-off (e.g. for subsequent removal) by appending 'exclude' to the command line. Input files with the final extension '.gz' are assumed to be gzip-compressed.  

Usage is: `filterGM.rb <in_GenMap.bed[.gz]> <cutoff> [exclude] > <out.bed>`.  

## indels2bed.rb  
The indels2bed.rb script identifies indels in the unfiltered all-sites, multi-sample VCF. It then generates a bed file that specifies a set window upstream/downstream around the indel for later exclusion in the pipeline. Input files with the final extension '.gz' are assumed to be gzip-compressed.  

Usage is: `indels2bed.rb <indels.vcf[.gz]> <bp_to_exclude_upstream/downstream> > <out.bed>`.  

Basic usage is: `calc_denovo_mutation_rate.rb [options] > <outfile>`. Help is available using `calc_denovo_mutation_rate.rb -h`.  

The following calc_denovo_mutation_rate.rb options are available:  
`-i, --input [FILE]`: Input VCF (Required). Files with the final extension '.gz' are assumed to be gzip-compressed.  
`-s, --sire [NAME]`: Sire's name in VCF (Required).  
`-d, --dam [NAME]`: Dam's name in VCF (Required).  

## kochDNp.rb  
The kochDNp.rb script calculates the Koch et al. DNp statistic [1] for candidate de novo mutations and filters them by a specified cutoff value. The methods in this script are accessed by calc_denovo_mutation_rate if the DNp filter is turned on. Alternatively, the kochDNp script can be executed directly for post-hoc analysis of previously identified candidate de novo mutations. Results are printed to STDOUT.  

Basic usage is: `kochDNp.rb [options]`. Help is available using `kochDNp.rb -h`.  

The following kochDNp.rb options are available:  
`-i, --input [FILE]`: Input file of candidate denovo mutations (Required). Input file can either be a VCF or a log previously generated using calc_denovo_mutation_rate or summarize_denovo. Files with the final extension '.gz' are assumed to be gzip-compressed.  
`-s, --sire [NAME]`: Sire's name in VCF (Required).  
`-d, --dam [NAME]`: Dam's name in VCF (Required).  
`-m, --mu [VALUE]`: Mutation rate (mu) for Koch DNp (Default = 1e-6)  
`-t, --theta [VALUE]`: Heterozygosity (theta) for Koch DNp (Default = 0.008)  
`-c, --cutoff [VALUE]`: Koch DNp cutoff (Default = 0.3)  

Program information options (Do not use other options with the following):  
`-v, --version`: Show program version.  
`-h, --help`: Show help.  

# logstats.sh  
The logstats.sh script extracts the VCFtools retained site counts from Nextflow's log file (.command.log) and checks the value for sanity. If VCFtools returns a non-sensical value (typically due to an overflow of sites causing VCFtools to report a negative site count), logstats.sh calculates these values using BCFtools [3] stats.  

## nextflow_split.rb  
The nextflow_split.rb script splits an all-sites VCF by chromosomes/contigs for parallelization of filtering and de novo mutation rate calculations using the RatesTools pipeline.  

Basic usage is: `nextflow_split.rb -i <in.VCF> -o <outdir>`. Help is available using `nextflow_split.rb -h`.  

Options available:  
`-i, --input [FILE]`: Input VCF (Required). Files with the final extension '.gz' are assumed to be gzip-compressed.  
`-o, --output [DIRECTORY]`: Output Directory (Default is current directory).  
`-W, --writecycles [VALUE]`: Number of variants to read before writing to disk (Default = 1000000).  
`-h, --help`: Show help (Do not use with other options).  

## plotDPGQ.R  
plotDPGQ.R summarizes extracted variant site DP and GQ values extracted using BCFtools [3] query. It also plots the distributions of these statistics.  

Usage is: `plotDPGQ.R <outfile_prefix>`.  

## RM2bed.rb  
RM2bed.rb converts a RepeatMasker [4] out file into a bed for later exclusion of repeat regions. Input files with the final extension '.gz' are assumed to be gzip-compressed.  

Usage: `RM2bed.rb <in_RM.txt[.gz]> > <out.bed>`.  

## simplify_bed.rb  
simplify_bed.rb merges overlapping entries in a bed file. It does not assume the input bed file is chromosome/contig- and coordinate-sorted, and therefore takes more computing resources and time to perform the merging than the [simplify_sorted_bed.rb](#simplify_sorted_bed.rb) script. Input files with the final extension '.gz' are assumed to be gzip-compressed.  

Usage: `simplify_bed.rb <in.bed[.gz]> > <out.bed>`.  

## simplify_sorted_bed.rb  
simplify_sorted_bed.rb merges overlapping entries in a chromosome/contig- and coordinate-sorted bed file. For an unsorted bed file, use the [simplify_sorted_bed.rb](#simplify_sorted_bed.rb) script. Input files with the final extension '.gz' are assumed to be gzip-compressed.  

Usage: `simplify_sorted_bed.rb <in.bed[.gz]> > <out.bed>`.  

## summarize_denovo.rb  
summarize_denovo.rb calculates the genomic DNM rate from a directory of per-chromosome calc_denovo_mutation_rate.rb logs.  

Usage: `summarize_denovo.rb <directory> > <out.txt>`.  

## References  
1. Koch, E.M., Schweizer, R.M., Schweizer, T.M., Stahler, D.R., Smith, D.W., Wayne, R.K., Novembre, J. (2019). De novo mutation rate estimation in wolves of known pedigree. *Mol Biol Evol*, __36__, 2536-2547, doi: [10.1093/molbev/msz159](https://academic.oup.com/mbe/article/36/11/2536/5531468?login=true).  
2. Pockrandt, C., Alzamel, M., Iliopoulos, C.S., Reinert, K. (2020) GenMap: ultra-fast computation of genome mappability. *Bioinformatics*, __36__, 3687–3692, doi: [10.1093/bioinformatics/btaa222](https://academic.oup.com/bioinformatics/article/36/12/3687/5815974?login=true).  
3. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H. (2021) Twelve years of SAMtools and BCFtools. *GigaScience*, __10__, giab008. DOI: [10.1093/gigascience/giab008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722).  
4. Smit, A.F.A., Hubley, R., Green, P. (2013-2015) *RepeatMasker Open-4.0*. (http://www.repeatmasker.org).  
