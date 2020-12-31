# RatesTools Ruby Scripts  

__Michael G. Campana & Ellie E. Armstrong, 2019-2020__  
Smithsonian Conservation Biology Institute  
Stanford University  

Here we document the usage and functions of the Ruby scripts included in the RatesTools package.  

## calc_denovo_mutation_rate.rb  
calc_denovo_mutation_rate.rb

Basic usage is: `calc_denovo_mutation_rate.rb [options]`. Help is available using `calc_denovo_mutation_rate.rb -h`.  

calc_denovo_mutation_rate options:
    -i, --input [FILE]               Input VCF
    -s, --sire [NAME]                Sire's name in VCF
    -d, --dam [NAME]                 Dam's name in VCF
        --parhom                     Require parents to be homozygous at DNM sites
        --minAD1                     Discard DNMs if parents have DNM alleles even if not called
        --minAF [VALUE]              Filter alleles by minimum frequency
    -w, --window [VALUE]             Sequence window length (bp) for bootstrapping (Default = 1000000)
    -S, --step [VALUE]               Window step (bp) for bootstrapping (Default = 1000000)
    -b, --bootstrap [VALUE]          Number of bootstrap replicates (Default = 0)
    -l, --minbootstraplength [VALUE] Minimum bootstrap window length (bp) to retain (Default = 1000000)
    -M, --minwindows [VALUE]         Minimum number of bootstrap windows to retain contig (Default = 10)
    -g, --gvcf                       Input is a gVCF (Default = false)
        --rng [VALUE]                Random number seed

## denovolib.rb  
This script provides a library of methods and classes used by the remaining Ruby scripts in the RatesTools pipeline.  

## filterGM.rb  

## indels2bed.rb  

## nextflow_split.rb  

## parallel_denovo.rb  
*WARNING: parallel_denovo.rb is deprecated in favor of ratestools.nf and nextflow_split.rb. It is included for reference.*  

parallel_denovo.rb splits a previously filtered, all-sites VCF by chromosome, parallelizes the calculation of the per-chromosome mutation rate using calc_denovo_mutation_rate, and the calculates the genomic mutation rate from the per-chromosome rates. It is written for the Univa Grid Engine (UGE)-based Smithsonian Institution High Performance Computing (SI/HPC) cluster 'Hydra' and will require manual revision for deployment on other systems.  

Basic usage is: `parallel_denovo.rb [options]`. Help is available using `parallel_denovo.rb -h`.  

All calc_denovo_mutation_rate.rb options are available in parallel_denovo.rb. See [calc_denovo_mutation_rate.rb](#calc_denovo_mutation_raterb) for details. Additionally, the following parallel_denovo.rb-specific options are available:  

parallel_denovo.rb Options:  
`-o, --output [DIRECTORY]`: Output Directory (Default is current directory)  
`-W, --writecycles [VALUE]`: Number of variants to read before writing to disk (Default = 1000000)  
 `--nosubmit`: Generate split VCFs and job files, but do not submit them. This allows the splitting to be performed as a job (rather than on the head node), even in systems that do not permit subjobs.  
 `-r, --restart`: Restart from previously split VCFs (Default = false). This permits revision of the calc_denovo_mutation_rate parameters without having to re-split the starting VCF file.  
`--submit`: Submit previously generated jobs and split VCFs (Implies -r).  

SI/HPC Options:  
`-q, --queue [VALUE]`: Qsub queue to use (Default = sThC.q)  
`-m, --memory [VALUE]`: Reserved memory (Default = 1G)
`-H, --himem`: Use high-memory queue (Default is false)  
 `-L, --lopri`: Use low priority queue (Default is false)  
 `-e, --email [VALUE]`: E-mail address to notify. 

## RM2bed.rb  

## simplify_bed.rb  

## simplify_sorted_bed.rb  

## summarize_denovo.rb  
