# RatesTools Ruby Scripts  

__Michael G. Campana & Ellie E. Armstrong, 2019-2020__  
Smithsonian Conservation Biology Institute  
Stanford University  

Here we document the usage and functions of the Ruby scripts included in the RatesTools package.  

## calc_denovo_mutation_rate.rb  

## denovolib.rb  
This script provides a library of methods and classes used by the remaining Ruby scripts in the RatesTools pipeline.  

## filterGM.rb  

## indels2bed.rb  

## nextflow_split.rb  

## parallel_denovo.rb  
*WARNING: parallel_denovo.rb is deprecated in favor of ratestools.nf and nextflow_split.rb. It is included for reference.*  

parallel_denovo splits a previously filtered, all-sites VCF by chromosome, parallelizes the calculation of the per-chromosome mutation rate using calc_denovo_mutation_rate, and the calculates the genomic mutation rate from the per-chromosome rates. It is written for the Univa Grid Engine (UGE)-based Smithsonian Institution High Performance Computing (SI/HPC) cluster 'Hydra' and will require manual revision for deployment on other systems.  

All calc_denovo_mutation_rate options are available in parallel_denovo. See [calc_denovo_mutation_rate.rb](#calc_denovo_mutation_raterb) for details. Additionally, the following parallel_denovo options are available:  

parallel_denovo Options:
    -o, --output [DIRECTORY]         Output Directory (Default is current directory)
    -W, --writecycles [VALUE]        Number of variants to read before writing to disk (Default = 1000000)
        --nosubmit                   Generate split VCFs and jobs, but do not submit them
    -r, --restart                    Restart from previously split VCFs (Default = false)
        --submit                     Submit previously generated jobs and split VCFs (Implies -r)

SI/HPC Options:
    -q, --queue [VALUE]              Qsub queue to use (Default = sThC.q
    -m, --memory [VALUE]             Reserved memory (Default = 1G)
    -H, --himem                      Use high-memory queue (Default is false)
    -L, --lopri                      Use low priority queue (Default is false)
    -e, --email [VALUE]              E-mail address to notify



## RM2bed.rb  

## simplify_bed.rb  

## simplify_sorted_bed.rb  

## summarize_denovo.rb  
