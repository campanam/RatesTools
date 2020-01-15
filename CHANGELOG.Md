# CHANGE LOG  
Michael G. Campana, 2019-2020  
Smithsonian Conservation Biology Institute  

## calc_denovo_mutation_rate  
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

## indels2bed  
### Version 0.1.0  
Initial script to convert VCF of indels into BED file for variant site exclusion  

## parallel_denovo  
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

## RM2bed  
### Version 0.2.0  
Outputs usage if no parameters specified  

### Version 0.1.0  
Initial script to convert [RepeatMasker](http://www.repeatmasker.org/) output to BED format  

## summarize_denovo  
### Version 0.3.0  
Unbootstrapped contigs are now excluded from calculations of bootstrap sequence lengths  

### Version 0.2.0  
Concatenates headerless VCF results of identified de novo mutations  

### Version 0.1.0  
Initial script to summarize parallelized de novo rate calculations  

## References
Campana, MG. 2018. BaitsTools: software for hybridization capture bait design. *Molecular Ecology Resources* __18__, 356-361. ( https://doi.org/10.1111/1755-0998.12721)  
Smit, AFA, Hubley, R & Green, P. RepeatMasker Open-4.0. 2013-2015 (http://www.repeatmasker.org).  
