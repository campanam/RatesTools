# RatesTools  

__Michael G. Campana and Ellie E. Armstrong, 2019-2021__  
Smithsonian Institution  
Stanford University  

Pipeline to calculate de novo mutation rates from parent-offspring trios  

This README provides basic details for installing, configuring and running the pipeline. Detailed documentation is available for the [Ruby scripts](doc/ruby_scripts.md) included in this package and for the [pipeline's operation](doc/pipeline_details.md).  

## Creative Commons 0 Waiver  
![image](https://user-images.githubusercontent.com/19614608/118704084-bf02f280-b7e4-11eb-8d59-0ce648313d9e.png)  
To the extent possible under law, the Smithsonian Institution and Stanford University have waived all copyright and related or neighboring rights to RatesTools; this work is published from the United States.  

## Citation  
We politely request that this work be cited as:  
Campana, M.G. & E.E. Armstrong. 2020. RatesTools: Pipeline to calculate de novo mutation rates from parent-offspring trios. Smithsonian Institution and Stanford University. (https://github.com/campanam/RatesTools).  

## Installation and Configuration  
### Install Nextflow and Ruby  
RatesTools requires [Nextflow](https://www.nextflow.io/) [1] v. 20.10.0, [Ruby](http://www.ruby-lang.org) v. 2.6.3 and [R](https://www.r-project.org/) [2] v. 4.0.2. Basic instructions for installing these languages are copied below. We recommend installing Ruby using the [Ruby Version Manager](https://rvm.io). See the official language documentation should you need help installing these languages.  

Install Nextflow: `curl -s https://get.nextflow.io | bash`  
Install R: Use the appropriate precompiled binary/installer available at the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org).  
Install the latest Ruby using Ruby Version Manager: `curl -sSL https://get.rvm.io | bash -s stable --ruby`  

### Install the RatesTools Scripts  
Clone the repository: `git clone https://github.com/campanam/RatesTools`  
Install the scripts: `cd RatesTools; make install` 

*By default, RatesTools scripts will be installed into the ~/ratestools directory. If you wish to change the default directory, specify the INSTALL parameter, e.g.:* `make INSTALL=/path/to/some/dir install`  

### Install the External Dependencies  
RatesTools requires the following external dependencies. See the documentation for these programs for their installation requirements. RatesTools requires the Genome Analysis Toolkit (GATK) [3] version 3.8-1 and Java 1.8 (due to GATK). Currently, RatesTools is not compatible with GATK 4 or other versions of Java. Otherwise, listed versions are those that have been tested and confirmed, but other versions may work. RatesTools can utilize modules to simplify deployment on computing clusters and limit dependency conflicts.  

* gzip  
* awk  
* sed  
* zcat  
* [BWA](http://bio-bwa.sourceforge.net/) [3] v. 0.7.17  
* [SAMtools](http://www.htslib.org/) [4,5] v. 1.9  
* [BCFtools](http://www.htslib.org/) [4,5] v. 1.9  
* [bgzip and tabix from HTSlib](http://www.htslib.org/) [5] v. 1.9  
* [Java](https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) v. 1.8  
* [Picard](https://broadinstitute.github.io/picard/) [6] v. 2.23.8  
* [Sambamba](https://lomereiter.github.io/sambamba/) [7] v. 0.7.1  
* [Genome Analysis Toolkit](https://github.com/broadgsa/gatk) v. 3.8-1  
* [VCFtools](https://vcftools.github.io/index.html) [8] v.0.1.16  
* [GenMap](https://github.com/cpockrandt/genmap) [9] v.1.2.0 with [SeqAn](https://github.com/seqan/seqan/tree/f548b50705be3f824a65a696943ea90a390564ce) [10] v. 2.4.1  
* [RepeatMasker](http://www.repeatmasker.org/) [11] v. 4.0.9  
* [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) [12] v. 2.0.1  
* [BEDTools](https://bedtools.readthedocs.io/en/latest/) [13] v. 2.28.0  

RatesTools requires the following R packages installed in your R environment:  
* [tidyverse] [14] v. 1.3.1 with [dplyr] v. 1.0.7 [15]
* [data.table] [16]
* [ggpubr] [17]


### Configure the Pipeline    
Assisted configuration of the RatesTools pipeline can be accomplished using the `configure.sh` bash script. The script copies the `nextflow.config` included with the repository and modifies the copy for the target system. The `configure.sh` script detects software installed on the local system and prompts the user to provide module files, paths to undetected files, and program options. The configuration file can also be manually edited using a text editor. However, please note that the `configure.sh` script requires an *unmodified* `nextflow.config` file to work.  

Please note that RatesTools automatically detects read pairs using globbing and the Nextflow Channel.fromFilePairs() method (https://www.nextflow.io/docs/latest/channel.html#fromfilepairs). The user will need to specify a globbing pattern corresponding to the data. RatesTools also assumes that the sample name (e.g. for the sire and dam) is the shared portion of the read pair file name, excluding text after the first difference or specified in the globbing pattern. It may be ideal to rename your reads to minimize the extraneous information in the read name (e.g. lane information). For instance, using a typical Illumina naming scheme:

Using the globbing pattern `*{R1,R2}_001.fastq.gz`:  
Read pairs `LION1_S01_L001_R1_001.fastq.gz` and `LION1_S01_L001_R2_001.fastq.gz` will be matched and have the sample name `LION1_S01_L001_`.  
Renaming the files to `LION1R1_001.fastq.gz` and `LION1R2_001.fastq.gz` will match these reads with the cleaner name `LION1`.  

### Platform-Specific Configuration  
Given the wide-variety of computing architectures and operating systems, we cannot provide specific optimized configurations for your computing system. The `nextflow.config` file includes two example [configuration profiles](https://www.nextflow.io/docs/latest/config.html#config-profiles): a default 'standard' profile for a local installation and a 'hydra' profile optimized for the SI/HPC Univa Grid Engine (UGE) computing cluster 'Hydra-5'. Please consult your computing staff to optimize the profile settings for your hardware. Additionally, by default, RatesTools stores process output in a directory named 'chkpnt' and limits the number of concurrent process forks to 1000. These settings can be changed by modifying the lines `process.storeDir = 'chkpnt'` and `process.maxForks = 1000` in the config file.  

## Running the Pipeline  
Enter `ratestools.nf -c <config_file>` to run the pipeline. Append `-resume` to restart a previous run or `-bg` to run RatesTools in the background. If you developed platform-specific configuration profiles, you can specify this using the `-profile <PROFILE>` option. See the Nextflow documentation for details. Final data are written to the specified output directory and its subdirectories.  

## References  
1. Di Tommaso, P., Chatzou, M., Floden, E.W., Prieto Barja, P., Palumbo, E., Notredame, C. (2017) Nextflow enables reproducible computational workflows. *Nat Biotechnol*, __35__, 316–319. DOI: [10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820).  
2. R Core Team (2020) *R: A language and environment for statistical computing.* R Foundation for Statistical Computing, Vienna, Austria. [URL](https://www.r-project.org/).  
3. McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., DePristo, M.A. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res*, __20__, 1297-1303. DOI: [10.1101/gr.107524.110](https://genome.cshlp.org/content/20/9/1297.abstract).  
4. Li, H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv*, [1303.3997v2](https://arxiv.org/abs/1303.3997).  
5. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., 1000 Genome Project Data Processing Subgroup (2009) The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25, 2078-2079. DOI: [10.1093/bioinformatics/btp352](https://academic.oup.com/bioinformatics/article/25/16/2078/204688).  
6. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H. (2021) Twelves years of SAMtools and BCFtools. *GigaScience*, __10__, giab008. DOI: [10.1093/gigascience/giab008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722).  
7. Broad Institute (2020). Picard v. 2.23.8 (https://broadinstitute.github.io/picard/).  
8. Tarasov, A., Vilella, A.J., Cuppen, E., Nijman, I.J., Prins, P. (2015) Sambamba: fast processing of NGS alignment formats. *Bioinformatics*, __31__, 2032–2034. DOI: [10.1093/bioinformatics/btv098](https://academic.oup.com/bioinformatics/article/31/12/2032/214758).  
9. Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., McVean, G., Durbin, R. (2011) The variant call format and VCFtools. *Bioinformatics*, __27__, 2156–2158. DOI: [10.1093/bioinformatics/btr330](https://academic.oup.com/bioinformatics/article/27/15/2156/402296).  
10. Pockrandt, C., Alzamel, M., Iliopoulos, C.S., Reinert, K. (2020) GenMap: ultra-fast computation of genome mappability. *Bioinformatics*, __36__, 3687–3692, doi: [10.1093/bioinformatics/btaa222](https://academic.oup.com/bioinformatics/article/36/12/3687/5815974).  
11. Reinert, K., Dadi, T.H., Ehrhardt, M., Hauswedell, H., Mehringer, S., Rahn, R., Kim. J., Pockrandt, C., Winkler, J., Siragusa, E., Urgese, G., Weese, D. (2017) The SeqAn C++ template library for efficient sequence analysis: A resource for programmers. *J Biotechnol*, __261__, 157-168. DOI: [10.1016/j.jbiotec.2017.07.017](https://www.sciencedirect.com/science/article/pii/S0168165617315420?via%3Dihub).  
12. Smit, A.F.A., Hubley, R., Green, P. (2013-2015) *RepeatMasker Open-4.0*. (http://www.repeatmasker.org).  
13. Flynn, J.M., Hubley, R., Goubert, C., Rosen, J. Clark,. A.G., Feschotte, C., Smit, A.F. (2020) RepeatModeler2 for automated genomic discovery of transposable element families. *Proc Natl Acad Sci U S A*, __117__, 9451-9457. DOI: [10.1073/pnas.1921046117](https://www.pnas.org/content/117/17/9451.short).  
14. Quinlan, A.R., Hall, I.M. (2010) BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, __26__, 841-842, doi: [10.1093/bioinformatics/btq0333](https://academic.oup.com/bioinformatics/article/26/6/841/244688).  
