# RatesTools  

<img align="right" src="NZP-20141024-032CPM_thumb.jpg">  

__Michael G. Campana and Ellie E. Armstrong, 2019-2024__  
Smithsonian Institution  
Stanford University  

Pipeline to calculate de novo mutation rates from parent-offspring trios  

This README provides basic details for installing, configuring and running the pipeline. Please note that as of version 1.0.0, RatesTools has upgraded to Nextflow DSL2. For the original DSL1 pipeline, please see versions <=0.5.16. Detailed documentation is available for the [Ruby and R scripts](doc/ruby_r_scripts.md) included in this package and for the [pipeline's operation](doc/pipeline_details.md). Test data are provided in the [Smithsonian Institution Figshare repository](https://dx.doi.org/10.25573/data.20250288) and a tutorial is available [here](doc/tutorial.md).  

## Table of Contents  
1. [Creative Commons 0 Waiver](#creative-commons-0-waiver)  
2. [Citation](#citation)  
3. [Conda-Assisted Installation](#conda-assisted-installation)  
4. [Manual Pipeline Installation](#manual-pipeline-installation)  
5. [Configure the Pipeline](#configure-the-pipeline)  
6. [Running the Pipeline](#running-the-pipeline)  
7. [References](#references)  

## Creative Commons 0 Waiver  
![image](https://user-images.githubusercontent.com/19614608/118704084-bf02f280-b7e4-11eb-8d59-0ce648313d9e.png)  
To the extent possible under law, the Smithsonian Institution and Stanford University have waived all copyright and related or neighboring rights to RatesTools; this work is published from the United States.  

## Citation  
<img align="left" src="Lion_Project.png" width="200">  

We politely request that this work be cited as:  

Armstrong, E.E., Campana, M.G. (2023) RatesTools: a Nextflow pipeline for detecting *de novo* germline mutations in pedigree sequence data. *Bioinformatics*, __39__, btac784. DOI: [10.1093/bioinformatics/btac784](https://doi.org/10.1093/bioinformatics/btac784).  

Preprint available on *bioRxiv*. DOI: [10.1101/2022.07.18.500472](https://doi.org/10.1101/2022.07.18.500472).  

For RatesTools versions >= 1.0, please also cite:  

Armstrong, E.E., Carey, S.B., Harkess, A., Zenato Lazzari, G., Solari, K.A., Maldonado, J.E., Fleischer, R.C., Aziz, N., Walsh, P., Koepfli, K.-P., Eizirik, E., Petrov, D.A., Campana, M.G. (2024) Parameterizing Pantherinae: de novo mutation rate estimates from *Panthera* and *Neofelis* pedigrees. *bioRxiv*, 2024.04.06.587788. DOI: [10.1101/2024.04.06.587788](https://doi.org/10.1101/2024.04.06.587788).

## Conda-Assisted Installation  
We provide a configuration profile "conda" in the default configuration file (`nextflow.config`) that installs all dependencies using [Conda](https://docs.conda.io/en/latest/). As of RatesTools 1.0.0, we recommend (and default to) the use of [Mamba](https://mamba.readthedocs.io/en/latest/) for environment construction. Using this profile, the user only needs to install [Nextflow](https://www.nextflow.io/) [1], Conda/Mamba and the RatesTools pipeline:  

Install Nextflow: `curl -s https://get.nextflow.io | bash`  
Install Conda (and/or Mamba): See installation instructions [here](https://conda.io/projects/conda/en/latest/user-guide/install/index.html) and [here](https://mamba.readthedocs.io/en/latest/installation.html#installation)  
Pull the current version of the RatesTools pipeline: `nextflow pull campanam/RatesTools -r main`  

## Manual Pipeline Installation  
We explicitly list software dependencies here as no installation system (e.g. via Conda or containerization) is universally supported across all computing architectures.  

### Install Nextflow, Ruby and R  
RatesTools requires [Nextflow](https://www.nextflow.io/) [1] v. >= 23.10.0, [Ruby](http://www.ruby-lang.org) v. >= 3.2.2, [R](https://www.r-project.org/) [2] v. 4.0.2 and [Bash](https://www.gnu.org/software/bash/) v. >= 4.2.46(2)-release. Basic instructions for installing these languages are copied below. We recommend installing Ruby using the [Ruby Version Manager](https://rvm.io). See the official language documentation should you need help installing these languages.  

Install Nextflow: `curl -s https://get.nextflow.io | bash`   
Install the latest Ruby using Ruby Version Manager: `curl -sSL https://get.rvm.io | bash -s stable --ruby`  
Install R: Use the appropriate precompiled binary/installer available at the [Comprehensive R Archive Network (CRAN)](https://cran.r-project.org).  

### Install the RatesTools Scripts  
Pull the current version of the pipeline: `nextflow pull campanam/RatesTools -r main`  
To specify another RatesTools release, replace `main` with the RatesTools release version (e.g. `v0.5.7`).  

### Install the External Dependencies  
RatesTools requires the following external dependencies. See the documentation for these programs for their installation requirements. RatesTools requires the Genome Analysis Toolkit (GATK) [3] v. 3.8-1 or v. >= 4.4.0.0 and Java v. 1.8 (GATK3) or v. 1.17 (GATK4). Currently, RatesTools is not compatible with other versions of Java. Otherwise, listed versions are those that have been tested and confirmed, but other versions may work. RatesTools can utilize [Environment Modules](http://modules.sourceforge.net/) modulefiles to simplify deployment on computing clusters and limit dependency conflicts (See the [tutorial](doc/tutorial.md)).  

* gzip  
* awk  
* sed  
* zcat  
* [BWA](http://bio-bwa.sourceforge.net/) [4] v. 0.7.17  
* [SAMtools](http://www.htslib.org/) [5,6] v. 1.18  
* [BCFtools](http://www.htslib.org/) [5,6] v. 1.18  
* [bgzip and tabix from HTSlib](http://www.htslib.org/) [6] v. 1.18  
* [Java](https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) v. 1.8 (GATK3) or v. 1.17 (GATK4)  
* [Picard](https://broadinstitute.github.io/picard/) [7] v. 2.23.8 (GATK3) or v. 3.1.07 (GATK4)  
* [Sambamba](https://lomereiter.github.io/sambamba/) [8] v. 0.8.2  
* [Genome Analysis Toolkit (GAKTK)](https://github.com/broadgsa/gatk) v. 3.8-1 or v. 4.4.0.0  
* [VCFtools](https://vcftools.github.io/index.html) [9] v.0.1.16  
* [GenMap](https://github.com/cpockrandt/genmap) [10] v.1.2.0 with [SeqAn](https://github.com/seqan/seqan/tree/f548b50705be3f824a65a696943ea90a390564ce) [11] v. 2.4.1  
* [RepeatMasker](http://www.repeatmasker.org/) [12] v. 4.1.5   
* [RepeatModeler](http://www.repeatmasker.org/RepeatModeler/) [13] v. 2.0.5  
* [BEDTools](https://bedtools.readthedocs.io/en/latest/) [14] v. 2.31.0  
* [WhatsHap](https://whatshap.readthedocs.io/en/latest/index.html) [15,16] v. 2.1  

RatesTools requires the following R packages installed in your R environment:  
* [tidyverse](https://www.tidyverse.org/) [17] v. 1.3.1 with [dplyr](https://CRAN.R-project.org/package=dplyr) [18] v. 1.0.7  and [ggplot2](https://ggplot2.tidyverse.org/) [19] v. 3.3.5.
* [data.table](https://rdatatable.gitlab.io/data.table/) [20] v. 1.14.2  
* [Hmisc](https://CRAN.R-project.org/package=Hmisc) [21] v.5.1-1  

### Conda-Assisted Installation of Java Dependencies  
To assist installation and execution of the Java dependencies, we provide built-in options to install GATK and Picard through Conda. See the [tutorial](doc/tutorial.md) for details.  

## Configure the Pipeline    
Assisted configuration of the RatesTools pipeline can be accomplished using the `configure.sh` bash script included with this repository. The script copies the `nextflow.config` included with this repository and modifies the copy for the target system. The `configure.sh` script detects software installed on the local system and prompts the user to provide modulefiles, paths to undetected files, and program options. The configuration file can also be manually edited using a text editor. However, please note that the `configure.sh` script requires an *unmodified* `nextflow.config` file to work.  

*NB: The most straightforward method to obtain the `configure.sh` and `nextflow.config` files is to clone this repository and move the files to a desired location:*  
Clone the repository: `git clone https://github.com/campanam/RatesTools`  
Move the files: `mv RatesTools/*config* /some/path/`  
Change to the specified directory: `cd /some/path`  
Execute the script: `bash configure.sh`  

### Specifying Library Information  
To specify sample and library information to RatesTools, provide a CSV with the following header and information:  

Sample,Library,Read1,Read2  
\<samp1\>,\<lib1\>,\<lib1.R1.fq.gz\>,\<lib1.R2.fq.gz\>   
\<samp2\>,\<lib2\>,\<lib2.R1.fq.gz\>,\<lib2.R2.fq.gz\>   
\<samp2\>,\<lib3\>,\<lib3.R1.fq.gz\>,\<lib3.R2.fq.gz\>   
...  

`Sample` designates the unique sample name. `Library` is the unique library name (multiple libraries can correspond to the same sample). `Read1` and `Read2` are the forward and reverse read files (FASTQ format) respectively.  

RatesTools assumes bidirectional sequencing for each library, but allows for multiple sequenced libraries per individual. RatesTools will merge the libraries by sample name assuming the libraries are independent. If an individual library has been sequenced multiple times, concatenate the reads from the library and treat as a single bidirectionally sequenced file.  

### Platform-Specific Configuration  
Given the wide-variety of computing architectures and operating systems, we cannot provide specific optimized configurations for your computing system. The `nextflow.config` file includes an example of a 'standard' [configuration profile](https://www.nextflow.io/docs/latest/config.html#config-profiles) for a local installation using modulefiles and a 'conda' configuration that installs all dependencies using Conda. Example configuration profiles for the analyses described in [Armstrong & Campana 2023](https://doi.org/10.1093/bioinformatics/btac784) are provided in the [Figshare repository](https://dx.doi.org/10.25573/data.20250288). Please consult your computing staff to optimize the profile settings for your hardware. We recommend storing configuration profiles in a system-wide central location for access by all users.  

## Running the Pipeline  
Enter `nextflow run campanam/RatesTools -r <version> -c <config_file>` to run the pipeline, where `version` is the installed RatesTools release. Append `-resume` to restart a previous run or `-bg` to run RatesTools in the background. If you developed platform-specific configuration profiles, you can specify this using the `-profile <PROFILE>` option. See the Nextflow documentation for details. Final data are written to the specified output directory and its subdirectories.  

## References  
1. Di Tommaso, P., Chatzou, M., Floden, E.W., Prieto Barja, P., Palumbo, E., Notredame, C. (2017) Nextflow enables reproducible computational workflows. *Nat Biotechnol*, __35__, 316–319. DOI: [10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820).  
2. R Core Team (2020) *R: A language and environment for statistical computing.* R Foundation for Statistical Computing, Vienna, Austria. (https://www.r-project.org/).  
3. McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., DePristo, M.A. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res*, __20__, 1297-1303. DOI: [10.1101/gr.107524.110](https://genome.cshlp.org/content/20/9/1297.abstract).  
4. Li, H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv*, [1303.3997v2](https://arxiv.org/abs/1303.3997).  
5. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., 1000 Genome Project Data Processing Subgroup (2009) The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25, 2078-2079. DOI: [10.1093/bioinformatics/btp352](https://academic.oup.com/bioinformatics/article/25/16/2078/204688).  
6. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H. (2021) Twelve years of SAMtools and BCFtools. *GigaScience*, __10__, giab008. DOI: [10.1093/gigascience/giab008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722).  
7. Broad Institute (2020). Picard v. 2.23.8 (https://broadinstitute.github.io/picard/).  
8. Tarasov, A., Vilella, A.J., Cuppen, E., Nijman, I.J., Prins, P. (2015) Sambamba: fast processing of NGS alignment formats. *Bioinformatics*, __31__, 2032–2034. DOI: [10.1093/bioinformatics/btv098](https://academic.oup.com/bioinformatics/article/31/12/2032/214758).  
9. Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., McVean, G., Durbin, R. (2011) The variant call format and VCFtools. *Bioinformatics*, __27__, 2156–2158. DOI: [10.1093/bioinformatics/btr330](https://academic.oup.com/bioinformatics/article/27/15/2156/402296).  
10. Pockrandt, C., Alzamel, M., Iliopoulos, C.S., Reinert, K. (2020) GenMap: ultra-fast computation of genome mappability. *Bioinformatics*, __36__, 3687–3692, DOI: [10.1093/bioinformatics/btaa222](https://academic.oup.com/bioinformatics/article/36/12/3687/5815974).  
11. Reinert, K., Dadi, T.H., Ehrhardt, M., Hauswedell, H., Mehringer, S., Rahn, R., Kim. J., Pockrandt, C., Winkler, J., Siragusa, E., Urgese, G., Weese, D. (2017) The SeqAn C++ template library for efficient sequence analysis: A resource for programmers. *J Biotechnol*, __261__, 157-168. DOI: [10.1016/j.jbiotec.2017.07.017](https://www.sciencedirect.com/science/article/pii/S0168165617315420?via%3Dihub).  
12. Smit, A.F.A., Hubley, R., Green, P. (2013-2015) *RepeatMasker Open-4.0*. (http://www.repeatmasker.org).  
13. Flynn, J.M., Hubley, R., Goubert, C., Rosen, J. Clark,. A.G., Feschotte, C., Smit, A.F. (2020) RepeatModeler2 for automated genomic discovery of transposable element families. *Proc Natl Acad Sci U S A*, __117__, 9451-9457. DOI: [10.1073/pnas.1921046117](https://www.pnas.org/content/117/17/9451.short).  
14. Quinlan, A.R., Hall, I.M. (2010) BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, __26__, 841-842, DOI: [10.1093/bioinformatics/btq0333](https://academic.oup.com/bioinformatics/article/26/6/841/244688).  
15. Martin, M., Patterson, M., Garg, S., Fischer, S.O., Pisanti, N., Klau, G.W., Schoenhuth, A., Marschall, T. (2016) WhatsHap: fast and accurate read-based phasing. *BioRxiv*, DOI: [10.1101/085050](https://www.biorxiv.org/content/10.1101/085050v2).  
16. Garg, S., Martin, M., Marschall, T. (2016) Read-based phasing of related individuals. *Bioinformatics*, __32__, i234-i242, DOI: [10.1093/bioinformatics/btw276](https://academic.oup.com/bioinformatics/article/32/12/i234/2288955).  
17. Wickham, H., Averick, M., Bryan, J., Chang, W., D'Agostino McGowan, L., François, R., Grolemund, G., Hayes, A., Henry, L., Hester, J., Kuhn, M., Pedersen, T.L., Miller, E., Bache, S.M., Müller, K., Ooms, J., Robinson, D., Seidel, D.P., Spinu, V., Takahashi, K., Vaughan, D., Wilke, C., Woo, K. Yutani, H. (2019). Welcome to the Tidyverse. *J Open Source Softw*, __4__, 1686. DOI: [10.21105/joss.01686](https://joss.theoj.org/papers/10.21105/joss.01686).  
18. Wickham, H., François, R., Henry, L., Müller, K. (2021) dplyr: a grammar of data manipulation. R package version 1.0.7 (https://dplyr.tidyverse.org/).  
19. Wickham, H. (2016) *ggplot2: Elegant Graphics for Data Analysis.* Springer-Verlag, New York, USA.  
20. Dowle, M., Srinivasan, A. (2021) data.table: extension of 'data.frame'. R package version 1.14.2. (https://r-datatable.com).  
21. Harrell, F.E., Jr. (2023) Hmisc: Harrell miscellaneous. R package version 5.1-1. (https://CRAN.R-project.org/package=Hmisc).

Image Credit: Conor Mallon. 2014. Smithsonian's National Zoo & Conservation Biology Institute. Smithsonian Institution. https://nationalzoo.si.edu/object/nzp_NZP-20141024-032CPM.  
