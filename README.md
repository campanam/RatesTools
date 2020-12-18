# RatesTools  

__Michael G. Campana & Ellie E. Armstrong, 2019-2020__  
Smithsonian Conservation Biology Institute  
Stanford University  

Workflow to calculate de novo mutation rates from parent-offspring trios  

## License  
This software is available under  

## Installation and Configuration  
### Install Nextflow and Ruby  
RatesTools requires [Nextflow](https://www.nextflow.io/) [1] and [Ruby](http://www.ruby-lang.org). Basic instructions for installing these languages are copied below. We recommend installing Ruby using the [Ruby Version Manager](https://rvm.io). See the official language documentation should you need help installing these languages.  

Install Nextflow: `curl -s https://get.nextflow.io | bash`  
Install the latest Ruby using Ruby Version Manager: `curl -sSL https://get.rvm.io | bash -s stable --ruby`  

### Install the RatesTools Scripts  
Clone the repository: `git clone https://github.com/campanam/RatesTools`  
Install the scripts: `cd RatesTools; make install` 

*By default, RatesTools scripts will be installed into the ~/ratestools directory. If you wish to change the default directory, specify the INSTALL parameter, e.g.:* `make INSTALL=/path/to/some/dir install`  

### Install the External Dependencies  
RatesTools requires the following external dependencies. See the documentation for these programs for their installation requirements. RatesTools requires the Genome Analysis Toolkit (GATK) [2] version 3.8-1 and Java 1.8 (due to GATK). Currently, RatesTools is not compatible with GATK 4 or other versions of Java. Otherwise, listed versions are those that have been tested and confirmed, but other versions may work. RatesTools can utilize modules to simplify deployment on computing clusters and limit dependency conflicts.  

* gzip  
* awk  
* Ruby v. 2.6.3  
* [Java](https://www.oracle.com/java/technologies/javase/javase-jdk8-downloads.html) v. 1.8  
* [BWA](http://bio-bwa.sourceforge.net/) [3] v. 0.7.17  
* [SAMtools](http://www.htslib.org/) [4] v. 1.9  

### Configure the Workflow    
Semi-automatic configuration of the RatesTools workflow can be accomplished using the `configure.sh` bash script. The script copies and modifies the `nextflow.config` included with the repository for the target system. The `nextflow.config` file can also be manually edited using a text editor.  

## Running the Workflow  
Enter `ratestools.nf -c <config_file>` to run the pipeline. Append `-resume` to restart a previous run or `-bg` to run RatesTools in the background. If you developed platform-specific configuration profiles, you can specify this using the `-profile <PROFILE>` option. See the Nextflow documentation for details.  

## References  
1. Di Tommaso, P., Chatzou, M., Floden, E.W., Prieto Barja, P., Palumbo, E., Notredame, C. (2017) Nextflow enables reproducible computational workflows. *Nat Biotechnol*, __35__, 316–319. DOI: [10.1038/nbt.3820](https://www.nature.com/articles/nbt.3820).  
2. McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., DePristo, M.A. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res*, __20__, 1297-1303. DOI: [](https://genome.cshlp.org/content/20/9/1297.abstract).  
3. Li, H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv*, [1303.3997v2](https://arxiv.org/abs/1303.3997).  
4. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., 1000 Genome Project Data Processing Subgroup (2009) The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25, 2078-2079. DOI: [10.1093/bioinformatics/btp352](https://academic.oup.com/bioinformatics/article/25/16/2078/204688).  
