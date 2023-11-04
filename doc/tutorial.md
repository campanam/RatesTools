# RatesTools Tutorial and Test Dataset

__Michael G. Campana and Ellie E. Armstrong, 2019-2023__  
Smithsonian Institution  
Stanford University  

Here we provide a brief tutorial for running RatesTools. This tutorial assumes that you have already installed the RatesTools pipeline and its dependencies (See [here](https://github.com/campanam/RatesTools#installation-and-configuration) for details). The dataset (available [here](https://dx.doi.org/10.25573/data.20250288)) can be used to test your RatesTools installation and configuration. Please note that results will vary slightly from the provided final example output files due to the random number generators used by the software dependencies. The wolf sequencing reads are subset from SRA accessions SRR1518530-SRR1518532 from [BioProject PRJNA255370](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA255370) [1] and SRR9095635-SRR9095636,SRR9095638 from [BioProject PRJNA543877](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA543877/) [2]. The reference genome chromosomes derive from the domestic dog genome assembly [canFam3.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000002285.3/) [3].  

## Table of Contents  
1. [File Setup](#file-setup)  
2. [Configure the Run](#configure-the-run)  
3. [Run the Analysis](#run-the-analysis)  
4. [References](#references)  

## File Setup  
1. Create a base directory for the run, e.g. `mkdir test`.  
2. Decompress the reference genome (`dog_test.fna.gz`) (`gunzip dog_test.fna.gz`). Place the reference genome and the list of required chromosomes (`dog_test_chr.csv`) within the base directory.  
3. Decompress the `RawData.tar.gz` archive (`tar xvfz RawData.tar.gz`). Place the directory `RawData` and the FASTQ sequence files within it into the base directory.  

## Configure the Run  
1. We provide two example configurations for running the test dataset: one using the Genome Analysis Toolkit (GATK) [4] v. 3.8-1 (`wolf_gatk3.config`) and the other for GATK v. >= 4.2.3.0 (`wolf_gatk4.config`). A custom configuration can be made by executing the `configure.sh` script and completing the prompts (See [here](https://github.com/campanam/RatesTools/blob/main/README.md#configure-the-pipeline)). Please note that, in the following steps, we describe how to configure the run manually using a text editor. The `configure.sh` script automates the setting of these values except for the final configuration of the executor profile, as this is highly variable across systems.  
2. Set the run parameters. These The parameters are located in the `params` directive. Take careful notice of the value enclosure of the strings (either `' '` or `" "`). Values enclosed in `' '` do not perform variable substitution, while those using `" "` do. Note the special variable $launchDir indicates the launch directory of the run.  
2a. `refseq` specifies the path to the reference sequence (in FASTA format).  
2b. `libraries` specifies a CSV giving sample information, library names, and their corresponding read files. See (README.md#specifying-library-information) for details.  
2c. `readDir` specifies the path to the sequence FASTQ files. Each library is expected to be represented by a file containing forward reads and another containing the corresponding reverse reads.  
2d. `bwa_alg` specifies the BWA [5] index algorithm. Setting this value to an empty string ("") instructs BWA to infer the indexing algorithm.  
2e. `markDuplicates` specifies whether duplicate marking should use Picard [7] (value "picard"), SAMtools [8] or Sambamba [9] (value "sambamba").  
2f. `picard_conda` determines whether to install Picard using Conda. Set the value to true to install using Conda.  
2g. `picard` specifies the path to the Picard jar file. This can be omitted if `picard_conda` is set to true.  
2h. `picard_java` specifies additional Java parameters for Picard operations. Enter the command-line Java options as a string.  
2i. `gatk_conda` determines whether to install GATK using Conda. Set the value to true to install using Conda.  
2j. `gatk` specifies the path to the GATK jar file. This can be omitted if `gatk_conda` is set to true.  
2k. `gatk_build` specifies whether the GATK jar is major version 3 or 4.  
2l. `gatk_java` specifies additional Java parameters for GATK operations. Enter the command-line Java options as a string.  
2m. `filter_bams` determines whether to filter BAM alignments (e.g. for secondary alignments and duplicates) before genotyping. Set the value to true to filter.  
2n. `region_filter` determines whether to remove low-quality regions from the called sites. Set to true to remove regions (recommended).
2o. `phase` determines whether to phase genotypes using WhatsHap trio-phasing.  
2p. `gm_tmpdir` specifies the path to a temporary directory for GenMap [10].  
2q. `gm_opts` specifies mapping parameters (as a string) for GenMap. Use cpus in the process configuration to set number of concurrent threads.  
2r. `rm_species` specifies the species library for RepeatMasker [11].  
2s. `rm_mask_opts` specifies parameters (as a string) for RepeatMasker. Use cpus in the process configuration to set number of concurrent threads.   
2t. `rm_model_opts` specifies parameters (as a string) for RepeatModeler [12]. Use cpus in the process configuration to set number of concurrent threads and omit the `-lib` parameter.   
2u. `indelpad` specifies the number of bases up- and downstream of an indel to remove.  
2v. `prefix` specifies the file prefix for output files.  
2w. `outdir` specifies the name of the output directory.  
2x. `dam` specifies the name of the dam. This must match the name derivied from file globbing of the reads.  
2y. `sire` specifies the name of the sire. This must match the name derivied from file globbing of the reads.  
2z. `vcftools_site_filters` specifies the site-specific filters using VCFtools [13] as a string. Setting this value to "NULL" bypasses this filter. See the [VCFtools](https://vcftools.github.io/) documentation for details. We recommend restricting to biallelic sites (either using VCFtools or GATK).  
2aa. `gatk_site_filters` specifies the site-specific filters using GATK as a string. Setting this value 'NULL' bypasses this filter. See the [GATK](https://gatk.broadinstitute.org/) documentation for details. Be sure to use GATK 3 syntax for GATK 3 runs and GATK 4 syntax for GATK 4 runs. We recommend restricting to biallelic sites (either using VCFtools or GATK).  
2ab. `chr_file` specifies the path to the list of chromosomes to retain. Set the value to "NULL" to ignore this filter.  
2ac. `min_contig_length` specifies the minimum length (in bp) of contigs to retain before applying site and region filters.
2ad. `min_filt_contig_length` specifies the minimum length (in bp) of contigs to retain after applying site and region filters.
2ae. `dnm_opts` specifies the options for calc_denovo_mutation_rate.rb as a string. See the [documentation](ruby_r_scripts.md#calc_denovo_mutation_raterb) for details.  
2af. `email` specifies an email address to send alerts regarding pipeline completion, termination and errors. Set to "NULL" to turn off email alerts.  

3. Update the module list. In the `modules` directive, there is a list of software. Enter the name of any modulefiles needed for each program. If no modules are needed, leave the value as an empty string.  

4. Configure the executor profiles for your system. If you are running locally, the standard local profile provided should be sufficient (but may need some adaptation depending on your hardware). The number of threads for parallelizable software (BWA, SAMtools, GATK, GenMap, RepeatMasker, RepeatModeler) can be controlled using the $task.cpus variable for each process. Executor profiles are passed to the pipeline using the config file specified at runtime using the `-c` option. If you are running on a cluster or a cloud service, consult the [Nextflow documentation](https://www.nextflow.io/docs/latest/config.html) and your system administrators to optimize the parameter profile.  

5. Place the final configuration file into the run base directory.  

## Run the Analysis  
1. Execute the analysis using the command: `nextflow run campanam/RatesTools -r <version> -c <config_file.config>`. You can use the `-bg` option to send the Nextflow process to the background (useful for cluster systems) and the `-profile <profile_name>` option to specify your custom parameter profile.  
2. The pipeline will generate the following files in the output directory:  
2a. `01_FinalBams`: The final BAM alignments used in the analysis.  
2b. `02_gVCFs`: Individual gVCFs for each sample.  
2c. `03_CombinedVCF`: Combined joint-genotyped all-sites VCF for all samples in the analysis.  
2d. `04_RepeatMasking`: Output of RepeatMasker and RepeatModeler.  
2e. `05_ExcludedRegions`: A BED file of all genomic regions containining indels, repetitive sequence and low-mappability sequence. These regions are excluded from the analysis.  
2f. `06_FilterChrVCFs`: All-sites VCF including only the chromosomes listed in the chromosome file specified in `params.chr_list`.  
2g. `07_TrioPhasedVCFs`: All-sites VCF after trio-phasing (optional output).  
2h. `08_SplitTrioVCFs`: VCFs of all parent-offspring trios after splitting the all-sample joint-genotyping VCF.  
2i. `09_gVCFs_GP_DQ`: Summary statistics and graphs of DP and GQ for all individuals.  
2j. `10_SplitVCFs`: Trio VCFs split by chromosome for parallelization.  
2k. `11_VCFtoolsSiteFilteredVCFs`: Trio chromosome VCFs filtered at the site level by VCFtools.  
2l. `12_GATKSiteFilteredVCFs`: Trio chromosome VCFs filtered at the site level by GATK.  
2m. `13_RegionFilteredVCFs`: Trio chromosome VCFs after removal of unreliable regions.  
2n. `14_SplitCalcDNMLogs`: De novo mutations (DNMs) and DNM rates calculated per chromosome for each trio.  
2o. `15_SummarizeDNMLogs`: Summarized DNMs and DNM rates for all chromosomes for each trio.  
2p. `16_SummaryStats`: Summary statistics of the total number of sites retained after each filtration step and the number of each mutation class. It also identifies overlapping candidate mutations between siblings and recalculates mutation rate estimates assuming these sites are erroneous.  
3. We provide example output of the `SummarizeDNMLogs` and `SummaryStats` files. A successful run should produce complete versions of these files with similar results.  

## References  
1. Fan, Z., Silva, P., Gronau, I., Wang, S., Serres Armero, A., Schweizer, R.M., Ramirez, O., Pollinger, J., Galaverni, M., Ortega Del-Vecchyo, D., Du, L., Zhang, W., Zhang, Z., Xing, J., Vilà, C., Marques-Bonet, T., Godinho, R., Yue, B., Wayne, R.K. (2016) Worldwide patterns of genomic variation and admixture in gray wolves. *Genome Res*, __26__, 163-173. DOI:[10.1101/gr.197517.115](https://genome.cshlp.org/content/26/2/163.short).  
2. Koch, E.M., Schweizer, R.M., Schweizer, T.M., Stahler, D.R., Smith, D.W., Wayne, R.K., Novembre, J. (2019). De novo mutation rate estimation in wolves of known pedigree. *Mol Biol Evol*, __36__, 2536-2547, doi: [10.1093/molbev/msz159](https://academic.oup.com/mbe/article/36/11/2536/5531468?login=true).  
3. Hoeppner, M.P., Lunquist, A., Pirun, M., Meadows, J.R.S., Zamani, N., Johnson, J., Sundström, G., Cook, A., FitzGerald, M.G., Swofford, R., Mauceli, E., Moghadam, B.T., Greka, A., Alföldi, J., Abouelleil, A., Aftuck, L., Bessette, D., Berlin, A., Brown, A., Gearin, G., Lui, A., Macdonald, J.P., Priest, M., Shea, T., Turner-Maier, J., Zimmer, A., Lander, E.S., di Palma, F., Lindblad-Toh, K., Grabherr, M.G. (2014) An improved canine genome and a comprehensive catalogue of coding genes and non-coding transcripts. *PLOS One*, __9__, e91172. DOI: [10.1371/journal.pone.0091172](https://doi.org/10.1371/journal.pone.0091172).  
4. McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., DePristo, M.A. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res*, __20__, 1297-1303. DOI: [10.1101/gr.107524.110](https://genome.cshlp.org/content/20/9/1297.abstract).  
5. Li, H. (2013) Aligning sequence reads, clone sequences and assembly contigs with BWA-MEM. *arXiv*, [1303.3997v2](https://arxiv.org/abs/1303.3997).  
6. Li, H., Handsaker, B., Wysoker, A., Fennell, T., Ruan, J., Homer, N., Marth, G., Abecasis, G., Durbin, R., 1000 Genome Project Data Processing Subgroup (2009) The Sequence Alignment/Map format and SAMtools. *Bioinformatics*, 25, 2078-2079. DOI: [10.1093/bioinformatics/btp352](https://academic.oup.com/bioinformatics/article/25/16/2078/204688).  
7. Broad Institute (2020). Picard v. 2.23.8 (https://broadinstitute.github.io/picard/).  
8. Danecek, P., Bonfield, J.K., Liddle, J., Marshall, J., Ohan, V., Pollard, M.O., Whitwham, A., Keane, T., McCarthy, S.A., Davies, R.M., Li, H. (2021) Twelves years of SAMtools and BCFtools. *GigaScience*, __10__, giab008. DOI: [10.1093/gigascience/giab008](https://academic.oup.com/gigascience/article/10/2/giab008/6137722).  
9. Tarasov, A., Vilella, A.J., Cuppen, E., Nijman, I.J., Prins, P. (2015) Sambamba: fast processing of NGS alignment formats. *Bioinformatics*, __31__, 2032–2034. DOI: [10.1093/bioinformatics/btv098](https://academic.oup.com/bioinformatics/article/31/12/2032/214758).  
10. Pockrandt, C., Alzamel, M., Iliopoulos, C.S., Reinert, K. (2020) GenMap: ultra-fast computation of genome mappability. *Bioinformatics*, __36__, 3687–3692, doi: [10.1093/bioinformatics/btaa222](https://academic.oup.com/bioinformatics/article/36/12/3687/5815974).  
11. Smit, A.F.A., Hubley, R., Green, P. (2013-2015) *RepeatMasker Open-4.0*. (http://www.repeatmasker.org).  
12. Flynn, J.M., Hubley, R., Goubert, C., Rosen, J. Clark,. A.G., Feschotte, C., Smit, A.F. (2020) RepeatModeler2 for automated genomic discovery of transposable element families. *Proc Natl Acad Sci U S A*, __117__, 9451-9457. DOI: [10.1073/pnas.1921046117](https://www.pnas.org/content/117/17/9451.short).  
13. Danecek, P., Auton, A., Abecasis, G., Albers, C.A., Banks, E., DePristo, M.A., Handsaker, R.E., Lunter, G., Marth, G.T., Sherry, S.T., McVean, G., Durbin, R. (2011) The variant call format and VCFtools. *Bioinformatics*, __27__, 2156–2158. DOI: [10.1093/bioinformatics/btr330](https://academic.oup.com/bioinformatics/article/27/15/2156/402296).  