# RatesTools Tutorial and Test Dataset

__Michael G. Campana and Ellie E. Armstrong, 2019-2022__  
Smithsonian Institution  
Stanford University  

Here we provide a brief tutorial for running RatesTools. This tutorial assumes that you have already installed the RatesTools pipeline and its dependencies (See [here](https://github.com/campanam/RatesTools#installation-and-configuration) for details). The dataset (available [here](###)) can be used to test your RatesTools installation and configuration. Please note that results will vary slightly from the provided final example output files due to the random number generators used by the software dependencies. The wolf sequencing reads are subset from SRA accessions SRR1518530-SRR1518532 from [BioProject](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA255370) [1] and SRR9095635-SRR9095636,SRR9095638 from [BioProject PRJNA543877](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA543877/) [2]. The reference genome chromosomes derive from the domestic dog genome assembly [canFam3.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_000002285.3/) [3].  

## File Setup  
1. Create a base directory for the run, e.g. `mkdir test`.  
2. Decompress the reference genome [dog_test.fna.gz](dog_test.fna.gz) (`gunzip dog_test.fna.gz`). Place the reference genome and the list of required chromosomes ([dog_test_chr.csv](dog_test_chr.csv)] within the base directory.  
3. Place the directory `RawData` within the base directory and the FASTQ sequence files within it into the base directory.  

## Configure the Run  
1. We provide two example configurations for running the test dataset: one using the Genome Analysis Toolkit (GATK) [4] v. 3.8-1 [wolf_test_gatk3.config](wolf_test_gatk3.config) and the other for GATK v. >= 4.2.3.0 [wolf_test_gatk4.config](wolf_test_gatk4.config). A custom configuration can be made by executing the `configure.sh` script and completing the prompts (See [here](https://github.com/campanam/RatesTools/blob/main/README.md#configure-the-pipeline)).  
2. Place the final configuration file into the run base directory.  

## References  
1. Fan, Z., Silva, P., Gronau, I., Wang, S., Serres Armero, A., Schweizer, R.M., Ramirez, O., Pollinger, J., Galaverni, M., Ortega Del-Vecchyo, D., Du, L., Zhang, W., Zhang, Z., Xing, J., Vilà, C., Marques-Bonet, T., Godinho, R., Yue, B., Wayne, R.K. (2016) Worldwide patterns of genomic variation and admixture in gray wolves. *Genome Res*, __26__, 163-173. DOI:[10.1101/gr.197517.115](https://genome.cshlp.org/content/26/2/163.short).  
2. Koch, E.M., Schweizer, R.M., Schweizer, T.M., Stahler, D.R., Smith, D.W., Wayne, R.K., Novembre, J. (2019). De novo mutation rate estimation in wolves of known pedigree. *Mol Biol Evol*, __36__, 2536-2547, doi: [10.1093/molbev/msz159](https://academic.oup.com/mbe/article/36/11/2536/5531468?login=true).  
3. Hoeppner, M.P., Lunquist, A., Pirun, M., Meadows, J.R.S., Zamani, N., Johnson, J., Sundström, G., Cook, A., FitzGerald, M.G., Swofford, R., Mauceli, E., Moghadam, B.T., Greka, A., Alföldi, J., Abouelleil, A., Aftuck, L., Bessette, D., Berlin, A., Brown, A., Gearin, G., Lui, A., Macdonald, J.P., Priest, M., Shea, T., Turner-Maier, J., Zimmer, A., Lander, E.S., di Palma, F., Lindblad-Toh, K., Grabherr, M.G. (2014) An improved canine genome and a comprehensive catalogue of coding genes and non-coding transcripts. *PLOS One*, __9__, e91172. DOI: [10.1371/journal.pone.0091172](https://doi.org/10.1371/journal.pone.0091172).  
4.  McKenna, A., Hanna, M., Banks, E., Sivachenko, A., Cibulskis, K., Kernytsky, A., Garimella, K., Altshuler, D., Gabriel, S., Daly, M., DePristo, M.A. (2010) The Genome Analysis Toolkit: a MapReduce framework for analyzing next-generation DNA sequencing data. *Genome Res*, __20__, 1297-1303. DOI: [10.1101/gr.107524.110](https://genome.cshlp.org/content/20/9/1297.abstract).  
