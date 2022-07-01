#! /bin/bash

#----------------------------------------------------------------------------------------
# Michael G. Campana and Ellie E. Armstrong, 2020-2022
# Smithsonian Institution and Stanford University

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to RatesTools;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.
 
# We politely request that this work be cited as:
# Campana, M.G. & E.E. Armstrong. 2020. RatesTools: Pipeline to calculate de novo
# mutation rates from parent-offspring trios. Smithsonian Institution and Stanford
# University. <https://github.com/campanam/RatesTools>.
#----------------------------------------------------------------------------------------

# Function to remove invalid Y/N responses
yes_no_answer () {
	read answer
	while [[ $answer != 'Y' && $answer != 'N' ]]; do
			echo "Enter 'Y' or 'N'"
			read answer
	done
}

# Function to get the path or module of non-jar binaries
get_path_module () {
	local binpath=`command -v $1`
	binmodule=''
	local answer=''
	if [ ! $binpath ]; then
		echo "Executable '$1' not found in PATH."
		echo "Use module file for $1? (Y/N)"
		yes_no_answer
		if [ $answer == 'Y' ]; then
			echo "Enter $1 module."
			read binmodule
			while [[ ! `module is-avail $binmodule` && $answer == 'Y' ]]; do
				echo "Module '$binmodule' not found. Re-enter? (Y/N)"
				yes_no_answer
				if [ $answer == 'Y' ]; then
					echo "Enter $1 module."
					read binmodule
				fi
			done
			binmodule=${binmodule//\//\\\/} # Escape all slashes in module name
			sed -i '' "s/$1 = \"\"/$1 = \"$binmodule\"/" $filename 
		fi
	else
		echo "$1 executable path: $binpath"
	fi
}

# Function to get path of jar binaries and single files
get_jar_path () {
	jar_path=`realpath $2`
	default_file=$2
	answer='Y'
	while [[ ! -f $jar_path && $answer == 'Y' ]]; do
		echo "$1 file not found. Re-enter? (Y/N)."
		yes_no_answer
		if [ $answer == 'Y' ]; then
			echo "Enter $1 path."
			read jar_path
			jar_path=`realpath $jar_path`
		fi
	done
	if [ -f $jar_path ]; then
	echo "$1 path: $jar_path"
		jar_path=${jar_path//\//\\\/} # Escape all slashes
		if [ $1 == 'Picard' ]; then
			stem='picard'
		elif [ $1 = 'GATK' ]; then
			stem='gatk'
		elif [ $1 = 'Refseq' ]; then
			stem='refseq'
			default_file=ref.fa
		elif [ $1 = 'Chromosome' ]; then
			stem='chr_file'
			default_file=chr.txt
		fi
		sed -i '' "s/$stem = \"\$baseDir\/$default_file\"/$stem = \"$jar_path\"/" $filename
	fi
}

echo 'Config file name?'
read filename
config_path=.
if [ ! -f nextflow.config ]; then
	echo "nextflow.config not found. Enter path to nextflow.config."
	read config_path
fi
cp $config_path/nextflow.config $filename
echo 'Enter path and file pattern for reads (See documentation).'
read reads
reads=${reads//\//\\\/} # Escape slashes
reads=${reads//\*/\\\*} # Escape asterisks
sed -i '' "s/reads = \"\$baseDir\/\*{R1,R2}_001.fastq\*/reads = \"$reads/" $filename
echo 'Sire name?'
read sire
sed -i '' "s/sire = \"SRR2\"/sire = \"$sire\"/" $filename
echo 'Dam name?'
read dam
sed -i '' "s/dam = \"SRR\"/dam = \"$dam\"/" $filename
echo 'Output directory name?'
read outdir
outdir=${outdir//\//\\\/} # Escape slashes
sed -i '' "s/outdir = \"test_results\"/outdir = \"$outdir\"/" $filename
echo "Output prefix?"
read prefix
sed -i '' "s/prefix = \"test\"/prefix = \"$prefix\"/" $filename
echo 'Reference sequence?'
read refseq
get_jar_path Refseq $refseq
echo 'Remove variants not assigned to specified chromosomes? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
	sed -i '' "s/chr_file = \"\$baseDir\/chr.txt\"/chr_file = \"NULL\"/" $filename
else
	echo "Enter list of retained chromosomes."
	read chr_file
	get_jar_path Chromosome $chr_file
fi
echo 'SAMtools configuration...'
get_path_module samtools
echo 'BWA configuration...'
get_path_module bwa
echo 'Use default BWA and SAMtools options? (20 threads, BWA auto-infers index algorithm)? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
	echo 'Enter number of BWA/SAMtools threads.'
	read bwa_threads
	sed -i '' "s/bwa_threads = 20/bwa_threads = $bwa_threads/" $filename
	echo 'Enter BWA index algorithm (bwtsw, is, rb2, auto-infer).'
	read bwa_alg
	while [[ $bwa_alg != 'bwtsw' && $bwa_alg != 'is' && $bwa_alg != 'rb2' && $bwa_alg != 'auto-infer' ]]; do
		echo 'Unknown algorithm. Re-enter BWA index algorithm (bwtsw, is, rb2, auto-infer).'
		read bwa_alg
	done
	if [ $bwa_alg != 'auto-infer' ]; then
		sed -i '' "s/bwa_alg = \"\"/bwa_alg = $bwa_alg/" $filename
	fi
fi
echo 'Sambamba configuration...'
get_path_module sambamba
echo 'Specify software for marking duplicates (sambamba or picard).'
read mkdup
while [[ $mkdup != 'sambamba' && $mkdup != 'picard' ]]; do
	echo 'Unknown software. Re-enter mark duplicates selection (sambamba or picard).'
	read mkdup
done
sed -i '' "s/markDuplicates = \"picard\"/markDuplicates = \"$mkdup\"/" $filename
echo 'gzip configuration...'
get_path_module gzip
echo 'bgzip configuration...'
get_path_module bgzip
echo 'tabix configuration...'
get_path_module tabix
echo 'GenMap configuration...'
get_path_module genmap
echo 'Use defaults for GenMap (8 threads, Temporary directory: /tmp)? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
	echo 'Enter number of GenMap threads.'
	read gm_threads
	sed -i '' "s/gm_threads = 8/gm_threads = $gm_threads/" $filename
	echo 'Enter path for GenMap temporary files.'
	read gm_tmpdir
	gm_tmpdir=${gm_tmpdir//\//\\\/} # Escape slashes
	sed -i '' "s/gm_tmpdir = \'\/tmp\'/gm_tmpdir = \"$gm_tmpdir\"/" $filename
fi
echo 'RepeatMasker configuration...'
get_path_module RepeatMasker
echo 'Specify RepeatMasker species.'
read species
sed -i '' "s/species = \"Felidae\"/species = \"$species\"/" $filename
echo 'RepeatModeler configuration...'
get_path_module RepeatModeler
echo 'Enter number of threads for RepeatMasker/RepeatModeler.'
read rm_pa
sed -i '' "s/rm_pa = 24/rm_pa = $rm_pa/" $filename
echo 'BEDTools configuration...'
get_path_module bedtools
echo 'BCFtools configuration...'
get_path_module bcftools
echo 'VCFtools configuration...'
get_path_module vcftools
echo 'awk configuration...'
get_path_module awk
echo 'R configuration...'
get_path_module R
echo 'Ruby configuration...'
get_path_module ruby
echo 'Testing for Picard...'
get_jar_path Picard picard.jar
echo 'Use default Java options for Picard? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
	echo 'Enter options to pass to Java.'
	read java_opts
	sed -i '' "s/picard_java = \"\"/picard_java = \"$java_opts\"/" $filename
fi
echo 'Testing for GATK...'
get_jar_path GATK GenomeAnalysisTK.jar
echo 'Enter GATK major version number (3 or 4).'
read $gatkver
while [[ $gatkver != '3' && $gatkver != '4' ]]; do
	echo 'Enter GATK major version number (3 or 4).'
	read $gatkver
done
if $gatkver != 3; then
	sed -i '' "s/gatk_build = 3/gatk_build = \$gatkver/" $filename
fi
echo 'Use default Java options for GATK? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
	echo 'Enter options to pass to Java.'
	read java_opts
	sed -i '' "s/gatk_java = \"\"/gatk_java = \"$java_opts\"/" $filename
fi
echo 'Enter number of nct threads for GATK.'
read gatk_nct
sed -i '' "s/gatk_nct = 16/gatk_nct = $gatk_nct/" $filename
echo 'Java configuration...'
get_path_module java
if `command -v java 2>&1 >/dev/null`; then # If java found, identify version
	java_version=`java -version 2>&1 >/dev/null | head -n1 | cut -d " " -f3`
	java_version=${java_version//\"/}
	echo Current Java version is $java_version.
	if [[ $java_version =~ ^1.8. || $java_version =~ ^8. ]]; then
		echo "Current Java environment is compatible with GATK."
	else
		echo "WARNING: GATK requires Java 1.8. Current Java environment is incompatible."
	fi
fi
echo 'Use default VCFtools site filters (--minDP 30 --minGQ 65 --maxDP 250 --max-missing 1 --min-alleles 1 --max-alleles 2)? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
	echo 'Enter filters to pass to VCFtools.'
	read site_filters
	sed -i '' "s/vcftools_site_filters = \"--minDP 30 --minGQ 65 --maxDP 250 --max-missing 1 --min-alleles 1 --max-alleles 2\"/vcftools_site_filters = \"$site_filters\"/" $filename
fi
echo 'Use default GATK site filters (QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || ReadPosRankSum < 15 || MQRankSum < -12.5)?(Y/N)'
yes_no_answer
if [[ $answer == 'Y' && $gatkver == '4' ]]; then
	sed -i '' "s/gatk_site_filters = \'--filterName \"filter\" --filterExpression/gatk_site_filters = \'--filter-name \"filter\" --filter-expression/" $filename
else if [[ $answer == 'N' && $gatkver == '3']]
	echo 'Enter filters to pass to GATK.'
	read gatk_site_filters
	sed -i '' "s/gatk_site_filters = \'--filterName \"filter\" --filterExpression \"QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || ReadPosRankSum < 15 || MQRankSum < -12.5\"\'/gatk_site_filters = \'--filterName \"filter\" --filterExpression \'$gatk_site_filters\'/" $filename
else if [[ $answer == 'N' && $gatkver == '4']]
	read gatk_site_filters
	sed -i '' "s/gatk_site_filters = \'--filterName \"filter\" --filterExpression \"QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || ReadPosRankSum < 15 || MQRankSum < -12.5\"\'/gatk_site_filters = \'--filter-name \"filter\" --filter-expression \'$gatk_site_filters\'/" $filename
fi
echo 'Use default calc_denovo_mutation_rate options (-b 100 -M 10 -w 100000 -l 100000 -S 50000 --parhom)? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
	echo 'Enter parameters to pass to calc_denovo_mutation_rate.'
	read dnm_opts
	sed -i '' "s/dnm_opts = \"-b 100 -M 10 -w 100000 -l 100000 -S 50000 --parhom\"/dnm_opts = \"$dnm_opts\"/" $filename
fi
