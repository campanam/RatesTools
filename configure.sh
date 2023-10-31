#! /bin/bash

# configure.sh script for RatesTools v1.0.0

#----------------------------------------------------------------------------------------
# Michael G. Campana and Ellie E. Armstrong, 2020-2023
# Smithsonian Institution and Stanford University

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to RatesTools;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.
 
# We politely request that this work be cited as:
# Armstrong, E.E. & M.G. Campana. 2023. RatesTools: a Nextflow pipeline for detecting
# de novo germline mutations in pedigree sequence data. Bioinformatics. 39: btac784.
# 10.1093/bioinformatics/btac784.
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
		sed -i '' "s/$stem = \"\$launchDir\/$default_file\"/$stem = \"$jar_path\"/" $filename
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
echo 'Sample information CSV?'
read csv
sed -i '' "s/libraries = \"\$launchDir\/data.csv/libraries = \"$csv/" $filename
echo 'Enter path for directory containing read files.'
read reads
reads=${reads//\//\\\/} # Escape slashes
sed -i '' "s/readDir = \"\$launchDir\/RawData\//readDir = \"$reads/" $filename
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
	sed -i '' "s/chr_file = \"\$launchDir\/chr.txt\"/chr_file = \"NULL\"/" $filename
else
	echo "Enter list of retained chromosomes."
	read chr_file
	get_jar_path Chromosome $chr_file
fi
echo 'SAMtools configuration...'
get_path_module samtools
echo 'BWA configuration...'
get_path_module bwa
echo 'Use default BWA index algorithm? (BWA auto-infers index algorithm)? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
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
echo 'Specify software for marking duplicates (sambamba, samtools or picard).'
read mkdup
while [[ $mkdup != 'sambamba' && $mkdup != 'picard' && $mkdup != 'samtools']]; do
	echo 'Unknown software. Re-enter mark duplicates selection (sambamba, samtools or picard).'
	read mkdup
done
sed -i '' "s/markDuplicates = \"picard\"/markDuplicates = \"$mkdup\"/" $filename
echo 'Remove filtered reads from final BAM alignments before genotyping (Y/N?)?'
yes_no_answer
if [ $answer == 'Y' ]; then sed -i 's/filter_bams = false/filter_bams = true/' $filename; fi
echo 'gzip configuration...'
get_path_module gzip
echo 'bgzip configuration...'
get_path_module bgzip
echo 'tabix configuration...'
get_path_module tabix
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
echo "Use Picard through Nextflow's Conda handling?"
yes_no_answer
if [ $answer == 'N' ]; then
	echo 'Testing for Picard...'
	get_jar_path Picard picard.jar
fi
echo 'Use default Java options for Picard? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
	echo 'Enter options to pass to Java.'
	read java_opts
	sed -i '' "s/picard_java = \"\"/picard_java = \"$java_opts\"/" $filename
else
	sed -i '' 's/picard_conda = false/picard_conda = true/' $filename
fi
echo "Use GATK through Nextflow's Conda handling?"
yes_no_answer
if [ $answer == 'N' ]; then
	echo 'Testing for GATK...'
	get_jar_path GATK GenomeAnalysisTK.jar
else
	sed -i '' 's/gatk_conda = false/gatk_conda = true/' $filename
fi
echo 'Enter GATK major version number (3 or 4).'
read gatkver
while [[ $gatkver != 3 && $gatkver != 4 ]]; do
	echo 'Enter GATK major version number (3 or 4).'
	read gatkver
done
if [ $gatkver != 3 ]; then
	sed -i '' "s/gatk_build = 3/gatk_build = \$gatkver/" $filename
fi
echo 'Use default Java options for GATK? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
	echo 'Enter options to pass to Java.'
	read java_opts
	sed -i '' "s/gatk_java = \"\"/gatk_java = \"$java_opts\"/" $filename
fi
echo 'Java configuration...'
get_path_module java
if `command -v java 2>&1 >/dev/null`; then # If java found, identify version
	java_version=`java -version 2>&1 >/dev/null | head -n1 | cut -d " " -f3`
	java_version=${java_version//\"/}
	echo Current Java version is $java_version.
	if [[ $java_version =~ ^1.8. || $java_version =~ ^8. || $java_version=~ ^1.17 || $java_version=~ ^17. ]]; then
		echo "Current Java environment is compatible with GATK."
	else
		echo "WARNING: GATK requires Java 1.8 or 1.17. Current Java environment is incompatible."
	fi
fi
echo 'Minimum length of contig before site filters?'
read minconlen
sed -i '' "s/min_contig_length = 1/min_contig_length = $minconlen/" $filename
echo 'Minimum length of contig after site filters?'
read minfiltconlen
sed -i '' "s/min_filt_contig_length = 1/min_filt_contig_length = $minfiltconlen/" $filename
echo 'Filter sites using VCFtools? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
	sed -i '' "s/vcftools_site_filters = \"--minDP 30 --minGQ 65 --maxDP 250 --max-missing 1 --min-alleles 1 --max-alleles 2\"/vcftools_site_filters = \"NULL\"/" $filename
else
	echo 'Use default VCFtools site filters (--minDP 30 --minGQ 65 --maxDP 250 --max-missing 1 --min-alleles 1 --max-alleles 2)? (Y/N)'
	yes_no_answer
	if [ $answer == 'N' ]; then
		echo 'Enter filters to pass to VCFtools.'
		read site_filters
		sed -i '' "s/vcftools_site_filters = \"--minDP 30 --minGQ 65 --maxDP 250 --max-missing 1 --min-alleles 1 --max-alleles 2\"/vcftools_site_filters = \"$site_filters\"/" $filename
	fi
fi
echo 'Filter sites using GATK? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
	sed -i '' "s/gatk_site_filters = \'--filterName \"filter\" --filterExpression \"QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || ReadPosRankSum < 15 || MQRankSum < -12.5\"\'/gatk_site_filters = \'NULL'/" $filename
else
	echo 'Use default GATK site filters (QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || ReadPosRankSum < 15 || MQRankSum < -12.5)? (Y/N)'
	yes_no_answer
	if [[ $answer == 'Y' && $gatkver == 4 ]]; then
		sed -i '' "s/gatk_site_filters = \'--filterName \"filter\" --filterExpression/gatk_site_filters = \'--filter-name \"filter\" --filter-expression/" $filename
	elif [[ $answer == 'N' && $gatkver == 3 ]]; then
		echo 'Enter filters to pass to GATK.'
		read gatk_site_filters
		sed -i '' "s/gatk_site_filters = \'--filterName \"filter\" --filterExpression \"QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || ReadPosRankSum < 15 || MQRankSum < -12.5\"\'/gatk_site_filters = \'--filterName \"filter\" --filterExpression \'$gatk_site_filters\'/" $filename
	elif [[ $answer == 'N' && $gatkver == 4 ]]; then
		echo 'Enter filters to pass to GATK.'
		read gatk_site_filters
		sed -i '' "s/gatk_site_filters = \'--filterName \"filter\" --filterExpression \"QUAL < 30.0 || QD < 2.0 || FS > 60.0 || MQ < 40.0 || SOR > 3.0 || ReadPosRankSum < 15 || MQRankSum < -12.5\"\'/gatk_site_filters = \'--filter-name \"filter\" --filter-expression \'$gatk_site_filters\'/" $filename
	fi
fi
echo 'Remove low-mappability and low-quality regions? (Y/N)'
yes_no_answer
if [ $answer == 'Y' ]; then
	echo 'GenMap configuration...'
	get_path_module genmap
	echo 'Use default temporary directory for GenMap (/tmp)? (Y/N)'
	yes_no_answer
	if [ $answer == 'N' ]; then
		echo 'Enter path for GenMap temporary files.'
		read gm_tmpdir
		gm_tmpdir=${gm_tmpdir//\//\\\/} # Escape slashes
		sed -i '' "s/gm_tmpdir = \'\/tmp\'/gm_tmpdir = \"$gm_tmpdir\"/" $filename
	fi
	echo 'Use default GenMap mapping options (-K 30 -E 2)? (Y/N)'
	yes_no_answer
	if [ $answer == 'N' ]; then
		echo 'Enter options to pass to GenMap map.'
		read gmopts
		sed -i "s/gm_opts = \'-K 30 -E 2\'/gm_opts = \'$gmopts\'/" $filename
	fi
	echo 'RepeatMasker configuration...'
	get_path_module RepeatMasker
	echo 'Specify RepeatMasker species.'
	read species
	sed -i '' "s/species = \"Felidae\"/species = \"$species\"/" $filename
	echo 'Use default RepeatMasker options (-gccalc -nolow -xsmall)? (Y/N)'
	yes_no_answer
	if [ $answer == 'N' ]; then
		echo 'Enter options to pass to RepeatMasker.'
		read rmopts
		sed -i "s/rm_model_opts = \'-gccalc -nolow -xsmall\'/rm_model_opts = \'rmopts\'/" $filename
	fi
	echo "Use RepeatModeler default options? (Y/N)'
	yes_no_answer
	if [ $answer == 'N' ]; then
		echo 'Enter options to pass to RepeatModeler.'
		read rmodelopts
		sed -i "s/rm_mask_opts = \'\'/rm_mask_opts = \'rmodelopts\'/" $filename
	fi
	echo 'RepeatModeler configuration...'
	get_path_module RepeatModeler
	echo 'Enter number of bases to remove on each side of an indel.'
	read indelpad
	sed -i '' "s/indelpad = 5/indelpad = $indelpad/" $filename
else
	sed -i 's/region_filter = true/region_filter = false/' $filename
fi
echo 'Use default calc_denovo_mutation_rate options (-b 100 -M 10 -w 100000 -l 100000 -S 50000 --parhom)? (Y/N)'
yes_no_answer
if [ $answer == 'N' ]; then
	echo 'Enter parameters to pass to calc_denovo_mutation_rate.'
	read dnm_opts
	sed -i '' "s/dnm_opts = \"-b 100 -M 10 -w 100000 -l 100000 -S 50000 --parhom\"/dnm_opts = \"$dnm_opts\"/" $filename
fi
echo 'Send emails regarding pipeline completion status and encountered errors? (Y/N)'
yes_no_answer
if [ $answer == 'Y' ]; then
	echo 'Enter email address.'
	read email
	sed -i '' "s/email = \"NULL\"/email = $email/" $filename
fi
