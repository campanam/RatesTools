#! /bin/bash

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
			binmodule=${binmodule//\//\\\/} # Escape all backslashes in module name
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
		jar_path=${jar_path//\//\\\/} # Escape all backslashes
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
echo 'Output directory name?'
read outdir
outdir=${outdir//\//\\\/} # Excape backslashes
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
echo 'Use default BWA and SAMtools options? (20 threads, BWA auto-infers index algorithm)?'
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
echo 'GenMap configuration...'
get_path_module genmap
echo 'Use defaults for GenMap (8 threads, Temporary directory: /tmp)?'
yes_no_answer
if [ $answer == 'N' ]; then
	echo 'Enter number of GenMap threads.'
	read gm_threads
	sed -i '' "s/gm_threads = 8/gm_threads = $gm_threads/" $filename
	echo 'Enter path for GenMap temporary files.'
	read gm_tmpdir
	gm_tmpdir=${gm_tmpdir//\//\\\/} # Excape backslashes
	sed -i '' "s/gm_tmpdir = \'\/tmp\'/gm_tmpdir = \"$gm_tmpdir\"/" $filename
fi
echo 'Ruby configuration...'
get_path_module ruby
echo 'RepeatMasker configuration...'
get_path_module RepeatMasker
echo 'RepeatModeler configuration...'
get_path_module RepeatModeler
echo 'VCFtools configuration...'
get_path_module vcftools
echo 'awk configuration...'
get_path_module awk
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
	if [[ $java_version =~ ^1.8. || $java_version =~ ^8. ]]; then
		echo "Current Java environment is compatible with GATK."
	else
		echo "WARNING: GATK requires Java 1.8. Current Java environment is incompatible."
	fi
fi