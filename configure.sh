#! /bin/bash

get_path_module () {
	local binpath=`command -v $1`
	binmodule=''
	local answer=''
	if [ ! $binpath ]; then
	echo "Executable '$1' not found in PATH."
	echo "Use module file for $1? (Y/N)"
	read answer
	while [[ $answer != 'Y' && $answer != 'N' ]]; do
		echo "Enter 'Y' or 'N'"
		read answer
	done
	if [ $answer == 'Y' ]; then
		echo "Enter $1 module."
		read binmodule
		while [[ ! `module is-avail $binmodule` && $answer == 'Y' ]]; do
			echo "Module '$binmodule' not found. Re-enter? (Y/N)"
			read answer
			while [[ $answer != 'Y' && $answer != 'N' ]]; do
				echo "Enter 'Y' or 'N'"
				read answer
			done
			if [ $answer == 'Y' ]; then
				echo "Enter $1 module."
				read binmodule
			fi
		done
		binmodule=${binmodule//\//\\\/} # Escape all backslashes in module name
		sed -i '' "s/$1 = \"\"/$1 = $binmodule/" $filename 
	fi
else
	echo "$1 executable path: "$binpath
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

echo 'SAMtools configuration...'
get_path_module samtools
echo 'BWA configuration...'
get_path_module bwa
echo 'Java configuration...' # Needs to test for 1.8
get_path_module java
echo 'Sambamba configuration...'
get_path_module sambamba
echo 'gzip configuration...'
get_path_module gzip
echo 'GenMap configuration...'
get_path_module genmap
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
