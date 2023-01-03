#!/usr/bin/env bash

#----------------------------------------------------------------------------------------
# logstats.sh 0.1.2
# Michael G. Campana and Ellie E. Armstrong, 2022-2023
# Smithsonian Institution and Stanford University

# CC0: To the extent possible under law, the Smithsonian Institution and Stanford 
# University have waived all copyright and related or neighboring rights to RatesTools;
# this work is published from the United States. You should have received a copy of the
# CC0 legal code along with this work. If not, see 
# <http://creativecommons.org/publicdomain/zero/1.0/>.
 
# We politely request that this work be cited as:
# Armstrong, E.E. & M.G. Campana. 2022. RatesTools: a Nextflow pipeline for detecting
# de novo germline mutations in pedigree sequence data. Bioinformatics. btac784.
# 10.1093/bioinformatics/btac784.
#----------------------------------------------------------------------------------------

logval=`tail -n2 $1 | head -n1`
if [[ $logval == 'File does not contain any sites' ]]; then
	bcftools stats <(gunzip -c $2) > tmp.txt
	bcftools stats <(gunzip -c $3) > tmp2.txt
	allval=`grep "number of records:" tmp.txt | cut -f 4`
	filtval=`grep "number of records:" tmp2.txt | cut -f 4`
else
	allval=`echo $logval | cut -f9 -d ' '`
	filtval=`echo $logval | cut -f4 -d ' '`
	if [ $allval -lt 0 ]; then
		bcftools stats <(gunzip -c $2) > tmp.txt
		allval=`grep "number of records:" tmp.txt | cut -f 4`
	fi
	if [ $filtval -lt 0 ]; then
		bcftools stats <(gunzip -c $3) > tmp2.txt
		filtval=`grep "number of records:" tmp2.txt | cut -f 4`
	fi
fi
# Mimic VCFtools output (Danecek et al. 2011. Bioinformatics. 27(15):2156-8. doi: 10.1093/bioinformatics/btr330).
echo 'After filtering, kept '$filtval' out of a possible '$allval' Sites'
if [[ $allval -ge $4 && $filtval -ge $5 ]]; then ln -s $3 ${3}_OK; done
