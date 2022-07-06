#----------------------------------------------------------------------------------------
# logstats.sh 0.1.0
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

logval=`tail -n2 .command.log | head -n1`
if [[ $logval == 'File does not contain any sites' ]]; then
	bcftools stats $0 > tmp.txt
	bcftools stats $1 > tmp2.txt
	allval=`grep "number of records:" tmp.txt | cut -f 4`
	filtval=`grep "number of records:" tmp2.txt | cut -f 4`
	echo 'After filtering, kept '$filtval' out of a possible '$allval' Sites'
else
	echo $logval
fi