#this should be applied after autosome scaffolds are extracted to look at the distributions of DP and GQ for variant call sites\
#requires bcftools

for file in *.vcf;do
	bn=`basename $file .vcf`
	#pull variants from gvcf files
	zcat ${bn}.g.vcf.gz | grep -v 'END' > ${bn}.variants.vcf
	#pull DP and GQ (depth and quality) from new variant files
	bcftools query -f "%CHROM %POS [ %DP] [ %GQ]\n" ${bn}.variants.vcf -o ${bn}.variants.vcf
done
