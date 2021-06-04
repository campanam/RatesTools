{\rtf1\ansi\ansicpg1252\cocoartf2513
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\margl1440\margr1440\vieww16540\viewh9580\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #this should be applied after autosome scaffolds are extracted to look at the distributions of DP and GQ for variant call sites\
#requires bcftools\
\
for file in *.vcf;do\
	bn=`basename $file .vcf`\
	#pull variants from gvcf files\
	zcat $\{bn\}.g.vcf.gz | grep -v 'END' > $\{bn\}.variants.vcf\
	#pull DP and GQ (depth and quality) from new variant files\
	bcftools query -f "%CHROM %POS [ %DP] [ %GQ]\\n" $\{bn\}.variants.vcf -o $\{bn\}.variants.vcf\
done}