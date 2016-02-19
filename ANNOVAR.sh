#!/bin/bash
# Annotate VCF files using ANNOVAR
# Written by David Coffey dcoffey@fhcrc.org
# Updated February 16, 2016

## Prerequisites (see Software_installation.sh)
# Download and install ANNOVAR
# Installed required perl databases
# Build cosmic76 database

## Variables
# export ANNOVAR=".../annovar/table_annovar.pl"
# export ANNOVAR_DATABASES=".../annovar/humandb/"
# export IGVTOOLS=".../IGVTools/2.3.26/igvtools"
# export SAMPLE="..."
# export LAST_SAMPLE="..."
# export VCF_FILE=".../$SAMPLE.filtered_calls.vcf"
# export ANNOVAR_DIRECTORY=".../ANNOVAR/"

START=`date +%s`
echo Begin ANNOVAR.sh for sample $SAMPLE on `date +"%B %d, %Y at %r"`

# Create and change to destination directory
mkdir $ANNOVAR_DIRECTORY
cd $ANNOVAR_DIRECTORY

# Annotate VCF file using multiple databases
perl $ANNOVAR $VCF_FILE $ANNOVAR_DATABASES \
-outfile $SAMPLE \
-build hg19 \
-protocol refGene,cytoBand,dgvMerged,ALL.sites.2015_08,exac03,dbnsfp30a,esp6500si_all,clinvar_20140929,snp138,nci60,cosmic70,cosmic76 \
-operation g,r,r,f,f,f,f,f,f,f,f,f \
-remove \
-vcfinput

# Clean up
rm $SAMPLE.avinput
mv $SAMPLE.hg19_multianno.txt $SAMPLE.filtered_calls_annotated.txt
mv $SAMPLE.hg19_multianno.vcf $SAMPLE.filtered_calls_annotated.vcf

# Index vcf file
$IGVTOOLS index $SAMPLE.filtered_calls_annotated.vcf

END=`date +%s`
MINUTES=$(((END-START)/60))
echo End ANNOVAR.sh for sample $SAMPLE.  The run time was $MINUTES minutes.

if [[ $SAMPLE = $LAST_SAMPLE ]]
then
	echo "The runtime was $MINUTES minutes" | mail -s "Finished ANNOVAR.sh for sample $LAST_SAMPLE" $EMAIL
fi
