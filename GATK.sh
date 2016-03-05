#!/bin/bash
# Perform GATK Variant Analysis of RNAseq Data using the Broad Institute best practices 3.5
# Adapted from https://www.broadinstitute.org/gatk/guide/article?id=3891
# Written by David Coffey dcoffey@fhcrc.org
# Updated March 4, 2016

## Prerequisites (see Software_installation.sh)
# Download and install GATK
# Download and install picard tools
# Download known variant files
# If creating an optional base recalibration plot, install R package "gsalib" and "gplots"
# If creating optional VCF statistical plots, install matplotlib (http://matplotlib.org/downloads.html)

## Variables
# export PICARD=".../picard.jar"
# export GATK=".../GenomeAnalysisTK.jar"
# export ANNOVAR=".../annovar/table_annovar.pl"
# export ANNOVAR_DATABASES=".../annovar/humandb/"
# export IGVTOOLS=".../igvtools.jar"
# export BAM_FILE=".../.bam"
# export SAMPLE="..."
# export VARIANT_DIRECTORY=".../GATK/$SAMPLE"
# export SCRATCH_DIRECTORY="../SCRATCH/$SAMPLE"
# export REFERENCE_FASTA_FILE=".../hg19.fasta"
# export MILLS_1000G=".../Mills_and_1000G_gold_standard.indels.hg19.vcf"
# export PHASE1_1000G=".../1000G_phase1.indels.hg19.vcf"
# export DBSNP_138=".../dbsnp_138.hg19.vcf"
# export LAST_SAMPLE="..."
# export EMAIL="..."

START=`date +%s`
echo Begin GATK.sh for sample $SAMPLE on `date +"%B %d, %Y at %r"`

mkdir -p $VARIANT_DIRECTORY
mkdir -p $SCRATCH_DIRECTORY

# Mark duplicates using picard tools
java -Djava.io.tmpdir=$SCRATCH_DIRECTORY/$SAMPLE.java.tmp -jar $PICARD MarkDuplicates \
INPUT=$BAM_FILE \
OUTPUT=$SCRATCH_DIRECTORY/$SAMPLE.dedup.bam \
CREATE_INDEX=true \
VALIDATION_STRINGENCY=SILENT \
METRICS_FILE=$SCRATCH_DIRECTORY/$SAMPLE.output.metrics

# Splits reads into exon segments and clip any sequences overhanging into the intronic regions
java -Djava.io.tmpdir=$SCRATCH_DIRECTORY/$SAMPLE.java.tmp -jar $GATK \
-T SplitNCigarReads \
-R $REFERENCE_FASTA_FILE \
-I $SCRATCH_DIRECTORY/$SAMPLE.dedup.bam \
-o $SCRATCH_DIRECTORY/$SAMPLE.split.bam \
-rf ReassignOneMappingQuality \
-RMQF 255 \
-RMQT 60 \
-U ALLOW_N_CIGAR_READS

# Identify regions in BAM file that need to be realigned (optional)
java -Djava.io.tmpdir=$SCRATCH_DIRECTORY/$SAMPLE.java.tmp -jar $GATK \
-T RealignerTargetCreator \
-R $REFERENCE_FASTA_FILE \
-I $SCRATCH_DIRECTORY/$SAMPLE.split.bam \
-known $MILLS_1000G \
-o $SCRATCH_DIRECTORY/$SAMPLE.realignment_targets.list \
-nct 1 \
-nt 24

# Perform realignment (optional)
java -Djava.io.tmpdir=$SCRATCH_DIRECTORY/$SAMPLE.java.tmp -jar $GATK \
-T IndelRealigner \
-R $REFERENCE_FASTA_FILE \
-I $SCRATCH_DIRECTORY/$SAMPLE.split.bam \
-known $MILLS_1000G \
-targetIntervals $SCRATCH_DIRECTORY/$SAMPLE.realignment_targets.list \
-o $SCRATCH_DIRECTORY/$SAMPLE.realigned.bam \
-nct 1 \
-nt 1

# Create base quality score recalibration table
java -Djava.io.tmpdir=$SCRATCH_DIRECTORY/$SAMPLE.java.tmp -jar $GATK \
-T BaseRecalibrator \
-R $REFERENCE_FASTA_FILE \
-I $SCRATCH_DIRECTORY/$SAMPLE.split.bam \
-knownSites $DBSNP_138 \
-knownSites $MILLS_1000G \
-knownSites $PHASE1_1000G \
-o $SCRATCH_DIRECTORY/$SAMPLE.recal_data.table \
-nct 8 \
-nt 1

# Do a second pass to analyze covariation remaining after recalibration
java -Djava.io.tmpdir=$SCRATCH_DIRECTORY/$SAMPLE.java.tmp -jar $GATK \
-T BaseRecalibrator \
-R $REFERENCE_FASTA_FILE \
-I $SCRATCH_DIRECTORY/$SAMPLE.realigned.bam \
-knownSites $DBSNP_138 \
-knownSites $MILLS_1000G \
-knownSites $PHASE1_1000G \
-BQSR $SCRATCH_DIRECTORY/$SAMPLE.recal_data.table \
-o $SCRATCH_DIRECTORY/$SAMPLE.post_recal_data.table \
-nct 8 \
-nt 1

# Generate before/after base quality score recalibration plots
java -Djava.io.tmpdir=$SCRATCH_DIRECTORY/$SAMPLE.java.tmp -jar $GATK \
-T AnalyzeCovariates \
-R $REFERENCE_FASTA_FILE \
-before $SCRATCH_DIRECTORY/$SAMPLE.recal_data.table \
-after $SCRATCH_DIRECTORY/$SAMPLE.post_recal_data.table \
-plots $VARIANT_DIRECTORY/$SAMPLE.recalibration_plots.pdf

# Recalibrate bam file using base quality score recalibration table
java -Djava.io.tmpdir=$SCRATCH_DIRECTORY/$SAMPLE.java.tmp -jar $GATK \
-T PrintReads \
-R $REFERENCE_FASTA_FILE \
-I $SCRATCH_DIRECTORY/$SAMPLE.realigned.bam \
-BQSR $SCRATCH_DIRECTORY/$SAMPLE.recal_data.table \
-o $SCRATCH_DIRECTORY/$SAMPLE.recal.bam \
-nct 8 \
-nt 1

# Call variants
java -Djava.io.tmpdir=$SCRATCH_DIRECTORY/$SAMPLE.java.tmp -jar $GATK \
-T HaplotypeCaller \
-R $REFERENCE_FASTA_FILE \
-I $SCRATCH_DIRECTORY/$SAMPLE.recal.bam \
-dontUseSoftClippedBases \
-stand_call_conf 20.0 \
-stand_emit_conf 20.0 \
-o $VARIANT_DIRECTORY/$SAMPLE.unfiltered_calls.vcf

# Filter clusters of at least 3 SNPs that are within a window of 35 bases and based on Fisher Strand values (FS > 30.0) and Qual By Depth values (QD < 2.0)
java -Djava.io.tmpdir=$SCRATCH_DIRECTORY/$SAMPLE.java.tmp -jar $GATK \
-T VariantFiltration \
-R $REFERENCE_FASTA_FILE \
-V $VARIANT_DIRECTORY/$SAMPLE.unfiltered_calls.vcf \
-window 35 \
-cluster 3 \
-filterName FS \
-filter "FS > 30.0" \
-filterName QD \
-filter "QD < 2.0" \
-o $VARIANT_DIRECTORY/$SAMPLE.filtered_calls.vcf 

# Create statistical plots from Filtered_calls.vcf file (optional)
module load samtools
bcftools stats $VARIANT_DIRECTORY/$SAMPLE.filtered_calls.vcf > $VARIANT_DIRECTORY/$SAMPLE.filtered_calls.vchk
plot-vcfstats $VARIANT_DIRECTORY/$SAMPLE.filtered_calls.vchk -p $VARIANT_DIRECTORY/$SAMPLE.filtered_calls_plots/

# Annotate VCF file using multiple databases
perl $ANNOVAR $VARIANT_DIRECTORY/$SAMPLE.filtered_calls.vcf $ANNOVAR_DATABASES \
-outfile $VARIANT_DIRECTORY/$SAMPLE \
-build hg19 \
-protocol refGene,cytoBand,dgvMerged,ALL.sites.2015_08,exac03,dbnsfp30a,esp6500si_all,clinvar_20140929,snp138,nci60,cosmic70,cosmic76 \
-operation g,r,r,f,f,f,f,f,f,f,f,f \
-remove \
-vcfinput

# Rename annotated vcf file
mv $VARIANT_DIRECTORY/$SAMPLE.hg19_multianno.txt $VARIANT_DIRECTORY/$SAMPLE.filtered_calls_annotated.txt
mv $VARIANT_DIRECTORY/$SAMPLE.hg19_multianno.vcf $VARIANT_DIRECTORY/$SAMPLE.filtered_calls_annotated.vcf

# Index vcf file
cd $SCRATCH_DIRECTORY
$IGVTOOLS index $VARIANT_DIRECTORY/$SAMPLE.filtered_calls_annotated.vcf

# Clean up
rm $VARIANT_DIRECTORY/$SAMPLE.avinput
rm -R $SCRATCH_DIRECTORY

END=`date +%s`
MINUTES=$(((END-START)/60))
echo End GATK.sh for sample $SAMPLE.  The run time was $MINUTES minutes.

if [[ $SAMPLE = $LAST_SAMPLE ]]
then
	echo "The runtime was $MINUTES minutes" | mail -s "Finished GATK.sh for sample $LAST_SAMPLE" $EMAIL
fi
