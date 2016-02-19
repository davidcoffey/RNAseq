#!/bin/bash
# Merge multiple STAR alignments for the same sample run on seperate lanes
# Written by David Coffey dcoffey@fhcrc.org
# Updated February 8, 2016

## Prerequisites (see Software_installation.sh)
# Download and install samtools
# Create STAR alignment for a multi-lane sample using STAR_alignment.sh script

## Variables
# export GENOME="Hg19"
# export SAMPLE=""
# export ALIGNMENT_DIRECTORY=".../STAR_alignment/$GENOME/$SAMPLE"
# export SAMTOOLS=".../samtools"
# export LAST_SAMPLE=""
# export EMAIL=""
# export AGGREGATE_READS_PER_GENE=".../Aggregate_ReadsPerGene.R"

START=`date +%s`
echo Begin Merge_alignments.sh for sample $SAMPLE on `date +"%B %d, %Y at %r"`

# Merge alignment files by lane
$SAMTOOLS merge $ALIGNMENT_DIRECTORY/$SAMPLE.Aligned.bam $ALIGNMENT_DIRECTORY/*Aligned.bam

# Merge chimeric files by lane
$SAMTOOLS merge $ALIGNMENT_DIRECTORY/$SAMPLE.Chimeric.out.bam $ALIGNMENT_DIRECTORY/*Chimeric.out.bam

# Merge SJ.out.tab
cat $ALIGNMENT_DIRECTORY/*SJ.out.tab >> $ALIGNMENT_DIRECTORY/$SAMPLE.SJ.out.tab

# Merge Chimeric.out.junction
cat $ALIGNMENT_DIRECTORY/*Chimeric.out.junction >> $ALIGNMENT_DIRECTORY/$SAMPLE.Chimeric.out.junction

# Merge Unmapped.out.mate
cat $ALIGNMENT_DIRECTORY/*Unmapped.out.mate1.gz >> $ALIGNMENT_DIRECTORY/$SAMPLE.Unmapped.out.mate1.gz
cat $ALIGNMENT_DIRECTORY/*Unmapped.out.mate2.gz >> $ALIGNMENT_DIRECTORY/$SAMPLE.Unmapped.out.mate2.gz

# Merge ReadsPerGene.out.tab
cat $ALIGNMENT_DIRECTORY/*ReadsPerGene.out.tab >> $ALIGNMENT_DIRECTORY/$SAMPLE.ReadsPerGene.out.tab
rm $ALIGNMENT_DIRECTORY/*L00*ReadsPerGene.out.tab
cd $ALIGNMENT_DIRECTORY
Rscript $AGGREGATE_READS_PER_GENE

# Clean up
rm $ALIGNMENT_DIRECTORY/*L00*Aligned.bam
rm $ALIGNMENT_DIRECTORY/*L00*Chimeric.out.bam
rm $ALIGNMENT_DIRECTORY/*L00*Chimeric.out.junction
rm $ALIGNMENT_DIRECTORY/*L00*SJ.out.tab
rm $ALIGNMENT_DIRECTORY/*L00*Unmapped.out.mate1.gz
rm $ALIGNMENT_DIRECTORY/*L00*Unmapped.out.mate2.gz

END=`date +%s`
MINUTES=$(((END-START)/60))
echo End Merge_alignments.sh for sample $SAMPLE.  The run time was $MINUTES minutes.

if [[ $SAMPLE = $LAST_SAMPLE ]]
then
	echo "The runtime was $MINUTES minutes" | mail -s "Finished Merge_alignments.sh for sample $LAST_SAMPLE" $EMAIL
fi
