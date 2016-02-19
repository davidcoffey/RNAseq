#!/bin/bash
# Align RNAseq to reference genome using STAR 2-pass method with chimeric detection
# Written by David Coffey dcoffey@fhcrc.org
# Updated February 11, 2016

## Prerequisites (see Software_installation.sh)
# Download and install STAR aligner
# Download and install samtools
# Create STAR genome using STAR_genome.sh script

## Variables
# export STAR=".../STAR_2.5.0a/bin/Linux_x86_64/STAR"
# export GENOME="Hg19"
# export GENOME_DIRECTORY=".../STAR_genomes/$GENOME"
# export SAMPLE="..."
# export READ1=".../R1.fastq.gz" # Fastq file must be gzipped
# export READ2=".../R2.fastq.gz" # Fastq file must be gzipped
# export PREFIX=".../STAR_alignment/$GENOME/$SAMPLE/$SAMPLE"
# export ALIGNMENT_DIRECTORY=".../STAR_alignment/$GENOME/$SAMPLE"
# export SAMTOOLS=".../samtools"
# export PICARD=".../picard.jar"
# export RGID="L001" # Read Group ID
# export RGPL="illumina" # Read group platform
# export RGPU="HiSeq_2500" # Read group platform unit
# export RGLB="TruSeqv2" # Read group library
# export LAST_SAMPLE=""
# export EMAIL=""

START=`date +%s`
echo Begin STAR_alignment.sh for sample $SAMPLE on `date +"%B %d, %Y at %r"`

# Align FASTQ file to reference genome using STAR 2-pass method
mkdir -p $ALIGNMENT_DIRECTORY

$STAR \
--genomeDir $GENOME_DIRECTORY \
--readFilesIn $READ1 $READ2 \
--readFilesCommand zcat \
--runThreadN 4 \
--outBAMsortingThreadN 2 \
--twopassMode Basic \
--chimOutType WithinBAM \
--chimSegmentMin 20 \
--chimJunctionOverhangMin 20 \
--outSAMstrandField intronMotif \
--outFilterIntronMotifs RemoveNoncanonical \
--outReadsUnmapped Fastx \
--quantMode GeneCounts \
--outSAMtype BAM SortedByCoordinate \
--outFileNamePrefix $PREFIX.

# Convert Chimeric.out.sam to Chimeric.out.bam
$SAMTOOLS view -bS $PREFIX.Chimeric.out.sam | samtools sort - $PREFIX.Chimeric.out
rm $PREFIX.Chimeric.out.sam

# Add read groups using picard tools
java -Djava.io.tmpdir=$PREFIX.java.tmp -jar $PICARD AddOrReplaceReadGroups \
INPUT=$PREFIX.Aligned.sortedByCoord.out.bam \
OUTPUT=$PREFIX.Aligned.bam \
SORT_ORDER=coordinate \
RGID=$RGID \
RGPL=$RGPL \
RGPU=$RGPU \
RGSM=$SAMPLE \
RGLB=$RGLB

# Clean up
rm -R $PREFIX._STARgenome
rm -R $PREFIX._STARpass1
rm -R $PREFIX._STARtmp
rm -R $PREFIX.java.tmp
gzip $PREFIX.Unmapped.out.mate1
gzip $PREFIX.Unmapped.out.mate2

if [[ -f $PREFIX.Aligned.bam ]]; then
    rm $PREFIX.Aligned.sortedByCoord.out.bam
fi

END=`date +%s`
MINUTES=$(((END-START)/60))
echo End STAR_alignment.sh for sample $SAMPLE.  The run time was $MINUTES minutes.

if [[ $SAMPLE = $LAST_SAMPLE ]]
then
	echo "The runtime was $MINUTES minutes" | mail -s "Finished STAR_alignment.sh for sample $LAST_SAMPLE" $EMAIL
fi