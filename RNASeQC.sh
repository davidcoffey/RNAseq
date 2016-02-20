#!/bin/bash
# Perform quality control on RNAseq alignments using RNA-SeQC
# Written by David Coffey dcoffey@fhcrc.org
# Updated February 10, 2016

## Prerequisites (see Software_installation.sh)
# Download and install STAR aligner
# Download and install RNA-SeQC
# Create STAR genome using STAR_genome.sh script
# Process RNAseq alignment using Process_alingment.sh script

## Variables
# export RNASEQC=".../RNA-SeQC_v1.1.8.jar"
# export SAMPLE_FILE="..."
# export GTF="..."
# export REFERENCE_FASTA_FILE=".../hg19.fasta"
# export GENOME="Hg19"
# export QC_DIRECTORY=".../Quality_control"
# export EMAIL="..."

START=`date +%s`
echo Begin RNASeQC.sh on `date +"%B %d, %Y at %r"`

mkdir -p $QC_DIRECTORY

java -jar $RNASEQC \
-n 1000 \
-s $SAMPLE_FILE \
-t $GTF \
-r $REFERENCE_FASTA_FILE \
-o $QC_DIRECTORY

END=`date +%s`
MINUTES=$(((END-START)/60))

echo End RNASeQC.sh.  The run time was $MINUTES minutes.

echo "The runtime was $MINUTES minutes" | mail -s "Finished RNASeQC.sh" $EMAIL
