#!/bin/bash
# Build STAR genome
# Written by David Coffey dcoffey@fhcrc.org
# Updated February 7, 2016

## Prerequisites (see Software_installation.sh)
# Download genome fasta file
# Download gene annotation gtf file
# Download and install STAR aligner

## Variables
# GENOME_DIRECTORY=".../STAR_genomes/Hg19"
# REFERENCE_FASTA_FILE=".../Hg19.fasta"
# REFERENCE_GTF_FILE=".../Hg19.gtf"
# STAR=".../STAR_2.5.0a/bin/Linux_x86_64/STAR"

START=`date +%s`
DATE=`date +"%B %d, %Y at %r"`

# Build Hg19-HHV4 STAR genome
$STAR \
--runMode genomeGenerate \
--genomeDir $GENOME_DIRECTORY \
--genomeFastaFiles $REFERENCE_FASTA_FILE \
--sjdbGTFfile $REFERENCE_GTF_FILE \
--sjdbOverhang 49 \
--runThreadN 12
# For small genomes include: --genomeSAindexNbases # (where # = min(14, log2(GenomeLength)/2 - 1)) 

END=`date +%s`
MINUTES=$(((END-START)/60))
echo STAR_genome.sh was executed on $DATE.
echo The run time was $MINUTES minutes.
