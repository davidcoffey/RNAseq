#!/bin/bash
# Download and extract files from the NCBI Sequence Read Archive http://www.ncbi.nlm.nih.gov/sra
# Written by David Coffey dcoffey@fhcrc.org
# Updated February 11, 2016

## Prerequisites (see Software_installation.sh)
# Download and install SRA toolkit

## Variables
# SRA_DIRECTORY="..."
# SAMPLES="..." # This is the sample "accession list" from the SRA download page
# SOURCE_DIRECTORY="ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/..."
# EMAIL="..."

START=`date +%s`
echo Begin Fetch_SRA.sh on `date +"%B %d, %Y at %r"`

# Create and change to destination directory
mkdir -p $SRA_DIRECTORY
cd $SRA_DIRECTORY

for S in ${SAMPLES}; do 
  SRA="${S}.sra"
  url="$SOURCE_DIRECTORY/${S}/${SRA}"
  if [ ! -e ${SRA} ]; then
  echo "Getting ${SRA}..."
  wget -nv ${url}
  else
    echo "Skipped ${SRA}..."
  fi
done

# Extract fastq files from SRA files
for S in ${SAMPLES}; do 
    SRA="${S}.sra"
	echo "Extracting ${SRA}..."
	fastq-dump --split-files -R ${SRA}
done

# Gzip fastq files
gzip $SRA_DIRECTORY/*

END=`date +%s`
MINUTES=$(((END-START)/60))
echo End Fetch_SRA.sh.  The run time was $MINUTES minutes.

echo "The runtime was $MINUTES minutes" | mail -s "Finished Fetch_SRA.sh" $EMAIL
