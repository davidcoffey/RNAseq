# RNAseq Pipeline v1.0
#### Updated February 19, 2016
These are a series of shell and R scripts used to process RNAseq files.  The typical work flow begins with downloading raw fastq files from the NCBI sequence read archive (http://www.ncbi.nlm.nih.gov/sra) or generating your own sequencing files.  This is followed by alignment to a reference genome.  Unaligned reads may then be aligned to alternative genomes such a pathogen genome.  After merging (for multilane samples) and processing, the resulting bam files can be run through a series of additional analyses such as GATK variant detection and STAR fusion gene detection.  Quality control analyses may also be performed on fastq files using FastQC and bam files using RNAseQC.  

Since the input of some scripts depend on the output of other scripts, they should be run in the following sequence:

1. Fetch_SRA.sh
2. FastQC.sh
3. Build_STAR_genome.sh
4. STAR_alignment.sh
5. Merge_alignment.sh (for multilane samples)
6. Process_alignment.sh
7. RNASeQC.sh
8. STAR_fusion.sh
9. GATK.sh


Unaligned fastq files may also be analyzed using MiXCR to extract antigen receptor sequences and HLA forest to predict HLA haplotypes.

* HLA_forest.sh
* MiXCR.sh
* Merge_MiXCR.R (for multilane samples)

All scripts were designed to be run on a single sample but can be looped for multiple files.  All the user needs to do is change the environmental variables and write a for loop as shown in the example below.

```
export VARIABLE="/path/to/software"

SAMPLES="001 001 003"

for S in ${SAMPLES}; do
    export SAMPLE=${S}
    script.sh
done
```

Instruction for software installation may be found in the [software installation file](https://github.com/davidcoffey/RNAseq/blob/master/Software_installation.md).