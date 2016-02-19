# RNAseq Pipeline

These are a series of shell and R scripts used to process RNAseq files.  The typical work flow begins with downloading raw fastq files from the NCBI sequence read archive (http://www.ncbi.nlm.nih.gov/sra) or generating your own sequencing files.  This is followed by alignemnt to a reference genome.  Unaligned reads may then be aligned to alternative genomes such a a viral genome.  After merging (for multilane samples) and processing, the resulting bam files can be run through a series of additional analyses such as GATK variant detection and STAR fusion gene detection.  Quality control analyses may also be performed on fastq files using FastQC and bam files using RNAseQC.  

Since the input of some scripts depend on the output of other scripts, they should be run in the following sequence:

1. Fetch_SRA.sh
2. FASTQC.sh
3. Build_hg19_STAR_genome.sh
4. Build_EBV_STAR_genome.sh (or other genome of your choice)
5. Build_EBV_STAR_genome.sh (if a gtf gene annotation file is not available)
6. STAR_alignment.sh
7. Merge_alignment.sh (may skip this step if all sample were run on a single lane)
8. Process_alignment.sh
9. RNASeQC.sh


The following scripts can be run in any order after alignment.

* STAR_fusion.sh
* GATK.sh
* ANNOVAR.sh
* MiXCR.sh (may use Merge_MiXCR.R for multilane samples)

Unaligned fastq files may also be analyzed using MiXCR to extract antigen receptor sequences and HLA forest to predict HLA haplotypes.

* HLA_forest.sh

All scripts were designed to be run on a single sample but can be looped for multiple files.  All the user needs to do is change the environment variables as shown in the example below.
```
export VARIABLE="/path/to/software"

SAMPLES="001 001 003"

for S in ${SAMPLES}; do
    export SAMPLE=${S}
    script.sh
done
```

Instruction for software installation may be found in the software installation file.