#!/bin/bash
# Software installation
# Written by David Coffey dcoffey@fhcrc.org
# Updated February 7, 2016

# Download Hg19 genome from UCSC
# http://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/
rsync -a -P rsync://hgdownload.cse.ucsc.edu/goldenPath/hg19/bigZips/hg19.2bit ./

# Download twoBitToFa and extract genome
http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/

# Extract Hg19 genome
chmod 700 twoBitToFa
./twoBitToFa hg19.2bit hg19.fasta

# Download Hg19 GFT file
# http://genome.ucsc.edu/cgi-bin/hgTables
#     clade: Mammal
#     genome: Human
#     assembly: Feb. 2009 (GRCh37/hg19)
#     group: Genes and Gene Predictions
#     track: UCSC Genes
#     table: knownGene
#     region: Select "genome" for the entire genome.
#     output format: GTF - gene transfer format
#     output file: enter a file name to save your results to a file, or leave blank to display results in the browser

# Picard tools
git clone https://github.com/broadinstitute/picard.git
cd picard/
ant clone-htsjdk
ant

# Samtools
http://www.htslib.org/download/
cd samtools-1.3
make
make prefix=/where/to/install install

# IGV tools
http://www.broadinstitute.org/software/igv/download

# STAR 2.5.0a 2015/11/7
git clone --recursive https://github.com/alexdobin/STAR.git
cd ./STAR/source
make STAR

# STAR Fusion v0.4.0
git clone https://github.com/STAR-Fusion/STAR-Fusion.git
cd ./STAR-Fusion
make

# Install perl and add path to perl modules to .bashrc
https://www.perl.org/get.html
export PERL5LIB=$PERL5LIB:~/perl5/lib/perl5

# Install Set::IntervalTree and DB_File perl modules from CPAN required by STAR fusion
perl -MCPAN -e shell
install Set::IntervalTree
install DB_File

# Install GATK
git clone https://github.com/broadgsa/gatk-protected.git
./install-CGAT-tools.sh --cgat-scripts
./install-CGAT-tools.sh --test

# Download known variant files used by GATK
ftp ftp.broadinstitute.org # username: gsapubftp-anonymous # password: <blank>
cd bundle/2.8/hg19
get dbsnp_138.hg19.vcf.gz
get Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
get 1000G_phase1.indels.hg19.sites.vcf.gz
bye

# Create fasta index and dictionary files required by GATK
samtools faidx .../hg19.fasta
picard CreateSequenceDictionary \
R=.../hg19.fasta \
O=.../hg19.dict

# Install R package "gsalib" and "gplots" for optional base recalibration plots in GATK
R
install.packages("gsalib")
install.packages("gplots")
q()

# Download ANNOVAR
# http://annovar.openbioinformatics.org/

# Install perl databases required for ANNOVAR
perl ~/Apps/annovar/annotate_variation.pl -downdb <database> -webfrom annovar ~/Apps/annovar/humandb/ -build hg19
# Databases: cytoBand, 1000g2015aug, exac03, dbnsfp30a, esp6500si_all, clinvar_20150330, snp138, dgvMerged, ljb23_sift, ljb23_metasvm, cosmic70, nci60

# Build Cosmic76 database for ANNOVAR
# Download comsic76 file CosmicCodingMuts.vcf (https://cancer.sanger.ac.uk/cosmic/download)
sftp "email"@sftp-cancer.sanger.ac.uk
get files/cosmic/grch37/cosmic/v76/VCF/CosmicCodingMuts.vcf
cat CosmicCodingMuts.vcf | grep -v "#" | cut -f 1,2,3,4,5 > CosmicCodingMuts.selected.columns.vcf
R
require(data.table)
cosmic76 = fread(input = "CosmicCodingMuts.selected.columns.vcf", header = FALSE, data.table = FALSE)
cosmic76$end = cosmic76$V2 + nchar(cosmic76$V4) - 1
cosmic76 = cosmic76[,c(1, 2, 6, 4, 5, 3)]
write.table(cosmic76, quote = FALSE, row.names = FALSE, col.names = FALSE, sep = "\t", file = "annovar/humandb/hg19_cosmic76.txt")
q()

# Install HLA forrest
git clone https://code.google.com/p/hlaforest/ 

# Install perl modules required for HLA forest
perl -MCPAN -e shell
install Math::Random

# Install SRA toolkit
# http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software

# Install RNA-SeQC
# http://www.broadinstitute.org/cancer/cga/rnaseqc_download

# Install MiXCR
# http://mixcr.milaboratory.com

# BED2GFT
git clone https://github.com/CGATOxford/cgat.git


