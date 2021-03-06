# Software installation

#### Download Hg19 genome file from Gencode
https://www.gencodegenes.org/human/release_19.html

#### Download Hg19 gtf file from Gencode
https://www.gencodegenes.org/human/release_19.html

#### Download and install picard tools
```
git clone https://github.com/broadinstitute/picard.git
cd picard/
ant clone-htsjdk
ant
```

#### Download and install  samtools
http://www.htslib.org/download/
```
cd samtools-1.3
make
make prefix=/where/to/install install
```

#### Download IGV tools
http://www.broadinstitute.org/software/igv/download

#### Download and install STAR 2.5.0a 2015/11/7
```
git clone --recursive https://github.com/alexdobin/STAR.git
cd ./STAR/source
make STAR
```

#### Download and install STAR Fusion v0.4.0
```
git clone https://github.com/STAR-Fusion/STAR-Fusion.git
cd ./STAR-Fusion
make
```

#### Download and install perl and add path to perl modules to .bashrc
https://www.perl.org/get.html
```
export PERL5LIB=$PERL5LIB:~/perl5/lib/perl5
```

#### Install Set::IntervalTree and DB_File perl modules from CPAN required by STAR fusion
```
perl -MCPAN -e shell
install Set::IntervalTree
install DB_File
```

#### Download and install GATK
```
git clone https://github.com/broadgsa/gatk-protected.git
./install-CGAT-tools.sh --cgat-scripts
./install-CGAT-tools.sh --test
```

#### Download known variant files used by GATK
```
ftp ftp.broadinstitute.org username: gsapubftp-anonymous password: <blank>
cd bundle/2.8/hg19
get dbsnp_138.hg19.vcf.gz
get Mills_and_1000G_gold_standard.indels.hg19.sites.vcf.gz
get 1000G_phase1.indels.hg19.sites.vcf.gz
bye
```
Depending on the GTF file you use, you may need to sort the vcf file and update the sequence dictionary using picard.
```
# Create reference dictionary file from fasta file
samtools faidx hg19.fasta 

# Update VCF sequence dictionary
picard UpdateVcfSequenceDictionary \
    I=Mills_and_1000G_gold_standard.indels.hg19.vcf \
    O=Mills_and_1000G_gold_standard.indels.hg19.updated.vcf \
    SEQUENCE_DICTIONARY=hg19.dict 
    
# Sort VCF using reference dictionary file
picard SortVcf \
    I=Mills_and_1000G_gold_standard.indels.hg19.updated.vcf \
    O=Mills_and_1000G_gold_standard.indels.hg19.updated.sorted.vcf \
    SEQUENCE_DICTIONARY=hg19.dict
```

#### Create fasta index and dictionary files required by GATK
```
samtools faidx .../hg19.fasta
picard CreateSequenceDictionary \
R=.../hg19.fasta \
O=.../hg19.dict
```

#### Install R package "gsalib" and "gplots" for optional base recalibration plots in GATK
```
R
install.packages("gsalib")
install.packages("gplots")
q()
```

#### Download ANNOVAR
http://annovar.openbioinformatics.org/

#### Install perl databases required for ANNOVAR
```
perl ~/Apps/annovar/annotate_variation.pl -downdb <database> -webfrom annovar ~/Apps/annovar/humandb/ -build hg19
Databases: cytoBand, 1000g2015aug, exac03, dbnsfp30a, esp6500si_all, clinvar_20150330, snp138, dgvMerged, ljb23_sift, ljb23_metasvm, cosmic70, nci60
```

#### Build Cosmic76 database for ANNOVAR
Download comsic76 file CosmicCodingMuts.vcf
https://cancer.sanger.ac.uk/cosmic/download
```
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
```

#### Download HLA forrest
```
git clone https://code.google.com/p/hlaforest/ 
```

#### Install perl modules required for HLA forest
```
perl -MCPAN -e shell
install Math::Random
```

#### Download SRA toolkit
http://www.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?cmd=show&f=software&m=software&s=software

#### Install RNA-SeQC
http://www.broadinstitute.org/cancer/cga/rnaseqc_download

#### Install MiXCR
http://mixcr.milaboratory.com

#### BED2GFT
```
git clone https://github.com/CGATOxford/cgat.git
```


