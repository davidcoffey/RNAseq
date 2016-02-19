# Aggregagate gene counts from merged STAR ReadsPerGene.out.tab files for the same sample run on seperate lanes
# Written by David Coffey dcoffey@fhcrc.org
# Updated February 9, 2016

## Prerequisites
# Install R
# Run Merge_alignments.sh to combine multiple STAR ReadsPerGene.out.tab files

# Import file from working directory
file.name <- grep("ReadsPerGene.out.tab", list.files(), value = TRUE)
ReadsPerGene <- read.delim(file = file.name, header = FALSE)

# Aggregagate gene counts from merged STAR ReadsPerGene.out.tab files
ReadsPerGene.aggregate <- aggregate(data = ReadsPerGene, .~V1, sum)

# Reorder rows so that N_ambiguous, N_multimapping, N_noFeature, and N_unmapped are listed on top
N_rows <- grep("N_", ReadsPerGene.aggregate$V1)
ReadsPerGene.aggregate <- ReadsPerGene.aggregate[c(N_rows, setdiff(1:nrow(ReadsPerGene.aggregate), N_rows)),]

# Overwrite existing file with new one
write.table(ReadsPerGene.aggregate, file = file.name, row.names = FALSE, col.names = FALSE, sep = "\t", quote = FALSE)
