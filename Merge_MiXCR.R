# Merge multiple MiXCR output files for the same sample run on seperate lanes
# Written by David Coffey dcoffey@fhcrc.org
# Updated February 17, 2016

## Prerequisites
# Install R
# Run MiCXR.sh on samples to be merged
# export SAMPLE and LOCUS environmental variables

# Load required packages

library(plyr)
library(data.table)

# Import environmental variables
sample = Sys.getenv("SAMPLE")
locus = Sys.getenv("LOCUS")

# Import files
file.names <- grep(sample, list.files(), value = TRUE)
file.info <- file.info(file.names)
file.names <- rownames(file.info)[file.info$size > 0]
if(length(file.names) > 0){
  file.list <- plyr::llply(file.names, data.table::fread, stringsAsFactors = FALSE, data.table = FALSE)
  merged.files <- plyr::ldply(file.list, data.frame)
  
  # Merge files and aggregate by count.  The alignemnt and gene names of the sequence with the high minimum quality score is kept.
  min.qual <- by(merged.files, merged.files$nucleotide, function(x) x[which.max(x$minQuality),])
  bind <- do.call("rbind", min.qual)
  bind$count <- NULL
  aggregated <- aggregate(data = merged.files, count ~ nucleotide, FUN = sum)
  new.file <- merge(bind, aggregated, by = "nucleotide")
  new.file$frequencyCount <- new.file$count/sum(new.file$count)
  
  # Save merged file
  write.table(new.file, file = paste(sample, locus, "tsv", sep ="."), row.names = FALSE, col.names = TRUE, sep = "\t", quote = FALSE)
}
