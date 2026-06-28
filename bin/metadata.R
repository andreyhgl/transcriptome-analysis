#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Make a metadata table to be used for all subsequent analysis
#
# Usage: Rscript bin/metadata.R --sample_info misc/sample_info.xlsx
# 
# > container: docker://ghcr.io/karlssonlaboratory/methylkit-env:6b7f121
#------------------------------------------------------------------------------


#---- Parse arguments ---------------------------------------------------------

# Expect two arguments to be passed, first mandatory
#testing:
#sample_info_file="misc/sample_info.xlsx";bad_samples_file="samples_to_remove.txt"
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(" > Script needs atleast one argument! Usage: metadata.R --sample_info sample_info.txt --bad_samples bad_samples.txt")
}

for (i in seq_along(args)) {
  if (args[i] == "--sample_info") sample_info_file <- args[i + 1]
  if (args[i] == "--bad_samples") bad_samples_file <- args[i + 1]
}

cat(paste(
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
  " > Samples info\t\t :", sample_info_file, "\n",
  #" > Bad samples?\t\t:", bad_samples_file,
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))


suppressPackageStartupMessages({
	library(openxlsx)
	library(gtools)
})


#---- metadata ----------------------------------------------------------------
sample_info <- read.xlsx(sample_info_file)[,1:7]

files <- list.files(
	path = "seqdata",
	pattern = "R1|R2",
	recursive = TRUE,
	full.names = TRUE
)

# Extract sample name, regex: /P36012_XXXX/
files <- data.frame(
	NGI.ID = sub(".*/(P36012_\\d+)/.*", "\\1", files),
	file = normalizePath(files)
)


#---- samplesheet -------------------------------------------------------------


dat <- merge(files, sample_info, by = "NGI.ID")
id <- unique(dat$sample)

samplesheet <- lapply(id, function(ID) {
  fastq <- dat[dat$sample %in% ID, ]

  data.frame(
    sample = paste0("sperm_", ID),
    fastq_1 = grep("R1", fastq$file, value = TRUE),
    fastq_2 = grep("R2", fastq$file, value = TRUE),

    #=========================================================================
    # From sequencing facility:
    # 1) Index 2 of all samples has been converted to reverse complement
    # prior to re-demultiplexing.
    #=========================================================================
    
    strandedness = "reverse"
  )
})
samplesheet <- Reduce(function(x, y) rbind(x, y), samplesheet)
samplesheet <- samplesheet[mixedorder(samplesheet$sample), ]


#---- Save files --------------------------------------------------------------
filename1 <- "metadata.csv"
write.csv(sample_info, file = filename1, row.names = FALSE, quote = FALSE)

filename2 <- "samplesheet.csv"
write.csv(samplesheet, file = filename2, row.names = FALSE, quote = FALSE)


#---- Done --------------------------------------------------------------------
cat(sprintf(
  " > Generated metadata for %s samples\n > Files saved: %s & %s\n",
  length(id), filename1, filename2
))