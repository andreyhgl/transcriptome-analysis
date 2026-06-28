#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# locate transcript quantification files (produced by salmon)
# import transcript index and metadata
# build dds-objects based on generations
#------------------------------------------------------------------------------

suppressPackageStartupMessages({
	library(DESeq2)
	library(tximport)
	library(gtools)
	library(scales)
})


#------------------------------------------------------------------------------


build_ddsObject <- function(meta, tx2gene, files) {
	# Analysis built on this:
	# https://bioconductor.org/packages/devel/bioc/vignettes/DESeq2/inst/doc/DESeq2.html
	metadata <- meta
	metadata$treatment <- factor(metadata$treatment)
	metadata$extraction_batch <- factor(metadata$extraction_batch)
	metadata$protocol <- factor(metadata$protocol)

	# ~~~~~~~~~~~~ Remove intercept, add "0 +" to the model ~~~~~~~~~~~~

	# USE TREATMENT AND DIET SOMEHOW TO MAKE THIS MORE SCALABLE
	# HARDCODED FOR NOW
	#grp <- factor(paste(metadata$treatment, metadata$diet, sep = "_"))
	#levels(grp) <- mixedsort(levels(grp))
	#metadata$grp <- grp
	#model <- as.formula(~ 0 + grp)
	#model <- as.formula(~ grp)
	#design <- model.matrix(model)

	#design <- ~ protocol + treatment
	design <- ~ extraction_batch + treatment

	cat(paste(
		"\n~~ processing files ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
		"Design:\t", paste0(design, collapse = ""), "\n",
		"Num. samples:\t", nrow(metadata),
		"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
	))

	files_sub <- files[names(files) %in% metadata$sample]
	txi <- tximport(files_sub, "salmon", tx2gene = tx2gene)
	
	# subset only samples found in the quant files import
	rows <- metadata$sample %in% colnames(txi$counts)
	metadata <- metadata[rows, ]

	dds <- DESeqDataSetFromTximport(
		txi,
		colData = metadata,
		design = design
	)
	#dds$grp <- grp

	# ~~ Pre-filtering steps ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

	# pre-filtering step, remove if
	# 1. < 10 counts in the smallest population group
	# 2. > 50% of values are 0

	# 1.
	#grpsize <- min(table(metadata$treatment, metadata$diet))
	out <- table(metadata$treatment)
	grpsize <- min(out[out > 0])
	keep <- rowSums(counts(dds) >= 10) >= grpsize

	# 2. Only keep rows (genes) which has < 50% of the samples w/ 0 counts
	manyzeros <- function(val) sum(val == 0) / length(val) < 0.5
	keep <- apply(counts(dds), 1, manyzeros) & keep

	# remove genes with few reads
	#keep2 <- rowSums(counts(dds)) > 20
	
	#stdev <- rowSds(counts(dds))
	#means <- rowMeans(counts(dds))
	
	#keep3 <- stdev < means

	#medians <- rowMedians(counts(dds))
	#keep4 <- stdev < medians
	# ---
	
	# decide which filtering step to use
	dds <- dds[keep, ]

	cat(paste0(
		"\n~~ Genes removed in filtering step ~~~~~~~~~~~~~~~\n\n",
		" > ", round(sum(!keep) / length(keep) * 100, 1),
		"% (", comma(sum(!keep)), " / ", comma(length(keep)), ") removed\n",
		" > ", round(sum(keep) / length(keep) * 100, 1), "% (",
		comma(sum(keep)), " / ", comma(length(keep)), ") kept\n",
		"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
	))

	# calculate diff. exp. genes
	# default test = "Wald"
	dds <- DESeq(dds, test = "Wald", quiet = TRUE)

	# likelihood ratio test requires: reduced = ~1
	#dds <- DESeq(dds, test = "LRT", reduced = ~1, quiet = TRUE)

	return(dds)
}


#------------------------------------------------------------------------------
meta <- read.csv("metadata.csv")

# import transcript index
tx2gene <- read.table("salmon.merged.tx2gene.tsv", header = TRUE)[, 1:2]
colnames(tx2gene) <- c("TXNAME", "GENEID")

# ~~ code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# locate all quantification files
files <- list.files(
	path = ".",
	pattern = "^quant.sf",
	recursive = TRUE,
	full.names = TRUE
)

# use realpath
files <- normalizePath(files)

# add sample id to the list, last folder contains sample id
names(files) <- sapply(files, function(path) basename(dirname(path)))
files <- files[mixedorder(names(files))]

# make sample names equal
meta$sample <- paste0("sperm_", meta$sample)

# build dds-objects
dds <- build_ddsObject(meta, tx2gene, files)


#---- Done --------------------------------------------------------------------


filename <- "DDS.Rds"

saveRDS(dds, file = filename)

cat(paste(
	"\n~~ DESeq2.R complete ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
	"Output:\t", filename,
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))