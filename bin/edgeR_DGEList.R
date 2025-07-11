#!/usr/bin/env Rscript

cat(paste(
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
	"Differential gene expression analysis with edgeR\n",
	"Generates a DGEList-object per generation",
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))

suppressPackageStartupMessages({
	library(tximport)
	library(edgeR)
	library(gtools)
})

# ~~ Functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

DE_analysis <- function(gen, meta, tx2gene, files){
	#gen="F0"
	metadata <- meta[meta$gen %in% gen, ]
	metadata$treatment <- factor(metadata$treatment)
	
	# quick-fix for second variable
	if (length(unique(metadata$diet)) > 1){
		metadata$diet <- factor(metadata$diet)
		mod <- as.formula(~ treatment * diet)
	} else {
		mod <- as.formula(~ treatment)
	}
	design <- model.matrix(mod, data = metadata)

	cat(
		"\n~~ processing files ~~~~~~~~~~~~~~~~~~~~~~~~~\n",
		"Generation:\t", gen, "\n",
		"Num. samples:\t", nrow(metadata), "\n",
		"Design:\t", paste(mod),
		"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
	)

	files_sub <- files[names(files) %in% metadata$id]
	txi <- tximport(files_sub, "salmon", tx2gene = tx2gene)

	# from tximport vignette ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
	# https://bioconductor.org/packages/devel/bioc/vignettes/tximport/inst/doc/tximport.html#edgeR
	cts <- txi$counts
	normMat <- txi$length

	# Obtaining per-observation scaling factors for length, adjusted to avoid
	# changing the magnitude of the counts.
	normMat <- normMat / exp( rowMeans(log(normMat)) )
	normCts <- cts / normMat

	# Computing effective library sizes from scaled counts, to account for
	# composition biases between samples.
	eff.lib <- calcNormFactors(normCts) * colSums(normCts)

	# Combining effective library sizes with the length factors, and calculating
	# offsets for a log-link GLM.
	normMat <- sweep(normMat, 2, eff.lib, "*")
	normMat <- log(normMat)

	# Creating a DGEList object for use in edgeR.
	y <- DGEList(cts)
	y <- scaleOffset(y, normMat)
	# y is now ready for estimate dispersion functions see edgeR User's Guide

	# filtering using the design information
	keep <- filterByExpr(y, design)
	y <- y[keep, , keep.lib.sizes = FALSE]
	y <- estimateDisp(y, design = design, robust = TRUE)

	return(y)
}

# ~~ Import datasets ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

meta <- read.csv("metadata.csv")

# import transcript index
tx2gene <- read.table("tx2gene.tsv", header = TRUE)[, 1:2]
colnames(tx2gene) <- c("TXNAME", "GENEID")

# ~~ variables ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

generations <- unique(meta$gen)
names(generations) <- generations

# ~~ code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# locate all coverage files
files <- list.files(
	path = ".",
	pattern = "^quant.sf",
	recursive = TRUE,
	full.names = TRUE
)

# add sample id to the list, last folder contains sample id
names(files) <- sapply(files, function(path) basename(dirname(path)))
files <- files[mixedorder(names(files))]

# build DGEList
DGE <- lapply(generations, \(gen){
	DE_analysis(gen, meta, tx2gene, files)
})

filename <- "DGEList.Rds"

saveRDS(DGE, file = filename)

cat(paste(
	"\n~~ edgeR_DGEList.R complete ~~~~~~~~~~~~~~~~~~~~~\n",
	"Output:\t", filename,
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))