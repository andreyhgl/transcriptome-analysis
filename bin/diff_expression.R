#!/usr/bin/env Rscript

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# build genexp table
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

suppressPackageStartupMessages({
	library(edgeR)
	library(gtools)
})


# ~~ functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

build_genexp_table <- function(gen){
	y <- DGE[[gen]]

	# The quasi-likelihood is recommended for bulk RNA-seq datasets due to 
	# uncertainty in dispersion estimation. (The likelihood ratio test is 
	# recommended for scRNA-seq dataset with no replicates). Uses a genewise 
	# negative binomial Generalized Linear Models (GLM)
	fit <- glmQLFit(y, design = y$design, robust = TRUE)

	# coefficients to test: DBP 10 & DBP 100
	coefficients <- colnames(y$design)[-1]
	names(coefficients) <- coefficients

	out <- lapply(coefficients, \(coef){
	  cat(paste(
	  	"Building table for:", gen, "&", coef, "\n"
	  ))

	  fit2 <- glmQLFTest(fit, coef = coef)
	  #summary(decideTests(fit2))
	  res <- topTags(fit2, n = Inf, adjust.method = "BH")
	  dat <- data.frame(res)
	  dat$gene_id <- rownames(dat)
	  dat$gen <- gen
	  dat$sign <- FALSE
	  dat$sign[dat$FDR <= 0.05 & abs(dat$logFC) > 0] <- TRUE
	  dat$type <- ifelse(dat$logFC > 0, "Up", "Down")
	  dat$coef <- coef

	  # extract gene related info
	  info <- subset(ens, gene_id %in% dat$gene_id)
	  info <- info[!duplicated(info$gene_id), ]

	  # combind and sort chromosomes
	  dat <- merge(dat, info, "gene_id")

	  # reorder columns
	  dat <- dat[, c(
	  	"gene_name",
	    "chr",
	    "start",
	    "end",
	    "strand",
	    "gene_info",
	    "gene_type",
	    "gene_type2",
	    "gene_id",
	    "sign",
	    "gen",
	    "coef",
	    "type",
	    "PValue",
	    "FDR",
	    "logCPM",
	    "logFC"
	  )]

	  dat <- dat[mixedorder(paste0(dat$chr, "_", dat$start)), ]
	  rownames(dat) <- NULL

	  return(dat)
	})

	return(out)
}

get_sign_gene_id <- function(dat){
	subset(dat, sign == TRUE)$gene_id
}

pull_genes <- function(gene){
	show(which(genes %in% gene))
	out <- lapply(generations, function(gen){
		out <- lapply(genexp_list[[gen]], function(x){
			subset(x, gene_id %in% gene & sign == TRUE)
		})
		Reduce(function(x,y) rbind(x,y), out)
	})
	Reduce(function(x,y) rbind(x,y), out)
}

# ~~ code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# Import datasets
DGE <- readRDS("DGEList.Rds")
ens <- read.csv("ensembl_dataset.csv.gz")

# Extract generatations
generations <- names(DGE)
names(generations) <- generations

# pull out pvalues for each generation and coefficient. edgeR::topTags()
genexp_list <- lapply(generations, build_genexp_table)
#genexp_list <- list_flatten(genexp_list)

#gene <- "ENSMUSG00000103377"
#a <- lapply(genexp_list, \(dat){
#	subset(dat, gene_id %in% gene)
#})
#list_rbind(a)

# find genes with significance in any generations or contrast
genes <- lapply(generations, function(gen){
	lapply(genexp_list[[gen]], get_sign_gene_id)
})
genes <- unique(unlist(genes))

# extract the results table for each significant gene in every contrast and gen
# save into 1 table
# gene_id can be duplicated if significant multiple times
genexp <- lapply(genes, pull_genes)
genexp <- Reduce(function(x,y) rbind(x,y), genexp)
rownames(genexp) <- NULL

filename <- "DGE_table.csv"

# save genexp table
write.csv(
	genexp, 
	file = filename, 
	quote = TRUE, 
	row.names = FALSE
)

cat(paste(
	"\n~~ diff-expression.R complete ~~~~~~~~~~~~\n",
	"Output:\t", filename,
	"\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
))