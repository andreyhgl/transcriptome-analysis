#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Make a ensembl table to be used for all subsequent analysis
#
# Usage: Rscript bin/ensembl.R --species --ensembl_version
# 
# > container: docker://ghcr.io/karlssonlaboratory/methylkit-env:6b7f121
#------------------------------------------------------------------------------


#---- Parse arguments ---------------------------------------------------------

# Expect two arguments to be passed
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 4) {
  stop(" > Script needs two argument! Usage: ensembl.R --species mouse --ensembl_version 115")
}

for (i in seq_along(args)) {
  if (args[i] == "--species") species <- args[i + 1]
  if (args[i] == "--ensembl_version") ensembl_version <- args[i + 1]
}


cat(paste(
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
  "Download the gene info from the ensembl database\n\n",
  " > only keep the annotated chromosomes\n",
  " > Species:\t\t", species, "\n",
  " > Biomart version:\t", ensembl_version,
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))


suppressPackageStartupMessages({
  library(biomaRt)
  library(gtools)
  library(data.table)
  library(scales)
})


#---- Functions ---------------------------------------------------------------


get_ensembl <- function(organism_dataset, ensembl_version){
  # listEnsembl() # list available datasets
  # ver <- listEnsembl()
  # ver <- unlist(strsplit(ver$version[1], " "))[3]
  # mart <- useEnsembl(biomart="genes") # download all gene lists
  # searchDatasets(mart=mart, pattern="mus") # identify house mouse dataset
  # mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
  # listAttributes(mart) # list of available attributes

  # get mart
  mart <- useEnsembl(
    biomart = "genes",
    dataset = organism_dataset,
    version = ensembl_version,
    verbose = TRUE
  ) 

  # columns to import
  cols <- c(
    "external_gene_name",
    "chromosome_name",
    "start_position",
    "end_position",
    "strand",
    "description",
    "gene_biotype",
    "ensembl_gene_id",
    "entrezgene_id"
  )

  # import & sort columns
  ens <- getBM(
    mart = mart, 
    attributes = cols,
    verbose = TRUE
  )
  ens <- ens[, cols]

  cat("\n > Building ensembl table...\n")

  # Custom naming, do not use
  colnames(ens) <- c(
    "gene_name",
    "chr",
    "start",
    "end",
    "strand",
    "gene_info",
    "gene_type",
    "ensembl_gene_id",
    "entrez_id"
  )

  # Only one entrez ID per ensembl ID, use the first
  ens <- ens[!duplicated(ens$ensembl_gene_id), ]

  ens$chr <- paste0("chr", ens$chr)
  ens$strand <- ifelse(ens$strand > 0, "+", "-")
  ens$size <- ens$end - ens$start
  
  # only save annotated chromosomes
  chrom <- grep("[.]", unique(ens$chr), value = TRUE)
  ens <- subset(ens, !chr %in% chrom)
  
  # sort chromosomes
  ens <- ens[mixedorder(paste0(ens$chr, "_", ens$start)), ]

  # remove un-needed info in gene_info column
  ens$gene_info <- sapply(ens$gene_info, function(x){gsub(" \\[.*\\]", "", x)})

  # reduce gene types ~~~~~~~~~~~~
  
  # collapse pseudogenes
  ens$gene_type2 <- ens$gene_type
  rows <- grep("pseudo", ens$gene_type2)
  ens$gene_type2[rows] <- "pseudogene"
  rows <- ens$gene_type %in% c(
    "snoRNA", "misc_RNA", "sRNA", "scaRNA", "snoRNA",
    "snRNA", "scRNA")
  ens$gene_type2[rows] <- "ncRNA"

  # 400+ IG / TR genes hid in protein coding
  #rows <- grep("_gene", ens$gene_type2)
  #ens$gene_type2[rows] <- "protein_coding"

  rows <- grep("_gene", ens$gene_type2)
  ens$gene_type2[rows] <- "protein_coding"

  ens <- ens[, c(
    "gene_name",
    "chr",
    "start",
    "end", 
    "strand",
    "size",
    "ensembl_gene_id",
    "gene_info", 
    "gene_type",
    "gene_type2"
  )]
  
  return(ens)
}


#---- Code --------------------------------------------------------------------


if (species == "mouse") organism_dataset <- "mmusculus_gene_ensembl"
if (species == "human") organism_dataset <- "hsapiens_gene_ensembl"

ens <- get_ensembl(organism_dataset, ensembl_version)

filename <- "ensembl_table.csv.gz"

fwrite(ens, file = filename)


#---- Done --------------------------------------------------------------------


cat(glue::glue("
 ~~ ensembl.R complete ~~~~~~~~~~~~~~~~~~~~~~~~
  > Output      : {filename}
  > num. genes  : {scales::comma(nrow(ens))}
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"))