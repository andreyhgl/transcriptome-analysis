#!/usr/bin/env Rscript

arg <- commandArgs(trailingOnly = TRUE)
ref_genome=arg[1]
#ref_genome="mouse"

cat(paste(
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n",
  "Download the gene info from the ensembl database\n\n",
  " > use version 113\n",
  " > only keep the annotated chromosomes\n",
  " > reference genome:", ref_genome,
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))

suppressPackageStartupMessages({
  library(biomaRt)
  library(gtools)
  library(readr)
  library(scales)
})


#oras://community.wave.seqera.io/library/r-biomartr_r-gtools_r-readr_r-scales:9732925029eac35c

if (ref_genome == "mouse") organism_dataset <- "mmusculus_gene_ensembl"
if (ref_genome == "human") organism_dataset <- "...."

get_ensembl <- function(organism_dataset = NULL){
  # listEnsembl() # list available datasets
  # mart <- useEnsembl(biomart="genes") # download all gene lists
  # searchDatasets(mart=mart, pattern="mus") # identify house mouse dataset
  # mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")
  # listAttributes(mart) # list of available attributes

  if (is.null(organism_dataset)) organism_dataset <- "hsapiens_gene_ensembl"

  # get mart
  mart <- useEnsembl(
    biomart = "genes",
    dataset = organism_dataset,
    #version = 113,
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
    "ensembl_gene_id"
  )

  # import & sort columns
  ens <- getBM(
    mart = mart, 
    attributes = cols,
    verbose = TRUE
  )
  ens <- ens[, cols]

  cat(" building ensembl table...\n")

  colnames(ens) <- c(
    "gene_name",
    "chr",
    "start",
    "end",
    "strand",
    "gene_info", 
    "gene_type",
    "gene_id"
  )

  ens$chr <- paste0("chr", ens$chr)
  ens$strand <- ifelse(ens$strand > 0, "+", "-")
  ens$size <- ens$end - ens$start

  # only save annotated chromosomes
  chrom <- grep("[.]", unique(ens$chr), value = TRUE)
  ens <- subset(ens, !chr %in% chrom)
  
  # sort chromosomes
  ens <- ens[mixedorder(paste0(ens$chr, "_", ens$start)), ]

  # remove un-needed info in gene_info column
  ens$gene_info <- sapply(ens$gene_info, function(x){
    gsub(" \\[.*\\]", "", x)
  })

  # reduce gene types
  
  # collapse pseudogenes
  ens$gene_type2 <- ens$gene_type
  rows <- grep("pseudo", ens$gene_type2)
  ens$gene_type2[rows] <- "pseudogene"
  rows <- ens$gene_type %in% c(
    "snoRNA", "misc_RNA", "sRNA", "scaRNA", "snoRNA",
    "snRNA", "scRNA")
  ens$gene_type2[rows] <- "ncRNA"

  # 400+ IG / TR genes hid in protein coding
  rows <- grep("_gene", ens$gene_type2)
  ens$gene_type2[rows] <- "protein_coding"

  ens <- ens[, c(
    "gene_name",
    "chr",
    "start",
    "end", 
    "strand",
    "size",
    "gene_id",
    "gene_info", 
    "gene_type",
    "gene_type2"
  )]
  
  return(ens)
}

ens <- get_ensembl(organism_dataset)

filename <- "ensembl_dataset.csv.gz"

write_csv(ens, file = filename)

cat(paste(
  "\n~~ ensembl_dataset.R complete ~~~~~~~~~~~~~~~~\n\n",
  "Output:\t", filename, "\n",
  "# genes:\t", comma(nrow(ens)),
  "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
))