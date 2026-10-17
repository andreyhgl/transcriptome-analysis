#!/usr/bin/env Rscript

#==============================================================================
# biomart_annotation.R
#
# Fetches a gene annotation table from Ensembl BioMart for one organism,
# pinned to an Ensembl release so the gene IDs match the salmon index. The
# columns are the ones deseq2_gene_plots.R reads.
#
# Input:
#   --organism  "human" or "mouse", or any species as "Mus musculus",
#               "mus_musculus", "mmusculus" or the BioMart dataset
#               "mmusculus_gene_ensembl"
#   --version   Ensembl release, e.g. 115 (0 = latest release with BioMart)
#
# Method:
#   Looks up the archive host of the release (e.g. 113 ->
#   oct2024.archive.ensembl.org), then biomaRt::useEnsembl() on that host and
#   getBM() for every gene. Since release 116, www.ensembl.org is the new
#   Ensembl site without BioMart, so biomaRt's own release lookup on www
#   fails. Needs internet access to ensembl.org.
#   The table is then tidied: genes on scaffolds and patches are dropped,
#   chromosomes get a "chr" prefix, rows are sorted by position, the
#   "[Source:...]" suffix is removed from descriptions, and gene_type2 groups
#   the biotypes (pseudogene, ncRNA, protein_coding incl. IG/TR segments).
#
# Usage:
# > biomart_annotation.R --organism mouse --version 113 --filename anno
#
# Container:
# > docker://quay.io/biocontainers/bioconductor-biomart:2.68.0--r45hdfd78af_0
#==============================================================================


#---- Parse arguments ---------------------------------------------------------

parse_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) %% 2 != 0) stop("Arguments must be --key value pairs")

  is_key <- seq_along(args) %% 2 == 1
  keys <- sub("^--", "", args[is_key])
  vals <- args[!is_key]

  unknown <- setdiff(keys, names(defaults))
  if (length(unknown) > 0) {
    stop("Unknown argument(s): --", paste(unknown, collapse = ", --"))
  }

  for (i in seq_along(keys)) {
    # methods::as() preserves the class based on the default value
    defaults[[keys[i]]] <- methods::as(vals[i], class(defaults[[keys[i]]]))
  }

  defaults
}

# Default
defaults <- list(
  filename = "gene_annotation",
  organism = "human",
  version = 115L
)

opt <- parse_args(defaults)

cat("\n > Arguments:\n---------------------------------\n")
cat(sprintf(" > %-12s : %s\n", names(opt), unlist(opt)), sep = "")
cat(strrep("-", 33), "\n")


#---- Environment -------------------------------------------------------------

suppressPackageStartupMessages({
  library(biomaRt)
})


#---- Settings ----------------------------------------------------------------

# Output column = BioMart attribute; names match deseq2_gene_plots.R
attributes <- c(
  gene_id = "ensembl_gene_id",
  gene_name = "external_gene_name",
  chr = "chromosome_name",
  start = "start_position",
  end = "end_position",
  strand = "strand",
  gene_type = "gene_biotype",
  description = "description"
)

# Small RNA biotypes collapsed into "ncRNA" in gene_type2
ncrna_types <- c("snoRNA", "misc_RNA", "sRNA", "scaRNA", "snRNA", "scRNA")

# Columns and their order in the output table
out_cols <- c(
  "gene_name", "chr", "start", "end", "strand", "size", "gene_id",
  "description", "gene_type", "gene_type2"
)

# Common names accepted by --organism, lower case
common_names <- c(
  human = "homo_sapiens",
  mouse = "mus_musculus"
)

version <- if (opt$version > 0) opt$version else NULL

# Keep biomaRt's cache in the task directory; HOME may not be writable
Sys.setenv(BIOMART_CACHE = file.path(getwd(), "biomart_cache"))

# Output files
out_table <- paste0(opt$filename, ".tsv")


#---- Functions ---------------------------------------------------------------

print_section <- function(title) {
  cat(sprintf("\n---- %s %s\n", title, strrep("-", 72 - nchar(title))))
}

# "mouse", "Mus musculus", "mus_musculus", "mmusculus" ->
# "mmusculus_gene_ensembl". Common names are looked up in `common_names`;
# for scientific names, Ensembl abbreviates every name but the last to its
# first letter, so "Canis lupus familiaris" -> "clfamiliaris_gene_ensembl".
ensembl_dataset <- function(organism) {
  x <- sub("_gene_ensembl$", "", tolower(trimws(organism)))
  if (x %in% names(common_names)) x <- common_names[[x]]

  parts <- strsplit(x, "[ _]+")[[1]]
  n <- length(parts)

  if (n > 1) {
    x <- paste0(paste(substr(parts[-n], 1, 1), collapse = ""), parts[n])
  }

  paste0(x, "_gene_ensembl")
}

# Release name as BioMart reports it, e.g. "Ensembl Genes 113"
ensembl_release <- function(version) {
  marts <- listEnsembl(version = version)
  marts$version[marts$biomart == "genes"]
}

fetch_annotation <- function(dataset, version, attributes) {
  mart <- useEnsembl(biomart = "genes", dataset = dataset, version = version)

  # No result cache: the table is saved by the pipeline anyway
  anno <- getBM(attributes = unname(attributes), mart = mart, useCache = FALSE)

  anno <- anno[, attributes]
  colnames(anno) <- names(attributes)

  anno[!duplicated(anno$gene_id), ]
}

# Genes on assembled chromosomes only: drops unplaced scaffolds
# ("JH584299.1") and alternative haplotypes / patches ("CHR_MG4151_PATCH")
keep_primary_chromosomes <- function(anno) {
  anno[!grepl("[.]|^CHR_", anno$chr), ]
}

# "chr" prefix, strand as "+" / "-", and gene length in bp. Ensembl
# coordinates are 1-based and inclusive, hence the + 1.
format_locations <- function(anno) {
  anno$chr <- paste0("chr", anno$chr)
  anno$strand <- ifelse(anno$strand > 0, "+", "-")
  anno$size <- anno$end - anno$start + 1
  anno
}

# Numbered chromosomes in numeric order, then the rest (MT, X, Y)
# alphabetically; by start within a chromosome
sort_by_position <- function(anno) {
  chr_num <- suppressWarnings(as.integer(sub("^chr", "", anno$chr)))
  anno[order(is.na(chr_num), chr_num, anno$chr, anno$start), ]
}

# Drops Ensembl's " [Source:MGI Symbol;Acc:MGI:...]" suffix
clean_description <- function(description) {
  sub(" \\[.*\\]", "", description)
}

# Groups biotypes: every pseudogene type -> "pseudogene", small RNAs ->
# "ncRNA", IG / TR gene segments (IG_V_gene, TR_C_gene, ...) ->
# "protein_coding". Other types are kept as they are.
simplify_gene_type <- function(gene_type) {
  type2 <- gene_type
  type2[grepl("pseudo", type2)] <- "pseudogene"
  type2[type2 %in% ncrna_types] <- "ncRNA"
  type2[grepl("_gene$", type2)] <- "protein_coding"
  type2
}

tidy_annotation <- function(anno) {
  anno <- keep_primary_chromosomes(anno) |>
    format_locations() |>
    sort_by_position()

  anno$description <- clean_description(anno$description)
  anno$gene_type2 <- simplify_gene_type(anno$gene_type)

  anno[, out_cols]
}


#---- Code --------------------------------------------------------------------

dataset <- ensembl_dataset(opt$organism)
release <- ensembl_release(version)

print_section("Ensembl BioMart")
cat(sprintf(" > %-12s : %s\n", "Organism", opt$organism))
cat(sprintf(" > %-12s : %s\n", "Dataset", dataset))
cat(sprintf(" > %-12s : %s\n", "Release", release))

if (is.null(version)) {
  cat(
    " > WARNING      : no version given, using the current release. Gene",
    "IDs may not match an older salmon index.\n"
  )
}

raw <- fetch_annotation(dataset, version, attributes)
if (nrow(raw) == 0) stop("BioMart returned no genes for ", dataset)

anno <- tidy_annotation(raw)

cat(sprintf(" > %-12s : %d\n", "Fetched", nrow(raw)))
cat(sprintf(" > %-12s : %d\n", "Off-primary", nrow(raw) - nrow(anno)))
cat(sprintf(" > %-12s : %d\n", "Genes", nrow(anno)))
cat(sprintf(
  " > %-12s : %d\n",
  "Without name",
  sum(is.na(anno$gene_name) | anno$gene_name == "")
))

write.table(
  anno,
  file = out_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


#---- Done --------------------------------------------------------------------

cat(sprintf(
  "\n > Script complete!\n%s\n > Outputs:\n   - %s\n",
  strrep("-", 33),
  out_table
))