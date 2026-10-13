#!/usr/bin/env Rscript

#==============================================================================
# deseq2_dds.R
#
# Locate salmon quantification files (quant.sf), import them with tximport
# and the tx2gene index, and build a DESeq2 dataset from the sample metadata.
# The design formula is built from the metadata column names in --design.
# Genes are pre-filtered on read counts and on the fraction of zero counts
# before DESeq() is run with the Wald test. The fitted object is saved as Rds.
#
# Input:
#   - metadata csv with a `sample` column and every column in --design
#   - tx2gene tsv from salmon (transcript id, gene id in the first 2 columns)
#   - quant.sf files under --quant_dir, one folder per sample
#     (the folder name is the sample id)
#
# Method:
#   DESeq2 vignette: https://bioconductor.org/packages/DESeq2
#   Design: --design "batch,id" gives ~ batch + id. The last column is the
#   variable of interest, as DESeq2 recommends, and defines the groups.
#   Gene filter: keep genes with >= --min_count reads in at least as many
#   samples as the smallest group, and < --max_zero_frac zero counts.
#
# Usage:
# > Rscript deseq2_dds.R --metadata metadata.csv --design id
# > Rscript deseq2_dds.R --design batch,id
#
# Container:
# > library://andreyhgl/singularity-r/rnaseq
#==============================================================================


#---- Parse arguments ---------------------------------------------------------

parse_args <- function(defaults, args = commandArgs(trailingOnly = TRUE)) {
  if (length(args) %% 2 != 0) stop("Arguments must be '--key value' pairs")

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
  filename = "dds",
  metadata = "metadata.csv",
  tx2gene = "salmon.merged.tx2gene.tsv",
  quant_dir = ".",
  design = "id",
  min_count = 10L,
  max_zero_frac = 0.5
)

opt <- parse_args(defaults)

cat("\n > Arguments:\n---------------------------------\n")
cat(sprintf(" > %-14s : %s\n", names(opt), unlist(opt)), sep = "")
cat(strrep("-", 33), "\n")


#---- Environment -------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(tximport)
  library(gtools)
  library(scales)
})


#---- Settings ----------------------------------------------------------------

# Design columns, separated by commas (or spaces / "+"): "batch,id"
design_cols <- strsplit(opt$design, "[,+~[:space:]]+")[[1]]
design_cols <- design_cols[design_cols != ""]
if (length(design_cols) == 0) stop("--design needs at least one column name")

design <- reformulate(design_cols)
group <- design_cols[length(design_cols)]

# Output files
out_rds <- paste0(opt$filename, ".Rds")


#---- Functions ---------------------------------------------------------------

print_section <- function(title) {
  cat("\n > ", title, "\n", strrep("-", 33), "\n", sep = "")
}

# Find all quant.sf files; the parent folder name is used as sample id
find_quant_files <- function(path) {
  files <- list.files(
    path = path,
    pattern = "^quant\\.sf$",
    recursive = TRUE,
    full.names = TRUE
  ) |>
    normalizePath()

  if (length(files) == 0) stop("No quant.sf files found in: ", path)

  names(files) <- basename(dirname(files))
  files[mixedorder(names(files))]
}

# Check required columns, harmonise sample ids and factor design variables
prepare_metadata <- function(meta, design_cols) {
  required <- unique(c("sample", design_cols))

  missing <- setdiff(required, colnames(meta))
  if (length(missing) > 0) {
    stop("Column(s) missing from metadata: ", paste(missing, collapse = ", "))
  }

  as_ordered_factor <- function(x) {
    lvls <- unique(as.character(x))
    ref <- lvls[tolower(lvls) %in% c("0", "control")]
    factor(x, levels = c(ref, setdiff(lvls, ref)))
  }

  meta[design_cols] <- lapply(meta[design_cols], as_ordered_factor)
  meta
}

# Import quant files for the metadata samples and build the DESeq2 dataset
build_dds <- function(meta, tx2gene, files, design) {
  files <- files[names(files) %in% meta$sample]
  if (length(files) == 0) stop("No quant files match the metadata samples")

  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

  # Order metadata rows as the count columns. With sample ids as rownames,
  # DESeq2 itself verifies that colData and counts line up.
  meta <- meta[match(colnames(txi$counts), meta$sample), ] |>
    droplevels()
  rownames(meta) <- meta$sample

  DESeqDataSetFromTximport(txi, colData = meta, design = design)
}

# Keep genes with >= min_count reads in at least as many samples as the
# smallest group, and with fewer than max_zero_frac zero counts
filter_genes <- function(dds, group, min_count, max_zero_frac) {
  counts <- counts(dds)

  group_sizes <- table(colData(dds)[[group]])
  min_group_size <- min(group_sizes[group_sizes > 0])

  keep_counts <- rowSums(counts >= min_count) >= min_group_size
  keep_zeros <- rowMeans(counts == 0) < max_zero_frac
  keep <- keep_counts & keep_zeros

  report_filtering(keep)
  dds[keep, ]
}

report_filtering <- function(keep) {
  n_genes <- length(keep)
  n <- c(removed = sum(!keep), kept = sum(keep))

  print_section("Gene filtering")
  msg <- sprintf(
    " > %.1f%% (%s / %s) %s\n",
    100 * n / n_genes,
    comma(n),
    comma(n_genes),
    names(n)
  )
  cat(msg, sep = "")
}


#---- Import datasets ---------------------------------------------------------

meta <- read.csv(opt$metadata) |>
  prepare_metadata(design_cols)

tx2gene <- read.table(opt$tx2gene, header = TRUE)[, 1:2] |>
  setNames(c("TXNAME", "GENEID"))

quant_files <- find_quant_files(opt$quant_dir)


#---- Code --------------------------------------------------------------------

dds <- build_dds(meta, tx2gene, quant_files, design)

print_section("DESeq2 dataset")
cat(sprintf(" > %-8s : %s\n", "Design", deparse(design)))
cat(sprintf(" > %-8s : %s\n", "Group", group))
cat(sprintf(" > %-8s : %s\n", "Levels", toString(levels(dds[[group]]))))
cat(sprintf(" > %-8s : %d of %d\n", "Samples", ncol(dds), nrow(meta)))

dds <- filter_genes(dds, group, opt$min_count, opt$max_zero_frac)
dds <- DESeq(dds, test = "Wald")
#dds_lrt <- DESeq(dds, test = "LRT", reduced = ~ 1)


#---- Done --------------------------------------------------------------------

saveRDS(dds, file = out_rds)

cat(sprintf(
  "\n > Script complete!\n%s\n > Outputs:\n   - %s\n",
  strrep("-", 33),
  out_rds
))