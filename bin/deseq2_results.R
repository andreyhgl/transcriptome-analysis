#!/usr/bin/env Rscript

#==============================================================================
# deseq2_results.R
#
# Extract Wald test results from a fitted DESeq2 dataset (deseq2_dds.R).
# The group is the last column of the design, as set in deseq2_dds.R. Every
# level of it is compared with the reference level (the first factor level,
# e.g. control). Log2 fold changes are also shrunk with apeglm. One table per
# contrast and a summary of significant genes are written.
#
# Input:
#   - Rds with a DESeqDataSet that has been run through DESeq()
#
# Method:
#   DESeq2 vignette: https://bioconductor.org/packages/DESeq2
#   results(name = coef) for the Wald test; lfcShrink(type = "apeglm") for
#   shrunken fold changes. p-values and padj come from results(), so
#   --lfc_threshold is applied in the test itself.
#
# Usage:
# > Rscript deseq2_results.R --dds dds.Rds
# > Rscript deseq2_results.R --dds dds.Rds --alpha 0.1 --lfc_threshold 1
#
# Container:
# > docker://ghcr.io/andreyhgl/deseq2:main
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
  filename = "dds",
  dds = "dds.Rds",
  alpha = 0.05,
  lfc_threshold = 0
)

opt <- parse_args(defaults)

cat("\n > Arguments:\n---------------------------------\n")
cat(sprintf(" > %-14s : %s\n", names(opt), unlist(opt)), sep = "")
cat(strrep("-", 33), "\n")


#---- Environment -------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(apeglm)
})


#---- Settings ----------------------------------------------------------------

# Output files
out_summary <- paste0(opt$filename, "_summary.tsv")
out_table <- function(contrast) paste0(opt$filename, "_", contrast, ".tsv")


#---- Functions ---------------------------------------------------------------

print_section <- function(title) {
  cat("\n > ", title, "\n", strrep("-", 33), "\n", sep = "")
}

# Make coefficients from the tested group levels
find_coefs <- function(dds, group) {
  if (!group %in% colnames(colData(dds))) {
    stop("Group '", group, "' is not a column in the DESeq2 colData")
  }

  coefs <- resultsNames(dds)
  coefs <- coefs[startsWith(coefs, paste0(group, "_")) & grepl("_vs_", coefs)]
  if (length(coefs) == 0) {
    stop("No coefficients for '", group, "' in: ", toString(resultsNames(dds)))
  }

  coefs
}

# Wald results with the raw and the apeglm-shrunken log2 fold change
contrast_table <- function(dds, coef, alpha, lfc_threshold) {
  res <- results(
    dds,
    name = coef,
    alpha = alpha,
    lfcThreshold = lfc_threshold
  )
  shrunk <- lfcShrink(dds, coef = coef, type = "apeglm", quiet = TRUE)

  tbl <- data.frame(gene_id = rownames(res), as.data.frame(res))
  tbl$log2FoldChange_shrunk <- shrunk[rownames(res), "log2FoldChange"]
  tbl[order(tbl$padj, tbl$pvalue), ]
}

count_significant <- function(tbl, contrast, alpha) {
  sig <- !is.na(tbl$padj) & tbl$padj < alpha

  data.frame(
    contrast = contrast,
    alpha = alpha,
    up = sum(sig & tbl$log2FoldChange > 0),
    down = sum(sig & tbl$log2FoldChange < 0),
    tested = sum(!is.na(tbl$padj))
  )
}


#---- Import datasets ---------------------------------------------------------

dds <- readRDS(opt$dds)


#---- Code --------------------------------------------------------------------

# Variable of interest: the last column of the design
group <- all.vars(design(dds)) |> tail(1)

coefs <- find_coefs(dds, group)
contrasts <- sub(paste0("^", group, "_"), "", coefs)
ref <- levels(dds[[group]])[1]

print_section("Contrasts")
cat(sprintf(" > %-8s : %s\n", "Design", deparse(design(dds))))
cat(sprintf(" > %-8s : %s\n", "Group", group))
cat(sprintf(" > %-8s : %s\n", "Ref", ref))
cat(sprintf("   - %s\n", contrasts), sep = "")

summary_list <- vector("list", length(coefs))
out_files <- character(length(coefs))

for (i in seq_along(coefs)) {
  tbl <- contrast_table(dds, coefs[i], opt$alpha, opt$lfc_threshold)

  out_files[i] <- out_table(contrasts[i])
  write.table(
    tbl,
    file = out_files[i],
    sep = "\t",
    quote = FALSE,
    row.names = FALSE
  )

  summary_list[[i]] <- count_significant(tbl, contrasts[i], opt$alpha)
}

summary_tbl <- do.call(rbind, summary_list)

print_section(sprintf("Significant genes (padj < %s)", opt$alpha))
print(summary_tbl[, c("contrast", "up", "down", "tested")], row.names = FALSE)

write.table(
  summary_tbl,
  file = out_summary,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


#---- Done --------------------------------------------------------------------

cat(sprintf(
  "\n > Script complete!\n%s\n > Outputs:\n%s",
  strrep("-", 33),
  paste0("   - ", c(out_files, out_summary), "\n", collapse = "")
))