#!/usr/bin/env Rscript

#==============================================================================
# sample_qc.R
#
# Quality check of the samples before any differential expression, without
# using the experimental design. Per-sample metrics, a PCA, the association
# of each principal component with the metadata, sample-to-sample
# correlations and confounded metadata columns are saved in one Rds for the
# Quarto report (sample_qc.qmd). A summary is printed to the log and samples
# that deviate from the rest are flagged.
#
# Input:
#   - metadata csv/tsv with one row per sample; --sample_col holds the id
#   - tx2gene tsv from salmon (transcript id and gene id in the first two
#     columns, no header; a third column is used as the gene name)
#   - quant.sf files under --quant_dir, one folder per sample (the folder
#     name is the sample id). The mapping rate is read from
#     aux_info/meta_info.json in the same folder, when it is there.
#
# Method:
#   Gene counts are imported with tximport. PCA as in DESeq2::plotPCA():
#   vst(blind = TRUE), then prcomp() on the --ntop most variable genes.
#   Each of the first PCs is regressed on every metadata column and
#   technical metric (R2 and F-test p-value of lm(PC ~ variable)).
#   Correlations are Spearman, on the vst values of the genes with at least
#   --min_count reads in half of the samples; samples are clustered on
#   1 - correlation with average linkage. Pairs of categorical columns are
#   compared with Cramer's V. A sample is flagged when a metric is more than
#   --mad median absolute deviations, and at least a minimum difference
#   (min_diff in Settings), from the median, on the bad side only (low
#   library size, mapping rate, detected genes or correlation; high share
#   of reads in the 10 most abundant genes).
#
# Usage:
# > Rscript sample_qc.R --metadata metadata.csv --tx2gene tx2gene.tsv \
#     --quant_dir star_salmon
# > Rscript sample_qc.R --columns id,season,year
#
# Container:
# > docker://ghcr.io/andreyhgl/deseq2-env:f070b01
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
  prefix = "sample_qc",
  metadata = "data/metadata.csv",
  tx2gene = "results/star_salmon/salmon.merged.tx2gene.tsv",
  quant_dir = "results/star_salmon/",
  sample_col = "sample",
  columns = "",
  ntop = 500L,
  min_count = 10L,
  mad = 3
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
  library(jsonlite)
})


#---- Settings ----------------------------------------------------------------

# Principal components kept per sample, and tested against the metadata
n_pcs <- 5

# A numeric column with more unique values than this is continuous; with
# fewer it is treated as groups (e.g. year or dose)
max_numeric_levels <- 8

# Top genes reported: per PC by loading, and by share of all reads
n_top_genes <- 10

# Points of each sample's count density curve
density_n <- 256

# Smallest distance from the median that is flagged, whatever the MAD.
# Keeps tight metrics (e.g. a mapping rate of 99.1 vs 99.5 %) from giving
# trivial outliers. Units are those of the tested values.
min_diff <- c(
  library_size = log10(2), # half the median (tested on log10)
  pct_mapped = 5,          # percentage points
  detected_genes = 2000,
  top10_pct = 5,           # percentage points
  median_cor = 0.02
)

# Output files
out_rds <- paste0(opt$prefix, ".Rds")
out_table <- paste0(opt$prefix, "_samples.tsv")


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

read_table <- function(path, ...) {
  sep <- if (grepl("\\.csv(\\.gz)?$", path)) "," else "\t"
  read.delim(path, sep = sep, ...)
}

# Transcript to gene table, with a gene name when there is a third column
read_tx2gene <- function(path) {
  tx2gene <- read_table(path, header = FALSE)
  if (ncol(tx2gene) < 2) stop("tx2gene needs at least 2 columns: ", path)

  # Drop a header row, e.g. "transcript_id  gene_id  gene_name"
  if (grepl("^(tx|transcript)", tx2gene[1, 1], ignore.case = TRUE)) {
    tx2gene <- tx2gene[-1, ]
  }

  if (ncol(tx2gene) < 3) tx2gene[[3]] <- tx2gene[[2]]
  tx2gene <- tx2gene[, 1:3]
  colnames(tx2gene) <- c("tx_id", "gene_id", "gene_name")
  tx2gene
}

# Metadata with a `sample` column that matches the quant folder names
read_metadata <- function(path, sample_col) {
  meta <- read_table(path)

  if (!sample_col %in% colnames(meta)) {
    stop("Column '", sample_col, "' is missing from the metadata: ", path)
  }

  meta$sample <- as.character(meta[[sample_col]])
  dup <- unique(meta$sample[duplicated(meta$sample)])
  if (length(dup) > 0) stop("Duplicated samples in metadata: ", toString(dup))

  meta
}

# Metadata columns to check: --columns, or every column that can separate
# samples (not constant, and not unique per sample unless it is numeric)
choose_columns <- function(meta, columns, sample_col) {
  if (columns != "") {
    cols <- strsplit(columns, "[,[:space:]]+")[[1]]
    cols <- cols[cols != ""]

    missing <- setdiff(cols, colnames(meta))
    if (length(missing) > 0) {
      stop("Column(s) missing from metadata: ", toString(missing))
    }
    return(cols)
  }

  cols <- setdiff(colnames(meta), c("sample", sample_col))
  keep <- vapply(cols, function(col) {
    x <- meta[[col]]
    n_unique <- length(unique(x[!is.na(x)]))
    n_unique >= 2 && (is.numeric(x) || n_unique < sum(!is.na(x)))
  }, logical(1))

  cols[keep]
}

# "continuous" for numeric columns with many values, else "categorical"
column_type <- function(x) {
  n_unique <- length(unique(x[!is.na(x)]))
  if (is.numeric(x) && n_unique > max_numeric_levels) {
    "continuous"
  } else {
    "categorical"
  }
}

# Factor with "0" / "control" first (the reference level), then the remaining
# levels in natural order (a2 before a10). NA stays NA.
as_ordered_factor <- function(x) {
  lvls <- unique(as.character(x[!is.na(x)]))
  lvls <- lvls[mixedorder(lvls)]
  ref <- lvls[tolower(lvls) %in% c("0", "control")]
  factor(as.character(x), levels = c(ref, setdiff(lvls, ref)))
}

# Reads processed and mapped by salmon, from aux_info/meta_info.json
salmon_stats <- function(files) {
  stats <- lapply(names(files), function(sample) {
    json <- file.path(dirname(files[[sample]]), "aux_info", "meta_info.json")
    info <- if (file.exists(json)) fromJSON(json) else list()

    num <- function(field) {
      if (is.null(info[[field]])) NA_real_ else as.numeric(info[[field]])
    }

    data.frame(
      sample = sample,
      n_processed = num("num_processed"),
      pct_mapped = num("percent_mapped")
    )
  })

  do.call(rbind, stats)
}

# Import quant files into a DESeq2 dataset. The design is ~ 1: QC must not
# depend on the experimental design.
build_dds <- function(meta, tx2gene, files) {
  cat(sprintf(" > Importing %d quant files...\n", length(files)))
  t0 <- Sys.time()

  txi <- suppressMessages(
    tximport(files, type = "salmon", tx2gene = tx2gene[, 1:2])
  )

  cat(sprintf(
    " > %d genes x %d samples (%.0f s)\n",
    nrow(txi$counts),
    ncol(txi$counts),
    as.numeric(difftime(Sys.time(), t0, units = "secs"))
  ))

  # Levels of samples without a quant file are dropped
  meta <- meta[match(colnames(txi$counts), meta$sample), , drop = FALSE] |>
    droplevels()
  rownames(meta) <- meta$sample

  cat(" > Building DESeq2 dataset: (design ~ 1)...\n")
  suppressMessages(
    DESeqDataSetFromTximport(txi, colData = meta, design = ~ 1)
  )
}

# Variance-stabilised values, blind to the design
run_vst <- function(dds) {
  # vst() estimates the dispersion trend on a 1000-gene subset
  vsd <- if (nrow(dds) < 1000) {
    varianceStabilizingTransformation(dds, blind = TRUE)
  } else {
    vst(dds, blind = TRUE)
  }
  assay(vsd)
}

# PCA of the ntop most variable genes (as plotPCA)
run_pca <- function(mat, ntop) {
  top <- order(MatrixGenerics::rowVars(mat), decreasing = TRUE)
  top <- head(top, min(ntop, nrow(mat)))
  pca <- prcomp(t(mat[top, , drop = FALSE]))

  list(
    x = pca$x,
    rotation = pca$rotation,
    var_pct = 100 * pca$sdev^2 / sum(pca$sdev^2),
    n_genes = length(top)
  )
}

# Share of a sample's reads in its n most abundant genes (%)
top_share <- function(cnt, n) {
  apply(cnt, 2, function(x) {
    100 * sum(head(sort(x, decreasing = TRUE), n)) / sum(x)
  })
}

# The genes with the largest mean share of the reads over the samples
top_genes_table <- function(cnt, gene_names, n) {
  share <- 100 * t(t(cnt) / colSums(cnt))
  top <- head(order(rowMeans(share), decreasing = TRUE), n)

  data.frame(
    gene_id = rownames(cnt)[top],
    gene_name = unname(gene_names[rownames(cnt)[top]]),
    mean_pct = round(rowMeans(share)[top], 2),
    max_pct = round(apply(share[top, , drop = FALSE], 1, max), 2),
    row.names = NULL
  )
}

# The genes with the largest absolute loading on each of the first PCs
loadings_table <- function(pca, gene_names, n_pcs, n) {
  tables <- lapply(seq_len(n_pcs), function(i) {
    load <- pca$rotation[, i]
    top <- head(order(abs(load), decreasing = TRUE), n)

    data.frame(
      pc = paste0("PC", i),
      gene_id = names(load)[top],
      gene_name = unname(gene_names[names(load)[top]]),
      loading = round(unname(load[top]), 3),
      row.names = NULL
    )
  })

  do.call(rbind, tables)
}

# Density curve of log2(CPM + 1) per sample, on a shared grid
density_table <- function(cnt) {
  lcpm <- log2(1e6 * t(t(cnt) / colSums(cnt)) + 1)
  upper <- max(lcpm)

  curves <- lapply(colnames(lcpm), function(sample) {
    d <- density(lcpm[, sample], from = 0, to = upper, n = density_n)
    data.frame(sample = sample, x = d$x, y = d$y)
  })

  do.call(rbind, curves)
}

# R2 and p-value of lm(pc ~ x); NA when x cannot be tested
pc_fit <- function(pc, x) {
  na <- c(r2 = NA_real_, p = NA_real_)

  ok <- !is.na(x) & !is.na(pc)
  pc <- pc[ok]
  x <- if (is.factor(x)) droplevels(x[ok]) else x[ok]

  n_unique <- length(unique(x))
  n_params <- if (is.factor(x)) n_unique else 2
  # Needs variation in x and residual degrees of freedom
  if (n_unique < 2 || length(pc) - n_params < 1) return(na)

  fit <- tryCatch(summary(lm(pc ~ x)), error = function(e) NULL)
  if (is.null(fit) || is.null(fit$fstatistic)) return(na)

  f <- fit$fstatistic
  c(
    r2 = fit$r.squared,
    p = unname(pf(f[1], f[2], f[3], lower.tail = FALSE))
  )
}

# One row per PC and variable; padj is BH over all tests
pc_association <- function(pcs, variables) {
  rows <- lapply(colnames(pcs), function(pc) {
    fits <- vapply(variables, function(x) pc_fit(pcs[, pc], x), numeric(2))

    data.frame(
      pc = pc,
      variable = names(variables),
      r2 = fits["r2", ],
      p = fits["p", ],
      row.names = NULL
    )
  })

  assoc <- do.call(rbind, rows)
  assoc$padj <- p.adjust(assoc$p, method = "BH")
  assoc
}

# Cramer's V of two categorical columns, and whether one is nested in the
# other (each level of one column occurs with a single level of the other)
compare_columns <- function(meta, a, b) {
  tab <- table(meta[[a]], meta[[b]], dnn = c(a, b))
  tab <- tab[rowSums(tab) > 0, colSums(tab) > 0, drop = FALSE]
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NULL)

  chi2 <- suppressWarnings(chisq.test(tab, correct = FALSE)$statistic)
  v <- sqrt(unname(chi2) / (sum(tab) * (min(dim(tab)) - 1)))

  a_in_b <- all(rowSums(tab > 0) == 1)
  b_in_a <- all(colSums(tab > 0) == 1)
  nested <- if (a_in_b && b_in_a) {
    "identical grouping"
  } else if (a_in_b) {
    paste(a, "within", b)
  } else if (b_in_a) {
    paste(b, "within", a)
  } else {
    ""
  }

  list(
    pair = data.frame(
      column_1 = a,
      column_2 = b,
      n = sum(tab),
      cramers_v = round(v, 3),
      nested = nested
    ),
    table = tab
  )
}

# TRUE where x is more than k MADs, and at least min_diff, from the median,
# on one side only
is_outlier <- function(x, k, side = c("low", "high"), min_diff = 0) {
  side <- match.arg(side)
  centre <- median(x, na.rm = TRUE)
  limit <- max(k * mad(x, na.rm = TRUE), min_diff, na.rm = TRUE)
  if (limit == 0) return(rep(FALSE, length(x)))

  dev <- if (side == "low") centre - x else x - centre
  out <- dev > limit
  !is.na(out) & out
}

# "low library size, low correlation" per sample; "" when nothing is flagged
flag_samples <- function(samples, k, min_diff) {
  test <- function(metric, side, x = samples[[metric]]) {
    is_outlier(x, k, side, min_diff[[metric]])
  }
  log_lib <- log10(samples$library_size)

  checks <- list(
    "low library size" = test("library_size", "low", log_lib),
    "low mapping rate" = test("pct_mapped", "low"),
    "few detected genes" = test("detected_genes", "low"),
    "high top-10 share" = test("top10_pct", "high"),
    "low correlation" = test("median_cor", "low")
  )

  flagged <- do.call(cbind, checks)
  apply(flagged, 1, function(row) toString(colnames(flagged)[row]))
}

# Number of samples flagged for one reason
flag_count <- function(reason) sum(grepl(reason, samples$flags, fixed = TRUE))

# Min, median and max of a metric, as text
metric_range <- function(x, digits = 0) {
  if (all(is.na(x))) return(rep("NA", 3))
  formatC(
    c(min(x, na.rm = TRUE), median(x, na.rm = TRUE), max(x, na.rm = TRUE)),
    format = "f",
    digits = digits,
    big.mark = ","
  )
}


#---- Import datasets ---------------------------------------------------------

meta <- read_metadata(opt$metadata, opt$sample_col)
tx2gene <- read_tx2gene(opt$tx2gene)
files <- find_quant_files(opt$quant_dir)


#---- Code --------------------------------------------------------------------

# Samples: only those with both metadata and a quant file are checked
meta_only <- setdiff(meta$sample, names(files))
quant_only <- setdiff(names(files), meta$sample)
files <- files[names(files) %in% meta$sample]

if (length(files) < 3) {
  stop(
    "Sample QC needs at least 3 samples with metadata and a quant.sf, found ",
    length(files), ". Check --sample_col: metadata has ",
    toString(head(meta$sample, 3)), "; quant folders are ",
    toString(head(c(names(files), quant_only), 3))
  )
}

cols <- choose_columns(meta, opt$columns, opt$sample_col)
col_types <- vapply(cols, function(col) column_type(meta[[col]]), "")
cat_cols <- cols[col_types == "categorical"]
meta[cat_cols] <- lapply(meta[cat_cols], as_ordered_factor)

print_section("Samples")
cat(sprintf(" > %-22s : %d\n", "With metadata and quant", length(files)))
cat(sprintf(" > %-22s : %d\n", "Metadata, no quant", length(meta_only)))
cat(sprintf(" > %-22s : %d\n", "Quant, no metadata", length(quant_only)))
if (length(meta_only) > 0) cat("   -", toString(meta_only), "\n")
if (length(quant_only) > 0) cat("   -", toString(quant_only), "\n")

# Counts
dds <- build_dds(meta, tx2gene, files)
meta <- as.data.frame(colData(dds))
n_quantified <- nrow(dds)

dds <- dds[rowSums(counts(dds)) > 0, ]
cnt <- counts(dds)
n_samples <- ncol(cnt)

gene_names <- setNames(tx2gene$gene_name, tx2gene$gene_id)
gene_names <- gene_names[!duplicated(names(gene_names))]

size_factor <- tryCatch(
  estimateSizeFactorsForMatrix(cnt),
  error = function(e) rep(NA_real_, n_samples)
)

# PCA and correlations on the variance-stabilised values
cat(" > Variance-stabilising the counts (vst, blind)...\n")
mat <- run_vst(dds)
pca <- run_pca(mat, opt$ntop)
n_pcs <- min(n_pcs, ncol(pca$x))
pcs <- pca$x[colnames(cnt), seq_len(n_pcs), drop = FALSE]

expressed <- rowSums(cnt >= opt$min_count) >= ceiling(n_samples / 2)
if (sum(expressed) < 10) {
  stop(
    "Fewer than 10 genes have ", opt$min_count,
    " reads in half of the samples; lower --min_count"
  )
}

cor_mat <- cor(mat[expressed, , drop = FALSE], method = "spearman")
median_cor <- vapply(seq_len(n_samples), function(i) {
  median(cor_mat[i, -i])
}, numeric(1))
clustering <- hclust(as.dist(1 - cor_mat), method = "average")

# One row per sample
salmon <- salmon_stats(files)
salmon <- salmon[match(colnames(cnt), salmon$sample), ]

samples <- data.frame(
  sample = colnames(cnt),
  meta[colnames(cnt), cols, drop = FALSE],
  n_processed = salmon$n_processed,
  pct_mapped = round(salmon$pct_mapped, 2),
  library_size = round(colSums(cnt)),
  detected_genes = colSums(cnt > 0),
  top10_pct = round(top_share(cnt, 10), 2),
  size_factor = round(size_factor, 3),
  median_cor = round(median_cor, 4),
  round(pcs, 3),
  row.names = NULL,
  check.names = FALSE
)
samples$flags <- flag_samples(samples, opt$mad, min_diff)

# Association of the PCs with the metadata and the technical metrics
technical <- list(
  "Library size (log10)" = log10(samples$library_size),
  "Mapping rate" = samples$pct_mapped,
  "Detected genes" = samples$detected_genes,
  "Top-10 share" = samples$top10_pct
)
variables <- c(as.list(samples[cols]), technical)
assoc <- pc_association(pcs, variables)
assoc$type <- ifelse(assoc$variable %in% cols, "metadata", "technical")

# Confounded metadata columns
pairs <- if (length(cat_cols) >= 2) {
  combn(cat_cols, 2, simplify = FALSE) |>
    lapply(function(p) compare_columns(samples, p[1], p[2]))
} else {
  list()
}
pairs <- pairs[!vapply(pairs, is.null, logical(1))]
names(pairs) <- vapply(pairs, function(p) {
  paste(p$pair$column_1, p$pair$column_2, sep = " x ")
}, character(1))

confounding <- do.call(rbind, lapply(pairs, `[[`, "pair"))
if (is.null(confounding)) {
  confounding <- data.frame(
    column_1 = character(0),
    column_2 = character(0),
    n = integer(0),
    cramers_v = numeric(0),
    nested = character(0)
  )
}
rownames(confounding) <- NULL

columns <- data.frame(
  column = cols,
  type = unname(col_types),
  levels = vapply(cols, function(col) {
    length(unique(samples[[col]][!is.na(samples[[col]])]))
  }, integer(1)),
  missing = vapply(cols, function(col) sum(is.na(samples[[col]])), integer(1)),
  row.names = NULL
)


#---- Summary -----------------------------------------------------------------

print_section("Metadata columns")
print(columns, row.names = FALSE)

print_section("Genes")
cat(sprintf(" > %-22s : %d\n", "Quantified", n_quantified))
cat(sprintf(" > %-22s : %d\n", "With reads", nrow(cnt)))
cat(sprintf(
  " > %-22s : %d (>= %d reads in half of the samples)\n",
  "Expressed",
  sum(expressed),
  opt$min_count
))

summary_tbl <- data.frame(
  metric = c(
    "Library size (reads)",
    "Mapping rate (%)",
    "Detected genes",
    "Top-10 share (%)",
    "Median correlation"
  ),
  rbind(
    metric_range(samples$library_size),
    metric_range(samples$pct_mapped, 1),
    metric_range(samples$detected_genes),
    metric_range(samples$top10_pct, 1),
    metric_range(samples$median_cor, 3)
  ),
  flagged = c(
    flag_count("low library size"),
    flag_count("low mapping rate"),
    flag_count("few detected genes"),
    flag_count("high top-10 share"),
    flag_count("low correlation")
  )
)
colnames(summary_tbl) <- c("metric", "min", "median", "max", "flagged")

print_section("Per-sample metrics")
print(summary_tbl, row.names = FALSE)

print_section(sprintf("PCA (%d most variable genes)", pca$n_genes))
cat(sprintf(
  " > %-5s : %5.1f %% of the variance\n",
  colnames(pcs),
  pca$var_pct[seq_len(n_pcs)]
), sep = "")

# The variable that explains most of each PC
print_section("Strongest association per PC (R2)")
for (pc in colnames(pcs)) {
  a <- assoc[assoc$pc == pc & !is.na(assoc$r2), ]
  if (nrow(a) == 0) next
  best <- a[which.max(a$r2), ]
  cat(sprintf(
    " > %-5s : %-22s R2 = %.2f, padj = %.2g (%s)\n",
    pc,
    best$variable,
    best$r2,
    best$padj,
    best$type
  ))
}

print_section("Metadata columns compared (Cramer's V)")
if (nrow(confounding) == 0) {
  cat(" > Fewer than two categorical columns\n")
} else {
  print(confounding, row.names = FALSE)
}

flagged <- samples[samples$flags != "", c("sample", "flags")]

print_section(sprintf(
  "Flagged samples (> %s MADs and min_diff from the median)",
  opt$mad
))
if (nrow(flagged) == 0) {
  cat(" > None\n")
} else {
  cat(sprintf(" > %-20s : %s\n", flagged$sample, flagged$flags), sep = "")
}


#---- Save files --------------------------------------------------------------

saveRDS(
  list(
    settings = list(
      columns = cols,
      sample_col = opt$sample_col,
      ntop = pca$n_genes,
      min_diff = min_diff,
      min_count = opt$min_count,
      mad = opt$mad,
      versions = c(
        DESeq2 = as.character(packageVersion("DESeq2")),
        tximport = as.character(packageVersion("tximport"))
      )
    ),
    samples = samples,
    columns = columns,
    checks = list(meta_only = meta_only, quant_only = quant_only),
    genes = c(
      quantified = n_quantified,
      with_reads = nrow(cnt),
      expressed = sum(expressed)
    ),
    summary = summary_tbl,
    density = density_table(cnt),
    top_genes = top_genes_table(cnt, gene_names, n_top_genes),
    pca = list(
      var_pct = pca$var_pct,
      loadings = loadings_table(pca, gene_names, min(2, n_pcs), n_top_genes)
    ),
    pc_assoc = assoc,
    cor = cor_mat,
    clustering = clustering,
    confounding = confounding,
    crosstabs = lapply(pairs, `[[`, "table")
  ),
  file = out_rds
)

write.table(
  samples,
  file = out_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


#---- Done --------------------------------------------------------------------

cat(sprintf(
  "\n > Script complete!\n%s\n > Outputs:\n%s",
  strrep("-", 33),
  paste0("   - ", c(out_rds, out_table), "\n", collapse = "")
))