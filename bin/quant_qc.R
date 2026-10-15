#!/usr/bin/env Rscript

#==============================================================================
# quant_qc.R
#
# Quality check of the salmon quantification files before differential
# expression. Gene counts are imported with tximport, variance-stabilised
# blind to the design, and summarised with a PCA of the most variable genes.
# One PCA panel is drawn per metadata column in --fill, combined into one
# figure with patchwork. A per-sample table reports library size, detected
# genes and the first principal components.
#
# Input:
#   - metadata csv with a `sample` column and every column in --fill
#   - tx2gene tsv from salmon (transcript id, gene id in the first 2 columns)
#   - quant.sf files under --quant_dir, one folder per sample
#     (the folder name is the sample id)
#
# Method:
#   vst(blind = TRUE) on genes with at least one read, then prcomp() on the
#   --ntop most variable genes, as in DESeq2::plotPCA().
#   Reference levels ("0", "control" or NA) are drawn in grey, the other
#   levels take the categorical palette in a fixed order. A column with more
#   than three other levels also gets one point shape per level, because
#   colour alone can't separate that many groups reliably in a scatter plot.
#
# Usage:
# > Rscript quant_qc.R --metadata metadata.csv --fill id,season,year
# > Rscript quant_qc.R --fill id --ntop 1000 --filename quant_qc
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
  filename = "quant_qc",
  metadata = "data/metadata.csv",
  tx2gene = "results/star_salmon/salmon.merged.tx2gene.tsv",
  quant_dir = "results/star_salmon/",
  fill = "id,season,year",
  ntop = 500L
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
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
})


#---- Settings ----------------------------------------------------------------

# Metadata columns to colour the PCA by, one panel each: "id,season,year"
fill_cols <- strsplit(opt$fill, "[,[:space:]]+")[[1]]
fill_cols <- fill_cols[fill_cols != ""]
if (length(fill_cols) == 0) stop("--fill needs at least one column name")

# Categorical palette, assigned in this fixed order and never recycled.
# Reference levels (control, 0, NA) are drawn in a neutral grey.
series_colours <- c(
  "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
  "#e87ba4", "#008300", "#6250d6", "#e34948"
)
ref_colour <- "#898781"

# Fillable point shapes: circle, triangle, square, diamond, down-triangle
point_shapes <- c(21, 24, 22, 23, 25)
max_colour_only <- 3

# Chart surface and ink
surface <- "#fcfcfb"
ink_primary <- "#0b0b0b"
ink_secondary <- "#52514e"
grid_colour <- "#e1e0d9"

# Figure size per panel (inches)
panel_width <- 4
fig_height <- 5

# Output files
out_fig <- paste0(opt$filename, "_pca.pdf")
out_table <- paste0(opt$filename, "_samples.tsv")


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

# Check required columns, harmonise sample ids and factor the fill columns
prepare_metadata <- function(meta, cols) {
  required <- unique(c("sample", cols))

  missing <- setdiff(required, colnames(meta))
  if (length(missing) > 0) {
    stop("Column(s) missing from metadata: ", paste(missing, collapse = ", "))
  }

  as_ordered_factor <- function(x) {
    lvls <- unique(as.character(x))
    ref <- lvls[tolower(lvls) %in% c("0", "control")]
    factor(x, levels = c(ref, setdiff(lvls, ref)))
  }

  meta[cols] <- lapply(meta[cols], as_ordered_factor)
  meta
}

# Import quant files for the metadata samples into a DESeq2 dataset. The
# design is ~ 1: QC must not depend on the experimental design.
build_dds <- function(meta, tx2gene, files) {
  files <- files[names(files) %in% meta$sample]
  if (length(files) == 0) stop("No quant files match the metadata samples")

  txi <- tximport(files, type = "salmon", tx2gene = tx2gene)

  meta <- meta[match(colnames(txi$counts), meta$sample), ] |> droplevels()
  rownames(meta) <- meta$sample

  DESeqDataSetFromTximport(txi, colData = meta, design = ~ 1)
}

# PCA of the vst values of the ntop most variable genes (as plotPCA)
run_pca <- function(dds, ntop) {
  dds <- dds[rowSums(counts(dds)) > 0, ]

  # vst() estimates the dispersion trend on a 1000-gene subset
  vsd <- if (nrow(dds) < 1000) {
    varianceStabilizingTransformation(dds, blind = TRUE)
  } else {
    vst(dds, blind = TRUE)
  }

  mat <- assay(vsd)
  rows <- order(MatrixGenerics::rowVars(mat), decreasing = TRUE)
  rows <- head(rows, min(ntop, nrow(mat)))
  pca <- prcomp(t(mat[rows, ]))

  list(
    x = pca$x,
    var_pct = 100 * pca$sdev^2 / sum(pca$sdev^2),
    n_genes = length(rows)
  )
}

# One row per sample: metadata, library size, detected genes and PCs
sample_table <- function(dds, pca, cols, n_pcs = 5) {
  cnts <- counts(dds)
  n_pcs <- min(n_pcs, ncol(pca$x))
  pcs <- round(pca$x[colnames(dds), seq_len(n_pcs), drop = FALSE], 3)

  data.frame(
    sample = colnames(dds),
    as.data.frame(colData(dds))[cols],
    library_size = colSums(cnts),
    detected_genes = colSums(cnts > 0),
    pcs,
    row.names = NULL,
    check.names = FALSE
  )
}

# Grey for reference levels (control / 0), then the palette in fixed order
level_colours <- function(x) {
  lvls <- levels(x)
  ref <- lvls[tolower(lvls) %in% c("0", "control")]
  other <- setdiff(lvls, ref)

  if (length(other) > length(series_colours)) {
    stop(
      "Too many levels to colour (max ", length(series_colours), "): ",
      toString(other)
    )
  }

  c(
    setNames(rep(ref_colour, length(ref)), ref),
    setNames(series_colours[seq_along(other)], other)
  )
}

theme_qc <- function() {
  theme_linedraw(base_size = 15) +
    theme(
      plot.background = element_rect(fill = surface, colour = NA),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(colour = grid_colour, linewidth = 0.3),
      plot.title = element_text(face = "bold", hjust = 0.5),
      legend.position = "bottom",
      legend.title = element_blank(),
    )
}

# PC1 vs PC2 with the points filled by one metadata column
plot_pca <- function(df, col, var_pct) {
  x <- df[[col]]
  colours <- level_colours(x)
  n_other <- sum(colours != ref_colour)
  legend <- guide_legend(nrow = 2)

  # Many groups: add a shape per level so they don't rely on colour alone
  use_shapes <- n_other > max_colour_only &&
    nlevels(x) <= length(point_shapes)

  p <- ggplot(df, aes(PC1, PC2, fill = .data[[col]])) +
    geom_text_repel(aes(label = sample))

  if (use_shapes) {
    shapes <- setNames(point_shapes[seq_len(nlevels(x))], levels(x))
    p <- p +
      geom_point(
        aes(shape = .data[[col]]),
        size = 3,
        colour = "black"
      ) +
      scale_shape_manual(values = shapes, na.value = 21, guide = legend)
  } else {
    p <- p + geom_point(shape = 21, size = 3, stroke = 0.6)
  }

  p +
    scale_fill_manual(values = colours, na.value = ref_colour, guide = legend) +
    labs(
      title = col,
      x = sprintf("PC1 (%.1f%%)", var_pct[1]),
      y = sprintf("PC2 (%.1f%%)", var_pct[2])
    ) +
    theme_qc()
}

report_samples <- function(qc) {
  lib <- qc$library_size

  print_section("Samples")
  cat(sprintf(" > %-15s : %d\n", "Samples", nrow(qc)))
  cat(sprintf(
    " > %-15s : %s / %s / %s (min / median / max)\n",
    "Library size",
    format(min(lib), big.mark = ","),
    format(median(lib), big.mark = ","),
    format(max(lib), big.mark = ",")
  ))
  cat(sprintf(
    " > %-15s : %s - %s\n",
    "Detected genes",
    format(min(qc$detected_genes), big.mark = ","),
    format(max(qc$detected_genes), big.mark = ",")
  ))
}


#---- Import datasets ---------------------------------------------------------

meta <- read.csv(opt$metadata) |>
  prepare_metadata(fill_cols)

tx2gene <- read.table(opt$tx2gene, header = TRUE)[, 1:2] |>
  setNames(c("TXNAME", "GENEID"))

quant_files <- find_quant_files(opt$quant_dir)


#---- Code --------------------------------------------------------------------

dds <- build_dds(meta, tx2gene, quant_files)
pca <- run_pca(dds, opt$ntop)
qc <- sample_table(dds, pca, fill_cols)

report_samples(qc)

print_section("PCA")
cat(sprintf(" > %-15s : %d\n", "Genes used", pca$n_genes))
cat(sprintf(
  " > PC%d variance    : %.1f%%\n",
  seq_len(min(5, length(pca$var_pct))),
  pca$var_pct[seq_len(min(5, length(pca$var_pct)))]
), sep = "")

plots <- lapply(fill_cols, function(col) plot_pca(qc, col, pca$var_pct))

fig <- wrap_plots(plots, nrow = 1) +
  plot_layout(axes = "collect") +
  plot_annotation(
    title = "PCA",
    subtitle = sprintf(
      "vst counts (blind to design), top %d variable genes, %d samples",
      pca$n_genes,
      nrow(qc)
    ),
    tag_levels = "A",
    theme = theme(
      plot.background = element_rect(fill = surface, colour = NA),
      plot.title = element_text(face = "bold")
    )
  )

ggsave(
  out_fig,
  fig,
  width = panel_width * length(plots),
  height = fig_height,
  bg = surface
)

write.table(
  qc,
  file = out_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


#---- Done --------------------------------------------------------------------

cat(sprintf(
  "\n > Script complete!\n%s\n > Outputs:\n   - %s\n   - %s\n",
  strrep("-", 33),
  out_fig,
  out_table
))