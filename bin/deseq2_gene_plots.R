#!/usr/bin/env Rscript

#==============================================================================
# deseq2_gene_plots.R
#
# Plot the normalised counts of the top differentially expressed genes, one
# PDF per gene (the Nextflow process merges them into one). All groups of the
# design's variable of interest share one plot, with significance brackets
# for every contrast in which the gene is significant. Gene name, location
# and description from an optional annotation table go in the title,
# subtitle and caption.
#
# Input:
#   - Rds with a DESeqDataSet that has been run through DESeq()
#   - results tables from deseq2_results.R in --results_dir
#     (<prefix>_<level>_vs_<reference>.tsv, one per contrast)
#   - optional annotation csv/tsv with a `gene_id` column, and any of
#     gene_name, chr, start, end, gene_type and description
#
# Method:
#   Genes: the --top_n genes with the lowest padj (< --alpha) per contrast,
#   or every significant gene with --top_n 0. Counts are DESeq2-normalised
#   and shown as log2(count + 1). Stars: padj < 0.05 *, < 0.01 **,
#   < 0.001 ***. Brackets sit in a panel above the counts, so they don't
#   change the counts' y axis.
#
# Usage:
# > Rscript deseq2_gene_plots.R --dds dds.Rds --results_dir results
# > Rscript deseq2_gene_plots.R --annotation genes.tsv --top_n 0 --cores 4
#
# Container:
# > docker://ghcr.io/andreyhgl/deseq2-env:main
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
  prefix = "dds",
  dds = "dds.Rds",
  results_dir = "results",
  annotation = "../annotation/ensembl_human_v115.tsv",
  sample_prefix = "sperm_",
  alpha = 0.05,
  top_n = 20L,
  include_label = FALSE,
  cores = 1L
)

opt <- parse_args(defaults)

cat("\n > Arguments:\n---------------------------------\n")
cat(sprintf(" > %-14s : %s\n", names(opt), unlist(opt)), sep = "")
cat(strrep("-", 33), "\n")


#---- Environment -------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(parallel)
})


#---- Settings ----------------------------------------------------------------

# Same colours as the PCA in quant_qc.R: grey for the reference level, then
# the categorical palette in a fixed order
series_colours <- c(
  "#2a78d6", "#eb6834", "#1baf7a", "#eda100",
  "#e87ba4", "#008300", "#6250d6", "#e34948"
)
ref_colour <- "#898781"

# Plot size (inches) and base font size
plot_width <- 5
plot_height <- 8
bracket_height <- 0.28
base_size <- 14

# Output files
out_dir <- paste0(opt$prefix, "_gene_plots")
out_table <- paste0(opt$prefix, "_gene_plots.tsv")


#---- Functions ---------------------------------------------------------------

print_section <- function(title) {
  cat("\n > ", title, "\n", strrep("-", 33), "\n", sep = "")
}

# Results tables for every level vs the reference, stacked with a `group`
# column naming the level
read_results <- function(dir, lvls) {
  ref <- lvls[1]
  files <- list.files(dir, pattern = "_vs_.*\\.tsv$", full.names = TRUE)

  tables <- lapply(lvls[-1], function(lvl) {
    file <- files[endsWith(files, paste0("_", lvl, "_vs_", ref, ".tsv"))]
    if (length(file) != 1) {
      stop(
        "Expected 1 results table for ", lvl, " vs ", ref, ", found ",
        length(file), " in: ", dir
      )
    }

    tbl <- read.delim(file)
    data.frame(gene_id = tbl$gene_id, group = lvl, padj = tbl$padj)
  })

  do.call(rbind, tables)
}

# Optional annotation, one row per gene_id; NULL when no file is given
read_annotation <- function(path) {
  if (path == "") return(NULL)

  sep <- if (grepl("\\.csv(\\.gz)?$", path)) "," else "\t"
  anno <- read.delim(path, sep = sep)

  if (!"gene_id" %in% colnames(anno)) {
    stop("Annotation needs a `gene_id` column: ", path)
  }

  anno[!duplicated(anno$gene_id), ]
}

# Top genes per contrast, one row per gene, ordered by best padj
select_genes <- function(sig, top_n) {
  if (nrow(sig) == 0) {
    return(data.frame(
      rank = integer(0),
      gene_id = character(0),
      best_padj = numeric(0),
      significant_in = character(0)
    ))
  }

  sig <- sig[order(sig$padj), ]
  top <- if (top_n > 0) {
    do.call(rbind, lapply(split(sig, sig$group), head, top_n))
  } else {
    sig
  }
  ids <- unique(top[order(top$padj), "gene_id"])

  # Report every significant contrast for the gene, not only the top ones
  best_padj <- tapply(sig$padj, sig$gene_id, min)
  contrasts <- tapply(sig$group, sig$gene_id, paste, collapse = ",")

  data.frame(
    rank = seq_along(ids),
    gene_id = ids,
    best_padj = signif(as.numeric(best_padj[ids]), 3),
    significant_in = as.character(contrasts[ids])
  )
}

# Annotation value for one gene, or NA when missing
anno_field <- function(anno, gene, col) {
  if (is.null(anno) || !col %in% colnames(anno)) return(NA_character_)

  val <- anno[[col]][match(gene, anno$gene_id)]
  if (is.na(val) || val == "") NA_character_ else as.character(val)
}

# Title, subtitle and caption for one gene from the annotation
gene_labels <- function(gene, anno) {
  field <- function(col) anno_field(anno, gene, col)

  name <- field("gene_name")
  title <- if (is.na(name)) gene else name

  loc <- c(field("chr"), field("start"), field("end"))
  subtitle <- if (anyNA(loc)) {
    if (is.na(name)) NULL else gene
  } else {
    sprintf(
      "%s:%s-%s",
      loc[1],
      scales::comma(as.numeric(loc[2])),
      scales::comma(as.numeric(loc[3]))
    )
  }

  # Drop Ensembl's "[Source:...]" suffix and start with a capital letter
  desc <- sub("\\s*\\[Source:.*\\]$", "", field("description"))
  if (!is.na(desc)) {
    desc <- paste0(toupper(substr(desc, 1, 1)), substr(desc, 2, nchar(desc)))
  }
  caption <- paste(na.omit(c(field("gene_type"), desc)), collapse = " | ")

  list(
    title = title,
    subtitle = subtitle,
    caption = if (caption == "") NULL else caption
  )
}

# Grey for the reference level, then the palette in fixed order
level_colours <- function(lvls) {
  other <- lvls[-1]
  if (length(other) > length(series_colours)) {
    stop("Too many groups to colour (max ", length(series_colours), ")")
  }

  setNames(c(ref_colour, series_colours[seq_along(other)]), lvls)
}

# Significance stars from padj; padj above 0.05 (with --alpha > 0.05) is
# printed as its value
padj_stars <- function(padj) {
  stars <- cut(
    padj,
    breaks = c(0, 0.001, 0.01, 0.05),
    labels = c("***", "**", "*"),
    right = FALSE
  )
  ifelse(is.na(stars), sprintf("%.2g", padj), as.character(stars))
}

# Significance brackets in their own panel, drawn above the counts so they
# don't stretch its y axis. Group i sits at x = i on the counts' discrete
# axis, so the same x range here (0.4 to n + 0.6, ggplot's default padding
# for discrete axes) lines the brackets up with the groups.
bracket_panel <- function(gene_sig, lvls) {
  n <- nrow(gene_sig)
  tip <- 0.3

  df <- data.frame(
    x = 1,
    xend = match(gene_sig$group, lvls),
    y = seq_len(n),
    label = padj_stars(gene_sig$padj)
  )

  ggplot(df) +
    geom_segment(aes(x = x, xend = xend, y = y, yend = y)) +
    geom_segment(aes(x = x, xend = x, y = y, yend = y - tip)) +
    geom_segment(aes(x = xend, xend = xend, y = y, yend = y - tip)) +
    geom_text(
      aes(x = (x + xend) / 2, y = y, label = label),
      vjust = -0.2,
      size = base_size / .pt
    ) +
    coord_cartesian(
      xlim = c(0.4, length(lvls) + 0.6),
      ylim = c(1 - tip, n + 0.8),
      expand = FALSE,
      clip = "off"
    ) +
    theme_void()
}

make_gene_plot <- function(gene, dds, group, sig, anno) {
  lvls <- levels(dds[[group]])

  data <- data.frame(
    group = dds[[group]],
    count = log2(counts(dds, normalized = TRUE)[gene, ] + 1),
    sample = sub(paste0("^", opt$sample_prefix), "", colnames(dds))
  )

  # One bracket per significant contrast, innermost (closest group) lowest
  gene_sig <- sig[sig$gene_id == gene, ]
  gene_sig <- gene_sig[order(match(gene_sig$group, lvls)), ]

  labels <- gene_labels(gene, anno)
  pos <- position_jitter(seed = 1337, width = 0.2)

  p <- ggplot(data, aes(group, count, fill = group, label = sample)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.shape = NA, width = 0.7) +
    geom_point(shape = 21, size = 2.5, position = pos) +
    scale_fill_manual(values = level_colours(lvls)) +
    scale_x_discrete(labels = function(x) gsub("_", "\n", x)) +
    scale_y_continuous(expand = expansion(mult = 0.05)) +
    labs(x = NULL, y = expression(log[2]("norm. counts + 1"))) +
    theme_classic(base_size) +
    theme(legend.position = "none", panel.grid.minor = element_blank())

  if (opt$include_label) {
    p <- p +
      geom_label_repel(
        size = 2.5,
        max.overlaps = 1,
        min.segment.length = 0.1,
        position = pos
      )
  }

  # Brackets panel on top: a fixed height per bracket, the counts get the rest
  fig <- if (nrow(gene_sig) > 0) {
    bracket_panel(gene_sig, lvls) / p +
      plot_layout(
        heights = unit(c(bracket_height * nrow(gene_sig), 1), c("in", "null"))
      )
  } else {
    wrap_plots(p)
  }

  fig +
    plot_annotation(
      title = labels$title,
      subtitle = labels$subtitle,
      caption = labels$caption,
      theme = theme(
        plot.title = element_text(
          size = base_size * 1.2,
          hjust = 0.5,
          face = "bold"
        ),
        plot.subtitle = element_text(size = base_size, hjust = 0.5),
        plot.caption = element_text(size = base_size * 0.8)
      )
    )
}

# Render and save the plot for row i of the gene table
save_gene_plot <- function(i, genes, dds, group, sig, anno) {
  gene <- genes$gene_id[i]
  p <- make_gene_plot(gene, dds, group, sig, anno)

  file <- file.path(out_dir, genes$file[i])
  ggsave(file, plot = p, width = plot_width, height = plot_height)
  file
}


#---- Import datasets ---------------------------------------------------------

dds <- readRDS(opt$dds)

# Variable of interest: the last column of the design
group <- all.vars(design(dds)) |> tail(1)

res <- read_results(opt$results_dir, levels(dds[[group]]))
anno <- read_annotation(opt$annotation)


#---- Code --------------------------------------------------------------------

sig <- res[!is.na(res$padj) & res$padj < opt$alpha, ]
genes <- select_genes(sig, opt$top_n)

if (!is.null(anno)) {
  genes$gene_name <- vapply(
    genes$gene_id,
    function(g) anno_field(anno, g, "gene_name"),
    character(1)
  )
}

# File names: rank and gene name (or id), with unsafe characters replaced
name <- if (is.null(genes$gene_name)) genes$gene_id else genes$gene_name
name <- ifelse(is.na(name), genes$gene_id, name)
genes$file <- sprintf(
  "%03d_%s.pdf",
  genes$rank,
  gsub("[^A-Za-z0-9._-]", "_", name)
)

print_section("Genes")
cat(sprintf(" > %-14s : %s\n", "Group", group))
cat(sprintf(" > %-14s : %s\n", "Significant", nrow(sig)))
cat(sprintf(" > %-14s : %d\n", "Genes to plot", nrow(genes)))

dir.create(out_dir, showWarnings = FALSE)

# No significant genes: one placeholder page, so the merged PDF still exists
if (nrow(genes) == 0) {
  msg <- sprintf("No genes with padj < %s in any contrast", opt$alpha)
  ggsave(
    file.path(out_dir, "000_no_genes.pdf"),
    plot = ggplot() + annotate("text", x = 0, y = 0, label = msg) +
      theme_void(),
    width = plot_width,
    height = plot_height
  )
}

# mclapply() returns errors as "try-error" objects instead of stopping
out <- mclapply(
  seq_len(nrow(genes)),
  save_gene_plot,
  genes = genes,
  dds = dds,
  group = group,
  sig = sig,
  anno = anno,
  mc.cores = opt$cores
)

failed <- vapply(out, inherits, logical(1), what = "try-error")
if (any(failed)) {
  stop(
    "Plotting failed for: ", toString(genes$gene_id[failed]), "\n",
    out[failed][[1]]
  )
}

write.table(
  genes,
  file = out_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)


#---- Done --------------------------------------------------------------------

cat(sprintf(
  "\n > Script complete!\n%s\n > Outputs:\n   - %s/ (%d plots)\n   - %s\n",
  strrep("-", 33),
  out_dir,
  length(out),
  out_table
))