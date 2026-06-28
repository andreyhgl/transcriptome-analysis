#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# get number of cores for process
# import genexp table
# extract the unique genes
# plot the gene counts per gene, include all generations in the same plot
# show significance indicators within the plots
# use log10 scale on y-axis, ylimits depend on number of sign. indicators
# add gene name and gene function to title and caption, respectively
#------------------------------------------------------------------------------


#---- Parse arguments ---------------------------------------------------------

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(" > Script needs atleast one argument! Usage: DESeq2_gene_plots.R --cores")
}

for (i in seq_along(args)) {
  if (args[i] == "--cores") ncores <- args[i + 1]
}

cat("\n > number of cores", ncores, "\n")

suppressPackageStartupMessages({
  library(DESeq2)
  
  library(ggplot2)
  library(ggsignif)
  library(ggrepel)
  library(cowplot)
  library(patchwork)
  library(RColorBrewer)
  
  library(scales)
  library(openxlsx)
  library(stringr)

  if (ncores > 1) library(parallel)
})


#---- Functions ---------------------------------------------------------------


make_gene_plot <- function(gene, dds, plot_size = 11) {
  #gene <- "ENSMUSG00000104328"

  # extract normalised and transformed counts
  data <- plotCounts(
    dds,
    gene = gene,
    intgroup = c("treatment", "protocol", "extraction_batch"),
    normalized = TRUE,
    transform = FALSE,
    returnData = TRUE
  )
  data$count <- log2(data$count + 1)
  data$sample <- sub("sperm_", "", rownames(data))
  rownames(data) <- NULL

  # extract gene and statistics info
  info <- genexp[genexp$ensembl_gene_id %in% gene, ]

  # calculate significance indicator values
  sign_info <- info[info$sign == "sign", ]
  
  # Parse the contrast into "treatment_100_vs_0" -> c("100", "0")
  comparisons_list <- lapply(sign_info$contrast, function(x) {
    parts <- sub("treatment_", "", x)
    groups <- strsplit(parts, "_vs_")[[1]]
    groups
  })

  # Map padj to significance stars (optional)
  annotations <- sapply(sign_info$padj, function(p) {
    if (p < 0.001) "***"
    else if (p < 0.01) "**"
    else if (p < 0.05) "*"
    else "ns"
  })


  subtitle <- paste0(
    info$chr[1], ":", 
    comma(info$start[1]), "-", comma(info$end[1])
  )
  caption <- paste(info$gene_type[1], "|", str_to_sentence(info$gene_info[1]))

  # set colors
  colors <- brewer.pal(8, "Set1")[c(3, 2, 4)]
  names(colors) <- c("0", "10", "100")

  pos <- position_jitter(seed = 1337, width = 0.2)

  # render plot
  ggplot(data, aes(treatment, count, fill = treatment, label = sample)) +
    stat_boxplot(geom = "errorbar", width = 0.2) +
    geom_boxplot(outlier.shape = NA, show.legend = FALSE, width = 0.5) +
    geom_label_repel(
      size = 2.5,
      show.legend = FALSE,
      #box.padding = 2,
      max.overlaps = 1,
      min.segment.length = 0.1,
      position = pos) +
    geom_jitter(shape = 21, size = 2.5, position = pos) +
    labs(
      y = expression(log[2]("norm. counts")),
      x = "DBP exposure [mg/kg]",
      title = info$gene_name,
      subtitle = subtitle,
      caption = caption) +
    theme_linedraw(plot_size) +
    theme(
      plot.title = element_text(hjust = 0.5),
      plot.subtitle = element_text(hjust = 0.5),
      #legend.position = "bottom",
      legend.position = "none") +
    background_grid() +
    geom_signif(
      comparison = comparisons_list,
      annotations = annotations,
      step_increase = 0.1, margin_top = 0.1) + 
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
    scale_fill_manual(values = colors)
}

save_gene_plot <- function(gene, width = 4, height = 5) {
  n <- grep(gene, genes)

  cat(" > Processing:", n, "\t", gene, "\n")

  p <- make_gene_plot(gene, dds = dds, plot_size = 11)

  filename <- paste0(n, ".pdf")

  ggsave(
    filename,
    plot = p,
    width = width,
    height = height
  )

  return(filename)
}


#---- Import datasets ---------------------------------------------------------


dds <- readRDS("DDS.Rds")
genexp <- read.csv("genexp_table.csv.gz")


#---- Code --------------------------------------------------------------------


genes <- unique(genexp[genexp$sign == "sign", "ensembl_gene_id"])
cat(" > number of unique genes:", comma(length(genes)), "\n")

# Parallel execution using mclapply (Linux/macOS only)
mclapply_output <- mclapply(
  genes,
  save_gene_plot,
  mc.cores = ncores
) 

# testing:
#writeLines(mclapply_output, "mclapply_output.txt")
print(mclapply_output)

#---- Done --------------------------------------------------------------------


cat(glue::glue("
 ~~ DESeq2_gene_plots.R complete ~~~~~~~~~~~~~~~~~~~~~~~~
  > Number plots rendered   : {length(mclapply_output)}
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"))