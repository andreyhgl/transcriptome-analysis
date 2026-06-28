#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# print quality control plots
# PCA
# cooks distance
# euclidean distance
#------------------------------------------------------------------------------


suppressPackageStartupMessages({
  #library(forcats)
  #library(scales)
  library(reshape2)

  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(cowplot)
  library(patchwork)
  library(ComplexHeatmap)
})


#---- functions ---------------------------------------------------------------

plot_PCA <- function(vsd, variable, ntop = 500, size = 14) {

  #print(variable)

  # Build the PCA data from the ntop most variable genes (row wise)
  rv <- rowVars(assay(vsd))
  select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(assay(vsd)[select, ]))
  var_perc <- round(pca$sdev^2 / sum(pca$sdev^2) * 100)

  data <- data.frame(
    PC1 = pca$x[, 1],
    PC2 = pca$x[, 2],
    colData(vsd)
  )
  data$sample <- sub("sperm_", "", data$sample)

  if (variable == "treatment") {
    title <- stringr::str_to_title(variable)
    fill <- "DBP [mg/kg]"
  } else if (variable == "extraction_batch") {
    title <- stringr::str_to_title(variable)
    fill <- "Batch"
  } else if (variable == "protocol") {
    title <- stringr::str_to_title(variable)
    fill <- "Protocol"
  } else if (variable == "Mreads") {
    title <- "Million sequencing reads per sample"
    fill <- variable
  }

  ggplot(data, aes(PC1, PC2, fill = .data[[variable]], label = sample)) +
    geom_label_repel(size = 4, show.legend = FALSE) +
    geom_point(size = 5, shape = 21) +
    theme_linedraw(size) +
    theme(plot.title = element_text(hjust = 0.5), legend.position = "bottom") +
    background_grid() +
    labs(
      title = title,
      fill = fill,
      x = paste("PC1:", var_perc[1], "% var."),
      y = paste("PC2:", var_perc[2], "% var.")
    )
}

plot_cooks_dist <- function(dds, title = NULL, size = 12) {
  # distance plot
  # Cook’s distance is a measure of how much a single sample is influencing the 
  # fitted coefficients for a gene, and a large value of Cook’s distance is 
  # intended to indicate an outlier count
  #par(mar=c(8,5,2,2))
  #boxplot(log10(assays(dds)[["cooks"]]), range=0, las=2, main = title)
  
  data <- log10(assays(dds)[["cooks"]])
  data <- melt(data, id.vars = rownames(data))
  
  ggplot(data, aes(Var2, value)) +
    stat_boxplot(geom = "errorbar", width = 0.3, coef = 10, lty = 1) +
    geom_boxplot(outlier.shape = NA, width = 0.5, fill = "lightblue", coef = 0) +
    theme_linedraw(size) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
    ) +
    background_grid() +
    labs(x = "", y = "", title = title)
}

plot_dist <- function(vsd, fontsize = 10) {
  # heatmap of sample-sample euclidean distance
  dists <- dist(t(assay(vsd)))
  mat <- as.matrix(dists)
  rownames(mat) <- vsd$treatment
  col <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
  
  hm <- pheatmap(
    mat, 
    clustering_distance_rows = dists,
    clustering_distance_cols = dists,
    col = col,
                 
    fontsize = fontsize,

    show_column_dend = FALSE, show_row_dend = FALSE,

    name = "Dist"
  )
  #grid.grabExpr(draw(hm, merge_legend = TRUE, heatmap_legend_side = "bottom"))
  grid.grabExpr(draw(hm))
}

plot_wrapper <- function(vsd, dds, title) {

  #cat("title in plot_wrapper ", title, "\n")

  #if (title == "Batch corrected") var <- "extraction_batch"
  #if (title == "Protocol corrected") var <- "protocol"
  #if (title == "Full model") var <- "protocol"

  #cat("var in plot_wrapper:", var, "\n")

  p1 <- plot_PCA(vsd, "treatment")
  p2 <- plot_PCA(vsd, "extraction_batch")
  p3 <- plot_PCA(vsd, "protocol")
  #p2 <- plot_PCA(vsd, var)
  #p3 <- plot_PCA(vsd, "Mreads")

  p4 <- plot_cooks_dist(dds)
  p5 <- plot_dist(vsd)

  # Add small distance on the left side with cowplot:plot_grid()
  # The distance plot is a grob-object unable to move for the annotation
  # automatically
  p5 <- plot_grid(NULL, p5, rel_widths = c(0.05, 1))

  layout <- "
    AAADDDDD
    AAADDDDD
    AAAEEEEE
    BBBEEEEE
    BBBEEEEE
    BBBEEEEE
    CCCEEEEE
    CCCEEEEE
    CCCEEEEE
  "

  wrap_plots(p1,p2,p3,p4,p5, design = layout) +
    plot_layout(axis_titles = "collect") +
    plot_annotation(
      title = title,
      tag_levels = "A",
      theme = theme(plot.title = element_text(size = 20, hjust = 0.5))
    )
}

#---- code --------------------------------------------------------------------

dds <- readRDS("DDS.Rds")
vsd <- varianceStabilizingTransformation(dds)


#==== Batch corrections =======================================================
# Remove extraction batch
dds_extraction_batch <- dds
design(dds_extraction_batch) <- ~ extraction_batch + treatment
vsd_extraction_batch <- varianceStabilizingTransformation(dds_extraction_batch)
vsd_batch_corr <- vsd_extraction_batch
vsd_mat <- assay(vsd_extraction_batch)
vsd_mat <- limma::removeBatchEffect(vsd_mat, batch = vsd$extraction_batch)
assay(vsd_batch_corr) <- vsd_mat

# Remove protocol batch
dds_protocol <- dds
design(dds_protocol) <- ~ protocol + treatment
vsd_protocol <- varianceStabilizingTransformation(dds_protocol)
vsd_protocol_corr <- vsd_protocol
vsd_mat <- assay(vsd_protocol_corr)
vsd_mat <- limma::removeBatchEffect(vsd_mat, batch = vsd$protocol)
assay(vsd_protocol_corr) <- vsd_mat
#==============================================================================


filename <- "QC_plots.pdf"
pdf(filename, width = 15, height = 15)

plot_wrapper(vsd, dds, "Full model")
plot_wrapper(vsd_batch_corr, dds, "Batch corrected")
plot_wrapper(vsd_protocol_corr, dds, "Protocol corrected")

dev.off()

cat(paste(
  "\n~~ DESeq2_QC_plots.R ~~~~~~~~~~~~~~~~~~~~~~~\n\n",
  filename, "generated",
  "\n\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n\n"
))