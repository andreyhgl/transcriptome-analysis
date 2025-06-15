#!/usr/bin/env Rscript

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# print quality control plots
# PCA
# library size
# euclidean distance
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

suppressPackageStartupMessages({
  library(edgeR)
  library(ggplot2)
  library(ggrepel)
  library(forcats)
  library(RColorBrewer)
  library(cowplot)
  library(ComplexHeatmap)
  library(scales)
})

# ~~ functions ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

plot_MDS <- function(y, title = NULL, gene_sel = "common", top = 500, show_id = TRUE){
  data <- plotMDS(y, gene.selection = gene_sel, top = top, plot = FALSE)
  var_perc <- round(data$var.explained * 100, 1)
  data <- data.frame(
    PC1 = data$x,
    PC2 = data$y,
    treatment = y$samples$group,
    id = rownames(y$samples)
  )

  title <- ifelse(gene_sel == "common", "Principal component", "Principal coordinate")
  title <- paste(title, "analysis")
  caption <- paste("Analysis based on", comma(top))

  ggplot(data, aes(PC1, PC2, fill = treatment, label = id)) +
    { if (show_id) geom_label_repel(size = 4, show.legend = FALSE, alpha = 0.5) } +
    geom_point(shape = 21, size = 3) +
    theme_bw(12) + 
    { if (gene_sel == "common") theme(legend.position = "none") } +
    { if (gene_sel == "pairwise") theme(legend.position = "bottom") } +
    labs(
      title = title, caption = caption,
      fill = "DBP exposure [mg/kg]",
      x = paste("PC1:", var_perc[1], "% var."),
      y = paste("PC2:", var_perc[2], "% var.")
    )
}
plot_libsize <- function(y, title = NULL){
  data <- data.frame(
    libsize = y$samples$lib.size * 1e-6,
    treatment = y$samples$group,
    id = fct_inorder(rownames(y$samples))
  )

  ggplot(data, aes(id, libsize, fill = treatment, label = id)) +
    geom_col(color = "black") +
    #geom_text(aes(y = 1), vjust = 0.5, angle = 90, color = "white", size = 4) + 
    theme_bw(12) + theme(
      legend.position = "bottom",
      axis.title.x = element_blank(),
      axis.text.x = element_text(angle = 90, vjust = 0.5)
    ) +
    labs(
      title = title,
      fill = "DBP [mg/kg]",
      y = "\nLib. size (Millions)"
    )
}
plot_dist <- function(y){
  dists <- dist(t(cpm(y, log = TRUE)))
  mat <- as.matrix(dists)
  rownames(mat) <- y$samples$group
  col <- colorRampPalette(rev(brewer.pal(9, "YlOrRd")))(255)

  hm <- pheatmap(
    mat, 
    clustering_distance_rows = dists,
    clustering_distance_cols = dists,
    color = col,

    fontsize = 12,

    show_column_dend = FALSE, show_row_dend = FALSE, name = "Dist",
    heatmap_legend_param = list(direction = "horizontal", title_position = "topcenter")
  )
  grid.grabExpr(draw(hm, heatmap_legend_side = "bottom", merge_legend = TRUE))
}
generate_plots <- function(y){
  p1 <- plot_MDS(y, gene_sel = "common")    # PCA
  p2 <- plot_MDS(y, gene_sel = "pairwise")  # PCoA
  p3 <- plot_libsize(y)
  p4 <- plot_dist(y)

  # cowplot
  PCs <- plot_grid(p1, p2, ncol = 1, labels = c("A", "B"), rel_heights = c(0.9,1))
  libsize <- plot_grid(p3, nrow = 1, labels = c("C"))
  
  left <- plot_grid(PCs, libsize, ncol = 1, rel_heights = c(1, 0.3), align = "v")
  
  right <- plot_grid(NULL, p4, NULL, nrow = 1, rel_widths = c(0.05, 1, 0.05), labels = c("D", NA, NA))
  plot_grid(left, right, nrow = 1, rel_widths = c(0.8, 1), align = "v")
}

# ~~ code ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

cat(paste0("\n", "Loading datasets...", "\n\n"))

DGE <- readRDS("DGEList.Rds")

filename <- "QC_plots.pdf"

pdf(filename, width = 15, height = 10)
for (i in seq_along(DGE)){
  cat(paste("Generating plots for", names(DGE)[i], "\n"))
  print(generate_plots(DGE[[i]]))
}
dev.off()

cat(paste(
  "\n~~ QC_plots.R complete ~~~~~~~~~~~~~~~~~~~~~\n",
  "Output:", filename,
  "\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n"
))