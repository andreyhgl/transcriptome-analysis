#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# Differential gene expression analysis with DESeq2
# Save significant genes as:
# > excel (supplementary): w/ sheets for all sign. and the common sign. genes
# > csv: all sign. genes 
#------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(DESeq2)
  library(gtools)
  library(stringr)
  library(scales)
  library(openxlsx)
  library(data.table)
})


#---- Functions ---------------------------------------------------------------


extract_DGE <- function(resname) {
  # DGE differentially expressed genes
  cat(" > Calculating DGE for", resname, "\n")
  
  res <- results(dds, name = resname)
  dat <- data.frame(res)
  dat$ensembl_gene_id <- rownames(res)
  dat$sign <- "NS"
  dat$sign[dat$padj <= 0.05 & abs(dat$log2FoldChange) > 0] <- "sign"
  dat$type <- ifelse(dat$log2FoldChange < 0, "Down", "Up")
  dat$contrast <- resname
  out <- merge(dat, ens, by = "ensembl_gene_id")
  out <- out[mixedorder(paste0(out$chr, out$start)), ]
  rownames(out) <- NULL

  print(table(out$sign))

  return(out)
}

find_common_genes <- function(DBP10, DBP100) {
  dat10 <- DBP10[DBP10$sign == "sign", ]
  dat100 <- DBP10[DBP100$sign == "sign", ]

  genes <- intersect(dat10$ensembl_gene_id, dat100$ensembl_gene_id)

  cat(" > Num. common genes:", length(genes), "\n")
  cat(sprintf(
    " > DBP10 %s%%\n > DBP100 %s%%\n",
    round(length(genes) / nrow(dat10) * 100, 2),
    round(length(genes) / nrow(dat100) * 100, 2)
  ))
  #print(table(DBP10$sign))
  #print(table(DBP100$sign))

  out <- ens[ens$ensembl_gene_id %in% genes, ]

  out$type <- dat10[dat10$ensembl_gene_id %in% genes, "type"]

  # check if gene expression type is in the same direction
  rows <- dat10[dat10$ensembl_gene_id %in% genes, "type"] == dat100[dat100$ensembl_gene_id %in% genes, "type"]

  # if not, paste = type 10 ; type 100
  if (any(!rows)) {
    out$type[!rows] <- paste0(
      dat10[dat10$ensembl_gene_id %in% genes, "type"],
      ";",
      dat100[dat100$ensembl_gene_id %in% genes, "type"]
    )
  }

  return(out)
}


#---- Import datasets ---------------------------------------------------------


dds <- readRDS("DDS.Rds")
ens <- read.csv("ensembl_table.csv.gz")


#---- Code --------------------------------------------------------------------


resnames <- grep("treatment", resultsNames(dds), value = TRUE)

DBP10 <- extract_DGE(resnames[1])
DBP100 <- extract_DGE(resnames[2])

common_genes <- find_common_genes(DBP10, DBP100)

sign_genes <- c(
  DBP10[DBP10$sign == "sign", "ensembl_gene_id"],
  DBP100[DBP100$sign == "sign", "ensembl_gene_id"]
)
sign_genes <- unique(sign_genes)

cat(" > Total number significant genes:", length(sign_genes), "\n")

genexp <- rbind(
  DBP10[DBP10$ensembl_gene_id %in% sign_genes, ],
  DBP100[DBP100$ensembl_gene_id %in% sign_genes, ]
)

genexp <- genexp[mixedorder(paste0(genexp$chr, "_", genexp$start)), ]


#---- Save excel file ---------------------------------------------------------


excel_file <- "significant_genes_S1.xlsx"
wb <- createWorkbook()
addWorksheet(wb, "all sign. genes")
addWorksheet(wb, "common sign. genes")
writeData(wb, "all sign. genes", genexp)
writeData(wb, "common sign. genes", common_genes)
saveWorkbook(wb, excel_file)


#---- Save csv ----------------------------------------------------------------


filename <- "genexp_table.csv.gz"
fwrite(genexp, file = filename)

cat(glue::glue("
 ~~ DESeq2_results.R complete ~~~~~~~~~~~~~~~~~
  > num. genes  : {scales::comma(nrow(genexp))}
  > Saved files :
                - {filename}
                - {excel_file}
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"))