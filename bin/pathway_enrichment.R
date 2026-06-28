#!/usr/bin/env Rscript

#------------------------------------------------------------------------------
# GENE ENRICHMENT ANALYSIS
#
# Extract GO results from GO, KEGG, Reactome databases
# rbind() resultsinto 1 data.frame per contrast
# Save as list, w/ structure: gen$contrast$type$database
# Contrasts: ~ batch + treatment = treatment100_vs_0
# Type: Up, Down, Both
# Database: GO, KEGG, Reactome
# Background genes (universe) = all tested genes
#
# Settings:
#   > 2 genes required to run GO analysis
#   qvalueCutoff = 0.2
# 
# Usage:
# > Rscript script.R --parameter1 file1
# 
# Container:
# > library://andreyhgl/singularity-r/rnaseq
#------------------------------------------------------------------------------


#---- Parse arguments ---------------------------------------------------------


args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 1) {
  stop(" > Script needs atleast one argument! Usage: script.R --parameter1 file1")
}

for (i in seq_along(args)) {
  if (args[i] == "--species") species <- args[i + 1]
}

cat(glue::glue("
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    > species         : {species}
  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"))


#---- Environment ------------------------------------------------------------


suppressPackageStartupMessages({
  library(DESeq2)
  library(gtools)
  library(openxlsx)

  # clusterProfiler
  # doi: 10.1089/omi.2011.0118
  # doi: 10.1016/j.xinn.2021.100141
  library(clusterProfiler)

  # reactome
  # doi: 10.1039/C5MB00663E
  library(ReactomePA)

  # https://bioconductor.org/packages//2.7/data/annotation/manuals/org.Mm.eg.db/man/org.Mm.eg.db.pdf
  if (species == 'mouse') library(org.Mm.eg.db)

  # https://bioconductor.org/packages//release/data/annotation/manuals/org.Hs.eg.db/man/org.Hs.eg.db.pdf
  if (species == 'human') library(org.Hs.eg.db)
})


#---- Functions ---------------------------------------------------------------


enrichment_analysis <- function(genelist, universe, database, species) {
  # database takes in "GO", "KEGG", "Reactome"
  # species takes in "human", "mouse"  

  if ( all(sapply(genelist, length) < 1) ) {
    cat("\t\t> No hits!\n")
    return(NULL)
  }

  cat(" > Running", database, "analysis...\n")

  if (database == "KEGG") {
    if (species == "human") org <- "hsa"
    if (species == "mouse") org <- "mmu"

    res <- clusterProfiler::enrichKEGG(
      gene = unname(genelist$entrez),
      organism = org,
      universe = unname(universe$entrez)
    )
    res <- data.frame(res)

    cat("    >", nrow(res), "hits\n")

    if (nrow(res) == 0) return(NULL)

    if (ncol(res) > 0) res <- subset(res, Count > 2)

    # map entrez ID => gene name
    res$geneID <- unlist(sapply(res$geneID, simplify = FALSE, \(gene){
      gene <- strsplit(gene, "/")[[1]]
      gene <- mapIds(org.Mm.eg.db, keys = gene, column = "SYMBOL" , keytype = "ENTREZID")
      paste(gene, collapse = "/")
    }))

    res$ensembl_gene_id <- sapply(
      res$geneID,
      simplify = FALSE,
      function(gene) {
        geneID <- strsplit(gene, "/")[[1]]
        ens[ens$gene_name %in% geneID, "ensembl_gene_id"]
    }) %>% unname

    # map to gene name
    res$gene_name <- sapply(
      res$geneID,
      simplify = FALSE,
      function(gene) {
        strsplit(gene, "/")[[1]]
    }) %>% unname

    res$ONTOLOGY <- database

  }

  if (database == "GO") {
    if (species == "human") org <- org.Hs.eg.db
    if (species == "mouse") org <- org.Mm.eg.db

    #==== Quick'n'dirty fix for ont = "ALL" error =============================

    onts <- c("BP", "CC", "MF")
    names(onts) <- onts

    # Run each ontology WITHOUT a q-value cutoff so nothing is pre-filtered out
    # Apply multiple testing correction manually
    res_list <- lapply(onts, function(o) {
      clusterProfiler::enrichGO(
        gene          = genelist$ensembl,
        OrgDb         = org,
        keyType       = "ENSEMBL",
        ont           = o,
        universe      = universe$ensembl,
        readable      = FALSE,
        pvalueCutoff  = 1,
        qvalueCutoff  = 1
      )
    })

    # Pool all terms into one data frame, tagging ontology
    res <- do.call(rbind, lapply(onts, function(o) {
      df <- as.data.frame(res_list[[o]])
      if (nrow(df) > 0) df$ONTOLOGY <- o
      df
    }))

    res$p.adjust <- p.adjust(res$pvalue, method = "BH")

    res <- res[res$p.adjust < 0.05, ]
    res <- res[order(res$p.adjust), ]

    cat("    >", nrow(res), "hits\n")

    if (nrow(res) == 0) return(NULL)
    if (ncol(res) > 0) res <- subset(res, Count > 2)

    # Combine into one table
    res$ensembl_gene_id <- sapply(
      res$geneID,
      simplify = FALSE,
      function(gene) {
        strsplit(gene, "/")[[1]]
    }) %>% unname

    res$gene_name <- sapply(
      res$geneID,
      simplify = FALSE,
      function(gene) {
        geneID <- strsplit(gene, "/")[[1]]
        ens[ens$ensembl_gene_id %in% geneID, "gene_name"]
    }) %>% unname
  }

  if (database == "Reactome") {
    res <- ReactomePA::enrichPathway(
      gene = genelist$entrez,
      organism = species,
      universe = universe$entrez,
      readable = TRUE
    )
    res <- data.frame(res)
    
    cat("    >", nrow(res), "hits\n")

    if (nrow(res) == 0) return(NULL)

    if (ncol(res) > 0) res <- subset(res, Count > 2)
    res$ONTOLOGY <- database

    res$gene_name <- sapply(
      res$geneID,
        simplify = FALSE,
        function(gene) strsplit(gene, "/")[[1]]
      )

    res$ensembl_gene_id <- sapply(
      res$gene_name,
      simplify = FALSE,
      function(gene) {
        rows <- match(gene, ens$gene_name)
        ens$ensembl_gene_id[rows]
    })
  }

  res$geneID <- NULL
  res$database <- database

  columns <- c(
    "database",
    "ONTOLOGY",
    "ID",
    "Description",
    "GeneRatio",
    "BgRatio",
    "pvalue",
    "p.adjust",
    #"qvalue",
    "Count",
    "gene_name",
    "ensembl_gene_id"
  )

  res <- res[, columns]
  #res <- res[order(-res$Count, res$qvalue), ]
  rownames(res) <- NULL

  return(res)
}

make_genelist <- function(genes, species) {
  if ( length(genes) < 2 ) return(NULL)
  
  if (species == "human") org <- org.Hs.eg.db
  if (species == "mouse") org <- org.Mm.eg.db

  ensembl <- na.omit(unique(genes))
  entrez <- mapIds(
    x = org,
    keys = ensembl,
    keytype="ENSEMBL",
    column = "ENTREZID"
  )
  entrez <- na.omit(entrez)
  
  return(list(ensembl = ensembl, entrez = entrez))
}


#---- Import datasets ---------------------------------------------------------


ens <- read.csv("ensembl_table.csv.gz")
genexp <- read.csv("genexp_table.csv.gz")
dds <- readRDS("DDS.Rds")


#---- Code --------------------------------------------------------------------


# Contrasts
contrasts <- unique(genexp$contrast)
names(contrasts) <- contrasts

cat(sprintf(
  " > Contrasts:\n%s\n",
  paste0("   - ", contrasts, collapse = "\n")
))

# Gene Regulation type
types <- list(
  Down = "Down",
  Up = "Up",
  Both = c("Up", "Down")
)

# Ontology databases
databases <- c("KEGG", "GO", "Reactome")
names(databases) <- databases


#---- Calculate pathway enrichment --------------------------------------------


pathway_enrichment <- lapply(contrasts, function(contr) {
  # Background genelist (universie) to test enrichment against
  universe <- make_genelist(rownames(dds), species)

  lapply(seq_along(types), function(index) {
    type <- types[index]

    rows <- genexp$contrast %in% contr & genexp$type %in% unlist(type) & genexp$sign == "sign"
    genes <- genexp$ensembl_gene_id[rows]
      
    cat(sprintf(
      "\n > Testing: %s\n\n",
      paste(
        names(contrasts)[contrasts %in% contr],
        names(types)[index],
        paste(sum(rows), "sign genes"),
        sep = " | ")
    ))

    # Genelist ofsignificant genes
    genelist <- make_genelist(genes, species)

    res <- lapply(databases, function(database) {
      enrichment_analysis(genelist, universe, database, species)
    })
    #Reduce(function(x,y) rbind(x,y), res)
    res <- do.call(rbind, res)

    if (!is.null(res)) {
      res$contrast <- contr
      res$type <- names(types)[index]
    }

    return(res)
  })
})

res <- unlist(pathway_enrichment, recursive = FALSE)

# Drop NULLs
res <- res[!sapply(res, is.null)]
res <- do.call(rbind, res)
rownames(res) <- NULL


#---- Save files --------------------------------------------------------------


# Fix gene_name and ensembl_gene_id columns, from list => normal
out <- res
out$gene_name <- sapply(out$gene_name, function(x) paste(x, collapse = "\n"))
out$ensembl_gene_id <- sapply(out$ensembl_gene_id, function(x) paste(x, collapse = "\n"))



# Excel
excel_file <- "pathway_enrichment_S2.xlsx"
wb <- createWorkbook()

addWorksheet(wb, "Both")
writeData(wb, "Both", out[out$type == "Both", ])

addWorksheet(wb, "Up")
writeData(wb, "Up", out[out$type == "Up", ])

addWorksheet(wb, "Down")
writeData(wb, "Down", out[out$type == "Down", ])

saveWorkbook(wb, excel_file)

# Rds
filename <- "pathway_enrichment.Rds"
saveRDS(res, file = filename)


#---- Done --------------------------------------------------------------------


cat(glue::glue("
 ~~ pathway_enrichment.R complete ~~~~~~~~~~~~~~~~~~~~~~~
  > Enrichment hits   : {filename}
  > Saved files       :
                      - {filename}
                      - {excel_file}
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

"))