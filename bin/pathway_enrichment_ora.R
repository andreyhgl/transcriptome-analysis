#!/usr/bin/env Rscript

#==============================================================================
# pathway_enrichment_ora.R
#
# Over-representation analysis of the significant genes of one DESeq2
# contrast in GO, KEGG and Reactome. Up, Down and Both (up and down together)
# are tested separately, against all genes in the results table as the
# background (universe). One table per contrast, in TSV and Excel; the Rds
# also holds the gene set sizes and settings, for the Quarto report.
#
# Input:
#   - results table for one contrast from deseq2_results.R, with gene_id
#     (Ensembl), log2FoldChange and padj
#   - optional annotation csv/tsv with gene_id and gene_name columns; genes
#     missing from it are named with the OrgDb symbol
#
# Method:
#   clusterProfiler (doi: 10.1089/omi.2011.0118, 10.1016/j.xinn.2021.100141)
#   and ReactomePA (doi: 10.1039/C5MB00663E). GO is run per ontology (BP, CC,
#   MF) and BH-adjusted over the pooled terms; KEGG and Reactome use their
#   own BH adjustment. Terms are kept with p.adjust < --alpha and at least
#   --min_count genes. Gene sets with fewer than --min_count genes are not
#   tested. KEGG is downloaded from the KEGG REST API at run time, so the
#   task needs internet access unless KEGG is left out of --databases.
#
# Usage:
# > Rscript pathway_enrichment_ora.R --results DDS_a_vs_control.tsv \
#     --species mouse --contrast a_vs_control --prefix DDS_a_vs_control
# > Rscript pathway_enrichment_ora.R --results x.tsv --databases GO,Reactome
#
# Container:
# > docker://ghcr.io/andreyhgl/enrichment-env:main
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
  prefix = "enrichment",
  results = "results.tsv",
  contrast = "",
  species = "human",
  annotation = "",
  databases = "GO,KEGG,Reactome",
  gene_alpha = 0.05,
  alpha = 0.05,
  min_count = 3L
)

opt <- parse_args(defaults)

cat("\n > Arguments:\n---------------------------------\n")
cat(sprintf(" > %-12s : %s\n", names(opt), unlist(opt)), sep = "")
cat(strrep("-", 33), "\n")


#---- Environment -------------------------------------------------------------

# OrgDb package, KEGG organism code and Reactome organism per species
species_db <- list(
  human = list(orgdb = "org.Hs.eg.db", kegg = "hsa", reactome = "human"),
  mouse = list(orgdb = "org.Mm.eg.db", kegg = "mmu", reactome = "mouse")
)

if (!opt$species %in% names(species_db)) {
  stop(
    "--species must be one of: ", toString(names(species_db)),
    ", not '", opt$species, "'"
  )
}
sp <- species_db[[opt$species]]

suppressPackageStartupMessages({
  library(AnnotationDbi)
  library(clusterProfiler)
  library(ReactomePA)
  library(openxlsx)
  library(sp$orgdb, character.only = TRUE)
})

orgdb <- get(sp$orgdb)


#---- Settings ----------------------------------------------------------------

# Gene regulation type: which directions of change go into each gene set
types <- list(
  Up = "Up",
  Down = "Down",
  Both = c("Up", "Down")
)

databases <- strsplit(opt$databases, ",")[[1]] |> trimws()
known_databases <- c("GO", "KEGG", "Reactome")
if (!all(databases %in% known_databases)) {
  stop(
    "--databases takes a comma-separated subset of: ",
    toString(known_databases)
  )
}
names(databases) <- databases

contrast <- if (opt$contrast == "") {
  sub("\\.tsv$", "", basename(opt$results))
} else {
  opt$contrast
}

# Columns of the output tables, in order
columns <- c(
  "contrast",
  "type",
  "database",
  "ONTOLOGY",
  "ID",
  "Description",
  "GeneRatio",
  "BgRatio",
  "pvalue",
  "p.adjust",
  "Count",
  "gene_name",
  "gene_id"
)

# Output files
out_rds <- paste0(opt$prefix, ".Rds")
out_table <- paste0(opt$prefix, ".tsv")
out_xlsx <- paste0(opt$prefix, ".xlsx")


#---- Functions ---------------------------------------------------------------

print_section <- function(title) {
  cat("\n > ", title, "\n", strrep("-", 33), "\n", sep = "")
}

# Ensembl IDs without a version suffix: ENSG00000141510.17 => ENSG00000141510
strip_version <- function(ids) sub("^(ENS[A-Z]*G\\d+)\\.\\d+$", "\\1", ids)

# Optional annotation as a named vector gene_id => gene_name; NULL without
read_annotation <- function(path) {
  if (path == "") return(NULL)

  sep <- if (grepl("\\.csv(\\.gz)?$", path)) "," else "\t"
  anno <- read.delim(path, sep = sep)

  if (!all(c("gene_id", "gene_name") %in% colnames(anno))) {
    stop("Annotation needs `gene_id` and `gene_name` columns: ", path)
  }

  anno$gene_id <- strip_version(anno$gene_id)
  anno <- anno[!duplicated(anno$gene_id), ]
  setNames(as.character(anno$gene_name), anno$gene_id)
}

# Ensembl IDs of a gene set with their Entrez IDs (NA when unmapped)
make_genelist <- function(ensembl, orgdb) {
  ensembl <- unique(ensembl[!is.na(ensembl)])
  known <- ensembl[ensembl %in% keys(orgdb, keytype = "ENSEMBL")]

  entrez <- if (length(known) > 0) {
    mapIds(
      orgdb,
      keys = known,
      keytype = "ENSEMBL",
      column = "ENTREZID",
      multiVals = "first"
    )
  } else {
    character(0)
  }

  data.frame(
    ensembl = ensembl,
    entrez = unname(entrez[ensembl]),
    row.names = NULL
  )
}

# Gene names for the universe: annotation first, then the OrgDb symbol, and
# the Ensembl ID when neither has one
make_name_map <- function(ensembl, orgdb, anno) {
  name_map <- setNames(ensembl, ensembl)
 
  known <- ensembl[ensembl %in% keys(orgdb, keytype = "ENSEMBL")]
  if (length(known) > 0) {
    symbol <- mapIds(
      orgdb,
      keys = known,
      keytype = "ENSEMBL",
      column = "SYMBOL",
      multiVals = "first"
    )
    symbol <- symbol[!is.na(symbol)]
    name_map[names(symbol)] <- symbol
  }
 
  if (!is.null(anno)) {
    hit <- intersect(ensembl, names(anno))
    hit <- hit[!is.na(anno[hit]) & anno[hit] != ""]
    name_map[hit] <- anno[hit]
  }
 
  name_map
}

# clusterProfiler result => data.frame tagged with the ontology
as_table <- function(res, ontology) {
  if (is.null(res)) return(data.frame())
 
  tbl <- as.data.frame(res)
  if (nrow(tbl) > 0) tbl$ONTOLOGY <- ontology
  tbl
}

# All terms are returned (cutoffs = 1) and filtered afterwards, so GO can be
# adjusted over the three ontologies together
run_go <- function(genelist, universe, orgdb) {
  tables <- lapply(c("BP", "CC", "MF"), function(ont) {
    res <- enrichGO(
      gene = genelist$ensembl,
      OrgDb = orgdb,
      keyType = "ENSEMBL",
      ont = ont,
      universe = universe$ensembl,
      readable = FALSE,
      pvalueCutoff = 1,
      qvalueCutoff = 1
    )
    as_table(res, ont)
  })

  res <- do.call(rbind, tables)
  if (nrow(res) == 0) return(res)

  res$p.adjust <- p.adjust(res$pvalue, method = "BH")
  res
}

run_kegg <- function(genelist, universe, organism) {
  res <- enrichKEGG(
    gene = unique(na.omit(genelist$entrez)),
    organism = organism,
    universe = unique(na.omit(universe$entrez)),
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  res <- as_table(res, "KEGG")

  # Older versions append the species, e.g. " - Mus musculus (house mouse)"
  if (nrow(res) > 0) {
    species_suffix <- " - [A-Z][a-z]+ [a-z]+ \\(.*\\)$"
    res$Description <- sub(species_suffix, "", res$Description)
  }
  res
}

run_reactome <- function(genelist, universe, organism) {
  res <- enrichPathway(
    gene = unique(na.omit(genelist$entrez)),
    organism = organism,
    universe = unique(na.omit(universe$entrez)),
    readable = FALSE,
    pvalueCutoff = 1,
    qvalueCutoff = 1
  )
  as_table(res, "Reactome")
}

# "id1/id2/id3" => vector of Ensembl IDs, mapping Entrez IDs back through the
# gene set so that only genes of the set are reported
split_genes <- function(gene_ids, genelist, keytype) {
  lapply(strsplit(gene_ids, "/", fixed = TRUE), function(ids) {
    if (keytype == "ENSEMBL") return(ids)
    unique(genelist$ensembl[genelist$entrez %in% ids])
  })
}

# Enriched terms of one gene set in one database (NULL when there are none),
# with the number of terms tested (terms with at least one gene of the set,
# which the BH adjustment runs over) and the number of hits
enrich <- function(database, genelist, universe) {
  cat("   > Running", database, "analysis...\n")
 
  res <- switch(
    database,
    GO = run_go(genelist, universe, orgdb),
    KEGG = run_kegg(genelist, universe, sp$kegg),
    Reactome = run_reactome(genelist, universe, sp$reactome)
  )
  tested <- nrow(res)
 
  if (nrow(res) > 0) {
    res <- res[res$p.adjust < opt$alpha & res$Count >= opt$min_count, ]
  }
  cat("    ", nrow(res), "hits of", tested, "tested terms\n")
 
  out <- list(terms = NULL, hits = nrow(res), tested = tested)
  if (nrow(res) == 0) return(out)
 
  keytype <- if (database == "GO") "ENSEMBL" else "ENTREZID"
  res$gene_id <- split_genes(res$geneID, genelist, keytype)
  res$gene_name <- lapply(res$gene_id, function(ids) unname(name_map[ids]))
  res$database <- database
  res$contrast <- contrast
 
  out$terms <- res[order(res$p.adjust), setdiff(columns, "type")]
  out
}

# Hits and tested terms per database for one type; NA when it was skipped
type_stats <- function(type, n_genes, hits = NA, tested = NA) {
  data.frame(
    type = type,
    genes = n_genes,
    database = databases,
    hits = hits,
    tested = tested,
    row.names = NULL
  )
}

# One row per type, one "hits / tested" column per database
summary_table <- function(stats) {
  tbl <- unique(stats[, c("type", "genes")])
 
  for (db in databases) {
    s <- stats[stats$database == db, ]
    s <- s[match(tbl$type, s$type), ]
    tbl[[db]] <- ifelse(
      is.na(s$tested),
      "skipped",
      sprintf("%d / %d", s$hits, s$tested)
    )
  }
 
  rownames(tbl) <- NULL
  tbl
}

# Output table with no rows, so that every output file is still written
empty_table <- function() {
  tbl <- data.frame(
    contrast = character(0),
    type = character(0),
    database = character(0),
    ONTOLOGY = character(0),
    ID = character(0),
    Description = character(0),
    GeneRatio = character(0),
    BgRatio = character(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    Count = integer(0)
  )
  tbl$gene_name <- list()
  tbl$gene_id <- list()
  tbl
}

# List columns => one string per term
collapse_genes <- function(tbl, sep) {
  tbl$gene_name <- vapply(tbl$gene_name, paste, character(1), collapse = sep)
  tbl$gene_id <- vapply(tbl$gene_id, paste, character(1), collapse = sep)
  tbl
}


#---- Import datasets ---------------------------------------------------------

deg <- read.delim(opt$results)
deg$gene_id <- strip_version(deg$gene_id)

anno <- read_annotation(opt$annotation)


#---- Code --------------------------------------------------------------------

# Background (universe):
# Tested genes from the results table w/ entrez ID
universe <- make_genelist(deg$gene_id, orgdb)
name_map <- make_name_map(universe$ensembl, orgdb, anno)

if (all(is.na(universe$entrez))) {
  stop(
    "None of the gene IDs are Ensembl IDs in ", sp$orgdb,
    ". Check --species."
  )
}

sig <- !is.na(deg$padj) & deg$padj < opt$gene_alpha

#----
# add lfc threshold here:
# if (opt$lfc != "") {
#    sig <- sig & abs(deg$log2FoldChange) > opt$lfc
# }
#----

direction <- ifelse(deg$log2FoldChange > 0, "Up", "Down")

print_section("Contrast")
cat(sprintf(" > %-12s : %s\n", "Contrast", contrast))
cat(sprintf(" > %-12s : %d genes\n", "Universe", nrow(universe)))
cat(sprintf(
  " > %-12s : %d with Entrez IDs\n",
  "",
  sum(!is.na(universe$entrez))
))
cat(sprintf(" > %-12s : %s\n", "Databases", toString(databases)))

enrichment <- lapply(names(types), function(type) {
  genes <- deg$gene_id[sig & direction %in% types[[type]]]
 
  print_section(sprintf(
    "%s | %s | %d sign genes",
    contrast,
    type,
    length(genes)
  ))
 
  # No term can reach --min_count genes in a smaller gene set
  if (length(genes) < opt$min_count) {
    cat("   > Too few genes, skipped\n")
    return(list(terms = NULL, stats = type_stats(type, length(genes))))
  }
 
  genelist <- make_genelist(genes, orgdb)
 
  res <- lapply(databases, enrich, genelist = genelist, universe = universe)
  tbl <- do.call(rbind, lapply(res, `[[`, "terms"))
  if (!is.null(tbl)) tbl$type <- type
 
  stats <- type_stats(
    type,
    length(genes),
    hits = vapply(res, `[[`, integer(1), "hits"),
    tested = vapply(res, `[[`, integer(1), "tested")
  )
 
  list(terms = tbl, stats = stats)
})

out <- do.call(rbind, lapply(enrichment, `[[`, "terms"))
if (is.null(out)) out <- empty_table()
 
out <- out[, columns]
rownames(out) <- NULL

stats <- do.call(rbind, lapply(enrichment, `[[`, "stats"))
 
print_section(sprintf(
  "Enriched terms (p.adjust < %s, Count >= %d): hits / tested",
  opt$alpha,
  opt$min_count
))
print(summary_table(stats), row.names = FALSE)


#---- Save files --------------------------------------------------------------

# Rds for the report: the terms (gene_name and gene_id as list columns), the
# gene set sizes, the background and the settings of this run
saveRDS(
  list(
    terms = out,
    stats = data.frame(contrast = contrast, stats),
    universe = data.frame(
      contrast = contrast,
      genes = nrow(universe),
      entrez = sum(!is.na(universe$entrez))
    ),
    settings = list(
      species = opt$species,
      orgdb = sp$orgdb,
      databases = unname(databases),
      gene_alpha = opt$gene_alpha,
      alpha = opt$alpha,
      min_count = opt$min_count,
      versions = c(
        clusterProfiler = as.character(packageVersion("clusterProfiler")),
        ReactomePA = as.character(packageVersion("ReactomePA")),
        orgdb = as.character(packageVersion(sp$orgdb))
      )
    )
  ),
  file = out_rds
)

write.table(
  collapse_genes(out, "\n"),
  file = out_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Excel: one sheet per type, one gene per line within a cell
wb <- createWorkbook()

for (type in c("Both", "Up", "Down")) {
  addWorksheet(wb, type)
  writeData(wb, type, collapse_genes(out[out$type == type, ], "\n"))
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)


#---- Done --------------------------------------------------------------------

cat(sprintf(
  "\n > Script complete!\n%s\n > Outputs:\n%s",
  strrep("-", 33),
  paste0("   - ", c(out_rds, out_table, out_xlsx), "\n", collapse = "")
))