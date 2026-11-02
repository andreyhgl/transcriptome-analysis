#!/usr/bin/env Rscript

#==============================================================================
# pathway_enrichment_gsea.R
#
# Gene set enrichment analysis (GSEA) of one DESeq2 contrast in GO, KEGG and
# Reactome. All genes in the results table are ranked by one score (the Wald
# statistic by default), with no significance cutoff. Terms enriched among
# up-regulated genes get a positive NES, terms among down-regulated genes a
# negative NES. One table per contrast, in TSV and Excel; the Rds also
# holds the hits per database, the ranked genes and the settings, for the
# Quarto report.
#
# Input:
#   - results table for one contrast from deseq2_results.R, with gene_id
#     (Ensembl) and the --rank_by column
#   - optional annotation csv/tsv with gene_id and gene_name columns; genes
#     missing from it are named with the OrgDb symbol
#
# Method:
#   clusterProfiler gseGO, gseKEGG and ReactomePA gsePathway, all with fgsea
#   (doi: 10.1101/060012). Genes with no score are left out. KEGG and
#   Reactome use Entrez IDs; when several genes map to one Entrez ID, the
#   gene with the largest absolute score is kept. GO is run per ontology
#   (BP, CC, MF) and BH-adjusted over the pooled terms; KEGG and Reactome use
#   their own BH adjustment. Terms are kept with p.adjust < --alpha. The
#   random seed is set before each database, so reruns give the same
#   p-values. With deseq2_lfc_threshold > 0, `stat` is the thresholded
#   statistic and is 0 for every gene within the threshold. KEGG is
#   downloaded from the KEGG REST API at run time, so the task needs internet
#   access unless KEGG is left out of --databases.
#
# Usage:
# > Rscript pathway_enrichment_gsea.R --results DDS_a_vs_control.tsv \
#     --species mouse --contrast a_vs_control --prefix DDS_a_vs_control_gsea
# > Rscript pathway_enrichment_gsea.R --results x.tsv \
#     --rank_by log2FoldChange_shrunk --databases GO,Reactome --cores 4
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
  prefix = "gsea",
  results = "results.tsv",
  contrast = "",
  species = "human",
  annotation = "",
  databases = "GO,KEGG,Reactome",
  rank_by = "stat",
  alpha = 0.05,
  min_size = 10L,
  max_size = 500L,
  seed = 1337L,
  cores = 1L
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
  library(BiocParallel)
  library(clusterProfiler)
  library(ReactomePA)
  library(openxlsx)
  library(sp$orgdb, character.only = TRUE)
})

orgdb <- get(sp$orgdb)

# fgsea runs on the registered backend; without this it uses every core on
# the node, not the cores given to the task
if (opt$cores > 1) {
  register(MulticoreParam(workers = opt$cores))
} else {
  register(SerialParam())
}


#---- Settings ----------------------------------------------------------------

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

# Columns of the output tables, in order. Count is the number of
# leading-edge (core enrichment) genes, listed in gene_name and gene_id
columns <- c(
  "contrast",
  "direction",
  "database",
  "ONTOLOGY",
  "ID",
  "Description",
  "setSize",
  "enrichmentScore",
  "NES",
  "pvalue",
  "p.adjust",
  "rank",
  "leading_edge",
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

# One OrgDb column for Ensembl IDs, NA when unmapped
map_ensembl <- function(ensembl, orgdb, column) {
  known <- ensembl[ensembl %in% keys(orgdb, keytype = "ENSEMBL")]
  if (length(known) == 0) return(rep(NA_character_, length(ensembl)))

  ids <- mapIds(
    orgdb,
    keys = known,
    keytype = "ENSEMBL",
    column = column,
    multiVals = "first"
  )
  unname(ids[ensembl])
}

# Gene names: annotation first, then the OrgDb symbol, and the Ensembl ID
# when neither has one
make_name_map <- function(ensembl, orgdb, anno) {
  name_map <- setNames(ensembl, ensembl)

  symbol <- map_ensembl(ensembl, orgdb, "SYMBOL")
  has_symbol <- !is.na(symbol)
  name_map[has_symbol] <- symbol[has_symbol]

  if (!is.null(anno)) {
    hit <- intersect(ensembl, names(anno))
    hit <- hit[!is.na(anno[hit]) & anno[hit] != ""]
    name_map[hit] <- anno[hit]
  }

  name_map
}

# Every gene with a score: ensembl, entrez and score, highest score first
make_ranking <- function(deg, rank_by, orgdb) {
  if (!rank_by %in% colnames(deg)) {
    numeric_cols <- colnames(deg)[vapply(deg, is.numeric, logical(1))]
    stop(
      "--rank_by '", rank_by, "' is not a column of the results table. ",
      "Numeric columns: ", toString(numeric_cols)
    )
  }

  ranking <- data.frame(ensembl = deg$gene_id, score = deg[[rank_by]])
  ranking <- ranking[!is.na(ranking$ensembl) & is.finite(ranking$score), ]
  ranking <- ranking[!duplicated(ranking$ensembl), ]
  ranking$entrez <- map_ensembl(ranking$ensembl, orgdb, "ENTREZID")

  ranking[order(ranking$score, decreasing = TRUE), ]
}

# One gene per Entrez ID: the one with the largest absolute score
dedupe_entrez <- function(ranking) {
  rnk <- ranking[!is.na(ranking$entrez), ]
  rnk <- rnk[order(abs(rnk$score), decreasing = TRUE), ]
  rnk <- rnk[!duplicated(rnk$entrez), ]
  rnk[order(rnk$score, decreasing = TRUE), ]
}

# Named, decreasing gene list for clusterProfiler, and the Ensembl gene
# behind each name (NULL for Ensembl, where they are the same)
gene_list <- function(ranking, keytype) {
  if (keytype == "ENSEMBL") {
    return(list(genes = setNames(ranking$score, ranking$ensembl)))
  }

  rnk <- dedupe_entrez(ranking)
  list(
    genes = setNames(rnk$score, rnk$entrez),
    id_map = setNames(rnk$ensembl, rnk$entrez)
  )
}

# clusterProfiler result => data.frame tagged with the ontology
as_table <- function(res, ontology) {
  if (is.null(res)) return(data.frame())

  tbl <- as.data.frame(res)
  if (nrow(tbl) > 0) tbl$ONTOLOGY <- ontology
  tbl
}

# All terms are returned (pvalueCutoff = 1) and filtered afterwards, so GO
# can be adjusted over the three ontologies together
run_go <- function(genes, orgdb) {
  tables <- lapply(c("BP", "CC", "MF"), function(ont) {
    set.seed(opt$seed)
    res <- gseGO(
      geneList = genes,
      OrgDb = orgdb,
      keyType = "ENSEMBL",
      ont = ont,
      minGSSize = opt$min_size,
      maxGSSize = opt$max_size,
      pvalueCutoff = 1,
      verbose = FALSE,
      seed = FALSE,
      by = "fgsea"
    )
    as_table(res, ont)
  })

  res <- do.call(rbind, tables)
  if (nrow(res) == 0) return(res)

  res$p.adjust <- p.adjust(res$pvalue, method = "BH")
  res
}

run_kegg <- function(genes, organism) {
  set.seed(opt$seed)
  res <- gseKEGG(
    geneList = genes,
    organism = organism,
    minGSSize = opt$min_size,
    maxGSSize = opt$max_size,
    pvalueCutoff = 1,
    verbose = FALSE,
    seed = FALSE,
    by = "fgsea"
  )
  res <- as_table(res, "KEGG")

  # Older versions append the species, e.g. " - Mus musculus (house mouse)"
  if (nrow(res) > 0) {
    species_suffix <- " - [A-Z][a-z]+ [a-z]+ \\(.*\\)$"
    res$Description <- sub(species_suffix, "", res$Description)
  }
  res
}

run_reactome <- function(genes, organism) {
  set.seed(opt$seed)
  res <- gsePathway(
    geneList = genes,
    organism = organism,
    minGSSize = opt$min_size,
    maxGSSize = opt$max_size,
    pvalueCutoff = 1,
    verbose = FALSE,
    seed = FALSE,
    by = "fgsea"
  )
  as_table(res, "Reactome")
}

# "id1/id2/id3" => vector of Ensembl IDs, in leading-edge order
split_genes <- function(core, id_map) {
  lapply(strsplit(core, "/", fixed = TRUE), function(ids) {
    if (is.null(id_map)) return(ids)
    unname(id_map[ids])
  })
}

# Enriched terms of the ranking in one database (NULL when there are none),
# with the number of ranked genes, the number of gene sets tested (sets
# within --min_size and --max_size, which the BH adjustment runs over) and
# the number of hits
enrich <- function(database, ranking) {
  cat("   > Running", database, "analysis...\n")

  keytype <- if (database == "GO") "ENSEMBL" else "ENTREZID"
  gl <- gene_list(ranking, keytype)
  genes <- gl$genes
  cat("     >", length(genes), "ranked genes\n")

  res <- switch(
    database,
    GO = run_go(genes, orgdb),
    KEGG = run_kegg(genes, sp$kegg),
    Reactome = run_reactome(genes, sp$reactome)
  )
  tested <- nrow(res)

  if (nrow(res) > 0) res <- res[res$p.adjust < opt$alpha, ]
  cat("     >", nrow(res), "hits of", tested, "tested gene sets\n")

  out <- list(
    terms = NULL,
    stats = data.frame(
      database = database,
      ranked = length(genes),
      hits = nrow(res),
      up = sum(res$NES > 0),
      down = sum(res$NES < 0),
      tested = tested
    )
  )
  if (nrow(res) == 0) return(out)

  res$gene_id <- split_genes(res$core_enrichment, gl$id_map)
  res$gene_name <- lapply(res$gene_id, function(ids) unname(name_map[ids]))
  res$Count <- lengths(res$gene_id)
  res$direction <- ifelse(res$NES > 0, "Up", "Down")
  res$database <- database
  res$contrast <- contrast

  out$terms <- res[order(res$p.adjust), columns]
  out
}

# Output table with no rows, so that every output file is still written
empty_table <- function() {
  tbl <- data.frame(
    contrast = character(0),
    direction = character(0),
    database = character(0),
    ONTOLOGY = character(0),
    ID = character(0),
    Description = character(0),
    setSize = integer(0),
    enrichmentScore = numeric(0),
    NES = numeric(0),
    pvalue = numeric(0),
    p.adjust = numeric(0),
    rank = integer(0),
    leading_edge = character(0),
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

ranking <- make_ranking(deg, opt$rank_by, orgdb)
name_map <- make_name_map(ranking$ensembl, orgdb, anno)

if (all(is.na(ranking$entrez))) {
  stop(
    "None of the gene IDs are Ensembl IDs in ", sp$orgdb,
    ". Check --species."
  )
}

print_section("Contrast")
cat(sprintf(" > %-12s : %s\n", "Contrast", contrast))
cat(sprintf(" > %-12s : %s\n", "Rank by", opt$rank_by))
cat(sprintf(" > %-12s : %d genes\n", "Ranked", nrow(ranking)))
cat(sprintf(
  " > %-12s : %d with Entrez IDs\n",
  "",
  sum(!is.na(ranking$entrez))
))
cat(sprintf(
  " > %-12s : %d up, %d down, %d zero\n",
  "Scores",
  sum(ranking$score > 0),
  sum(ranking$score < 0),
  sum(ranking$score == 0)
))
cat(sprintf(" > %-12s : %s\n", "Databases", toString(databases)))

print_section(sprintf("%s | GSEA", contrast))

enrichment <- lapply(databases, enrich, ranking = ranking)

out <- do.call(rbind, lapply(enrichment, `[[`, "terms"))
if (is.null(out)) out <- empty_table()

out <- out[, columns]
rownames(out) <- NULL

stats <- do.call(rbind, lapply(enrichment, `[[`, "stats"))
rownames(stats) <- NULL

print_section(sprintf(
  "Enriched gene sets (p.adjust < %s): hits / tested",
  opt$alpha
))
print(
  data.frame(
    database = stats$database,
    ranked = stats$ranked,
    up = stats$up,
    down = stats$down,
    hits = sprintf("%d / %d", stats$hits, stats$tested)
  ),
  row.names = FALSE
)


#---- Save files --------------------------------------------------------------

# Rds for the report, in the same layout as the ORA Rds: the terms
# (gene_name and gene_id as list columns), the hits per database, the
# ranked genes and the settings of this run
saveRDS(
  list(
    terms = out,
    stats = data.frame(contrast = contrast, stats),
    universe = data.frame(
      contrast = contrast,
      genes = nrow(ranking),
      entrez = sum(!is.na(ranking$entrez))
    ),
    settings = list(
      method = "GSEA",
      species = opt$species,
      orgdb = sp$orgdb,
      databases = unname(databases),
      rank_by = opt$rank_by,
      alpha = opt$alpha,
      min_size = opt$min_size,
      max_size = opt$max_size,
      seed = opt$seed,
      versions = c(
        clusterProfiler = as.character(packageVersion("clusterProfiler")),
        ReactomePA = as.character(packageVersion("ReactomePA")),
        fgsea = as.character(packageVersion("fgsea")),
        orgdb = as.character(packageVersion(sp$orgdb))
      )
    )
  ),
  file = out_rds
)

write.table(
  collapse_genes(out, "/"),
  file = out_table,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

# Excel: one sheet per direction, one gene per line within a cell
wb <- createWorkbook()

for (direction in c("Up", "Down")) {
  addWorksheet(wb, direction)
  writeData(
    wb,
    direction,
    collapse_genes(out[out$direction == direction, ], "\n")
  )
}

saveWorkbook(wb, out_xlsx, overwrite = TRUE)


#---- Done --------------------------------------------------------------------

cat(sprintf(
  "\n > Script complete!\n%s\n > Outputs:\n%s",
  strrep("-", 33),
  paste0("   - ", c(out_rds, out_table, out_xlsx), "\n", collapse = "")
))