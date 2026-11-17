#!/usr/bin/env Rscript

#==============================================================================
# pathway_enrichment_go_themes.R
#
# Group redundant GO terms into themes. The enriched GO terms of every
# contrast (pathway_enrichment_ora.R Rds files) are pooled, so a theme is the
# same set of terms in every contrast and the contrasts can be compared by
# theme. Per ontology, terms are clustered by semantic similarity, and each
# theme is named after its term with the smallest adjusted p-value.
#
# Input:
#   - pathway_enrichment_ora.R Rds files in --enrichment_dir, each a
#     list(terms, stats, universe, settings); only GO terms are used
#
# Method:
#   GOSemSim (doi: 10.1093/bioinformatics/btq064) termSim() between all
#   pairs of terms of one ontology, with --measure (default Rel, as rrvgo).
#   Complete-linkage clustering of 1 - similarity, cut at --threshold, as
#   rrvgo::reduceSimMatrix() (doi: 10.17912/micropub.biology.000811): every
#   pair of terms in a theme has a similarity >= 1 - threshold. The score of
#   a term is -log10 of its smallest p.adjust over all contrasts and types.
#   Terms missing from GO.db (obsolete IDs) get similarity 0, so they form a
#   theme of their own.
#
# Usage:
# > Rscript go_themes.R --enrichment_dir enrichment
# > Rscript go_themes.R --enrichment_dir . --threshold 0.9 --measure Wang
#
# Container:
# > docker://ghcr.io/andreyhgl/enrichment-env:c237546
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
  prefix = "go_themes",
  enrichment_dir = "analysis/pathway_enrichment/ora",
  orgdb = "",
  measure = "Wang", # use "Rel"
  threshold = 0.7
)

opt <- parse_args(defaults)

cat("\n > Arguments:\n---------------------------------\n")
cat(sprintf(" > %-14s : %s\n", names(opt), unlist(opt)), sep = "")
cat(strrep("-", 33), "\n")


#---- Environment -------------------------------------------------------------

suppressPackageStartupMessages({
  library(GOSemSim)
})


#---- Settings ----------------------------------------------------------------

known_measures <- c("Rel", "Resnik", "Lin", "Jiang", "Wang")
if (!opt$measure %in% known_measures) {
  stop("--measure must be one of: ", toString(known_measures))
}
if (opt$threshold <= 0 || opt$threshold >= 1) {
  stop("--threshold must be between 0 and 1, not ", opt$threshold)
}

# Information-content measures need the GO annotation counts of the OrgDb
needs_ic <- opt$measure != "Wang"

ontologies <- c("BP", "MF", "CC")

# Output files
out_rds <- paste0(opt$prefix, ".Rds")
out_table <- paste0(opt$prefix, ".tsv")


#---- Functions ---------------------------------------------------------------

print_section <- function(title) {
  cat("\n > ", title, "\n", strrep("-", 33), "\n", sep = "")
}

package_version <- function(pkg) {
  tryCatch(as.character(packageVersion(pkg)), error = function(e) NA)
}

# Every pathway_enrichment_ora.R Rds in a folder; other Rds files are skipped
read_enrichment <- function(dir) {
  files <- list.files(dir, pattern = "\\.Rds$", full.names = TRUE)
  rds <- lapply(files, readRDS)

  is_enr <- vapply(rds, function(x) {
    is.list(x) && all(c("terms", "settings") %in% names(x))
  }, logical(1))
  if (!any(is_enr)) stop("No enrichment Rds files in: ", dir)

  rds[is_enr]
}

# Enriched GO terms of all contrasts, one row per term with its best padj
best_go_terms <- function(enr) {
  cols <- c("database", "ONTOLOGY", "ID", "Description", "p.adjust")
  terms <- do.call(rbind, lapply(enr, function(e) e$terms[, cols]))
  terms <- terms[terms$database == "GO", ]
  if (nrow(terms) == 0) return(terms[, -1])

  key <- paste(terms$ONTOLOGY, terms$ID)
  best <- terms[order(terms$p.adjust), ]
  best <- best[!duplicated(paste(best$ONTOLOGY, best$ID)), ]
  best$n_hits <- as.integer(table(key)[paste(best$ONTOLOGY, best$ID)])

  data.frame(
    ONTOLOGY = best$ONTOLOGY,
    ID = best$ID,
    Description = best$Description,
    best_padj = best$p.adjust,
    score = -log10(best$p.adjust),
    n_hits = best$n_hits
  )
}

# Pairwise similarity of GO terms; terms unknown to GO.db get 0, not NA
similarity_matrix <- function(ids, sem) {
  sim <- termSim(ids, ids, semData = sem, method = opt$measure)

  full <- matrix(0, length(ids), length(ids), dimnames = list(ids, ids))
  r <- intersect(ids, rownames(sim))
  c <- intersect(ids, colnames(sim))
  full[r, c] <- sim[r, c]
  full[is.na(full)] <- 0

  full <- pmax(full, t(full))
  diag(full) <- 1
  full
}

# Themes of one ontology: complete-linkage clusters of 1 - similarity,
# named after the member with the highest score
cluster_terms <- function(terms, ont) {
  terms <- terms[order(-terms$score, -terms$n_hits), ]

  if (nrow(terms) == 1) {
    cluster <- 1L
    sim <- matrix(1, 1, 1, dimnames = list(terms$ID, terms$ID))
  } else {
    sem <- godata(
      annoDb = opt$orgdb,
      ont = ont,
      computeIC = needs_ic
    )
    sim <- similarity_matrix(terms$ID, sem)
    tree <- hclust(as.dist(1 - sim), method = "complete")
    cluster <- cutree(tree, h = opt$threshold)[terms$ID]
  }

  # terms are sorted by score, so the first of each cluster is its best
  rep_id <- tapply(terms$ID, cluster, `[`, 1)[as.character(cluster)]

  # theme ids in order of the representative's score: BP-001, BP-002, ...
  rank <- match(rep_id, unique(rep_id))
  width <- nchar(max(rank))

  data.frame(
    ONTOLOGY = ont,
    theme_id = sprintf(paste0("%s-%0", width, "d"), ont, rank),
    theme = terms$Description[match(rep_id, terms$ID)],
    theme_term = unname(rep_id),
    theme_size = as.integer(table(rank)[as.character(rank)]),
    ID = terms$ID,
    Description = terms$Description,
    best_padj = terms$best_padj,
    similarity = round(sim[cbind(terms$ID, unname(rep_id))], 3),
    n_hits = terms$n_hits
  )
}


#---- Import datasets ---------------------------------------------------------

enr <- read_enrichment(opt$enrichment_dir)


#---- Code --------------------------------------------------------------------

if (opt$orgdb == "") opt$orgdb <- enr[[1]]$settings$orgdb
suppressPackageStartupMessages(library(opt$orgdb, character.only = TRUE))

go <- best_go_terms(enr)

print_section("Settings")
cat(sprintf(" > %-12s : %d\n", "Rds files", length(enr)))
cat(sprintf(" > %-12s : %s\n", "OrgDb", opt$orgdb))
cat(sprintf(" > %-12s : %s\n", "Measure", opt$measure))
cat(sprintf(
  " > %-12s : %s (similarity >= %s within a theme)\n",
  "Threshold",
  opt$threshold,
  1 - opt$threshold
))

print_section("GO themes")

themes_list <- lapply(ontologies, function(ont) {
  terms <- go[go$ONTOLOGY == ont, ]
  if (nrow(terms) == 0) {
    cat(sprintf(" > %s : no enriched terms\n", ont))
    return(NULL)
  }

  res <- cluster_terms(terms, ont)
  cat(sprintf(
    " > %s : %d terms in %d themes\n",
    ont,
    nrow(res),
    length(unique(res$theme_id))
  ))
  res
})

themes <- do.call(rbind, themes_list)
if (is.null(themes)) {
  themes <- data.frame(
    ONTOLOGY = character(), theme_id = character(), theme = character(),
    theme_term = character(), theme_size = integer(), ID = character(),
    Description = character(), best_padj = numeric(),
    similarity = numeric(), n_hits = integer()
  )
}
rownames(themes) <- NULL

settings <- list(
  orgdb = opt$orgdb,
  measure = opt$measure,
  threshold = opt$threshold,
  versions = c(
    GOSemSim = package_version("GOSemSim"),
    GO.db = package_version("GO.db")
  )
)


#---- Save files --------------------------------------------------------------

# Rds for the report: list(themes, settings)
saveRDS(list(themes = themes, settings = settings), file = out_rds)

write.table(
  themes,
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