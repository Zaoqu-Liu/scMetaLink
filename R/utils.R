#' @title scMetaLink Utility Functions
#' @description Internal utility functions for scMetaLink package
#' @name utils
#' @keywords internal

#' Check Package Dependencies
#'
#' @param packages Character vector of package names
#' @param type Character. "error" to stop, "warning" to warn
#' @keywords internal
.check_dependencies <- function(packages, type = "error") {
  missing <- packages[!sapply(packages, requireNamespace, quietly = TRUE)]

  if (length(missing) > 0) {
    msg <- sprintf("Missing packages: %s. Please install them.",
                   paste(missing, collapse = ", "))
    if (type == "error") {
      stop(msg, call. = FALSE)
    } else {
      warning(msg, call. = FALSE)
    }
  }

  invisible(TRUE)
}

#' Safe Log Transform
#'
#' @param x Numeric vector
#' @param base Log base
#' @param offset Offset to avoid log(0)
#' @keywords internal
.safe_log <- function(x, base = 2, offset = 1) {
  log(x + offset, base = base)
}

#' Scale Matrix to 0-1
#'
#' @param mat Matrix
#' @param by_row Logical. Scale by row
#' @keywords internal
.scale_01 <- function(mat, by_row = TRUE) {
  if (by_row) {
    t(apply(mat, 1, function(x) {
      rng <- range(x, na.rm = TRUE)
      if (rng[1] == rng[2]) return(rep(0.5, length(x)))
      (x - rng[1]) / (rng[2] - rng[1])
    }))
  } else {
    apply(mat, 2, function(x) {
      rng <- range(x, na.rm = TRUE)
      if (rng[1] == rng[2]) return(rep(0.5, length(x)))
      (x - rng[1]) / (rng[2] - rng[1])
    })
  }
}

#' Z-score Normalization
#'
#' @param mat Matrix
#' @param by_row Logical. Normalize by row
#' @keywords internal
.zscore <- function(mat, by_row = TRUE) {
  if (by_row) {
    t(apply(mat, 1, function(x) {
      s <- sd(x, na.rm = TRUE)
      if (s == 0) return(rep(0, length(x)))
      (x - mean(x, na.rm = TRUE)) / s
    }))
  } else {
    apply(mat, 2, function(x) {
      s <- sd(x, na.rm = TRUE)
      if (s == 0) return(rep(0, length(x)))
      (x - mean(x, na.rm = TRUE)) / s
    })
  }
}

#' Get Color Palette
#'
#' @param n Number of colors
#' @param palette Character. Palette name
#' @keywords internal
.get_colors <- function(n, palette = "default") {
  if (palette == "default") {
    if (n <= 8) {
      RColorBrewer::brewer.pal(max(3, n), "Set2")[1:n]
    } else if (n <= 12) {
      RColorBrewer::brewer.pal(n, "Set3")
    } else {
      viridis::viridis(n)
    }
  } else if (palette == "paired") {
    if (n <= 12) {
      RColorBrewer::brewer.pal(max(3, n), "Paired")[1:n]
    } else {
      c(RColorBrewer::brewer.pal(12, "Paired"),
        viridis::viridis(n - 12))
    }
  } else {
    viridis::viridis(n, option = palette)
  }
}

#' Print Progress
#'
#' @param current Current iteration
#' @param total Total iterations
#' @param prefix Prefix message
#' @keywords internal
.print_progress <- function(current, total, prefix = "Progress") {
  pct <- round(current / total * 100)
  cat(sprintf("\r%s: %d/%d (%d%%)", prefix, current, total, pct))
  if (current == total) cat("\n")
}

#' Validate Input Data
#'
#' @param expr_data Expression matrix
#' @param cell_meta Cell metadata
#' @keywords internal
.validate_input <- function(expr_data, cell_meta) {
  errors <- character()

  # Check expression data
  if (!is.matrix(expr_data) && !inherits(expr_data, "dgCMatrix")) {
    errors <- c(errors, "expression_data must be a matrix or sparse dgCMatrix")
  }

  if (is.null(rownames(expr_data))) {
    errors <- c(errors, "expression_data must have row names (gene symbols)")
  }

  if (is.null(colnames(expr_data))) {
    errors <- c(errors, "expression_data must have column names (cell IDs)")
  }

  # Check cell_meta
  if (!is.data.frame(cell_meta)) {
    errors <- c(errors, "cell_meta must be a data.frame")
  }

  if (is.null(rownames(cell_meta))) {
    errors <- c(errors, "cell_meta must have row names (cell IDs)")
  }

  # Check overlap
  if (length(intersect(colnames(expr_data), rownames(cell_meta))) == 0) {
    errors <- c(errors, "No matching cell IDs between expression_data and cell_meta")
  }

  if (length(errors) > 0) {
    stop(paste(errors, collapse = "\n"), call. = FALSE)
  }

  invisible(TRUE)
}

#' Export Results to CSV
#'
#' @param object scMetaLink object
#' @param output_dir Character. Output directory
#' @param prefix Character. File prefix
#'
#' @return Invisibly returns file paths
#' @export
exportResults <- function(object, output_dir = ".", prefix = "scMetaLink") {

  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  files <- character()

  # Export production scores
  if (!is.null(object@production_scores)) {
    f <- file.path(output_dir, paste0(prefix, "_production_scores.csv"))
    write.csv(object@production_scores, f)
    files <- c(files, f)
  }

  # Export sensing scores
  if (!is.null(object@sensing_scores)) {
    f <- file.path(output_dir, paste0(prefix, "_sensing_scores.csv"))
    write.csv(object@sensing_scores, f)
    files <- c(files, f)
  }

  # Export significant interactions
  if (nrow(object@significant_interactions) > 0) {
    f <- file.path(output_dir, paste0(prefix, "_significant_interactions.csv"))
    write.csv(object@significant_interactions, f, row.names = FALSE)
    files <- c(files, f)
  }

  # Export pathway aggregated
  if (nrow(object@pathway_aggregated) > 0) {
    f <- file.path(output_dir, paste0(prefix, "_pathway_aggregated.csv"))
    write.csv(object@pathway_aggregated, f, row.names = FALSE)
    files <- c(files, f)
  }

  message(sprintf("Exported %d files to %s", length(files), output_dir))
  invisible(files)
}

#' Save scMetaLink Object
#'
#' @param object scMetaLink object
#' @param file Character. File path (RDS format)
#'
#' @return Invisibly returns file path
#' @export
saveScMetaLink <- function(object, file) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (!grepl("\\.rds$", file, ignore.case = TRUE)) {
    file <- paste0(file, ".rds")
  }
  saveRDS(object, file)
  message(sprintf("Saved scMetaLink object to %s", file))
  invisible(file)
}

#' Load scMetaLink Object
#'
#' @param file Character. File path
#'
#' @return scMetaLink object
#' @export
loadScMetaLink <- function(file) {
  if (!file.exists(file)) {
    stop("File not found: ", file)
  }
  obj <- readRDS(file)
  if (!inherits(obj, "scMetaLink")) {
    warning("Loaded object is not a scMetaLink object")
  }
  obj
}

#' Get Database Information
#'
#' @return data.frame with database statistics
#' @export
getDatabaseInfo <- function() {
  db <- .load_metalinksdb()

  info <- data.frame(
    Component = c(
      "Total metabolites",
      "Total proteins/genes",
      "Total interactions",
      "Ligand-receptor interactions",
      "Produce-degrade interactions",
      "Signaling metabolites (with receptors)",
      "Metabolic metabolites (with enzymes)",
      "Pathways",
      "Extracellular metabolites",
      "Disease associations"
    ),
    Count = c(
      nrow(db$metabolites),
      nrow(db$proteins),
      nrow(db$edges),
      sum(db$edges$type == "lr"),
      sum(db$edges$type == "pd"),
      length(unique(db$edges$hmdb[db$edges$type == "lr"])),
      length(unique(db$edges$hmdb[db$edges$type == "pd"])),
      length(unique(db$pathway$pathway)),
      length(unique(db$cell_location$hmdb[db$cell_location$cell_location == "Extracellular"])),
      nrow(db$disease)
    ),
    stringsAsFactors = FALSE
  )

  info
}

#' Search Metabolite in Database
#'
#' @param query Character. Search query (name or HMDB ID)
#' @param exact Logical. Exact match only
#'
#' @return data.frame with matching metabolites
#' @export
searchMetabolite <- function(query, exact = FALSE) {
  if (missing(query) || is.null(query) || query == "") {
    stop("query must be provided")
  }
  
  db <- .load_metalinksdb()

  if (exact) {
    matches <- db$metabolites[db$metabolites$hmdb == query |
                              db$metabolites$metabolite == query, ]
  } else {
    matches <- db$metabolites[grepl(query, db$metabolites$hmdb, ignore.case = TRUE) |
                              grepl(query, db$metabolites$metabolite, ignore.case = TRUE), ]
  }

  if (nrow(matches) == 0) {
    message("No metabolites found matching: ", query)
    return(data.frame())
  }

  # Add interaction counts
  matches$n_receptors <- sapply(matches$hmdb, function(h) {
    sum(db$edges$hmdb == h & db$edges$type == "lr")
  })
  matches$n_enzymes <- sapply(matches$hmdb, function(h) {
    sum(db$edges$hmdb == h & db$edges$type == "pd")
  })

  matches
}

#' Search Gene in Database
#'
#' @param query Character. Gene symbol
#'
#' @return data.frame with gene information and associated metabolites
#' @export
searchGene <- function(query) {
  if (missing(query) || is.null(query) || query == "") {
    stop("query must be provided")
  }
  
  db <- .load_metalinksdb()

  # Find gene
  gene_info <- db$proteins[grepl(query, db$proteins$gene_symbol, ignore.case = TRUE), ]

  if (nrow(gene_info) == 0) {
    message("Gene not found: ", query)
    return(data.frame())
  }

  # Get associated metabolites
  gene_edges <- db$edges[db$edges$uniprot %in% gene_info$uniprot, ]
  gene_edges <- merge(gene_edges, gene_info[, c("uniprot", "gene_symbol", "protein_type")],
                      by = "uniprot")
  gene_edges <- merge(gene_edges, db$metabolites[, c("hmdb", "metabolite")],
                      by = "hmdb")

  cols_to_return <- c("gene_symbol", "protein_type", "metabolite", "hmdb", "type", "mor", "combined_score")
  cols_to_return <- cols_to_return[cols_to_return %in% names(gene_edges)]
  gene_edges[, cols_to_return]
}

#' Format Number with Commas
#' @keywords internal
.format_number <- function(x) {
  format(x, big.mark = ",", scientific = FALSE)
}

#' Check if Running in Interactive Mode
#' @keywords internal
.is_interactive <- function() {
  interactive()
}
