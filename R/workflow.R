#' @title Run Complete scMetaLink Analysis
#' @description One-step function to run the complete metabolite-mediated cell
#'   communication analysis pipeline
#'
#' @param expression_data A matrix or sparse matrix of normalized expression
#' @param cell_meta A data.frame with cell metadata
#' @param cell_type_column Character. Column name for cell type annotations
#' @param method Character. Expression scoring method: "mean", "proportion", "combined"
#' @param min_cells Integer. Minimum cells per cell type
#' @param min_pct Numeric. Minimum percentage of expressing cells (0-1)
#' @param n_permutations Integer. Number of permutations for significance testing
#' @param pvalue_threshold Numeric. P-value threshold for significance
#' @param n_cores Integer. Number of cores for parallel computing
#' @param verbose Logical. Print progress messages
#'
#' @return A scMetaLink object with all analysis completed
#' @export
#'
#' @details
#' This function runs the complete scMetaLink pipeline:
#' \enumerate{
#'   \item Create scMetaLink object from expression data
#'   \item Infer metabolite production potential (MPP)
#'   \item Infer metabolite sensing capability (MSC)
#'   \item Compute cell-cell communication scores
#'   \item Perform permutation-based significance testing
#'   \item Filter significant interactions
#'   \item Aggregate by pathway
#' }
#'
#' @examples
#' \donttest{
#' # Run complete analysis
#' result <- runScMetaLink(
#'   expression_data = expr_matrix,
#'   cell_meta = cell_metadata,
#'   cell_type_column = "cell_type",
#'   n_permutations = 1000
#' )
#'
#' # View significant interactions
#' sig <- getSignificantInteractions(result)
#' head(sig)
#' }
runScMetaLink <- function(expression_data,
                          cell_meta,
                          cell_type_column = "cell_type",
                          method = "combined",
                          min_cells = 10,
                          min_pct = 0.1,
                          n_permutations = 1000,
                          pvalue_threshold = 0.05,
                          n_cores = 1,
                          verbose = TRUE) {
  # Parameter validation
  if (!method %in% c("mean", "proportion", "combined")) {
    stop("method must be one of: 'mean', 'proportion', 'combined'")
  }
  if (min_pct < 0 || min_pct > 1) {
    stop("min_pct must be between 0 and 1")
  }
  if (pvalue_threshold <= 0 || pvalue_threshold > 1) {
    stop("pvalue_threshold must be between 0 and 1")
  }

  # Step 1: Create object
  if (verbose) message("=== Step 1/6: Creating scMetaLink object ===")
  object <- createScMetaLink(
    expression_data = expression_data,
    cell_meta = cell_meta,
    cell_type_column = cell_type_column,
    min_cells = min_cells
  )

  # Step 2: Infer production
  if (verbose) message("\n=== Step 2/6: Inferring metabolite production ===")
  object <- inferProduction(
    object = object,
    method = method,
    min_pct = min_pct,
    verbose = verbose
  )

  # Step 3: Infer sensing
  if (verbose) message("\n=== Step 3/6: Inferring metabolite sensing ===")
  object <- inferSensing(
    object = object,
    method = method,
    min_pct = min_pct,
    verbose = verbose
  )

  # Step 4: Compute communication
  if (verbose) message("\n=== Step 4/6: Computing cell communication ===")
  object <- computeCommunication(
    object = object,
    n_permutations = n_permutations,
    n_cores = n_cores,
    verbose = verbose
  )

  # Step 5: Filter significant interactions
  if (verbose) message("\n=== Step 5/6: Filtering significant interactions ===")
  object <- filterSignificantInteractions(
    object = object,
    pvalue_threshold = pvalue_threshold
  )

  # Step 6: Aggregate by pathway (only if significant interactions exist)
  if (nrow(object@significant_interactions) > 0) {
    if (verbose) message("\n=== Step 6/6: Aggregating by pathway ===")
    object <- aggregateByPathway(object)
  } else {
    if (verbose) message("\n=== Step 6/6: Skipping pathway aggregation (no significant interactions) ===")
  }

  if (verbose) {
    message("\n=== Analysis Complete ===")
    message(sprintf("Significant interactions: %d", nrow(object@significant_interactions)))
    if (nrow(object@pathway_aggregated) > 0) {
      message(sprintf("Pathways involved: %d", length(unique(object@pathway_aggregated$pathway))))
    }
  }

  object
}

#' @title Run scMetaLink from Seurat Object
#' @description Convenience function to run analysis from Seurat object
#'
#' @param seurat_obj A Seurat object
#' @param cell_type_column Character. Column in meta.data for cell type
#' @param assay Character. Assay to use
#' @param slot Character. Slot to use ("data" or "counts")
#' @param ... Additional arguments passed to runScMetaLink
#'
#' @return A scMetaLink object
#' @export
runScMetaLinkSeurat <- function(seurat_obj,
                                cell_type_column = "cell_type",
                                assay = "RNA",
                                slot = "data",
                                ...) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required")
  }

  # Extract data
  if (slot == "data") {
    expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = "data")
  } else {
    expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = "counts")
  }

  cell_meta <- seurat_obj@meta.data

  runScMetaLink(
    expression_data = expr_data,
    cell_meta = cell_meta,
    cell_type_column = cell_type_column,
    ...
  )
}

#' @title Compare Two Conditions
#' @description Compare metabolite-mediated communication between two conditions
#'
#' @param object1 scMetaLink object for condition 1
#' @param object2 scMetaLink object for condition 2
#' @param condition_names Character vector of length 2 for condition names
#' @param method Character. Comparison method: "difference", "ratio", or "log2fc"
#'
#' @return data.frame with differential communication results
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare tumor vs normal (requires two scMetaLink objects)
#' diff <- compareCommunication(tumor_obj, normal_obj,
#'   condition_names = c("Tumor", "Normal")
#' )
#' }
compareCommunication <- function(object1,
                                 object2,
                                 condition_names = c("Condition1", "Condition2"),
                                 method = "log2fc") {
  # Parameter validation
  if (!inherits(object1, "scMetaLink") || !inherits(object2, "scMetaLink")) {
    stop("Both objects must be scMetaLink objects")
  }
  if (length(condition_names) != 2) {
    stop("condition_names must be a character vector of length 2")
  }
  if (!method %in% c("difference", "ratio", "log2fc")) {
    stop("method must be one of: 'difference', 'ratio', 'log2fc'")
  }

  sig1 <- object1@significant_interactions
  sig2 <- object2@significant_interactions

  if (nrow(sig1) == 0 && nrow(sig2) == 0) {
    stop("Both objects have no significant interactions")
  }

  # Use a separator that is unlikely to appear in cell type names
  sep <- ":::"

  # Create keys using safe separator
  sig1$key <- paste(sig1$sender, sig1$receiver, sig1$metabolite_id, sep = sep)
  sig2$key <- paste(sig2$sender, sig2$receiver, sig2$metabolite_id, sep = sep)

  # Get all unique interactions
  all_keys <- unique(c(sig1$key, sig2$key))

  # Build comparison table
  comparison <- data.frame(
    key = all_keys,
    stringsAsFactors = FALSE
  )

  # Add scores from both conditions
  comparison$score1 <- sig1$communication_score[match(comparison$key, sig1$key)]
  comparison$score2 <- sig2$communication_score[match(comparison$key, sig2$key)]

  # Replace NA with 0
  comparison$score1[is.na(comparison$score1)] <- 0
  comparison$score2[is.na(comparison$score2)] <- 0

  # Calculate difference/ratio
  if (method == "difference") {
    comparison$change <- comparison$score1 - comparison$score2
  } else if (method == "ratio") {
    comparison$change <- (comparison$score1 + 0.01) / (comparison$score2 + 0.01)
  } else if (method == "log2fc") {
    comparison$change <- log2((comparison$score1 + 0.01) / (comparison$score2 + 0.01))
  }

  # Parse key back to components using the same separator
  key_parts <- strsplit(comparison$key, sep, fixed = TRUE)
  comparison$sender <- sapply(key_parts, `[`, 1)
  comparison$receiver <- sapply(key_parts, `[`, 2)
  comparison$metabolite_id <- sapply(key_parts, `[`, 3)

  # Add metabolite names
  db <- .load_metalinksdb()
  comparison <- merge(comparison, db$metabolites[, c("hmdb", "metabolite")],
    by.x = "metabolite_id", by.y = "hmdb", all.x = TRUE
  )

  # Rename columns
  names(comparison)[names(comparison) == "score1"] <- paste0("score_", condition_names[1])
  names(comparison)[names(comparison) == "score2"] <- paste0("score_", condition_names[2])

  # Sort by absolute change
  comparison <- comparison[order(-abs(comparison$change)), ]
  comparison <- comparison[, c(
    "sender", "receiver", "metabolite_id", "metabolite",
    paste0("score_", condition_names), "change"
  )]
  rownames(comparison) <- NULL

  comparison
}

#' @title Identify Cell Type Specific Metabolites
#' @description Find metabolites specifically produced or sensed by cell types
#'
#' @param object scMetaLink object
#' @param type Character. "production" or "sensing"
#' @param specificity_threshold Numeric. Z-score threshold for specificity
#'
#' @return data.frame with cell type specific metabolites
#' @export
identifyCellTypeSpecificMetabolites <- function(object,
                                                type = "production",
                                                specificity_threshold = 1.5) {
  if (!type %in% c("production", "sensing")) {
    stop("type must be 'production' or 'sensing'")
  }
  if (specificity_threshold <= 0) {
    stop("specificity_threshold must be positive")
  }

  if (type == "production") {
    scores <- object@production_scores
  } else {
    scores <- object@sensing_scores
  }

  if (is.null(scores)) {
    stop(sprintf("%s scores not calculated", type))
  }

  # Calculate z-scores
  z_scores <- t(apply(scores, 1, function(x) {
    s <- sd(x, na.rm = TRUE)
    if (s == 0) {
      return(rep(0, length(x)))
    }
    (x - mean(x, na.rm = TRUE)) / s
  }))

  # Find specific metabolites
  results <- data.frame()

  for (ct in colnames(z_scores)) {
    ct_specific <- rownames(z_scores)[z_scores[, ct] > specificity_threshold]

    if (length(ct_specific) > 0) {
      for (met in ct_specific) {
        results <- rbind(results, data.frame(
          cell_type = ct,
          metabolite_id = met,
          z_score = z_scores[met, ct],
          raw_score = scores[met, ct],
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  if (nrow(results) == 0) {
    message("No cell type specific metabolites found")
    return(data.frame())
  }

  # Add metabolite names
  db <- object@database
  results <- merge(results, db$metabolites[, c("hmdb", "metabolite", "metabolite_subclass")],
    by.x = "metabolite_id", by.y = "hmdb", all.x = TRUE
  )

  results <- results[order(-results$z_score), ]
  rownames(results) <- NULL

  results
}

#' @title Get Summary Statistics
#' @description Get summary statistics of scMetaLink analysis
#'
#' @param object scMetaLink object
#'
#' @return list with summary statistics
#' @export
getSummaryStats <- function(object) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  stats <- list()

  # Basic info
  stats$n_genes <- nrow(object@expression_data)
  stats$n_cells <- ncol(object@expression_data)
  stats$n_cell_types <- length(unique(object@cell_meta[[object@cell_type_column]]))
  stats$cell_types <- unique(object@cell_meta[[object@cell_type_column]])

  # Production/sensing
  if (!is.null(object@production_scores)) {
    stats$n_metabolites_production <- nrow(object@production_scores)
  }
  if (!is.null(object@sensing_scores)) {
    stats$n_metabolites_sensing <- nrow(object@sensing_scores)
  }

  # Communication
  if (!is.null(object@communication_scores)) {
    stats$n_metabolites_communication <- dim(object@communication_scores)[3]
    stats$total_communication_pairs <- sum(object@communication_scores > 0)
  }

  # Significant interactions
  if (nrow(object@significant_interactions) > 0) {
    stats$n_significant <- nrow(object@significant_interactions)
    stats$n_unique_metabolites <- length(unique(object@significant_interactions$metabolite_id))
    stats$n_unique_sender_receiver <- nrow(unique(object@significant_interactions[, c("sender", "receiver")]))

    # Top metabolites
    met_counts <- table(object@significant_interactions$metabolite_name)
    stats$top_metabolites <- names(sort(met_counts, decreasing = TRUE))[1:min(10, length(met_counts))]

    # Top cell type pairs
    pair_counts <- table(paste(object@significant_interactions$sender,
      object@significant_interactions$receiver,
      sep = " -> "
    ))
    stats$top_pairs <- names(sort(pair_counts, decreasing = TRUE))[1:min(10, length(pair_counts))]
  }

  # Pathway
  if (nrow(object@pathway_aggregated) > 0) {
    stats$n_pathways <- length(unique(object@pathway_aggregated$pathway))

    # Top pathways
    pw_scores <- aggregate(communication_score ~ pathway, data = object@pathway_aggregated, FUN = sum)
    pw_scores <- pw_scores[order(-pw_scores$communication_score), ]
    stats$top_pathways <- pw_scores$pathway[1:min(10, nrow(pw_scores))]
  }

  stats
}
