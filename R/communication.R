#' @title Compute Metabolite-Mediated Cell Communication
#' @description Calculate communication scores between all cell type pairs
#'   mediated by metabolites. Communication strength represents the potential
#'   for signal transmission from sender to receiver cells via metabolites.
#'
#' @param object A scMetaLink object with production and sensing scores
#' @param method Character. Communication score method: "geometric" (default), "product", "harmonic"
#' @param min_production Numeric. Minimum production score threshold (0-1)
#' @param min_sensing Numeric. Minimum sensing score threshold (0-1)
#' @param population.size Logical. Whether to weight communication by cell type
#'   population sizes. When TRUE, communication strength is scaled by the relative
#'   abundance of sender and receiver cell types, reflecting the biological reality
#'   that larger populations contribute more to overall tissue signaling.
#' @param n_permutations Integer. Number of permutations for significance testing (0 to skip)
#' @param n_cores Integer. Number of cores for parallel computing
#' @param seed Integer. Random seed for reproducibility
#' @param verbose Logical. Print progress messages
#'
#' @details
#' The communication score combines production and sensing capabilities:
#' \itemize{
#'   \item geometric: sqrt(production * sensing) - balanced measure
#'   \item product: production * sensing - emphasizes strong bilateral signals
#'   \item harmonic: 2*production*sensing/(production+sensing) - penalizes imbalance
#' }
#'
#' When \code{population.size = TRUE}, scores are multiplied by:
#' sqrt(n_sender/n_total * n_receiver/n_total), which accounts for the
#' relative contribution of each cell type to tissue-level communication.
#'
#' @return Updated scMetaLink object with communication_scores and p-values
#' @export
computeCommunication <- function(object,
                                 method = "geometric",
                                 min_production = 0.1,
                                 min_sensing = 0.1,
                                 population.size = FALSE,
                                 n_permutations = 100,
                                 n_cores = 1,
                                 seed = 42,
                                 verbose = TRUE) {
  if (!inherits(object, "scMetaLink")) stop("object must be a scMetaLink object")
  if (is.null(object@production_scores)) stop("Run inferProduction() first")
  if (is.null(object@sensing_scores)) stop("Run inferSensing() first")

  # Parameter validation
  if (!method %in% c("geometric", "product", "harmonic")) {
    stop("method must be one of: 'geometric', 'product', 'harmonic'")
  }
  if (min_production < 0 || min_production > 1) {
    stop("min_production must be between 0 and 1")
  }
  if (min_sensing < 0 || min_sensing > 1) {
    stop("min_sensing must be between 0 and 1")
  }
  if (n_permutations < 0) {
    stop("n_permutations must be non-negative")
  }

  set.seed(seed)

  prod_scores <- object@production_scores
  sens_scores <- object@sensing_scores
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column

  # Get common metabolites
  common_mets <- intersect(rownames(prod_scores), rownames(sens_scores))
  if (length(common_mets) == 0) stop("No common metabolites between production and sensing")

  prod_scores <- prod_scores[common_mets, , drop = FALSE]
  sens_scores <- sens_scores[common_mets, , drop = FALSE]

  cell_types <- colnames(prod_scores)
  n_types <- length(cell_types)
  n_mets <- length(common_mets)

  if (verbose) message(sprintf("Computing communication for %d metabolites, %d cell types...", n_mets, n_types))

  # Calculate population size weights if requested
  pop_weights <- NULL
  if (population.size) {
    if (verbose) message("  Calculating population size weights...")
    cell_counts <- table(cell_meta[[cell_type_col]])
    total_cells <- sum(cell_counts)

    # Create weight matrix: sqrt(p_sender * p_receiver)
    pop_fractions <- as.numeric(cell_counts[cell_types]) / total_cells
    pop_weights <- sqrt(outer(pop_fractions, pop_fractions))
    dimnames(pop_weights) <- list(cell_types, cell_types)
  }

  # Calculate communication scores using vectorized approach
  if (verbose) message("  Calculating communication scores...")

  comm_scores <- .compute_communication_scores(
    prod_scores = prod_scores,
    sens_scores = sens_scores,
    method = method,
    min_production = min_production,
    min_sensing = min_sensing,
    pop_weights = pop_weights
  )

  dimnames(comm_scores) <- list(sender = cell_types, receiver = cell_types, metabolite = common_mets)

  object@communication_scores <- comm_scores

  # Permutation test
  if (n_permutations > 0) {
    if (verbose) message(sprintf("  Running %d permutations (n_cores=%d)...", n_permutations, n_cores))

    pvalues <- .run_permutation_test_consistent(
      object = object,
      comm_scores = comm_scores,
      common_mets = common_mets,
      n_permutations = n_permutations,
      n_cores = n_cores,
      method = method,
      min_production = min_production,
      min_sensing = min_sensing,
      population.size = population.size,
      verbose = verbose
    )

    object@communication_pvalues <- pvalues
  }

  object@parameters$communication <- list(
    method = method,
    min_production = min_production,
    min_sensing = min_sensing,
    population.size = population.size,
    n_permutations = n_permutations,
    seed = seed
  )

  if (verbose) message("Done!")
  object
}

#' Compute Communication Scores (Vectorized)
#'
#' @param prod_scores Production score matrix (metabolites x cell_types)
#' @param sens_scores Sensing score matrix (metabolites x cell_types)
#' @param method Communication method
#' @param min_production Minimum production threshold
#' @param min_sensing Minimum sensing threshold
#' @param pop_weights Optional population weight matrix (cell_types x cell_types)
#' @return 3D array of communication scores (sender x receiver x metabolite)
#' @keywords internal
.compute_communication_scores <- function(prod_scores, sens_scores, method,
                                          min_production, min_sensing,
                                          pop_weights = NULL) {
  n_mets <- nrow(prod_scores)
  n_types <- ncol(prod_scores)

  comm_scores <- array(0, dim = c(n_types, n_types, n_mets))

  for (i in seq_len(n_mets)) {
    p <- prod_scores[i, ]
    s <- sens_scores[i, ]
    p[p < min_production] <- 0
    s[s < min_sensing] <- 0

    if (method == "geometric") {
      score_mat <- sqrt(outer(p, s))
    } else if (method == "product") {
      score_mat <- outer(p, s)
    } else if (method == "harmonic") {
      ps <- outer(p, s)
      sum_ps <- outer(p, rep(1, n_types)) + outer(rep(1, n_types), s)
      score_mat <- 2 * ps / (sum_ps + 1e-10)
    }

    # Apply population size weights if provided
    if (!is.null(pop_weights)) {
      score_mat <- score_mat * pop_weights
    }

    comm_scores[, , i] <- score_mat
  }

  comm_scores
}

#' Permutation Test - Consistent with inferProduction/inferSensing
#' @description This function replicates EXACTLY the same computation as
#'   inferProduction and inferSensing to ensure statistical validity
#' @keywords internal
.run_permutation_test_consistent <- function(object, comm_scores, common_mets,
                                             n_permutations, n_cores, method,
                                             min_production, min_sensing,
                                             population.size = FALSE, verbose) {
  expr_data <- object@expression_data
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column
  db <- object@database
  params <- object@parameters

  cell_types <- dimnames(comm_scores)[[1]]
  n_types <- length(cell_types)
  n_mets <- length(common_mets)
  genes <- rownames(expr_data)

  # ============================================================
  # PRODUCTION: Replicate EXACTLY the same logic as inferProduction
  # ============================================================
  prod_params <- params$production
  consider_degradation <- prod_params$consider_degradation %||% TRUE
  consider_secretion <- prod_params$consider_secretion %||% TRUE
  normalize_production <- prod_params$normalize %||% TRUE
  min_expression <- prod_params$min_expression %||% 0
  scoring_method <- prod_params$method %||% "combined"
  mean_method_prod <- prod_params$mean_method %||% "arithmetic"

  # Get production enzymes (same as inferProduction)
  prod_enzymes <- .get_production_enzymes(db)
  prod_enzymes <- prod_enzymes[!is.na(prod_enzymes$gene_symbol), ]
  prod_enzymes <- prod_enzymes[prod_enzymes$gene_symbol %in% genes, ]
  prod_enzymes <- prod_enzymes[prod_enzymes$hmdb %in% common_mets, ]

  # Build production mapping
  prod_map <- Matrix::sparseMatrix(
    i = match(prod_enzymes$hmdb, common_mets),
    j = match(prod_enzymes$gene_symbol, genes),
    x = 1,
    dims = c(length(common_mets), length(genes)),
    dimnames = list(common_mets, genes)
  )
  prod_enzyme_counts <- Matrix::rowSums(prod_map)
  prod_enzyme_counts[prod_enzyme_counts == 0] <- 1

  # Degradation mapping (same as inferProduction)
  deg_map <- NULL
  deg_enzyme_counts <- NULL
  if (consider_degradation) {
    deg_enzymes <- .get_degradation_enzymes(db)
    deg_enzymes <- deg_enzymes[!is.na(deg_enzymes$gene_symbol), ]
    deg_enzymes <- deg_enzymes[deg_enzymes$gene_symbol %in% genes, ]
    deg_enzymes <- deg_enzymes[deg_enzymes$hmdb %in% common_mets, ]

    if (nrow(deg_enzymes) > 0) {
      deg_map <- Matrix::sparseMatrix(
        i = match(deg_enzymes$hmdb, common_mets),
        j = match(deg_enzymes$gene_symbol, genes),
        x = 1,
        dims = c(length(common_mets), length(genes)),
        dimnames = list(common_mets, genes)
      )
      deg_enzyme_counts <- Matrix::rowSums(deg_map)
      deg_enzyme_counts[deg_enzyme_counts == 0] <- 1
    }
  }

  # Secretion weights (same as inferProduction)
  secretion_weights <- NULL
  if (consider_secretion) {
    extra_mets <- .get_extracellular_metabolites(db)
    secretion_weights <- ifelse(common_mets %in% extra_mets, 1.0, 0.5)
  }

  # ============================================================
  # SENSING: Replicate EXACTLY the same logic as inferSensing
  # ============================================================
  sens_params <- params$sensing
  weight_by_affinity <- sens_params$weight_by_affinity %||% TRUE
  include_transporters <- sens_params$include_transporters %||% TRUE
  normalize_sensing <- sens_params$normalize %||% TRUE
  use_hill <- sens_params$use_hill %||% FALSE
  hill_n <- sens_params$hill_n %||% 1
  hill_Kh <- sens_params$hill_Kh %||% 0.5

  # Get sensing proteins (same as inferSensing)
  sensing_proteins <- .get_sensing_proteins(db)
  sensing_proteins <- sensing_proteins[!is.na(sensing_proteins$gene_symbol), ]

  # Add transporters if requested
  if (include_transporters) {
    in_trans <- .get_transporters(db, direction = "in")
    in_trans <- in_trans[!is.na(in_trans$gene_symbol), ]
    if (nrow(in_trans) > 0) {
      sensing_proteins <- rbind(sensing_proteins, in_trans)
    }
  }

  # Filter to available genes
  sensing_proteins <- sensing_proteins[sensing_proteins$gene_symbol %in% genes, ]
  sensing_proteins <- sensing_proteins[sensing_proteins$hmdb %in% common_mets, ]

  # Remove duplicates, keep highest score
  sensing_proteins <- sensing_proteins[order(-sensing_proteins$combined_score), ]
  sensing_proteins <- sensing_proteins[!duplicated(paste(sensing_proteins$hmdb, sensing_proteins$gene_symbol)), ]

  # Build sensing mapping with affinity weights
  if (weight_by_affinity) {
    sens_weights <- sensing_proteins$combined_score
    sens_weights[is.na(sens_weights)] <- 500
    sens_weights <- sens_weights / 1000
  } else {
    sens_weights <- rep(1, nrow(sensing_proteins))
  }

  sens_map <- Matrix::sparseMatrix(
    i = match(sensing_proteins$hmdb, common_mets),
    j = match(sensing_proteins$gene_symbol, genes),
    x = sens_weights,
    dims = c(length(common_mets), length(genes)),
    dimnames = list(common_mets, genes)
  )
  sens_receptor_counts <- Matrix::rowSums(sens_map > 0)
  sens_receptor_counts[sens_receptor_counts == 0] <- 1

  # ============================================================
  # POPULATION SIZE WEIGHTS
  # ============================================================
  pop_weights <- NULL
  if (population.size) {
    cell_counts <- table(cell_meta[[cell_type_col]])
    total_cells <- sum(cell_counts)
    pop_fractions <- as.numeric(cell_counts[cell_types]) / total_cells
    pop_weights <- sqrt(outer(pop_fractions, pop_fractions))
    dimnames(pop_weights) <- list(cell_types, cell_types)
  }

  # ============================================================
  # PERMUTATION FUNCTION
  # ============================================================
  run_single_perm <- function(perm_id) {
    # Shuffle cell type labels
    shuffled_labels <- sample(cell_meta[[cell_type_col]])
    shuffled_meta <- cell_meta
    shuffled_meta[[cell_type_col]] <- shuffled_labels

    # Calculate expression profiles with appropriate mean_method
    perm_profiles <- .calculate_celltype_expression_fast(
      expr_data = expr_data,
      cell_meta = shuffled_meta,
      cell_type_col = cell_type_col,
      min_expression = min_expression,
      mean_method = mean_method_prod
    )

    # Apply same scoring method as original
    if (scoring_method == "combined") {
      perm_gene_scores <- perm_profiles$mean_expr * perm_profiles$pct_expr
    } else if (scoring_method == "mean") {
      perm_gene_scores <- perm_profiles$mean_expr
    } else {
      perm_gene_scores <- perm_profiles$pct_expr
    }

    # Ensure column order matches cell_types
    perm_gene_scores <- perm_gene_scores[, cell_types, drop = FALSE]

    # ---- PRODUCTION SCORES ----
    perm_prod <- as.matrix(prod_map %*% perm_gene_scores) / prod_enzyme_counts

    # Degradation adjustment
    if (!is.null(deg_map)) {
      deg_scores <- as.matrix(deg_map %*% perm_gene_scores) / deg_enzyme_counts
      perm_prod <- perm_prod - 0.5 * deg_scores
      perm_prod[perm_prod < 0] <- 0
    }

    # Secretion weights
    if (!is.null(secretion_weights)) {
      perm_prod <- perm_prod * secretion_weights
    }

    # Normalize
    if (normalize_production) {
      perm_prod <- .normalize_scores_permutation(perm_prod)
    }

    # ---- SENSING SCORES ----
    perm_sens <- as.matrix(sens_map %*% perm_gene_scores) / sens_receptor_counts

    # Apply Hill transformation if used in original
    if (use_hill) {
      sens_min <- min(perm_sens, na.rm = TRUE)
      sens_max <- max(perm_sens, na.rm = TRUE)
      if (sens_max > sens_min) {
        sens_norm <- (perm_sens - sens_min) / (sens_max - sens_min)
        perm_sens <- .hill_transform(sens_norm, n = hill_n, Kh = hill_Kh)
      }
    }

    # Normalize
    if (normalize_sensing) {
      perm_sens <- .normalize_scores_permutation(perm_sens)
    }

    # ---- COMMUNICATION SCORES ----
    perm_comm <- .compute_communication_scores(
      prod_scores = perm_prod,
      sens_scores = perm_sens,
      method = method,
      min_production = min_production,
      min_sensing = min_sensing,
      pop_weights = pop_weights
    )

    # Return comparison result
    perm_comm >= comm_scores
  }

  # ============================================================
  # RUN PERMUTATIONS
  # ============================================================
  if (n_cores > 1 && requireNamespace("future.apply", quietly = TRUE)) {
    # Parallel execution
    future::plan(future::multisession, workers = n_cores)
    on.exit(future::plan(future::sequential), add = TRUE)

    if (verbose) message("    Using parallel processing...")

    results <- future.apply::future_lapply(
      seq_len(n_permutations),
      run_single_perm,
      future.seed = TRUE
    )

    count_matrix <- Reduce(`+`, results)
  } else {
    # Sequential execution with progress bar
    count_matrix <- array(0, dim = dim(comm_scores))

    if (verbose) pb <- txtProgressBar(min = 0, max = n_permutations, style = 3)

    for (perm in seq_len(n_permutations)) {
      count_matrix <- count_matrix + run_single_perm(perm)
      if (verbose) setTxtProgressBar(pb, perm)
    }

    if (verbose) close(pb)
  }

  # Calculate p-values
  pvalues <- (count_matrix + 1) / (n_permutations + 1)
  dimnames(pvalues) <- dimnames(comm_scores)
  pvalues
}

#' Normalize Scores for Permutation (same logic as inferProduction/inferSensing)
#' @keywords internal
.normalize_scores_permutation <- function(scores) {
  if (nrow(scores) == 0) {
    return(scores)
  }

  # Z-score per metabolite (row)
  row_means <- rowMeans(scores, na.rm = TRUE)
  row_sds <- apply(scores, 1, sd, na.rm = TRUE)
  row_sds[row_sds == 0] <- 1
  scores <- (scores - row_means) / row_sds

  # Scale to 0-1
  min_val <- min(scores, na.rm = TRUE)
  max_val <- max(scores, na.rm = TRUE)
  if (max_val > min_val) {
    scores <- (scores - min_val) / (max_val - min_val)
  }

  scores
}

#' Null-coalescing operator
#' @keywords internal
`%||%` <- function(a, b) if (is.null(a)) b else a

#' Filter Significant Interactions
#' @param object scMetaLink object
#' @param pvalue_threshold Numeric. P-value threshold (default: 0.05)
#' @param adjust_method Character. Multiple testing correction method: "BH", "bonferroni",
#'   "holm", "none", or "metabolite_stratified" (recommended, performs BH within each metabolite)
#' @param min_score Numeric. Minimum communication score (default: 0)
#' @return Updated scMetaLink object
#' @export
#'
#' @details
#' The \code{adjust_method} parameter controls how p-values are corrected for multiple testing:
#' \itemize{
#'   \item "metabolite_stratified" (recommended): Performs BH correction within each metabolite,
#'     then combines results. This is less conservative than global correction and more
#'     biologically meaningful since metabolites are independent biological signals.
#'   \item "BH", "bonferroni", "holm": Standard global correction methods. Can be very
#'     conservative when testing many interactions.
#'   \item "none": No correction, uses raw p-values. Use with caution.
#' }
filterSignificantInteractions <- function(object,
                                          pvalue_threshold = 0.05,
                                          adjust_method = "metabolite_stratified",
                                          min_score = 0) {
  if (is.null(object@communication_scores)) {
    stop("Communication scores not calculated. Run computeCommunication() first.")
  }

  # Parameter validation
  if (pvalue_threshold <= 0 || pvalue_threshold > 1) {
    stop("pvalue_threshold must be between 0 and 1")
  }

  valid_methods <- c(p.adjust.methods, "metabolite_stratified")
  if (!adjust_method %in% valid_methods) {
    stop("Invalid adjust_method. Use 'metabolite_stratified', 'BH', 'bonferroni', 'holm', or 'none'.")
  }

  comm_scores <- object@communication_scores
  pvalues <- object@communication_pvalues

  cell_types <- dimnames(comm_scores)[[1]]
  metabolites <- dimnames(comm_scores)[[3]]

  # If no permutation test was done, just filter by score
  if (is.null(pvalues)) {
    idx <- which(comm_scores > min_score, arr.ind = TRUE)

    if (nrow(idx) == 0) {
      message("No interactions passed the minimum score threshold")
      object@significant_interactions <- data.frame()
      return(object)
    }

    results <- data.frame(
      sender = cell_types[idx[, 1]],
      receiver = cell_types[idx[, 2]],
      metabolite = metabolites[idx[, 3]],
      communication_score = comm_scores[idx],
      stringsAsFactors = FALSE
    )
  } else if (adjust_method == "metabolite_stratified") {
    # Stratified BH correction: correct within each metabolite
    results_list <- list()

    for (met in metabolites) {
      met_idx <- which(metabolites == met)
      met_pvals <- pvalues[, , met_idx]
      met_scores <- comm_scores[, , met_idx]

      # Find interactions above min_score for this metabolite
      valid_idx <- which(met_scores > min_score)

      if (length(valid_idx) > 0) {
        valid_pvals <- met_pvals[valid_idx]
        valid_scores <- met_scores[valid_idx]

        # BH correction within this metabolite
        adj_pvals <- p.adjust(valid_pvals, method = "BH")

        # Find significant interactions
        sig_mask <- adj_pvals < pvalue_threshold

        if (any(sig_mask)) {
          sig_idx <- valid_idx[sig_mask]
          idx_2d <- arrayInd(sig_idx, dim(met_pvals))

          results_list[[met]] <- data.frame(
            sender = cell_types[idx_2d[, 1]],
            receiver = cell_types[idx_2d[, 2]],
            metabolite = met,
            communication_score = valid_scores[sig_mask],
            pvalue = valid_pvals[sig_mask],
            pvalue_adjusted = adj_pvals[sig_mask],
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (length(results_list) > 0) {
      results <- do.call(rbind, results_list)
      rownames(results) <- NULL
    } else {
      results <- data.frame()
    }
  } else {
    # Standard global correction
    idx <- which(comm_scores > min_score, arr.ind = TRUE)

    if (nrow(idx) == 0) {
      message("No interactions passed the minimum score threshold")
      object@significant_interactions <- data.frame()
      return(object)
    }

    results <- data.frame(
      sender = cell_types[idx[, 1]],
      receiver = cell_types[idx[, 2]],
      metabolite = metabolites[idx[, 3]],
      communication_score = comm_scores[idx],
      stringsAsFactors = FALSE
    )

    results$pvalue <- pvalues[idx]

    if (adjust_method == "none") {
      results$pvalue_adjusted <- results$pvalue
    } else {
      results$pvalue_adjusted <- p.adjust(results$pvalue, method = adjust_method)
    }

    results <- results[results$pvalue_adjusted < pvalue_threshold, ]
  }

  if (nrow(results) == 0) {
    message("No significant interactions found")
    object@significant_interactions <- data.frame()
    return(object)
  }

  # Add metabolite names
  db <- object@database
  results <- merge(results, db$metabolites[, c("hmdb", "metabolite")],
    by.x = "metabolite", by.y = "hmdb", all.x = TRUE
  )
  names(results)[names(results) == "metabolite.y"] <- "metabolite_name"
  names(results)[names(results) == "metabolite"] <- "metabolite_id"

  results <- results[order(-results$communication_score), ]
  rownames(results) <- NULL

  message(sprintf("Found %d significant interactions", nrow(results)))

  object@significant_interactions <- results
  object
}

#' Summarize Communication by Cell Type Pairs
#' @param object scMetaLink object
#' @param aggregate_method Character. "sum", "mean", or "count"
#' @return data.frame with summarized communication
#' @export
summarizeCommunicationPairs <- function(object, aggregate_method = "sum") {
  if (nrow(object@significant_interactions) == 0) {
    stop("No significant interactions. Run filterSignificantInteractions() first.")
  }

  if (!aggregate_method %in% c("sum", "mean", "count")) {
    stop("aggregate_method must be one of: 'sum', 'mean', 'count'")
  }

  sig <- object@significant_interactions

  if (aggregate_method == "sum") {
    result <- aggregate(communication_score ~ sender + receiver, data = sig, FUN = sum)
  } else if (aggregate_method == "mean") {
    result <- aggregate(communication_score ~ sender + receiver, data = sig, FUN = mean)
  } else {
    result <- aggregate(communication_score ~ sender + receiver, data = sig, FUN = length)
    names(result)[3] <- "n_interactions"
  }

  result[order(-result[[3]]), ]
}

#' Get Communication Summary Matrix
#' @param object scMetaLink object
#' @param aggregate_method Character. How to aggregate metabolites
#' @return Matrix of aggregated communication scores
#' @export
getCommunicationMatrix <- function(object, aggregate_method = "sum") {
  summary_df <- summarizeCommunicationPairs(object, aggregate_method)

  cell_types <- unique(c(summary_df$sender, summary_df$receiver))
  mat <- matrix(0,
    nrow = length(cell_types), ncol = length(cell_types),
    dimnames = list(cell_types, cell_types)
  )

  for (i in seq_len(nrow(summary_df))) {
    mat[summary_df$sender[i], summary_df$receiver[i]] <- summary_df[[3]][i]
  }

  mat
}
