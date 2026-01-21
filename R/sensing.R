#' @title Infer Metabolite Sensing Capability
#' @description Infer metabolite sensing capability for each cell type based on
#'   receptor and transporter expression. Sensing capability reflects a cell type's
#'   capacity to detect and respond to extracellular metabolites.
#'
#' @param object A scMetaLink object
#' @param method Character. Scoring method: "mean", "proportion", or "combined"
#' @param mean_method Character. Method for calculating mean expression:
#'   "arithmetic" (standard mean) or "trimean" (more robust to outliers and dropout).
#' @param min_expression Numeric. Minimum expression threshold
#' @param min_pct Numeric. Minimum percentage of expressing cells (0-1)
#' @param weight_by_affinity Logical. Weight by receptor-metabolite affinity score
#' @param include_transporters Logical. Include uptake transporters in sensing
#' @param use_hill Logical. Apply Hill function transformation to model receptor
#'   binding saturation kinetics. When TRUE, high expression levels show diminishing
#'   returns, reflecting biological receptor saturation.
#' @param hill_n Numeric. Hill coefficient (cooperativity). Default 1 (no cooperativity).
#'   Values > 1 indicate positive cooperativity.
#' @param hill_Kh Numeric. Half-maximal response threshold (0-1 scale after normalization).
#'   Default 0.5. Lower values mean saturation occurs at lower expression levels.
#' @param normalize Logical. Normalize scores across cell types
#' @param verbose Logical. Print progress messages
#'
#' @details
#' The Hill function transformation models receptor-ligand binding dynamics:
#' \deqn{P = \frac{E^n}{K_h^n + E^n}}
#' where E is expression, n is the Hill coefficient, and Kh is the half-maximal
#' threshold. This reflects the biological reality that receptor signaling
#' saturates at high ligand/receptor concentrations.
#'
#' @return Updated scMetaLink object with sensing_scores slot filled
#' @export
inferSensing <- function(object,
                         method = "combined",
                         mean_method = "arithmetic",
                         min_expression = 0,
                         min_pct = 0.1,
                         weight_by_affinity = TRUE,
                         include_transporters = TRUE,
                         use_hill = FALSE,
                         hill_n = 1,
                         hill_Kh = 0.5,
                         normalize = TRUE,
                         verbose = TRUE) {
  # Input validation
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (!method %in% c("mean", "proportion", "combined")) {
    stop("method must be one of: 'mean', 'proportion', 'combined'")
  }
  if (!mean_method %in% c("arithmetic", "trimean")) {
    stop("mean_method must be one of: 'arithmetic', 'trimean'")
  }
  if (min_pct < 0 || min_pct > 1) {
    stop("min_pct must be between 0 and 1")
  }
  if (min_expression < 0) {
    stop("min_expression must be non-negative")
  }
  if (hill_n <= 0) {
    stop("hill_n must be positive")
  }
  if (hill_Kh <= 0 || hill_Kh >= 1) {
    stop("hill_Kh must be between 0 and 1 (exclusive)")
  }

  db <- object@database
  expr_data <- object@expression_data
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column
  cell_types <- unique(cell_meta[[cell_type_col]])

  if (verbose) message(sprintf("Inferring sensing capability for %d cell types...", length(cell_types)))

  # Step 1: Calculate cell type expression profiles
  if (verbose) message("  Calculating cell type expression profiles...")
  expr_profiles <- .calculate_celltype_expression_fast(
    expr_data = expr_data,
    cell_meta = cell_meta,
    cell_type_col = cell_type_col,
    min_expression = min_expression,
    mean_method = mean_method
  )

  # Combined score
  if (method == "combined") {
    gene_scores <- expr_profiles$mean_expr * expr_profiles$pct_expr
  } else if (method == "mean") {
    gene_scores <- expr_profiles$mean_expr
  } else {
    gene_scores <- expr_profiles$pct_expr
  }

  # Step 2: Get sensing proteins
  if (verbose) message("  Building metabolite-receptor mapping...")

  sensing_proteins <- .get_sensing_proteins(db)
  sensing_proteins <- sensing_proteins[!is.na(sensing_proteins$gene_symbol), ]

  if (include_transporters) {
    in_trans <- .get_transporters(db, direction = "in")
    in_trans <- in_trans[!is.na(in_trans$gene_symbol), ]
    if (nrow(in_trans) > 0) {
      sensing_proteins <- rbind(sensing_proteins, in_trans)
    }
  }

  # Filter to available genes
  genes <- rownames(gene_scores)
  sensing_proteins <- sensing_proteins[sensing_proteins$gene_symbol %in% genes, ]

  if (nrow(sensing_proteins) == 0) {
    stop("No sensing proteins found in the expression data")
  }

  if (verbose) message(sprintf("  Found %d receptor-metabolite pairs", nrow(sensing_proteins)))

  # Remove duplicates, keep highest score
  sensing_proteins <- sensing_proteins[order(-sensing_proteins$combined_score), ]
  sensing_proteins <- sensing_proteins[!duplicated(paste(sensing_proteins$hmdb, sensing_proteins$gene_symbol)), ]

  metabolites <- unique(sensing_proteins$hmdb)

  # Step 3: Create sparse mapping matrix with affinity weights
  if (weight_by_affinity) {
    # Normalize scores to 0-1
    weights <- sensing_proteins$combined_score
    weights[is.na(weights)] <- 500 # Default for missing scores
    weights <- weights / 1000
  } else {
    weights <- rep(1, nrow(sensing_proteins))
  }

  met_gene_map <- Matrix::sparseMatrix(
    i = match(sensing_proteins$hmdb, metabolites),
    j = match(sensing_proteins$gene_symbol, genes),
    x = weights,
    dims = c(length(metabolites), length(genes)),
    dimnames = list(metabolites, genes)
  )

  # Count receptors per metabolite
  receptor_counts <- Matrix::rowSums(met_gene_map > 0)
  receptor_counts[receptor_counts == 0] <- 1

  # Step 4: Matrix multiplication for sensing scores
  if (verbose) message("  Computing sensing scores (matrix multiplication)...")

  sensing_scores <- as.matrix(met_gene_map %*% gene_scores) / receptor_counts

  # Step 5: Apply Hill function transformation (if requested)
  if (use_hill) {
    if (verbose) message("  Applying Hill function transformation...")

    # First normalize to 0-1 range for Hill function
    # Then apply Hill: P = E^n / (Kh^n + E^n)
    sensing_min <- min(sensing_scores, na.rm = TRUE)
    sensing_max <- max(sensing_scores, na.rm = TRUE)

    if (sensing_max > sensing_min) {
      # Scale to 0-1
      sensing_norm <- (sensing_scores - sensing_min) / (sensing_max - sensing_min)

      # Apply Hill function
      sensing_scores <- .hill_transform(sensing_norm, n = hill_n, Kh = hill_Kh)
    }
  }

  # Step 6: Normalize
  if (normalize && nrow(sensing_scores) > 0) {
    if (verbose) message("  Normalizing scores...")

    row_means <- rowMeans(sensing_scores, na.rm = TRUE)
    row_sds <- apply(sensing_scores, 1, sd, na.rm = TRUE)
    row_sds[row_sds == 0] <- 1
    sensing_scores <- (sensing_scores - row_means) / row_sds

    min_val <- min(sensing_scores, na.rm = TRUE)
    max_val <- max(sensing_scores, na.rm = TRUE)
    if (max_val > min_val) {
      sensing_scores <- (sensing_scores - min_val) / (max_val - min_val)
    }
  }

  object@sensing_scores <- sensing_scores
  object@parameters$sensing <- list(
    method = method,
    mean_method = mean_method,
    min_expression = min_expression,
    min_pct = min_pct,
    weight_by_affinity = weight_by_affinity,
    include_transporters = include_transporters,
    use_hill = use_hill,
    hill_n = hill_n,
    hill_Kh = hill_Kh,
    normalize = normalize
  )

  if (verbose) message(sprintf("  Computed sensing scores for %d metabolites", nrow(sensing_scores)))
  if (verbose) message("Done!")
  object
}

#' Get Top Sensing Cell Types for a Metabolite
#' @param object scMetaLink object
#' @param metabolite Character. Metabolite HMDB ID or name
#' @param top_n Integer. Number of top cell types to return
#' @return data.frame with top sensing cell types
#' @export
getTopSensors <- function(object, metabolite, top_n = 5) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (is.null(object@sensing_scores)) {
    stop("Sensing scores not calculated. Run inferSensing() first.")
  }
  if (top_n < 1) {
    stop("top_n must be at least 1")
  }

  if (metabolite %in% rownames(object@sensing_scores)) {
    met_id <- metabolite
  } else {
    db <- object@database
    met_match <- db$metabolites$hmdb[grepl(metabolite, db$metabolites$metabolite, ignore.case = TRUE)]
    if (length(met_match) == 0) stop("Metabolite not found")
    met_id <- met_match[1]
    if (!met_id %in% rownames(object@sensing_scores)) {
      stop("Metabolite not found in sensing scores")
    }
  }

  scores <- object@sensing_scores[met_id, ]
  scores <- sort(scores, decreasing = TRUE)

  n_return <- min(top_n, length(scores))

  data.frame(
    cell_type = names(scores)[1:n_return],
    sensing_score = as.numeric(scores[1:n_return]),
    rank = 1:n_return,
    stringsAsFactors = FALSE
  )
}

#' Hill Function Transformation
#'
#' @description Apply Hill function to model receptor-ligand binding saturation.
#'   The Hill function captures the sigmoidal relationship between expression
#'   and biological response, reflecting receptor saturation at high concentrations.
#'
#' @param x Numeric matrix or vector of expression values (should be 0-1 scaled)
#' @param n Hill coefficient (cooperativity parameter). n=1 gives standard
#'   Michaelis-Menten kinetics; n>1 indicates positive cooperativity.
#' @param Kh Half-maximal response threshold. The value of x at which response = 0.5.
#'
#' @return Transformed values in the same structure as input
#' @keywords internal
.hill_transform <- function(x, n = 1, Kh = 0.5) {
  # Hill function: P = x^n / (Kh^n + x^n)
  # Handles edge cases to avoid numerical issues
  x_n <- x^n
  Kh_n <- Kh^n

  result <- x_n / (Kh_n + x_n)

  # Handle edge cases
  result[is.na(result)] <- 0
  result[is.infinite(result)] <- 1

  result
}

#' Get Receptors for a Metabolite
#' @param metabolite Character. Metabolite HMDB ID or name
#' @param include_transporters Logical. Include transporters
#' @return data.frame with receptor information
#' @export
getMetaboliteReceptors <- function(metabolite, include_transporters = TRUE) {
  db <- .load_metalinksdb()

  if (!metabolite %in% db$metabolites$hmdb) {
    met_match <- db$metabolites$hmdb[grepl(metabolite, db$metabolites$metabolite, ignore.case = TRUE)]
    if (length(met_match) == 0) stop("Metabolite not found")
    metabolite <- met_match[1]
  }

  sensing <- .get_sensing_proteins(db)
  receptors <- sensing[sensing$hmdb == metabolite, ]

  if (include_transporters) {
    trans_in <- .get_transporters(db, direction = "in")
    trans <- trans_in[trans_in$hmdb == metabolite, ]
    if (nrow(trans) > 0) receptors <- rbind(receptors, trans)
  }

  if (nrow(receptors) == 0) {
    message("No receptors found for this metabolite")
    return(data.frame())
  }

  receptors <- receptors[!duplicated(receptors$gene_symbol), ]
  receptors <- receptors[order(-receptors$combined_score), ]

  # Select columns to return
  cols <- c("gene_symbol", "protein_type", "combined_score", "metabolite")
  cols <- cols[cols %in% names(receptors)]
  receptors[, cols]
}
