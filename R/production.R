#' @title Infer Metabolite Production Potential
#' @description Infer metabolite production potential for each cell type based on
#'   enzyme expression patterns. Production potential reflects a cell type's capacity
#'   to synthesize and secrete metabolites for intercellular communication.
#'
#' @param object A scMetaLink object
#' @param method Character. Scoring method: "mean", "proportion", or "combined"
#' @param mean_method Character. Method for calculating mean expression:
#'   "arithmetic" (standard mean) or "trimean" (more robust to outliers and dropout).
#'   Trimean is recommended for single-cell data with high dropout rates.
#' @param min_expression Numeric. Minimum expression threshold
#' @param min_pct Numeric. Minimum percentage of expressing cells (0-1)
#' @param consider_degradation Logical. Subtract degradation enzyme expression
#' @param consider_secretion Logical. Weight by secretion potential
#' @param normalize Logical. Normalize scores across cell types
#' @param verbose Logical. Print progress messages
#'
#' @return Updated scMetaLink object with production_scores slot filled
#' @export
inferProduction <- function(object,
                            method = "combined",
                            mean_method = "arithmetic",
                            min_expression = 0,
                            min_pct = 0.1,
                            consider_degradation = TRUE,
                            consider_secretion = TRUE,
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

  db <- object@database
  expr_data <- object@expression_data
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column
  cell_types <- unique(cell_meta[[cell_type_col]])

  if (verbose) message(sprintf("Inferring production potential for %d cell types...", length(cell_types)))

  # Step 1: Calculate cell type expression profiles (vectorized)
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

  # Step 2: Build metabolite-gene mapping matrix
  if (verbose) message("  Building metabolite-gene mapping...")
  
  prod_enzymes <- .get_production_enzymes(db)
  prod_enzymes <- prod_enzymes[!is.na(prod_enzymes$gene_symbol), ]
  prod_enzymes <- prod_enzymes[prod_enzymes$gene_symbol %in% rownames(gene_scores), ]
  
  if (nrow(prod_enzymes) == 0) {
    stop("No production enzymes found in the expression data")
  }
  
  if (verbose) message(sprintf("  Found %d production enzyme-metabolite pairs", nrow(prod_enzymes)))
  
  metabolites <- unique(prod_enzymes$hmdb)
  genes <- rownames(gene_scores)
  
  # Create sparse mapping matrix (metabolites x genes)
  met_gene_map <- Matrix::sparseMatrix(
    i = match(prod_enzymes$hmdb, metabolites),
    j = match(prod_enzymes$gene_symbol, genes),
    x = 1,
    dims = c(length(metabolites), length(genes)),
    dimnames = list(metabolites, genes)
  )
  
  # Count enzymes per metabolite for averaging
  enzyme_counts <- Matrix::rowSums(met_gene_map)
  enzyme_counts[enzyme_counts == 0] <- 1  # Avoid division by zero
  
  # Step 3: Matrix multiplication for production scores
  if (verbose) message("  Computing production scores (matrix multiplication)...")
  
  production_scores <- as.matrix(met_gene_map %*% gene_scores) / enzyme_counts
  
  # Step 4: Subtract degradation (if requested)
  if (consider_degradation) {
    if (verbose) message("  Adjusting for degradation...")
    
    deg_enzymes <- .get_degradation_enzymes(db)
    deg_enzymes <- deg_enzymes[!is.na(deg_enzymes$gene_symbol), ]
    deg_enzymes <- deg_enzymes[deg_enzymes$gene_symbol %in% genes, ]
    deg_enzymes <- deg_enzymes[deg_enzymes$hmdb %in% metabolites, ]
    
    if (nrow(deg_enzymes) > 0) {
      deg_map <- Matrix::sparseMatrix(
        i = match(deg_enzymes$hmdb, metabolites),
        j = match(deg_enzymes$gene_symbol, genes),
        x = 1,
        dims = c(length(metabolites), length(genes)),
        dimnames = list(metabolites, genes)
      )
      deg_counts <- Matrix::rowSums(deg_map)
      deg_counts[deg_counts == 0] <- 1
      
      deg_scores <- as.matrix(deg_map %*% gene_scores) / deg_counts
      production_scores <- production_scores - 0.5 * deg_scores
      production_scores[production_scores < 0] <- 0
    }
  }
  
  # Step 5: Weight by secretion potential
  if (consider_secretion) {
    if (verbose) message("  Applying secretion potential weights...")
    
    extra_mets <- .get_extracellular_metabolites(db)
    
    # Secretion weight: 1.0 for extracellular, 0.5 for others
    secretion_weights <- ifelse(metabolites %in% extra_mets, 1.0, 0.5)
    production_scores <- production_scores * secretion_weights
  }
  
  # Step 6: Normalize
  if (normalize && nrow(production_scores) > 0) {
    if (verbose) message("  Normalizing scores...")
    
    # Z-score per metabolite
    row_means <- rowMeans(production_scores, na.rm = TRUE)
    row_sds <- apply(production_scores, 1, sd, na.rm = TRUE)
    row_sds[row_sds == 0] <- 1
    production_scores <- (production_scores - row_means) / row_sds
    
    # Scale to 0-1
    min_val <- min(production_scores, na.rm = TRUE)
    max_val <- max(production_scores, na.rm = TRUE)
    if (max_val > min_val) {
      production_scores <- (production_scores - min_val) / (max_val - min_val)
    }
  }
  
  object@production_scores <- production_scores
  object@parameters$production <- list(
    method = method,
    mean_method = mean_method,
    min_expression = min_expression,
    min_pct = min_pct,
    consider_degradation = consider_degradation,
    consider_secretion = consider_secretion,
    normalize = normalize
  )
  
  if (verbose) message(sprintf("  Computed production scores for %d metabolites", nrow(production_scores)))
  if (verbose) message("Done!")
  object
}

#' Calculate Cell Type Expression Profiles
#'
#' @description Efficiently calculate mean expression and percentage of expressing
#'   cells for each gene across cell types using vectorized operations.
#'
#' @param expr_data Expression matrix (genes x cells)
#' @param cell_meta Cell metadata data.frame
#' @param cell_type_col Column name for cell type annotation
#' @param min_expression Minimum expression threshold
#' @param mean_method Method for calculating mean: "arithmetic" or "trimean"
#'
#' @details
#' The trimean is defined as: (Q1 + 2*Q2 + Q3) / 4, where Q1, Q2, Q3 are the
#' 25th, 50th, and 75th percentiles respectively. This provides a robust
#' estimate of central tendency that is less sensitive to outliers and the
#' high dropout rates common in single-cell RNA-seq data.
#'
#' @return List with mean_expr and pct_expr matrices
#' @keywords internal
.calculate_celltype_expression_fast <- function(expr_data, cell_meta, cell_type_col,
                                                 min_expression = 0,
                                                 mean_method = "arithmetic") {

  cell_types <- unique(cell_meta[[cell_type_col]])
  n_genes <- nrow(expr_data)
  n_types <- length(cell_types)

  mean_expr <- matrix(0, nrow = n_genes, ncol = n_types,
                      dimnames = list(rownames(expr_data), cell_types))
  pct_expr <- mean_expr

  for (ct in cell_types) {
    cells <- rownames(cell_meta)[cell_meta[[cell_type_col]] == ct]
    ct_expr <- expr_data[, cells, drop = FALSE]

    # Calculate percentage expressing
    if (inherits(ct_expr, "dgCMatrix")) {
      pct_expr[, ct] <- Matrix::rowMeans(ct_expr > min_expression)
    } else {
      pct_expr[, ct] <- rowMeans(ct_expr > min_expression, na.rm = TRUE)
    }

    # Calculate mean expression
    if (mean_method == "arithmetic") {
      if (inherits(ct_expr, "dgCMatrix")) {
        mean_expr[, ct] <- Matrix::rowMeans(ct_expr)
      } else {
        mean_expr[, ct] <- rowMeans(ct_expr, na.rm = TRUE)
      }
    } else if (mean_method == "trimean") {
      # Trimean: (Q1 + 2*Q2 + Q3) / 4
      # More robust to outliers and dropout in single-cell data
      mean_expr[, ct] <- .row_trimean(ct_expr)
    }
  }

  list(mean_expr = mean_expr, pct_expr = pct_expr)
}

#' Calculate Row-wise Trimean
#'
#' @description Calculate the trimean for each row of a matrix.
#'   Trimean = (Q1 + 2*Q2 + Q3) / 4, providing a robust central tendency estimate.
#'
#' @param mat Matrix (genes x cells), can be sparse or dense
#' @return Numeric vector of trimean values for each row
#' @keywords internal
.row_trimean <- function(mat) {

  if (inherits(mat, "dgCMatrix")) {
    mat <- as.matrix(mat)
  }

  # Vectorized quantile calculation
  # apply is used here as matrixStats::rowQuantiles requires additional dependency
  q <- apply(mat, 1, function(x) {
    quantile(x, probs = c(0.25, 0.5, 0.75), na.rm = TRUE, names = FALSE)
  })

  # q is 3 x n_genes matrix
  # Trimean = (Q1 + 2*Q2 + Q3) / 4
  (q[1, ] + 2 * q[2, ] + q[3, ]) / 4
}

#' Get Top Producing Cell Types for a Metabolite
#' @param object scMetaLink object
#' @param metabolite Character. Metabolite HMDB ID or name
#' @param top_n Integer. Number of top cell types to return
#' @return data.frame with top producing cell types
#' @export
getTopProducers <- function(object, metabolite, top_n = 5) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }
  if (is.null(object@production_scores)) {
    stop("Production scores not calculated. Run inferProduction() first.")
  }
  if (top_n < 1) {
    stop("top_n must be at least 1")
  }
  
  if (metabolite %in% rownames(object@production_scores)) {
    met_id <- metabolite
  } else {
    db <- object@database
    met_match <- db$metabolites$hmdb[grepl(metabolite, db$metabolites$metabolite, ignore.case = TRUE)]
    if (length(met_match) == 0) stop("Metabolite not found")
    met_id <- met_match[1]
    if (!met_id %in% rownames(object@production_scores)) {
      stop("Metabolite not found in production scores")
    }
  }
  
  scores <- object@production_scores[met_id, ]
  scores <- sort(scores, decreasing = TRUE)
  
  n_return <- min(top_n, length(scores))
  
  data.frame(
    cell_type = names(scores)[1:n_return],
    production_score = as.numeric(scores[1:n_return]),
    rank = 1:n_return,
    stringsAsFactors = FALSE
  )
}
