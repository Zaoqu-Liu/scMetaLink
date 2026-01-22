# ============================================================================
# Lactate Signaling Analysis Functions
# ============================================================================
#
# Exported functions for lactate-mediated cell communication analysis.
# This module implements analysis of both direct (HCAR1) and indirect
# (proton-sensing GPCRs) lactate signaling pathways.
#
# Scientific basis:
# - Lactate (pKa=3.86) fully dissociates at physiological pH 7.4
# - Direct signaling: Lactate binds HCAR1/GPR81
# - Indirect signaling: H+ from lactate dissociation activates GPR4/65/68/132
#
# References:
# - PLOS Biology 2024 (HCAR1 cryo-EM structure)
# - Nature Metabolism 2024 (GPR81 function)
# - Cell Press iScience 2019 (MCT transport direction)
# - Reactome R-HSA-444731 (proton-sensing GPCRs)
#
# ============================================================================


#' @title Infer Lactate-Mediated Cell Communication
#'
#' @description
#' Infers both direct and indirect lactate signaling in single-cell data.
#'
#' \strong{Direct signaling}: Lactate binds to HCAR1 (GPR81), the only confirmed
#' lactate GPCR (validated by cryo-EM structure, PLOS Biology 2024).
#'
#' \strong{Indirect signaling}: Lactate dissociation (pKa=3.86) produces H+ ions
#' that activate proton-sensing GPCRs (GPR4, GPR65, GPR68, GPR132).
#'
#' @param object A scMetaLink object with expression data
#' @param include_direct Logical. Include direct lactate-HCAR1 signaling. Default TRUE.
#' @param include_indirect Logical. Include indirect lactate-H+-GPCR signaling. Default TRUE.
#' @param method Character. Scoring method: "combined" (recommended), "mean", or "proportion".
#'   Default "combined".
#' @param comm_method Character. Communication score method: "geometric" (default),
#'   "product", or "harmonic".
#' @param min_pct Numeric. Minimum percentage of expressing cells (0-1). Default 0.1.
#' @param min_production Numeric. Minimum production score threshold (0-1). Default 0.
#' @param min_sensing Numeric. Minimum sensing score threshold (0-1). Default 0.
#' @param consider_uptake Logical. Include MCT uptake transporters in direct sensing.
#'   Default TRUE.
#' @param consider_degradation Logical. Subtract degradation enzyme expression from
#'   production scores. Default TRUE.
#' @param normalize Logical. Normalize scores across cell types. Default TRUE.
#' @param n_permutations Integer. Number of permutations for significance testing.
#'   Set to 0 to skip. Default 100.
#' @param seed Integer. Random seed for reproducibility. Default 42.
#' @param verbose Logical. Print progress messages. Default TRUE.
#'
#' @return Updated scMetaLink object with lactate_signaling results stored in
#'   the parameters slot, containing:
#'   \item{production}{Lactate production scores per cell type}
#'   \item{direct_sensing}{Direct sensing scores (HCAR1)}
#'   \item{indirect_sensing}{Indirect sensing scores (proton GPCRs)}
#'   \item{direct_communication}{Communication matrix via direct pathway}
#'   \item{indirect_communication}{Communication matrix via indirect pathway}
#'   \item{combined_communication}{Sum of direct and indirect communication}
#'   \item{pvalues}{Permutation-based p-values (if n_permutations > 0)}
#'   \item{gene_contributions}{Expression contribution of each gene}
#'   \item{parameters}{Analysis parameters used}
#'
#' @export
#'
#' @references
#' \itemize{
#'   \item HCAR1 structure: PLOS Biology (2024)
#'   \item GPR81 function: Nature Metabolism (2024)
#'   \item Lactate signaling: Signal Transduction and Targeted Therapy (2024)
#'   \item Proton-sensing GPCRs: Reactome R-HSA-444731
#' }
#'
#' @examples
#' \donttest{
#' # Load example data
#' data(crc_example)
#'
#' # Create scMetaLink object
#' obj <- createScMetaLink(crc_expr, crc_meta, "cell_type")
#'
#' # Run lactate signaling analysis
#' obj <- inferLactateSignaling(obj)
#'
#' # Access results
#' lactate_results <- obj@parameters$lactate_signaling
#' head(lactate_results$direct_communication)
#'
#' # Analyze only indirect (proton) signaling
#' obj <- inferLactateSignaling(obj, include_direct = FALSE)
#' }
inferLactateSignaling <- function(object,
                                   include_direct = TRUE,
                                   include_indirect = TRUE,
                                   method = "combined",
                                   comm_method = "geometric",
                                   min_pct = 0.1,
                                   min_production = 0,
                                   min_sensing = 0,
                                   consider_uptake = TRUE,
                                   consider_degradation = TRUE,
                                   normalize = TRUE,
                                   n_permutations = 100,
                                   seed = 42,
                                   verbose = TRUE) {
  # Input validation
 if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  if (!include_direct && !include_indirect) {
    stop("At least one of include_direct or include_indirect must be TRUE")
  }

  if (!method %in% c("combined", "mean", "proportion")) {
    stop("method must be one of: 'combined', 'mean', 'proportion'")
  }

  if (!comm_method %in% c("geometric", "product", "harmonic")) {
    stop("comm_method must be one of: 'geometric', 'product', 'harmonic'")
  }

  if (min_pct < 0 || min_pct > 1) {
    stop("min_pct must be between 0 and 1")
  }

  if (n_permutations < 0) {
    stop("n_permutations must be non-negative")
  }

  set.seed(seed)

  # Get data from object
  expr_data <- object@expression_data
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column
  cell_types <- unique(cell_meta[[cell_type_col]])

  if (verbose) {
    message("Inferring lactate-mediated cell communication...")
    message(sprintf("  Cell types: %d", length(cell_types)))
    message(sprintf("  Direct pathway: %s", include_direct))
    message(sprintf("  Indirect pathway: %s", include_indirect))
  }

  # Get gene sets
  gene_sets <- .get_lactate_gene_sets()

  # Get all genes needed
  all_genes <- unique(c(
    gene_sets$production$synthesis,
    gene_sets$production$export,
    gene_sets$degradation$enzymes,
    gene_sets$direct_sensing$receptor,
    gene_sets$indirect_sensing$proton_receptors,
    gene_sets$uptake$import
  ))

  # Check gene availability
  available_genes <- intersect(all_genes, rownames(expr_data))
  if (verbose) {
    message(sprintf("  Available lactate genes: %d/%d", length(available_genes), length(all_genes)))
  }

  if (length(available_genes) == 0) {
    stop("No lactate-related genes found in expression data")
  }

  # Step 1: Calculate gene expression scores
  if (verbose) message("  Calculating gene expression scores...")

  gene_scores <- .calc_lactate_gene_scores(
    expr_data = expr_data,
    cell_meta = cell_meta,
    cell_type_col = cell_type_col,
    genes = all_genes,
    min_expression = 0,
    method = method
  )

  if (is.null(gene_scores) || ncol(gene_scores) == 0) {
    stop("Failed to calculate gene scores")
  }

  # Step 2: Calculate production scores
  if (verbose) message("  Calculating lactate production scores...")

  production_scores <- .calc_lactate_production(
    gene_scores = gene_scores,
    gene_sets = gene_sets,
    consider_export = TRUE,
    consider_degradation = consider_degradation
  )

  # Step 3: Calculate sensing scores
  direct_sensing <- NULL
  indirect_sensing <- NULL

  if (include_direct) {
    if (verbose) message("  Calculating direct sensing scores (HCAR1)...")
    direct_sensing <- .calc_direct_sensing(
      gene_scores = gene_scores,
      gene_sets = gene_sets,
      include_uptake = consider_uptake
    )
  }

  if (include_indirect) {
    if (verbose) message("  Calculating indirect sensing scores (GPR4/65/68/132)...")
    indirect_sensing <- .calc_proton_sensing(
      gene_scores = gene_scores,
      gene_sets = gene_sets,
      use_weights = TRUE
    )
  }

  # Normalize scores if requested
  if (normalize) {
    if (verbose) message("  Normalizing scores...")
    production_scores <- .normalize_lactate_scores(production_scores)
    if (!is.null(direct_sensing)) {
      direct_sensing <- .normalize_lactate_scores(direct_sensing)
    }
    if (!is.null(indirect_sensing)) {
      indirect_sensing <- .normalize_lactate_scores(indirect_sensing)
    }
  }

  # Step 4: Calculate communication scores
  direct_comm <- NULL
  indirect_comm <- NULL
  combined_comm <- NULL

  if (include_direct && !is.null(direct_sensing)) {
    if (verbose) message("  Calculating direct communication scores...")
    direct_comm <- .calc_lactate_communication(
      production_scores = production_scores,
      sensing_scores = direct_sensing,
      method = comm_method,
      min_production = min_production,
      min_sensing = min_sensing
    )
  }

  if (include_indirect && !is.null(indirect_sensing)) {
    if (verbose) message("  Calculating indirect communication scores...")
    indirect_comm <- .calc_lactate_communication(
      production_scores = production_scores,
      sensing_scores = indirect_sensing,
      method = comm_method,
      min_production = min_production,
      min_sensing = min_sensing
    )
  }

  # Calculate combined communication
  if (!is.null(direct_comm) && !is.null(indirect_comm)) {
    combined_comm <- direct_comm + indirect_comm
  } else if (!is.null(direct_comm)) {
    combined_comm <- direct_comm
  } else if (!is.null(indirect_comm)) {
    combined_comm <- indirect_comm
  }

  # Step 5: Permutation testing
  pvalues <- list()

  if (n_permutations > 0) {
    if (include_direct && !is.null(direct_comm)) {
      if (verbose) message(sprintf("  Running permutation test for direct pathway (%d permutations)...", n_permutations))
      pvalues$direct <- .run_lactate_permutation(
        object = object,
        comm_scores = direct_comm,
        pathway = "direct",
        n_permutations = n_permutations,
        method = method,
        comm_method = comm_method,
        min_production = min_production,
        min_sensing = min_sensing,
        include_uptake = consider_uptake,
        verbose = verbose
      )
    }

    if (include_indirect && !is.null(indirect_comm)) {
      if (verbose) message(sprintf("  Running permutation test for indirect pathway (%d permutations)...", n_permutations))
      pvalues$indirect <- .run_lactate_permutation(
        object = object,
        comm_scores = indirect_comm,
        pathway = "indirect",
        n_permutations = n_permutations,
        method = method,
        comm_method = comm_method,
        min_production = min_production,
        min_sensing = min_sensing,
        use_weights = TRUE,
        verbose = verbose
      )
    }
  }

  # Get gene contributions
  gene_contributions <- .get_gene_contributions(gene_scores, gene_sets)

  # Store results
  lactate_results <- list(
    production = production_scores,
    direct_sensing = direct_sensing,
    indirect_sensing = indirect_sensing,
    direct_communication = direct_comm,
    indirect_communication = indirect_comm,
    combined_communication = combined_comm,
    pvalues = pvalues,
    gene_contributions = gene_contributions,
    parameters = list(
      include_direct = include_direct,
      include_indirect = include_indirect,
      method = method,
      comm_method = comm_method,
      min_pct = min_pct,
      min_production = min_production,
      min_sensing = min_sensing,
      consider_uptake = consider_uptake,
      consider_degradation = consider_degradation,
      normalize = normalize,
      n_permutations = n_permutations,
      seed = seed,
      timestamp = Sys.time()
    )
  )

  object@parameters$lactate_signaling <- lactate_results

  if (verbose) {
    message("Done!")
    message(sprintf("  Production scores for %d cell types", length(production_scores)))
    if (!is.null(direct_comm)) {
      message(sprintf("  Direct communication: %d x %d matrix", nrow(direct_comm), ncol(direct_comm)))
    }
    if (!is.null(indirect_comm)) {
      message(sprintf("  Indirect communication: %d x %d matrix", nrow(indirect_comm), ncol(indirect_comm)))
    }
  }

  object
}


#' @title Get Lactate Signaling Gene Sets
#'
#' @description
#' Returns curated, literature-validated gene sets involved in lactate signaling pathways.
#' Useful for pathway analysis, visualization, and validation.
#'
#' \strong{Gene selection criteria}:
#' \itemize{
#'   \item Synthesis: LDHA family (excludes LDHB which prefers reverse reaction)
#'   \item Export: MCT4 (SLC16A3) as primary exporter, plus other MCTs and BSG chaperone
#'   \item Direct sensing: HCAR1 only (the sole confirmed lactate GPCR)
#'   \item Indirect sensing: Classic proton-sensing GPCRs (GPR4, GPR65, GPR68, GPR132)
#'   \item Uptake: MCT1 (SLC16A1) as primary importer, plus BSG chaperone
#' }
#'
#' @param category Character. Gene category to return:
#'   \itemize{
#'     \item "all" (default): All gene sets
#'     \item "production": Synthesis and export genes
#'     \item "degradation": Lactate consumption enzymes
#'     \item "direct_sensing": HCAR1 receptor
#'     \item "indirect_sensing": Proton-sensing GPCRs
#'     \item "uptake": Import transporters
#'   }
#'
#' @return A named list of gene vectors. If category is "all", returns the complete
#'   nested list structure. Otherwise, returns the specific category.
#'
#' @export
#'
#' @examples
#' # Get all gene sets
#' genes <- getLactateGenes()
#' names(genes)
#'
#' # Get synthesis enzymes
#' genes$production$synthesis
#' # [1] "LDHA" "LDHC" "LDHAL6A" "LDHAL6B"
#'
#' # Get proton-sensing receptors
#' genes$indirect_sensing$proton_receptors
#' # [1] "GPR4" "GPR65" "GPR68" "GPR132"
#'
#' # Get only production genes
#' prod_genes <- getLactateGenes("production")
#'
#' # Get direct sensing receptor
#' getLactateGenes("direct_sensing")$receptor
#' # [1] "HCAR1"
getLactateGenes <- function(category = "all") {
  valid_categories <- c("all", "production", "degradation",
                        "direct_sensing", "indirect_sensing", "uptake")

  if (!category %in% valid_categories) {
    stop(paste("category must be one of:", paste(valid_categories, collapse = ", ")))
  }

  gene_sets <- .get_lactate_gene_sets()

  if (category == "all") {
    return(gene_sets)
  } else {
    return(gene_sets[[category]])
  }
}


#' @title Get Top Lactate Producers
#'
#' @description
#' Returns cell types ranked by lactate production potential.
#'
#' @param object scMetaLink object with lactate signaling results
#' @param top_n Integer. Number of top cell types to return. Default 5.
#'
#' @return data.frame with cell types, production scores, and ranks
#' @export
#'
#' @examples
#' \donttest{
#' data(crc_example)
#' obj <- createScMetaLink(crc_expr, crc_meta, "cell_type")
#' obj <- inferLactateSignaling(obj)
#' getTopLactateProducers(obj)
#' }
getTopLactateProducers <- function(object, top_n = 5) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  if (is.null(object@parameters$lactate_signaling)) {
    stop("Lactate signaling not calculated. Run inferLactateSignaling() first.")
  }

  production <- object@parameters$lactate_signaling$production
  production <- sort(production, decreasing = TRUE)

  n_return <- min(top_n, length(production))

  data.frame(
    cell_type = names(production)[1:n_return],
    production_score = as.numeric(production[1:n_return]),
    rank = 1:n_return,
    stringsAsFactors = FALSE
  )
}


#' @title Get Top Lactate Sensors
#'
#' @description
#' Returns cell types ranked by lactate sensing capability.
#'
#' @param object scMetaLink object with lactate signaling results
#' @param pathway Character. "direct" (HCAR1), "indirect" (proton GPCRs), or "both".
#'   Default "both".
#' @param top_n Integer. Number of top cell types to return. Default 5.
#'
#' @return data.frame with cell types, sensing scores, and ranks
#' @export
#'
#' @examples
#' \donttest{
#' data(crc_example)
#' obj <- createScMetaLink(crc_expr, crc_meta, "cell_type")
#' obj <- inferLactateSignaling(obj)
#'
#' # Get top sensors for indirect pathway
#' getTopLactateSensors(obj, pathway = "indirect")
#' }
getTopLactateSensors <- function(object, pathway = "both", top_n = 5) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  if (is.null(object@parameters$lactate_signaling)) {
    stop("Lactate signaling not calculated. Run inferLactateSignaling() first.")
  }

  if (!pathway %in% c("direct", "indirect", "both")) {
    stop("pathway must be 'direct', 'indirect', or 'both'")
  }

  lactate_res <- object@parameters$lactate_signaling

  if (pathway == "direct") {
    if (is.null(lactate_res$direct_sensing)) {
      stop("Direct sensing not calculated. Run inferLactateSignaling() with include_direct=TRUE")
    }
    sensing <- lactate_res$direct_sensing
    pathway_label <- "direct"
  } else if (pathway == "indirect") {
    if (is.null(lactate_res$indirect_sensing)) {
      stop("Indirect sensing not calculated. Run inferLactateSignaling() with include_indirect=TRUE")
    }
    sensing <- lactate_res$indirect_sensing
    pathway_label <- "indirect"
  } else {
    # Combine both pathways
    direct <- lactate_res$direct_sensing
    indirect <- lactate_res$indirect_sensing

    if (is.null(direct) && is.null(indirect)) {
      stop("No sensing scores available")
    } else if (is.null(direct)) {
      sensing <- indirect
      pathway_label <- "indirect"
    } else if (is.null(indirect)) {
      sensing <- direct
      pathway_label <- "direct"
    } else {
      sensing <- direct + indirect
      pathway_label <- "combined"
    }
  }

  sensing <- sort(sensing, decreasing = TRUE)
  n_return <- min(top_n, length(sensing))

  data.frame(
    cell_type = names(sensing)[1:n_return],
    sensing_score = as.numeric(sensing[1:n_return]),
    pathway = pathway_label,
    rank = 1:n_return,
    stringsAsFactors = FALSE
  )
}


#' @title Get Lactate Signaling Summary
#'
#' @description
#' Returns a summary of lactate-mediated communication between cell types.
#'
#' @param object scMetaLink object with lactate signaling results
#' @param pathway Character. "direct", "indirect", or "combined". Default "combined".
#' @param top_n Integer. Number of top interactions to return. Default 10.
#'
#' @return data.frame summarizing top cell-cell communications
#' @export
#'
#' @examples
#' \donttest{
#' data(crc_example)
#' obj <- createScMetaLink(crc_expr, crc_meta, "cell_type")
#' obj <- inferLactateSignaling(obj)
#' getLactateSignalingSummary(obj)
#' }
getLactateSignalingSummary <- function(object, pathway = "combined", top_n = 10) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  if (is.null(object@parameters$lactate_signaling)) {
    stop("Lactate signaling not calculated. Run inferLactateSignaling() first.")
  }

  lactate_res <- object@parameters$lactate_signaling

  # Select communication matrix
  if (pathway == "direct") {
    comm <- lactate_res$direct_communication
    pvals <- lactate_res$pvalues$direct
  } else if (pathway == "indirect") {
    comm <- lactate_res$indirect_communication
    pvals <- lactate_res$pvalues$indirect
  } else {
    comm <- lactate_res$combined_communication
    pvals <- NULL  # Combined doesn't have separate p-values
  }

  if (is.null(comm)) {
    stop(paste("No communication scores for pathway:", pathway))
  }

  # Convert to data frame
  results <- data.frame(
    sender = rep(rownames(comm), ncol(comm)),
    receiver = rep(colnames(comm), each = nrow(comm)),
    communication_score = as.vector(comm),
    stringsAsFactors = FALSE
  )

  # Add p-values if available
  if (!is.null(pvals)) {
    results$pvalue <- as.vector(pvals)
  }

  # Filter non-zero and sort
  results <- results[results$communication_score > 0, ]
  results <- results[order(-results$communication_score), ]

  # Return top N
  n_return <- min(top_n, nrow(results))
  if (n_return > 0) {
    results <- results[1:n_return, ]
    results$rank <- 1:n_return
  }

  rownames(results) <- NULL
  results
}


#' @title Check Lactate Gene Availability
#'
#' @description
#' Checks which lactate signaling genes are available in the expression data.
#' Useful for quality control before running analysis.
#'
#' @param object scMetaLink object
#'
#' @return data.frame showing gene availability by category
#' @export
#'
#' @examples
#' \donttest{
#' data(crc_example)
#' obj <- createScMetaLink(crc_expr, crc_meta, "cell_type")
#' checkLactateGenes(obj)
#' }
checkLactateGenes <- function(object) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  gene_sets <- .get_lactate_gene_sets()
  expr_genes <- rownames(object@expression_data)

  results <- data.frame(
    category = character(),
    subcategory = character(),
    gene = character(),
    available = logical(),
    stringsAsFactors = FALSE
  )

  # Production - synthesis
  for (g in gene_sets$production$synthesis) {
    results <- rbind(results, data.frame(
      category = "production",
      subcategory = "synthesis",
      gene = g,
      available = g %in% expr_genes,
      stringsAsFactors = FALSE
    ))
  }

  # Production - export
  for (g in gene_sets$production$export) {
    results <- rbind(results, data.frame(
      category = "production",
      subcategory = "export",
      gene = g,
      available = g %in% expr_genes,
      stringsAsFactors = FALSE
    ))
  }

  # Degradation
  for (g in gene_sets$degradation$enzymes) {
    results <- rbind(results, data.frame(
      category = "degradation",
      subcategory = "enzymes",
      gene = g,
      available = g %in% expr_genes,
      stringsAsFactors = FALSE
    ))
  }

  # Direct sensing
  for (g in gene_sets$direct_sensing$receptor) {
    results <- rbind(results, data.frame(
      category = "direct_sensing",
      subcategory = "receptor",
      gene = g,
      available = g %in% expr_genes,
      stringsAsFactors = FALSE
    ))
  }

  # Indirect sensing
  for (g in gene_sets$indirect_sensing$proton_receptors) {
    results <- rbind(results, data.frame(
      category = "indirect_sensing",
      subcategory = "proton_receptors",
      gene = g,
      available = g %in% expr_genes,
      stringsAsFactors = FALSE
    ))
  }

  # Uptake
  for (g in gene_sets$uptake$import) {
    results <- rbind(results, data.frame(
      category = "uptake",
      subcategory = "import",
      gene = g,
      available = g %in% expr_genes,
      stringsAsFactors = FALSE
    ))
  }

  # Add summary
  cat("Lactate Gene Availability Summary:\n")
  cat("==================================\n")
  summary_df <- aggregate(available ~ category, data = results, FUN = function(x) {
    sprintf("%d/%d (%.0f%%)", sum(x), length(x), 100 * sum(x) / length(x))
  })
  print(summary_df, row.names = FALSE)
  cat("\n")

  invisible(results)
}


# =============================================================================
# Spatial Lactate Signaling Analysis
# =============================================================================

#' @title Infer Spatial Lactate Signaling
#'
#' @description
#' Infers lactate signaling in spatial transcriptomics data with distance-weighted
#' communication scores. Supports both direct (HCAR1) and indirect (proton-sensing
#' GPCRs) pathways with spatial context.
#'
#' @param object A spatial scMetaLink object (created with createScMetaLinkFromSpatial)
#' @param max_distance Numeric. Maximum communication distance in micrometers.
#'   Spot pairs beyond this distance are considered non-interacting. Default 200 um.
#' @param distance_decay Character. Distance decay function:
#'   \itemize{
#'     \item "gaussian": Gaussian decay exp(-d^2/(2*sigma^2)) (default)
#'     \item "exponential": Exponential decay exp(-d/sigma)
#'     \item "linear": Linear decay max(0, 1 - d/max_distance)
#'     \item "none": No distance weighting, use cell type level aggregation
#'   }
#' @param sigma Numeric. Decay parameter for gaussian/exponential (in um). Default 50 um.
#'   Literature suggests lactate has medium-range diffusion (~50-80 um).
#' @param include_direct Logical. Include direct lactate-HCAR1 signaling. Default TRUE.
#' @param include_indirect Logical. Include indirect lactate-H+-GPCR signaling. Default TRUE.
#' @param method Character. Scoring method: "combined", "mean", or "proportion". Default "combined".
#' @param comm_method Character. Communication score method: "geometric", "product", "harmonic".
#'   Default "geometric".
#' @param aggregate_by Character. Aggregation level:
#'   \itemize{
#'     \item "cell_type": Aggregate by cell type (default)
#'     \item "spot": Return spot-level scores (memory intensive)
#'     \item "both": Return both levels
#'   }
#' @param min_production Numeric. Minimum production score threshold. Default 0.
#' @param min_sensing Numeric. Minimum sensing score threshold. Default 0.
#' @param normalize Logical. Normalize scores. Default TRUE.
#' @param verbose Logical. Print progress messages. Default TRUE.
#'
#' @return Updated scMetaLink object with spatial_lactate_signaling results
#'   stored in the parameters slot, containing:
#'   \item{spot_production}{Spot-level production scores}
#'   \item{spot_direct_sensing}{Spot-level direct sensing scores}
#'   \item{spot_indirect_sensing}{Spot-level indirect sensing scores}
#'   \item{celltype_communication}{Cell type level communication (if aggregate_by includes "cell_type")}
#'   \item{spot_communication}{Spot level communication (if aggregate_by includes "spot")}
#'   \item{parameters}{Analysis parameters}
#'
#' @export
#'
#' @examples
#' \donttest{
#' data(st_expr)
#' data(st_meta)
#' data(st_scalefactors)
#'
#' # Create spatial object
#' obj <- createScMetaLinkFromSpatial(
#'   expression_data = st_expr,
#'   spatial_coords = st_meta[, c("x", "y")],
#'   cell_meta = st_meta,
#'   cell_type_column = "cell_type",
#'   scale_factors = st_scalefactors
#' )
#'
#' # Run spatial lactate analysis
#' obj <- inferSpatialLactateSignaling(obj)
#'
#' # Access results
#' spatial_results <- obj@parameters$spatial_lactate_signaling
#' }
inferSpatialLactateSignaling <- function(object,
                                          max_distance = 200,
                                          distance_decay = "gaussian",
                                          sigma = 50,
                                          include_direct = TRUE,
                                          include_indirect = TRUE,
                                          method = "combined",
                                          comm_method = "geometric",
                                          aggregate_by = "cell_type",
                                          min_production = 0,
                                          min_sensing = 0,
                                          normalize = TRUE,
                                          verbose = TRUE) {
  # Input validation
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  # Check if spatial data is available
  if (is.null(object@parameters$spatial_coords) || !isTRUE(object@parameters$is_spatial)) {
    stop("Object does not contain spatial information. Use createScMetaLinkFromSpatial() to create a spatial object.")
  }

  if (!include_direct && !include_indirect) {
    stop("At least one of include_direct or include_indirect must be TRUE")
  }

  if (!distance_decay %in% c("gaussian", "exponential", "linear", "none")) {
    stop("distance_decay must be one of: 'gaussian', 'exponential', 'linear', 'none'")
  }

  if (!aggregate_by %in% c("cell_type", "spot", "both")) {
    stop("aggregate_by must be one of: 'cell_type', 'spot', 'both'")
  }

  if (max_distance <= 0) {
    stop("max_distance must be positive")
  }

  if (sigma <= 0) {
    stop("sigma must be positive")
  }

  # Get data
  expr_data <- object@expression_data
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column
  spatial_coords <- object@parameters$spatial_coords
  scale_factors <- object@parameters$scale_factors

  # Convert coordinates to micrometers if scale_factors available
  pixels_per_um <- scale_factors$pixels_per_um
  if (is.null(pixels_per_um)) pixels_per_um <- 1

  coords_um <- data.frame(
    x = spatial_coords$x / pixels_per_um,
    y = spatial_coords$y / pixels_per_um
  )
  rownames(coords_um) <- rownames(spatial_coords)

  spots <- colnames(expr_data)
  n_spots <- length(spots)
  cell_types <- unique(cell_meta[[cell_type_col]])

  if (verbose) {
    message("Inferring spatial lactate signaling...")
    message(sprintf("  Spots: %d", n_spots))
    message(sprintf("  Cell types: %d", length(cell_types)))
    message(sprintf("  Distance decay: %s (sigma=%d um, max=%d um)", distance_decay, sigma, max_distance))
  }

  # Get gene sets
  gene_sets <- .get_lactate_gene_sets()

  # Get all genes needed
  all_genes <- unique(c(
    gene_sets$production$synthesis,
    gene_sets$production$export,
    gene_sets$degradation$enzymes,
    gene_sets$direct_sensing$receptor,
    gene_sets$indirect_sensing$proton_receptors,
    gene_sets$uptake$import
  ))

  available_genes <- intersect(all_genes, rownames(expr_data))
  if (verbose) {
    message(sprintf("  Available lactate genes: %d/%d", length(available_genes), length(all_genes)))
  }

  if (length(available_genes) == 0) {
    stop("No lactate-related genes found in expression data")
  }

  # Step 1: Calculate spot-level gene scores
  if (verbose) message("  Calculating spot-level expression scores...")

  # For spot-level, we use the expression directly
  spot_expr <- expr_data[available_genes, spots, drop = FALSE]

  # Normalize to 0-1 scale per gene
  spot_scores <- t(apply(spot_expr, 1, function(x) {
    if (max(x) > min(x)) {
      (x - min(x)) / (max(x) - min(x))
    } else {
      rep(0, length(x))
    }
  }))
  colnames(spot_scores) <- spots

  # Step 2: Calculate spot-level production scores
  if (verbose) message("  Calculating spot-level production scores...")

  synthesis_genes <- intersect(gene_sets$production$synthesis, rownames(spot_scores))
  export_genes <- intersect(gene_sets$production$export, rownames(spot_scores))
  degradation_genes <- intersect(gene_sets$degradation$enzymes, rownames(spot_scores))

  spot_production <- rep(0, n_spots)
  names(spot_production) <- spots

  for (i in seq_len(n_spots)) {
    spot <- spots[i]
    synth_score <- if (length(synthesis_genes) > 0) mean(spot_scores[synthesis_genes, spot]) else 0
    export_score <- if (length(export_genes) > 0) mean(spot_scores[export_genes, spot]) else 1
    degrad_score <- if (length(degradation_genes) > 0) mean(spot_scores[degradation_genes, spot]) else 0

    spot_production[spot] <- synth_score * export_score - 0.3 * degrad_score
  }
  spot_production[spot_production < 0] <- 0

  # Step 3: Calculate spot-level sensing scores
  spot_direct_sensing <- NULL
  spot_indirect_sensing <- NULL

  if (include_direct) {
    if (verbose) message("  Calculating spot-level direct sensing scores...")
    receptor_genes <- intersect(gene_sets$direct_sensing$receptor, rownames(spot_scores))
    uptake_genes <- intersect(gene_sets$uptake$import, rownames(spot_scores))

    spot_direct_sensing <- rep(0, n_spots)
    names(spot_direct_sensing) <- spots

    for (i in seq_len(n_spots)) {
      spot <- spots[i]
      receptor_score <- if (length(receptor_genes) > 0) mean(spot_scores[receptor_genes, spot]) else 0
      uptake_score <- if (length(uptake_genes) > 0) mean(spot_scores[uptake_genes, spot]) else 0
      spot_direct_sensing[spot] <- receptor_score + 0.5 * uptake_score
    }
  }

  if (include_indirect) {
    if (verbose) message("  Calculating spot-level indirect sensing scores...")
    proton_genes <- intersect(gene_sets$indirect_sensing$proton_receptors, rownames(spot_scores))
    weights <- gene_sets$indirect_sensing$weights

    spot_indirect_sensing <- rep(0, n_spots)
    names(spot_indirect_sensing) <- spots

    for (i in seq_len(n_spots)) {
      spot <- spots[i]
      if (length(proton_genes) > 0) {
        expr_vals <- spot_scores[proton_genes, spot]
        gene_weights <- weights[proton_genes]
        gene_weights[is.na(gene_weights)] <- 1.0
        spot_indirect_sensing[spot] <- sum(expr_vals * gene_weights) / sum(gene_weights)
      }
    }
  }

  # Normalize if requested
  if (normalize) {
    spot_production <- .normalize_lactate_scores(spot_production)
    if (!is.null(spot_direct_sensing)) {
      spot_direct_sensing <- .normalize_lactate_scores(spot_direct_sensing)
    }
    if (!is.null(spot_indirect_sensing)) {
      spot_indirect_sensing <- .normalize_lactate_scores(spot_indirect_sensing)
    }
  }

  # Step 4: Calculate distance matrix and weights
  if (verbose) message("  Computing spatial distance weights...")

  dist_matrix <- as.matrix(dist(coords_um))

  # Calculate distance weights
  if (distance_decay == "gaussian") {
    weight_matrix <- exp(-dist_matrix^2 / (2 * sigma^2))
  } else if (distance_decay == "exponential") {
    weight_matrix <- exp(-dist_matrix / sigma)
  } else if (distance_decay == "linear") {
    weight_matrix <- pmax(0, 1 - dist_matrix / max_distance)
  } else {
    weight_matrix <- matrix(1, n_spots, n_spots)
  }

  # Apply distance threshold
  weight_matrix[dist_matrix > max_distance] <- 0
  diag(weight_matrix) <- 0  # No self-communication

  # Step 5: Calculate spatial communication
  results <- list(
    spot_production = spot_production,
    spot_direct_sensing = spot_direct_sensing,
    spot_indirect_sensing = spot_indirect_sensing
  )

  # Cell type level aggregation
  if (aggregate_by %in% c("cell_type", "both")) {
    if (verbose) message("  Aggregating by cell type...")

    # Calculate cell type level scores
    ct_production <- tapply(spot_production, cell_meta[[cell_type_col]], mean)
    ct_direct <- if (!is.null(spot_direct_sensing)) {
      tapply(spot_direct_sensing, cell_meta[[cell_type_col]], mean)
    } else NULL
    ct_indirect <- if (!is.null(spot_indirect_sensing)) {
      tapply(spot_indirect_sensing, cell_meta[[cell_type_col]], mean)
    } else NULL

    # Calculate communication matrices
    direct_comm <- NULL
    indirect_comm <- NULL

    if (include_direct && !is.null(ct_direct)) {
      direct_comm <- .calc_lactate_communication(
        production_scores = ct_production,
        sensing_scores = ct_direct,
        method = comm_method,
        min_production = min_production,
        min_sensing = min_sensing
      )
    }

    if (include_indirect && !is.null(ct_indirect)) {
      indirect_comm <- .calc_lactate_communication(
        production_scores = ct_production,
        sensing_scores = ct_indirect,
        method = comm_method,
        min_production = min_production,
        min_sensing = min_sensing
      )
    }

    combined_comm <- NULL
    if (!is.null(direct_comm) && !is.null(indirect_comm)) {
      combined_comm <- direct_comm + indirect_comm
    } else if (!is.null(direct_comm)) {
      combined_comm <- direct_comm
    } else if (!is.null(indirect_comm)) {
      combined_comm <- indirect_comm
    }

    results$celltype_production <- ct_production
    results$celltype_direct_sensing <- ct_direct
    results$celltype_indirect_sensing <- ct_indirect
    results$celltype_direct_communication <- direct_comm
    results$celltype_indirect_communication <- indirect_comm
    results$celltype_combined_communication <- combined_comm
  }

  # Spot level communication (memory intensive)
  if (aggregate_by %in% c("spot", "both")) {
    if (verbose) message("  Computing spot-level communication (this may take a while)...")

    # For spot level, we compute weighted communication considering spatial proximity
    spot_direct_comm <- NULL
    spot_indirect_comm <- NULL

    if (include_direct && !is.null(spot_direct_sensing)) {
      spot_direct_comm <- matrix(0, n_spots, n_spots, dimnames = list(spots, spots))
      for (i in seq_len(n_spots)) {
        for (j in seq_len(n_spots)) {
          if (weight_matrix[i, j] > 0) {
            p <- spot_production[i]
            s <- spot_direct_sensing[j]
            if (p >= min_production && s >= min_sensing) {
              if (comm_method == "geometric") {
                spot_direct_comm[i, j] <- sqrt(p * s) * weight_matrix[i, j]
              } else if (comm_method == "product") {
                spot_direct_comm[i, j] <- p * s * weight_matrix[i, j]
              } else {
                spot_direct_comm[i, j] <- 2 * p * s / (p + s + 1e-10) * weight_matrix[i, j]
              }
            }
          }
        }
      }
    }

    if (include_indirect && !is.null(spot_indirect_sensing)) {
      spot_indirect_comm <- matrix(0, n_spots, n_spots, dimnames = list(spots, spots))
      for (i in seq_len(n_spots)) {
        for (j in seq_len(n_spots)) {
          if (weight_matrix[i, j] > 0) {
            p <- spot_production[i]
            s <- spot_indirect_sensing[j]
            if (p >= min_production && s >= min_sensing) {
              if (comm_method == "geometric") {
                spot_indirect_comm[i, j] <- sqrt(p * s) * weight_matrix[i, j]
              } else if (comm_method == "product") {
                spot_indirect_comm[i, j] <- p * s * weight_matrix[i, j]
              } else {
                spot_indirect_comm[i, j] <- 2 * p * s / (p + s + 1e-10) * weight_matrix[i, j]
              }
            }
          }
        }
      }
    }

    results$spot_direct_communication <- spot_direct_comm
    results$spot_indirect_communication <- spot_indirect_comm
    results$distance_weights <- weight_matrix
  }

  # Store parameters
  results$parameters <- list(
    max_distance = max_distance,
    distance_decay = distance_decay,
    sigma = sigma,
    include_direct = include_direct,
    include_indirect = include_indirect,
    method = method,
    comm_method = comm_method,
    aggregate_by = aggregate_by,
    min_production = min_production,
    min_sensing = min_sensing,
    normalize = normalize,
    pixels_per_um = pixels_per_um,
    timestamp = Sys.time()
  )

  object@parameters$spatial_lactate_signaling <- results

  if (verbose) {
    message("Done!")
    message(sprintf("  Spot production scores: %d spots", length(spot_production)))
    if (aggregate_by %in% c("cell_type", "both") && !is.null(results$celltype_combined_communication)) {
      message(sprintf("  Cell type communication: %d x %d",
                      nrow(results$celltype_combined_communication),
                      ncol(results$celltype_combined_communication)))
    }
  }

  object
}


#' @title Get Spatial Lactate Hotspots
#'
#' @description
#' Identifies spatial hotspots of lactate production and sensing.
#'
#' @param object scMetaLink object with spatial lactate signaling results
#' @param type Character. Type of hotspot: "production", "direct_sensing",
#'   "indirect_sensing", or "all". Default "all".
#' @param top_n Integer. Number of top spots to return. Default 20.
#'
#' @return data.frame with spot IDs, coordinates, scores, and cell types
#' @export
getSpatialLactateHotspots <- function(object, type = "all", top_n = 20) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  if (is.null(object@parameters$spatial_lactate_signaling)) {
    stop("Spatial lactate signaling not calculated. Run inferSpatialLactateSignaling() first.")
  }

  spatial_res <- object@parameters$spatial_lactate_signaling
  spatial_coords <- object@parameters$spatial_coords
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column

  results <- data.frame(
    spot = character(),
    x = numeric(),
    y = numeric(),
    cell_type = character(),
    score_type = character(),
    score = numeric(),
    stringsAsFactors = FALSE
  )

  if (type %in% c("production", "all")) {
    prod <- spatial_res$spot_production
    prod_sorted <- sort(prod, decreasing = TRUE)
    n <- min(top_n, length(prod_sorted))
    top_spots <- names(prod_sorted)[1:n]

    for (spot in top_spots) {
      results <- rbind(results, data.frame(
        spot = spot,
        x = spatial_coords[spot, "x"],
        y = spatial_coords[spot, "y"],
        cell_type = cell_meta[spot, cell_type_col],
        score_type = "production",
        score = prod[spot],
        stringsAsFactors = FALSE
      ))
    }
  }

  if (type %in% c("direct_sensing", "all") && !is.null(spatial_res$spot_direct_sensing)) {
    sens <- spatial_res$spot_direct_sensing
    sens_sorted <- sort(sens, decreasing = TRUE)
    n <- min(top_n, length(sens_sorted))
    top_spots <- names(sens_sorted)[1:n]

    for (spot in top_spots) {
      results <- rbind(results, data.frame(
        spot = spot,
        x = spatial_coords[spot, "x"],
        y = spatial_coords[spot, "y"],
        cell_type = cell_meta[spot, cell_type_col],
        score_type = "direct_sensing",
        score = sens[spot],
        stringsAsFactors = FALSE
      ))
    }
  }

  if (type %in% c("indirect_sensing", "all") && !is.null(spatial_res$spot_indirect_sensing)) {
    sens <- spatial_res$spot_indirect_sensing
    sens_sorted <- sort(sens, decreasing = TRUE)
    n <- min(top_n, length(sens_sorted))
    top_spots <- names(sens_sorted)[1:n]

    for (spot in top_spots) {
      results <- rbind(results, data.frame(
        spot = spot,
        x = spatial_coords[spot, "x"],
        y = spatial_coords[spot, "y"],
        cell_type = cell_meta[spot, cell_type_col],
        score_type = "indirect_sensing",
        score = sens[spot],
        stringsAsFactors = FALSE
      ))
    }
  }

  rownames(results) <- NULL
  results
}


# =============================================================================
# Lactate Signaling Visualization Functions
# =============================================================================

#' @title Plot Lactate Signaling Heatmap
#'
#' @description
#' Creates a heatmap visualization of lactate-mediated cell communication.
#' Shows sender cell types on rows and receiver cell types on columns.
#'
#' @param object scMetaLink object with lactate signaling results
#' @param pathway Character. Which pathway to plot: "direct", "indirect", or "combined".
#'   Default "combined".
#' @param cluster_rows Logical. Cluster rows. Default TRUE.
#' @param cluster_cols Logical. Cluster columns. Default TRUE.
#' @param show_values Logical. Show score values in cells. Default FALSE.
#' @param colors Character vector. Color palette. Default viridis.
#' @param title Character. Plot title. Default auto-generated.
#'
#' @return A ggplot2 object or ComplexHeatmap object (if available)
#' @export
#'
#' @examples
#' \donttest{
#' data(crc_example)
#' obj <- createScMetaLink(crc_expr, crc_meta, "cell_type")
#' obj <- inferLactateSignaling(obj)
#' plotLactateSignaling(obj, pathway = "indirect")
#' }
plotLactateSignaling <- function(object,
                                  pathway = "combined",
                                  cluster_rows = TRUE,
                                  cluster_cols = TRUE,
                                  show_values = FALSE,
                                  colors = NULL,
                                  title = NULL) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  if (is.null(object@parameters$lactate_signaling)) {
    stop("Lactate signaling not calculated. Run inferLactateSignaling() first.")
  }

  if (!pathway %in% c("direct", "indirect", "combined")) {
    stop("pathway must be 'direct', 'indirect', or 'combined'")
  }

  lactate_res <- object@parameters$lactate_signaling

  # Get communication matrix
  if (pathway == "direct") {
    mat <- lactate_res$direct_communication
    if (is.null(title)) title <- "Lactate Direct Signaling (HCAR1)"
  } else if (pathway == "indirect") {
    mat <- lactate_res$indirect_communication
    if (is.null(title)) title <- "Lactate Indirect Signaling (H+ -> GPCRs)"
  } else {
    mat <- lactate_res$combined_communication
    if (is.null(title)) title <- "Lactate Combined Signaling"
  }

  if (is.null(mat)) {
    stop(paste("No communication data for pathway:", pathway))
  }

  # Try to use ComplexHeatmap if available
  if (requireNamespace("ComplexHeatmap", quietly = TRUE) &&
      requireNamespace("circlize", quietly = TRUE)) {

    if (is.null(colors)) {
      if (requireNamespace("viridis", quietly = TRUE)) {
        colors <- viridis::viridis(100)
      } else {
        colors <- grDevices::colorRampPalette(c("white", "blue", "red"))(100)
      }
    }

    col_fun <- circlize::colorRamp2(
      seq(min(mat, na.rm = TRUE), max(mat, na.rm = TRUE), length.out = length(colors)),
      colors
    )

    ht <- ComplexHeatmap::Heatmap(
      mat,
      name = "Score",
      col = col_fun,
      column_title = title,
      row_title = "Sender (Producer)",
      column_title_side = "top",
      row_names_side = "left",
      column_names_side = "bottom",
      cluster_rows = cluster_rows,
      cluster_columns = cluster_cols,
      show_row_dend = cluster_rows,
      show_column_dend = cluster_cols,
      row_names_gp = grid::gpar(fontsize = 10),
      column_names_gp = grid::gpar(fontsize = 10),
      column_names_rot = 45,
      cell_fun = if (show_values) {
        function(j, i, x, y, width, height, fill) {
          grid::grid.text(sprintf("%.2f", mat[i, j]), x, y,
            gp = grid::gpar(fontsize = 8)
          )
        }
      } else {
        NULL
      },
      heatmap_legend_param = list(
        title = "Communication\nScore",
        title_position = "topcenter"
      )
    )

    return(ht)

  } else {
    # Fallback to base R heatmap
    if (is.null(colors)) {
      colors <- grDevices::colorRampPalette(c("white", "blue", "red"))(100)
    }

    # Order by clustering if requested
    if (cluster_rows) {
      row_order <- hclust(dist(mat))$order
      mat <- mat[row_order, ]
    }
    if (cluster_cols) {
      col_order <- hclust(dist(t(mat)))$order
      mat <- mat[, col_order]
    }

    graphics::par(mar = c(8, 8, 4, 2))
    graphics::image(t(mat)[, rev(seq_len(nrow(mat))), drop = FALSE], col = colors,
                    axes = FALSE, main = title)
    graphics::axis(1, at = seq(0, 1, length.out = ncol(mat)),
                   labels = colnames(mat), las = 2, cex.axis = 0.8)
    graphics::axis(2, at = seq(0, 1, length.out = nrow(mat)),
                   labels = rev(rownames(mat)), las = 2, cex.axis = 0.8)
    graphics::mtext("Receiver (Sensor)", side = 1, line = 6)
    graphics::mtext("Sender (Producer)", side = 2, line = 6)

    invisible(mat)
  }
}


#' @title Plot Lactate Pathway Comparison
#'
#' @description
#' Creates a side-by-side comparison of direct and indirect lactate signaling pathways.
#'
#' @param object scMetaLink object with lactate signaling results
#' @param show_production Logical. Also show production scores. Default TRUE.
#'
#' @return A ggplot2 object or base R plot
#' @export
plotLactatePathwayComparison <- function(object, show_production = TRUE) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  if (is.null(object@parameters$lactate_signaling)) {
    stop("Lactate signaling not calculated. Run inferLactateSignaling() first.")
  }

  lactate_res <- object@parameters$lactate_signaling

  # Check if both pathways are available
  has_direct <- !is.null(lactate_res$direct_communication)
  has_indirect <- !is.null(lactate_res$indirect_communication)

  if (!has_direct && !has_indirect) {
    stop("No communication data available")
  }

  # Prepare data
  cell_types <- names(lactate_res$production)

  plot_data <- data.frame(
    cell_type = cell_types,
    production = lactate_res$production[cell_types],
    stringsAsFactors = FALSE
  )

  if (has_direct) {
    plot_data$direct_sensing <- lactate_res$direct_sensing[cell_types]
  }

  if (has_indirect) {
    plot_data$indirect_sensing <- lactate_res$indirect_sensing[cell_types]
  }

  # Try ggplot2 if available
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    # Reshape for plotting
    plot_long <- data.frame(
      cell_type = character(),
      score_type = character(),
      score = numeric(),
      stringsAsFactors = FALSE
    )

    if (show_production) {
      plot_long <- rbind(plot_long, data.frame(
        cell_type = cell_types,
        score_type = "Production",
        score = plot_data$production,
        stringsAsFactors = FALSE
      ))
    }

    if (has_direct) {
      plot_long <- rbind(plot_long, data.frame(
        cell_type = cell_types,
        score_type = "Direct Sensing (HCAR1)",
        score = plot_data$direct_sensing,
        stringsAsFactors = FALSE
      ))
    }

    if (has_indirect) {
      plot_long <- rbind(plot_long, data.frame(
        cell_type = cell_types,
        score_type = "Indirect Sensing (GPCRs)",
        score = plot_data$indirect_sensing,
        stringsAsFactors = FALSE
      ))
    }

    # Order cell types by production
    ct_order <- names(sort(lactate_res$production, decreasing = TRUE))
    plot_long$cell_type <- factor(plot_long$cell_type, levels = ct_order)

    p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = .data$cell_type, y = .data$score, fill = .data$score_type)) +
      ggplot2::geom_bar(stat = "identity", position = "dodge") +
      ggplot2::theme_minimal() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(
        x = "Cell Type",
        y = "Score",
        fill = "Score Type",
        title = "Lactate Signaling: Production and Sensing Comparison"
      ) +
      ggplot2::scale_fill_brewer(palette = "Set2")

    return(p)

  } else {
    # Base R fallback
    n_types <- length(cell_types)
    n_scores <- 1 + has_direct + has_indirect

    graphics::par(mar = c(10, 4, 4, 8), xpd = TRUE)

    # Create empty plot
    plot(NULL, xlim = c(0.5, n_types + 0.5), ylim = c(0, 1.2),
         xlab = "", ylab = "Score", main = "Lactate Signaling Comparison",
         xaxt = "n")

    # Bar width
    bar_width <- 0.8 / n_scores

    # Plot production bars
    col_idx <- 1
    if (show_production) {
      for (i in seq_along(cell_types)) {
        graphics::rect(i - 0.4 + (col_idx - 1) * bar_width, 0,
                      i - 0.4 + col_idx * bar_width, plot_data$production[i],
                      col = "steelblue")
      }
      col_idx <- col_idx + 1
    }

    # Plot direct sensing
    if (has_direct) {
      for (i in seq_along(cell_types)) {
        graphics::rect(i - 0.4 + (col_idx - 1) * bar_width, 0,
                      i - 0.4 + col_idx * bar_width, plot_data$direct_sensing[i],
                      col = "coral")
      }
      col_idx <- col_idx + 1
    }

    # Plot indirect sensing
    if (has_indirect) {
      for (i in seq_along(cell_types)) {
        graphics::rect(i - 0.4 + (col_idx - 1) * bar_width, 0,
                      i - 0.4 + col_idx * bar_width, plot_data$indirect_sensing[i],
                      col = "forestgreen")
      }
    }

    # Add x-axis labels
    graphics::axis(1, at = seq_along(cell_types), labels = cell_types, las = 2, cex.axis = 0.8)

    # Add legend
    legend_labels <- c()
    legend_colors <- c()
    if (show_production) {
      legend_labels <- c(legend_labels, "Production")
      legend_colors <- c(legend_colors, "steelblue")
    }
    if (has_direct) {
      legend_labels <- c(legend_labels, "Direct (HCAR1)")
      legend_colors <- c(legend_colors, "coral")
    }
    if (has_indirect) {
      legend_labels <- c(legend_labels, "Indirect (GPCRs)")
      legend_colors <- c(legend_colors, "forestgreen")
    }

    graphics::legend("topright", inset = c(-0.2, 0),
                     legend = legend_labels, fill = legend_colors,
                     title = "Score Type", cex = 0.8)

    invisible(plot_data)
  }
}


#' @title Plot Spatial Lactate Signaling
#'
#' @description
#' Visualizes spatial distribution of lactate production and sensing scores.
#'
#' @param object scMetaLink object with spatial lactate signaling results
#' @param type Character. What to plot: "production", "direct_sensing",
#'   "indirect_sensing", or "all". Default "production".
#' @param point_size Numeric. Size of spots. Default 2.
#' @param title Character. Plot title. Default auto-generated.
#'
#' @return A ggplot2 object or base R plot
#' @export
plotSpatialLactate <- function(object,
                                type = "production",
                                point_size = 2,
                                title = NULL) {
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  if (is.null(object@parameters$spatial_lactate_signaling)) {
    stop("Spatial lactate signaling not calculated. Run inferSpatialLactateSignaling() first.")
  }

  spatial_res <- object@parameters$spatial_lactate_signaling
  spatial_coords <- object@parameters$spatial_coords

  if (!type %in% c("production", "direct_sensing", "indirect_sensing", "all")) {
    stop("type must be 'production', 'direct_sensing', 'indirect_sensing', or 'all'")
  }

  # Prepare plot data
  spots <- names(spatial_res$spot_production)
  plot_data <- data.frame(
    spot = spots,
    x = spatial_coords[spots, "x"],
    y = spatial_coords[spots, "y"],
    production = spatial_res$spot_production[spots],
    stringsAsFactors = FALSE
  )

  if (!is.null(spatial_res$spot_direct_sensing)) {
    plot_data$direct_sensing <- spatial_res$spot_direct_sensing[spots]
  }

  if (!is.null(spatial_res$spot_indirect_sensing)) {
    plot_data$indirect_sensing <- spatial_res$spot_indirect_sensing[spots]
  }

  # Use ggplot2 if available
  if (requireNamespace("ggplot2", quietly = TRUE)) {

    if (type == "all") {
      # Create faceted plot
      plot_long <- data.frame(
        spot = character(),
        x = numeric(),
        y = numeric(),
        score_type = character(),
        score = numeric(),
        stringsAsFactors = FALSE
      )

      plot_long <- rbind(plot_long, data.frame(
        spot = spots,
        x = plot_data$x,
        y = plot_data$y,
        score_type = "Production",
        score = plot_data$production,
        stringsAsFactors = FALSE
      ))

      if (!is.null(spatial_res$spot_direct_sensing)) {
        plot_long <- rbind(plot_long, data.frame(
          spot = spots,
          x = plot_data$x,
          y = plot_data$y,
          score_type = "Direct Sensing",
          score = plot_data$direct_sensing,
          stringsAsFactors = FALSE
        ))
      }

      if (!is.null(spatial_res$spot_indirect_sensing)) {
        plot_long <- rbind(plot_long, data.frame(
          spot = spots,
          x = plot_data$x,
          y = plot_data$y,
          score_type = "Indirect Sensing",
          score = plot_data$indirect_sensing,
          stringsAsFactors = FALSE
        ))
      }

      p <- ggplot2::ggplot(plot_long, ggplot2::aes(x = x, y = y, color = score)) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::facet_wrap(~score_type) +
        ggplot2::theme_minimal() +
        ggplot2::coord_fixed() +
        ggplot2::labs(
          title = if (is.null(title)) "Spatial Lactate Signaling" else title,
          color = "Score"
        )

    } else {
      # Single plot
      if (type == "production") {
        score_vals <- plot_data$production
        auto_title <- "Spatial Lactate Production"
      } else if (type == "direct_sensing") {
        if (is.null(plot_data$direct_sensing)) {
          stop("Direct sensing scores not available")
        }
        score_vals <- plot_data$direct_sensing
        auto_title <- "Spatial Direct Lactate Sensing (HCAR1)"
      } else {
        if (is.null(plot_data$indirect_sensing)) {
          stop("Indirect sensing scores not available")
        }
        score_vals <- plot_data$indirect_sensing
        auto_title <- "Spatial Indirect Lactate Sensing (GPCRs)"
      }

      plot_data$score <- score_vals

      p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = x, y = y, color = score)) +
        ggplot2::geom_point(size = point_size) +
        ggplot2::scale_color_viridis_c() +
        ggplot2::theme_minimal() +
        ggplot2::coord_fixed() +
        ggplot2::labs(
          title = if (is.null(title)) auto_title else title,
          color = "Score"
        )
    }

    return(p)

  } else {
    # Base R fallback
    if (type != "production" && type != "all") {
      if (type == "direct_sensing" && is.null(plot_data$direct_sensing)) {
        stop("Direct sensing scores not available")
      }
      if (type == "indirect_sensing" && is.null(plot_data$indirect_sensing)) {
        stop("Indirect sensing scores not available")
      }
    }

    score_vals <- if (type == "production") {
      plot_data$production
    } else if (type == "direct_sensing") {
      plot_data$direct_sensing
    } else {
      plot_data$indirect_sensing
    }

    colors <- grDevices::colorRampPalette(c("blue", "yellow", "red"))(100)
    score_scaled <- (score_vals - min(score_vals)) / (max(score_vals) - min(score_vals))
    point_colors <- colors[pmax(1, ceiling(score_scaled * 100))]

    graphics::par(mar = c(4, 4, 4, 6))
    plot(plot_data$x, plot_data$y, col = point_colors, pch = 19, cex = point_size,
         xlab = "X", ylab = "Y",
         main = if (is.null(title)) paste("Spatial Lactate", type) else title,
         asp = 1)

    invisible(plot_data)
  }
}
