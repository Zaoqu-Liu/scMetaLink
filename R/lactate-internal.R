# ============================================================================
# Lactate Signaling Internal Functions
# ============================================================================
#
# Internal helper functions for lactate-mediated cell communication analysis.
# These functions implement the core calculations for both direct (HCAR1) and
# indirect (proton-sensing GPCRs) lactate signaling pathways.
#
# Gene selection based on:
# - PLOS Biology 2024 (HCAR1 structure)
# - Nature Metabolism 2024 (GPR81 function)
# - Cell Press iScience 2019 (MCT direction)
# - Reactome pathway database (R-HSA-444731)
#
# ============================================================================

#' Get Curated Lactate Gene Sets
#'
#' @description
#' Returns literature-validated gene sets for lactate signaling analysis.
#' Gene selection is based on rigorous scientific review of 2024 literature.
#'
#' @return A named list containing gene sets for each category
#' @keywords internal
.get_lactate_gene_sets <- function() {
  list(
    # ================================================================
    # PRODUCTION: Lactate synthesis and export
    # ================================================================
    production = list(
      # Core synthesis enzymes (Pyruvate -> Lactate direction)
      # LDHA preferentially catalyzes pyruvate-to-lactate (Warburg effect)
      synthesis = c(
        "LDHA",
        "LDHC",
        "LDHAL6A",
        "LDHAL6B"
        # NOTE: LDHB excluded - preferentially catalyzes reverse reaction
      ),

      # Export transporters (lactate efflux)
      # MCT4 is the PRIMARY exporter in glycolytic/tumor cells
      export = c(
        "SLC16A3",
        "SLC16A1",
        "SLC16A7",
        "SLC16A8",
        "AQP9",
        "BSG"
      )
    ),

    # ================================================================
    # DEGRADATION: Lactate consumption
    # ================================================================
    degradation = list(
      # LDHB preferentially catalyzes lactate-to-pyruvate
      enzymes = c(
        "LDHB",
        "LDHD"
      )
    ),

    # ================================================================
    # DIRECT SENSING: Lactate receptor
    # ================================================================
    direct_sensing = list(
      # HCAR1 is the ONLY confirmed lactate GPCR
      # Confirmed by cryo-EM structure (PLOS Biology 2024)
      receptor = c(
        "HCAR1"
      )
    ),

    # ================================================================
    # INDIRECT SENSING: Proton-sensing GPCRs
    # ================================================================
    indirect_sensing = list(
      # Classic proton-sensing GPCRs (Reactome R-HSA-444731)
      # Sense pH via histidine protonation on extracellular domain
      proton_receptors = c(
        "GPR4",
        "GPR65",
        "GPR68",
        "GPR132"
      ),
      # Weights reflecting pH sensitivity strength
      # GPR132 has weak proton sensitivity, primarily senses oxidized fatty acids
      weights = c(
        GPR4 = 1.0,
        GPR65 = 1.0,
        GPR68 = 1.0,
        GPR132 = 0.5
      )
    ),

    # ================================================================
    # UPTAKE: Lactate import
    # ================================================================
    uptake = list(
      # MCT1 is the PRIMARY importer in oxidative cells
      import = c(
        "SLC16A1",
        "SLC16A7",
        "BSG"
      )
    )
  )
}


#' Calculate Cell Type Expression Profiles for Lactate Genes
#'
#' @description
#' Calculates mean expression and percentage of expressing cells for lactate-related
#' genes across cell types.
#'
#' @param expr_data Expression matrix (genes x cells)
#' @param cell_meta Cell metadata data.frame
#' @param cell_type_col Column name for cell type annotation
#' @param genes Character vector of genes to analyze
#' @param min_expression Minimum expression threshold
#' @param method Scoring method: "combined", "mean", or "proportion"
#'
#' @return Matrix of gene scores (genes x cell_types)
#' @keywords internal
.calc_lactate_gene_scores <- function(expr_data, cell_meta, cell_type_col,
                                       genes, min_expression = 0,
                                       method = "combined") {
  # Filter to available genes
  available_genes <- intersect(genes, rownames(expr_data))

  if (length(available_genes) == 0) {
    return(NULL)
  }

  cell_types <- unique(cell_meta[[cell_type_col]])
  n_genes <- length(available_genes)
  n_types <- length(cell_types)

  mean_expr <- matrix(0,
    nrow = n_genes, ncol = n_types,
    dimnames = list(available_genes, cell_types)
  )
  pct_expr <- mean_expr

  for (ct in cell_types) {
    cells <- rownames(cell_meta)[cell_meta[[cell_type_col]] == ct]
    if (length(cells) == 0) next

    ct_expr <- expr_data[available_genes, cells, drop = FALSE]

    # Calculate percentage expressing
    if (inherits(ct_expr, "dgCMatrix")) {
      pct_expr[, ct] <- Matrix::rowMeans(ct_expr > min_expression)
      mean_expr[, ct] <- Matrix::rowMeans(ct_expr)
    } else {
      pct_expr[, ct] <- rowMeans(ct_expr > min_expression, na.rm = TRUE)
      mean_expr[, ct] <- rowMeans(ct_expr, na.rm = TRUE)
    }
  }

  # Calculate final scores based on method
  if (method == "combined") {
    gene_scores <- mean_expr * pct_expr
  } else if (method == "mean") {
    gene_scores <- mean_expr
  } else {
    gene_scores <- pct_expr
  }

  gene_scores
}


#' Calculate Lactate Production Score
#'
#' @description
#' Calculates lactate production potential for each cell type based on:
#' Production = mean(Synthesis) * mean(Export) - alpha * mean(Degradation)
#'
#' @param gene_scores Matrix of gene scores (genes x cell_types)
#' @param gene_sets Lactate gene sets from .get_lactate_gene_sets()
#' @param consider_export Logical. Weight by export potential
#' @param consider_degradation Logical. Subtract degradation
#' @param degradation_weight Numeric. Weight for degradation subtraction (default 0.3)
#'
#' @return Named numeric vector of production scores per cell type
#' @keywords internal
.calc_lactate_production <- function(gene_scores,
                                      gene_sets,
                                      consider_export = TRUE,
                                      consider_degradation = TRUE,
                                      degradation_weight = 0.3) {
  cell_types <- colnames(gene_scores)
  production_scores <- rep(0, length(cell_types))
  names(production_scores) <- cell_types

  # Get available genes
  synthesis_genes <- intersect(gene_sets$production$synthesis, rownames(gene_scores))
  export_genes <- intersect(gene_sets$production$export, rownames(gene_scores))
  degradation_genes <- intersect(gene_sets$degradation$enzymes, rownames(gene_scores))

  for (ct in cell_types) {
    # Synthesis score
    if (length(synthesis_genes) > 0) {
      synthesis_score <- mean(gene_scores[synthesis_genes, ct], na.rm = TRUE)
    } else {
      synthesis_score <- 0
    }

    # Export score
    if (consider_export && length(export_genes) > 0) {
      export_score <- mean(gene_scores[export_genes, ct], na.rm = TRUE)
      production_scores[ct] <- synthesis_score * export_score
    } else {
      production_scores[ct] <- synthesis_score
    }

    # Subtract degradation
    if (consider_degradation && length(degradation_genes) > 0) {
      degradation_score <- mean(gene_scores[degradation_genes, ct], na.rm = TRUE)
      production_scores[ct] <- production_scores[ct] - degradation_weight * degradation_score
    }
  }

  # Ensure non-negative

  production_scores[production_scores < 0] <- 0

  production_scores
}


#' Calculate Direct Lactate Sensing Score (HCAR1)
#'
#' @description
#' Calculates direct lactate sensing capability based on HCAR1 expression,
#' optionally weighted by uptake transporter expression.
#'
#' DirectSensing = HCAR1_expression + beta * mean(Uptake_transporters)
#'
#' @param gene_scores Matrix of gene scores (genes x cell_types)
#' @param gene_sets Lactate gene sets from .get_lactate_gene_sets()
#' @param include_uptake Logical. Include MCT uptake transporters
#' @param uptake_weight Numeric. Weight for uptake contribution (default 0.5)
#'
#' @return Named numeric vector of direct sensing scores per cell type
#' @keywords internal
.calc_direct_sensing <- function(gene_scores,
                                  gene_sets,
                                  include_uptake = TRUE,
                                  uptake_weight = 0.5) {
  cell_types <- colnames(gene_scores)
  sensing_scores <- rep(0, length(cell_types))
  names(sensing_scores) <- cell_types

  # Get available genes
  receptor_genes <- intersect(gene_sets$direct_sensing$receptor, rownames(gene_scores))
  uptake_genes <- intersect(gene_sets$uptake$import, rownames(gene_scores))

  for (ct in cell_types) {
    # HCAR1 expression
    if (length(receptor_genes) > 0) {
      receptor_score <- mean(gene_scores[receptor_genes, ct], na.rm = TRUE)
    } else {
      receptor_score <- 0
    }

    sensing_scores[ct] <- receptor_score

    # Add uptake contribution
    if (include_uptake && length(uptake_genes) > 0) {
      uptake_score <- mean(gene_scores[uptake_genes, ct], na.rm = TRUE)
      sensing_scores[ct] <- sensing_scores[ct] + uptake_weight * uptake_score
    }
  }

  sensing_scores
}


#' Calculate Indirect Proton Sensing Score (GPR4/65/68/132)
#'
#' @description
#' Calculates indirect lactate sensing capability through proton-sensing GPCRs.
#' Uses weighted mean based on pH sensitivity strength of each receptor.
#'
#' @param gene_scores Matrix of gene scores (genes x cell_types)
#' @param gene_sets Lactate gene sets from .get_lactate_gene_sets()
#' @param use_weights Logical. Use receptor-specific weights
#'
#' @return Named numeric vector of indirect sensing scores per cell type
#' @keywords internal
.calc_proton_sensing <- function(gene_scores,
                                  gene_sets,
                                  use_weights = TRUE) {
  cell_types <- colnames(gene_scores)
  sensing_scores <- rep(0, length(cell_types))
  names(sensing_scores) <- cell_types

  # Get available genes
  receptor_genes <- intersect(gene_sets$indirect_sensing$proton_receptors,
    rownames(gene_scores))

  if (length(receptor_genes) == 0) {
    return(sensing_scores)
  }

  # Get weights
  if (use_weights) {
    weights <- gene_sets$indirect_sensing$weights[receptor_genes]
    # Handle missing weights
    weights[is.na(weights)] <- 1.0
  } else {
    weights <- rep(1.0, length(receptor_genes))
    names(weights) <- receptor_genes
  }

  for (ct in cell_types) {
    # Weighted mean of proton receptor expression
    expr_values <- gene_scores[receptor_genes, ct]
    sensing_scores[ct] <- sum(expr_values * weights, na.rm = TRUE) / sum(weights)
  }

  sensing_scores
}


#' Calculate Lactate Communication Scores
#'
#' @description
#' Combines production and sensing scores into cell-cell communication scores.
#'
#' @param production_scores Named numeric vector of production scores
#' @param sensing_scores Named numeric vector of sensing scores
#' @param method Character. Communication method: "geometric", "product", "harmonic"
#' @param min_production Numeric. Minimum production threshold
#' @param min_sensing Numeric. Minimum sensing threshold
#'
#' @return Matrix of communication scores (sender x receiver)
#' @keywords internal
.calc_lactate_communication <- function(production_scores,
                                         sensing_scores,
                                         method = "geometric",
                                         min_production = 0,
                                         min_sensing = 0) {
  senders <- names(production_scores)
  receivers <- names(sensing_scores)

  comm_scores <- matrix(0,
    nrow = length(senders),
    ncol = length(receivers),
    dimnames = list(sender = senders, receiver = receivers)
  )

  # Apply thresholds
  prod <- production_scores
  sens <- sensing_scores
  prod[prod < min_production] <- 0
  sens[sens < min_sensing] <- 0

  # Calculate communication scores
  for (s in senders) {
    for (r in receivers) {
      p <- prod[s]
      se <- sens[r]

      if (method == "geometric") {
        comm_scores[s, r] <- sqrt(p * se)
      } else if (method == "product") {
        comm_scores[s, r] <- p * se
      } else if (method == "harmonic") {
        if (p + se > 0) {
          comm_scores[s, r] <- 2 * p * se / (p + se)
        }
      }
    }
  }

  comm_scores
}


#' Normalize Scores
#'
#' @description
#' Z-score normalization followed by min-max scaling to 0-1.
#'
#' @param scores Numeric vector or matrix of scores
#'
#' @return Normalized scores in same structure as input
#' @keywords internal
.normalize_lactate_scores <- function(scores) {
  if (is.matrix(scores)) {
    # Normalize per row (metabolite)
    row_means <- rowMeans(scores, na.rm = TRUE)
    row_sds <- apply(scores, 1, sd, na.rm = TRUE)
    row_sds[row_sds == 0] <- 1
    scores <- (scores - row_means) / row_sds
  } else {
    # Normalize vector
    scores <- (scores - mean(scores, na.rm = TRUE)) / sd(scores, na.rm = TRUE)
  }

  # Scale to 0-1
  min_val <- min(scores, na.rm = TRUE)
  max_val <- max(scores, na.rm = TRUE)
  if (max_val > min_val) {
    scores <- (scores - min_val) / (max_val - min_val)
  }

  scores
}


#' Run Permutation Test for Lactate Communication
#'
#' @description
#' Performs permutation testing by shuffling cell type labels.
#'
#' @param object scMetaLink object
#' @param comm_scores Communication score matrix
#' @param pathway Character. "direct" or "indirect"
#' @param n_permutations Number of permutations
#' @param method Scoring method
#' @param min_production Minimum production threshold
#' @param min_sensing Minimum sensing threshold
#' @param include_uptake Include uptake in direct sensing
#' @param use_weights Use weights in indirect sensing
#' @param verbose Print progress
#'
#' @return Matrix of p-values
#' @keywords internal
.run_lactate_permutation <- function(object,
                                      comm_scores,
                                      pathway = "direct",
                                      n_permutations = 100,
                                      method = "combined",
                                      comm_method = "geometric",
                                      min_production = 0,
                                      min_sensing = 0,
                                      include_uptake = TRUE,
                                      use_weights = TRUE,
                                      verbose = TRUE) {
  expr_data <- object@expression_data
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column
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

  count_matrix <- matrix(0,
    nrow = nrow(comm_scores),
    ncol = ncol(comm_scores),
    dimnames = dimnames(comm_scores)
  )

  if (verbose) pb <- txtProgressBar(min = 0, max = n_permutations, style = 3)

  for (perm in seq_len(n_permutations)) {
    # Shuffle cell type labels
    shuffled_labels <- sample(cell_meta[[cell_type_col]])
    shuffled_meta <- cell_meta
    shuffled_meta[[cell_type_col]] <- shuffled_labels

    # Calculate gene scores
    perm_gene_scores <- .calc_lactate_gene_scores(
      expr_data = expr_data,
      cell_meta = shuffled_meta,
      cell_type_col = cell_type_col,
      genes = all_genes,
      method = method
    )

    if (is.null(perm_gene_scores)) next

    # Calculate production
    perm_prod <- .calc_lactate_production(
      gene_scores = perm_gene_scores,
      gene_sets = gene_sets
    )

    # Calculate sensing based on pathway
    if (pathway == "direct") {
      perm_sens <- .calc_direct_sensing(
        gene_scores = perm_gene_scores,
        gene_sets = gene_sets,
        include_uptake = include_uptake
      )
    } else {
      perm_sens <- .calc_proton_sensing(
        gene_scores = perm_gene_scores,
        gene_sets = gene_sets,
        use_weights = use_weights
      )
    }

    # Calculate communication
    perm_comm <- .calc_lactate_communication(
      production_scores = perm_prod,
      sensing_scores = perm_sens,
      method = comm_method,
      min_production = min_production,
      min_sensing = min_sensing
    )

    # Count how often permuted >= observed
    count_matrix <- count_matrix + (perm_comm >= comm_scores)

    if (verbose) setTxtProgressBar(pb, perm)
  }

  if (verbose) close(pb)

  # Calculate p-values
  pvalues <- (count_matrix + 1) / (n_permutations + 1)

  pvalues
}


#' Get Gene Contributions for Lactate Signaling
#'
#' @description
#' Returns the contribution of each gene to production and sensing scores.
#'
#' @param gene_scores Matrix of gene scores
#' @param gene_sets Lactate gene sets
#'
#' @return List with production and sensing gene contributions
#' @keywords internal
.get_gene_contributions <- function(gene_scores, gene_sets) {
  # Production genes
  prod_genes <- c(
    gene_sets$production$synthesis,
    gene_sets$production$export
  )
  prod_available <- intersect(prod_genes, rownames(gene_scores))

  if (length(prod_available) > 0) {
    prod_contrib <- data.frame(
      gene = prod_available,
      category = ifelse(prod_available %in% gene_sets$production$synthesis,
        "synthesis", "export"),
      t(gene_scores[prod_available, , drop = FALSE]),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  } else {
    prod_contrib <- data.frame()
  }

  # Sensing genes
  sens_genes <- c(
    gene_sets$direct_sensing$receptor,
    gene_sets$indirect_sensing$proton_receptors,
    gene_sets$uptake$import
  )
  sens_available <- intersect(sens_genes, rownames(gene_scores))

  if (length(sens_available) > 0) {
    sens_contrib <- data.frame(
      gene = sens_available,
      category = ifelse(sens_available %in% gene_sets$direct_sensing$receptor,
        "direct_receptor",
        ifelse(sens_available %in% gene_sets$indirect_sensing$proton_receptors,
          "proton_receptor", "uptake")),
      t(gene_scores[sens_available, , drop = FALSE]),
      stringsAsFactors = FALSE,
      row.names = NULL
    )
  } else {
    sens_contrib <- data.frame()
  }

  list(
    production = prod_contrib,
    sensing = sens_contrib
  )
}
