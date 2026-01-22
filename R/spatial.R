#' @title scMetaLink Spatial Transcriptomics Module
#' @description Functions for analyzing metabolite-mediated cell communication
#'   in spatial transcriptomics data.
#' @name spatial
#' @keywords internal
NULL

# =============================================================================
# Object Creation
# =============================================================================

#' @title Create scMetaLink Object from Spatial Data
#' @description Initialize a scMetaLink object from spatial transcriptomics data
#'   with spatial coordinate information.
#'
#' @param expression_data A matrix or dgCMatrix of normalized expression values (genes x spots)
#' @param spatial_coords A data.frame or matrix with spatial coordinates (spots x 2).
#'   Must have row names matching column names of expression_data.
#' @param cell_meta A data.frame containing spot metadata (e.g., cell type from deconvolution)
#' @param cell_type_column Character. Column name in cell_meta containing cell type annotations
#' @param scale_factors List. Optional scale factors for coordinate conversion.
#'   Should contain 'pixels_per_um' for distance calculations.
#'   **IMPORTANT**: For Visium data, coordinates are in pixels, not micrometers.
#'   Without correct scale_factors, distance-based parameters will be wrong.
#' @param min_cells Integer. Minimum number of spots per cell type (default: 5)
#'
#' @details
#' This function creates a scMetaLink object with spatial information stored
#' in additional slots. The spatial coordinates are used for distance-weighted
#' communication analysis.
#'
#' **Important: Cell Type Annotation for Visium Data**
#'
#' For 10x Visium data, each spot (55 micrometer diameter) typically contains 1-10 cells
#' of potentially different types. Therefore, cell type annotations should
#' ideally come from deconvolution methods such as:
#' \itemize{
#'   \item RCTD (spacexr package)
#'   \item cell2location
#'   \item SPOTlight
#'   \item CARD
#' }
#'
#' If using dominant cell type assignment per spot, be aware that this is a
#' simplification and may miss important cell type heterogeneity within spots.
#'
#' **Important: Coordinate Units**
#'
#' For 10x Visium data, coordinates are typically in pixels. You MUST provide
#' scale_factors with 'pixels_per_um' to convert to micrometers for biologically
#' meaningful distance calculations. Without this, sigma=50 will be interpreted
#' as 50 pixels instead of 50 micrometers.
#'
#' @return A scMetaLink object with spatial information
#' @export
#'
#' @examples
#' \donttest{
#' data(st_colon)
#'
#' obj <- createScMetaLinkFromSpatial(
#'   expression_data = st_expr,
#'   spatial_coords = st_meta[, c("x", "y")],
#'   cell_meta = st_meta,
#'   cell_type_column = "cell_type",
#'   scale_factors = st_scalefactors
#' )
#' }
createScMetaLinkFromSpatial <- function(expression_data,
                                        spatial_coords,
                                        cell_meta,
                                        cell_type_column = "cell_type",
                                        scale_factors = NULL,
                                        min_cells = 5) {
  # Input validation
  if (!is.matrix(expression_data) && !inherits(expression_data, "Matrix")) {
    stop("expression_data must be a matrix or sparse Matrix")
  }

  # Convert to dgCMatrix if needed (more efficient for operations)
  if (inherits(expression_data, "dgTMatrix")) {
    expression_data <- methods::as(expression_data, "CsparseMatrix")
  }

  if (!is.data.frame(cell_meta)) {
    stop("cell_meta must be a data.frame")
  }

  if (!cell_type_column %in% colnames(cell_meta)) {
    stop(paste("Column", cell_type_column, "not found in cell_meta"))
  }

  # Validate spatial coordinates
  if (is.null(spatial_coords)) {
    stop("spatial_coords must be provided for spatial analysis")
  }

  if (!is.data.frame(spatial_coords) && !is.matrix(spatial_coords)) {
    stop("spatial_coords must be a data.frame or matrix")
  }

  if (ncol(spatial_coords) < 2) {
    stop("spatial_coords must have at least 2 columns (x, y coordinates)")
  }

  # Ensure row names
  if (is.null(rownames(expression_data))) {
    stop("expression_data must have row names (gene symbols)")
  }

  if (is.null(colnames(expression_data))) {
    stop("expression_data must have column names (spot IDs)")
  }

  if (is.null(rownames(spatial_coords))) {
    stop("spatial_coords must have row names (spot IDs)")
  }

  if (is.null(rownames(cell_meta))) {
    stop("cell_meta must have row names (spot IDs)")
  }

  # Match spots across all data
  common_spots <- Reduce(intersect, list(
    colnames(expression_data),
    rownames(cell_meta),
    rownames(spatial_coords)
  ))

  if (length(common_spots) == 0) {
    stop("No matching spot IDs between expression_data, cell_meta, and spatial_coords")
  }

  # Filter to common spots
  expression_data <- expression_data[, common_spots, drop = FALSE]
  cell_meta <- cell_meta[common_spots, , drop = FALSE]
  spatial_coords <- as.data.frame(spatial_coords[common_spots, 1:2, drop = FALSE])
  colnames(spatial_coords) <- c("x", "y")

  # Filter cell types with minimum spots
  cell_type_counts <- table(cell_meta[[cell_type_column]])
  valid_types <- names(cell_type_counts)[cell_type_counts >= min_cells]

  if (length(valid_types) == 0) {
    stop(paste("No cell types with at least", min_cells, "spots"))
  }

  if (length(valid_types) < length(cell_type_counts)) {
    removed_types <- setdiff(names(cell_type_counts), valid_types)
    message(sprintf("Removed %d cell types with < %d spots: %s",
                    length(removed_types), min_cells,
                    paste(removed_types, collapse = ", ")))
  }

  keep_spots <- cell_meta[[cell_type_column]] %in% valid_types
  expression_data <- expression_data[, keep_spots, drop = FALSE]
  cell_meta <- cell_meta[keep_spots, , drop = FALSE]
  spatial_coords <- spatial_coords[keep_spots, , drop = FALSE]

  message(sprintf(
    "Created spatial scMetaLink object with %d genes, %d spots, %d cell types",
    nrow(expression_data), ncol(expression_data), length(valid_types)
  ))

  # Load database
  db <- .load_metalinksdb()

  # Process scale factors with intelligent detection
  if (is.null(scale_factors)) {
    # Detect if coordinates appear to be in pixels (large values)
    coord_range <- max(c(max(spatial_coords$x) - min(spatial_coords$x),
                         max(spatial_coords$y) - min(spatial_coords$y)))

    if (coord_range > 1000) {
      warning(paste0(
        "Coordinates appear to be in pixels (range: ", round(coord_range), "). ",
        "For Visium data, please provide scale_factors with 'pixels_per_um'. ",
        "Without correct scaling, distance parameters (sigma, threshold) will be ",
        "interpreted as pixels instead of micrometers, leading to incorrect results. ",
        "Assuming coordinates are in micrometers for now (pixels_per_um = 1)."
      ))
    } else {
      message("  No scale_factors provided, assuming coordinates in micrometers")
    }
    scale_factors <- list(pixels_per_um = 1)
  }

  # Create object using standard method then add spatial info
ObjectWithSpatialInfo <- new("scMetaLink",
    expression_data = expression_data,
    cell_meta = cell_meta,
    cell_type_column = cell_type_column,
    parameters = list(
      min_cells = min_cells,
      created_at = Sys.time(),
      is_spatial = TRUE,
      scale_factors = scale_factors
    ),
    database = db
  )

  # Store spatial information in parameters
  ObjectWithSpatialInfo@parameters$spatial_coords <- spatial_coords
  ObjectWithSpatialInfo@parameters$is_spatial <- TRUE

  ObjectWithSpatialInfo
}

#' @title Create scMetaLink from Seurat Spatial Object
#' @description Initialize a scMetaLink object from a Seurat object with spatial data
#'
#' @param seurat_obj A Seurat object with spatial assay and coordinates
#' @param cell_type_column Character. Column name in meta.data for cell type
#' @param assay Character. Assay to use (default: "Spatial")
#' @param slot Character. Slot to use (default: "data" for normalized data)
#' @param image Character. Name of the spatial image to use (default: first available)
#' @param min_cells Integer. Minimum spots per cell type
#'
#' @return A scMetaLink object with spatial information
#' @export
createScMetaLinkFromSeuratSpatial <- function(seurat_obj,
                                               cell_type_column = "cell_type",
                                               assay = "Spatial",
                                               slot = "data",
                                               image = NULL,
                                               min_cells = 5) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required. Please install it.")
  }

  # Extract expression data
  if (slot == "data") {
    expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = "data")
  } else {
    expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = "counts")
  }

  # Extract metadata
  cell_meta <- seurat_obj@meta.data

  # Extract spatial coordinates
  if (is.null(image)) {
    image <- names(seurat_obj@images)[1]
  }

  if (is.null(image) || !image %in% names(seurat_obj@images)) {
    stop("No spatial image found in Seurat object")
  }

  spatial_coords <- Seurat::GetTissueCoordinates(seurat_obj, image = image)

  # Get scale factors if available
  scale_factors <- list(pixels_per_um = 1)
  if (!is.null(seurat_obj@images[[image]]@scale.factors)) {
    sf <- seurat_obj@images[[image]]@scale.factors
    # Standard Visium spot is 55 um
    if (!is.null(sf$spot)) {
      scale_factors$pixels_per_um <- sf$spot / 55
      scale_factors$spot_diameter_um <- 55
    }
  }

  # Create scMetaLink object
  createScMetaLinkFromSpatial(
    expression_data = expr_data,
    spatial_coords = spatial_coords,
    cell_meta = cell_meta,
    cell_type_column = cell_type_column,
    scale_factors = scale_factors,
    min_cells = min_cells
  )
}

# =============================================================================
# Spatial Communication Analysis
# =============================================================================

#' @title Compute Spatial Communication
#' @description Calculate spatially-weighted metabolite-mediated communication.
#'   This function incorporates spatial distance between spots to weight
#'   communication strength.
#'
#' @param object A scMetaLink object with spatial information and production/sensing scores
#' @param method Character. Spatial weighting method:
#'   \itemize{
#'     \item "knn": K-nearest neighbors only (recommended for Visium, most honest
#'       given the resolution limitations)
#'     \item "gaussian": Gaussian decay exp(-d^2/2*sigma^2)
#'     \item "exponential": Exponential decay exp(-d/lambda)
#'     \item "linear": Linear decay max(0, 1 - d/d_max)
#'     \item "threshold": Binary cutoff (1 if d <= threshold, 0 otherwise)
#'   }
#' @param k_neighbors Integer. Number of nearest neighbors for knn method. Default: 6.
#'   For Visium hexagonal grid, 6 corresponds to immediate neighbors.
#'   **Note**: The weight matrix is symmetrized, so actual neighbor count per spot
#'   may exceed k (typically k to 2k).
#' @param symmetric Logical. Whether to symmetrize the KNN weight matrix (default: TRUE).
#'   If TRUE, when spot A is a neighbor of B, B is also considered a neighbor of A.
#'   This ensures bidirectional communication potential.
#' @param distance_threshold Numeric. Maximum distance for communication in micrometers.
#'   Spot pairs beyond this distance are considered non-interacting. Default: 150 um.
#'   Note: Most metabolites have effective diffusion ranges of 10-200 um in tissue.
#'   For Visium data (100 um spot spacing), values >200 um are rarely meaningful.
#' @param sigma Numeric. Sigma parameter for Gaussian decay (in um). Default: 50 um.
#'   Represents the characteristic decay distance. Literature suggests:
#'   \itemize{
#'     \item Fast-turnover metabolites (adenosine, ATP): 20-30 um
#'     \item Medium-range metabolites (lactate): 50-80 um
#'     \item Stable metabolites (amino acids): 80-120 um
#'   }
#' @param lambda Numeric. Lambda parameter for exponential decay (in um). Default: 50 um.
#' @param comm_method Character. Communication score method: "geometric", "product", "harmonic"
#' @param min_production Numeric. Minimum production score threshold (0-1).
#'   Cell types with production score below this are considered non-producers.
#' @param min_sensing Numeric. Minimum sensing score threshold (0-1).
#'   Cell types with sensing score below this are considered non-sensors.
#' @param analysis_level Character. Level of analysis:
#'   \itemize{
#'     \item "region": Aggregate by cell type (default, recommended).
#'       Uses the full inferProduction/inferSensing logic.
#'     \item "spot": Compute spot-to-spot communication.
#'       **WARNING**: This is a simplified implementation for exploratory analysis.
#'       It does NOT include: degradation adjustment, secretion weighting,
#'       Hill transformation, or full normalization. Use for hotspot identification only.
#'   }
#' @param n_permutations Integer. Number of permutations for significance testing.
#'   Default: 1000. For publication, recommend >= 1000 permutations.
#'   The permutation test shuffles cell type labels while preserving spatial
#'   structure, then **recalculates production and sensing scores** from scratch
#'   (consistent with non-spatial version), which is the scientifically correct
#'   null model for testing whether communication patterns are associated with
#'   cell type identity.
#' @param n_cores Integer. Number of cores for parallel computing
#' @param seed Integer. Random seed for reproducibility
#' @param verbose Logical. Print progress messages
#'
#' @details
#' **Important Notes on Spatial Resolution:**
#'
#' For Visium data (55 um spots, ~100 um spacing), each spot contains 1-10 cells.
#' This means:
#' \itemize{
#'   \item Cell type annotations should ideally come from deconvolution methods
#'     (e.g., RCTD, cell2location, SPOTlight)
#'   \item The "knn" method (k=6) is recommended as it only considers immediate
#'     neighbors, which is most honest given the resolution limitations
#'   \item Distance-weighted methods (gaussian, exponential) may provide
#'     "false precision" when spot spacing is similar to metabolite diffusion distance
#' }
#'
#' **Permutation Test Design:**
#'
#' The permutation test shuffles cell type labels (not spatial positions) and
#' **recalculates production/sensing scores from expression data**. This is
#' consistent with the non-spatial version and tests whether the observed
#' communication pattern is associated with cell type identity beyond what
#' would be expected by chance.
#'
#' **Spatial communication is modeled as:**
#' \deqn{C_{spatial}(i -> j, m) = MPP(m, i) * MSC(m, j) * w(d_{ij})}
#'
#' where w(d) is the spatial weight function:
#' \itemize{
#'   \item Gaussian: w(d) = exp(-d^2 / 2*sigma^2)
#'   \item Exponential: w(d) = exp(-d / lambda)
#'   \item Linear: w(d) = max(0, 1 - d/d_max)
#'   \item Threshold: w(d) = I(d <= d_threshold)
#' }
#'
#' @return Updated scMetaLink object with spatial_communication slot
#' @export
#'
#' @examples
#' \donttest{
#' data(st_colon)
#'
#' # Create object and run analysis
#' obj <- createScMetaLinkFromSpatial(st_expr, st_meta[,c("x","y")],
#'                                    st_meta, "cell_type", st_scalefactors)
#' obj <- inferProduction(obj)
#' obj <- inferSensing(obj)
#'
#' # Compute spatial communication (knn method recommended for Visium)
#' obj <- computeSpatialCommunication(obj, method = "knn", k_neighbors = 6)
#'
#' # Alternative: Gaussian decay with conservative parameters
#' obj <- computeSpatialCommunication(obj, method = "gaussian",
#'                                    sigma = 50,  # 50 um decay
#'                                    distance_threshold = 150)  # max 150 um
#' }
computeSpatialCommunication <- function(object,
                                        method = "knn",
                                        k_neighbors = 6,
                                        symmetric = TRUE,
                                        distance_threshold = 150,
                                        sigma = 50,
                                        lambda = 50,
                                        comm_method = "geometric",
                                        min_production = 0.1,
                                        min_sensing = 0.1,
                                        analysis_level = "region",
                                        n_permutations = 1000,
                                        n_cores = 1,
                                        seed = 42,
                                        verbose = TRUE) {
  # Validation
  if (!inherits(object, "scMetaLink")) {
    stop("object must be a scMetaLink object")
  }

  if (!isTRUE(object@parameters$is_spatial)) {
    stop("Object does not contain spatial information. Use createScMetaLinkFromSpatial()")
  }

  if (is.null(object@production_scores)) {
    stop("Run inferProduction() first")
  }

  if (is.null(object@sensing_scores)) {
    stop("Run inferSensing() first")
  }

  if (!method %in% c("knn", "gaussian", "exponential", "linear", "threshold")) {
    stop("method must be one of: 'knn', 'gaussian', 'exponential', 'linear', 'threshold'")
  }

  if (!comm_method %in% c("geometric", "product", "harmonic")) {
    stop("comm_method must be one of: 'geometric', 'product', 'harmonic'")
  }

  if (!analysis_level %in% c("region", "spot")) {
    stop("analysis_level must be 'region' or 'spot'")
  }

  # Scientific validity warnings
  if (distance_threshold > 300 && verbose) {
    warning(paste0(
      "distance_threshold = ", distance_threshold, " um is larger than the effective ",
      "diffusion range of most metabolites (typically 10-200 um). ",
      "Consider using a smaller threshold for biological relevance."
    ))
  }

  if (sigma > 100 && method == "gaussian" && verbose) {
    warning(paste0(
      "sigma = ", sigma, " um may be too large. Literature suggests: ",
      "adenosine/ATP: 20-30 um, lactate: 50-80 um, amino acids: 80-120 um."
    ))
  }

  if (n_permutations > 0 && n_permutations < 1000 && verbose) {
    message(paste0(
      "  Note: n_permutations = ", n_permutations, ". ",
      "For publication-quality results, recommend >= 1000 permutations."
    ))
  }

  set.seed(seed)

  # Get spatial coordinates
  spatial_coords <- object@parameters$spatial_coords
  scale_factors <- object@parameters$scale_factors
  pixels_per_um <- scale_factors$pixels_per_um %||% 1

  # Convert coordinates to micrometers
  coords_um <- spatial_coords / pixels_per_um
  n_spots <- nrow(coords_um)

  if (verbose) {
    message(sprintf("Computing spatial communication (%s level)...", analysis_level))
    message(sprintf("  Method: %s, threshold: %.0f um, %d spots", method, distance_threshold, n_spots))
  }

  # Memory warning for large datasets
  if (n_spots > 10000 && method != "knn" && verbose) {
    estimated_memory_gb <- n_spots^2 * 8 / 1e9
    warning(paste0(
      "Dataset has ", n_spots, " spots. Distance matrix will require ~",
      round(estimated_memory_gb, 1), " GB memory. ",
      "Consider using method='knn' for large datasets to reduce memory usage."
    ))
  }

  # Compute distance matrix
  if (verbose) message("  Computing distance matrix...")
  dist_matrix <- as.matrix(dist(coords_um))

  # Compute spatial weights
  if (verbose) message("  Computing spatial weights...")
  weight_matrix <- .compute_spatial_weights(
    dist_matrix = dist_matrix,
    method = method,
    threshold = distance_threshold,
    sigma = sigma,
    lambda = lambda,
    k = k_neighbors,
    symmetric = symmetric
  )

  # Get production and sensing scores
  prod_scores <- object@production_scores
  sens_scores <- object@sensing_scores
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column

  # Get common metabolites
  common_mets <- intersect(rownames(prod_scores), rownames(sens_scores))
  prod_scores <- prod_scores[common_mets, , drop = FALSE]
  sens_scores <- sens_scores[common_mets, , drop = FALSE]

  cell_types <- colnames(prod_scores)
  n_types <- length(cell_types)
  n_mets <- length(common_mets)

  if (verbose) message(sprintf("  Metabolites: %d, Cell types: %d", n_mets, n_types))

  # Compute communication based on analysis level
  if (analysis_level == "region") {
    # Region-level: aggregate by cell type with spatial weighting
    result <- .compute_spatial_communication_region(
      prod_scores = prod_scores,
      sens_scores = sens_scores,
      weight_matrix = weight_matrix,
      cell_meta = cell_meta,
      cell_type_col = cell_type_col,
      cell_types = cell_types,
      common_mets = common_mets,
      comm_method = comm_method,
      min_production = min_production,
      min_sensing = min_sensing,
      verbose = verbose
    )

    object@communication_scores <- result$comm_scores

    # Run spatial permutation test (FIXED: now recalculates production/sensing)
    if (n_permutations > 0) {
      if (verbose) message(sprintf("  Running %d spatial permutations...", n_permutations))

      pvalues <- .run_spatial_permutation_consistent(
        object = object,
        comm_scores = result$comm_scores,
        weight_matrix = weight_matrix,
        common_mets = common_mets,
        n_permutations = n_permutations,
        n_cores = n_cores,
        comm_method = comm_method,
        min_production = min_production,
        min_sensing = min_sensing,
        verbose = verbose
      )

      object@communication_pvalues <- pvalues
    }

  } else {
    # Spot-level analysis
    if (verbose) {
      message("  Computing spot-level communication...")
      message("  WARNING: Spot-level analysis is simplified and does not include")
      message("           degradation adjustment, secretion weighting, or Hill transformation.")
      message("           Use for exploratory analysis and hotspot identification only.")
    }

    result <- .compute_spatial_communication_spot(
      object = object,
      weight_matrix = weight_matrix,
      common_mets = common_mets,
      comm_method = comm_method,
      min_production = min_production,
      min_sensing = min_sensing,
      verbose = verbose
    )

    # Store in parameters since spot-level is a different structure
    object@parameters$spot_communication <- result
  }

  # Store parameters
  object@parameters$spatial_communication <- list(
    method = method,
    k_neighbors = k_neighbors,
    symmetric = symmetric,
    distance_threshold = distance_threshold,
    sigma = sigma,
    lambda = lambda,
    comm_method = comm_method,
    analysis_level = analysis_level,
    n_permutations = n_permutations,
    seed = seed
  )

  if (verbose) message("Done!")
  object
}

# =============================================================================
# Internal Functions
# =============================================================================

#' Compute Spatial Weight Matrix
#' @param dist_matrix Distance matrix (spots x spots)
#' @param method Spatial weighting method
#' @param threshold Distance threshold in micrometers
#' @param sigma Sigma parameter for Gaussian decay
#' @param lambda Lambda parameter for exponential decay
#' @param k Number of nearest neighbors for knn method
#' @param symmetric Whether to symmetrize the KNN weight matrix
#' @keywords internal
.compute_spatial_weights <- function(dist_matrix, method, threshold, sigma, lambda,
                                     k = 6, symmetric = TRUE) {
  n <- nrow(dist_matrix)
  weight_matrix <- matrix(0, n, n)

  if (method == "knn") {
    # K-nearest neighbors: most honest for low-resolution data like Visium
    # For each spot, only consider k nearest neighbors
    for (i in seq_len(n)) {
      dists <- dist_matrix[i, ]
      dists[i] <- Inf  # Exclude self
      nn_idx <- order(dists)[1:min(k, n - 1)]
      weight_matrix[i, nn_idx] <- 1
    }
    # Optionally make symmetric (if i is neighbor of j, j is also neighbor of i)
    if (symmetric) {
      weight_matrix <- pmax(weight_matrix, t(weight_matrix))
    }
  } else if (method == "gaussian") {
    # Gaussian decay: exp(-d^2 / 2*sigma^2)
    weight_matrix <- exp(-dist_matrix^2 / (2 * sigma^2))
  } else if (method == "exponential") {
    # Exponential decay: exp(-d / lambda)
    weight_matrix <- exp(-dist_matrix / lambda)
  } else if (method == "linear") {
    # Linear decay: max(0, 1 - d/threshold)
    weight_matrix <- pmax(0, 1 - dist_matrix / threshold)
  } else if (method == "threshold") {
    # Binary threshold
    weight_matrix <- (dist_matrix <= threshold) * 1.0
  }

  # Apply distance threshold (zero out very distant pairs) for non-knn methods
  if (method != "knn") {
    weight_matrix[dist_matrix > threshold] <- 0
  }

  # Zero diagonal (no self-communication for spots)
  diag(weight_matrix) <- 0

  weight_matrix
}

#' Compute Region-Level Spatial Communication
#' @keywords internal
.compute_spatial_communication_region <- function(prod_scores, sens_scores, weight_matrix,
                                                  cell_meta, cell_type_col, cell_types,
                                                  common_mets, comm_method,
                                                  min_production, min_sensing, verbose) {
  n_types <- length(cell_types)
  n_mets <- length(common_mets)

  # Initialize output array
  comm_scores <- array(0, dim = c(n_types, n_types, n_mets))
  dimnames(comm_scores) <- list(sender = cell_types, receiver = cell_types, metabolite = common_mets)

  # Get spot indices for each cell type
  type_spots <- lapply(cell_types, function(ct) {
    which(cell_meta[[cell_type_col]] == ct)
  })
  names(type_spots) <- cell_types

  # Compute spatial-weighted communication for each metabolite
  for (m_idx in seq_len(n_mets)) {
    met <- common_mets[m_idx]
    p <- prod_scores[met, ]
    s <- sens_scores[met, ]

    # Apply thresholds (absolute thresholds, consistent with non-spatial version)
    p[p < min_production] <- 0
    s[s < min_sensing] <- 0

    for (i in seq_len(n_types)) {
      sender_type <- cell_types[i]
      sender_spots <- type_spots[[sender_type]]
      sender_prod <- p[sender_type]

      if (sender_prod == 0 || length(sender_spots) == 0) next

      for (j in seq_len(n_types)) {
        receiver_type <- cell_types[j]
        receiver_spots <- type_spots[[receiver_type]]
        receiver_sens <- s[receiver_type]

        if (receiver_sens == 0 || length(receiver_spots) == 0) next

        # Get spatial weights between sender and receiver spots
        spatial_weights <- weight_matrix[sender_spots, receiver_spots, drop = FALSE]

        # Aggregate spatial weight (mean of non-zero weights)
        nonzero_weights <- spatial_weights[spatial_weights > 0]
        if (length(nonzero_weights) == 0) next
        mean_weight <- mean(nonzero_weights)

        # Compute communication score
        if (comm_method == "geometric") {
          score <- sqrt(sender_prod * receiver_sens) * mean_weight
        } else if (comm_method == "product") {
          score <- sender_prod * receiver_sens * mean_weight
        } else if (comm_method == "harmonic") {
          score <- 2 * sender_prod * receiver_sens / (sender_prod + receiver_sens + 1e-10) * mean_weight
        }

        comm_scores[i, j, m_idx] <- score
      }
    }
  }

  list(comm_scores = comm_scores)
}

#' Compute Spot-Level Spatial Communication (Simplified)
#'
#' @description Simplified spot-level communication for exploratory analysis.
#'   This does NOT include: degradation adjustment, secretion weighting,
#'   Hill transformation, or full normalization.
#'
#' @keywords internal
.compute_spatial_communication_spot <- function(object, weight_matrix, common_mets,
                                                comm_method, min_production, min_sensing,
                                                verbose) {
  expr_data <- object@expression_data
  db <- object@database
  spots <- colnames(expr_data)
  n_spots <- length(spots)
  n_mets <- length(common_mets)

  # Get production enzymes
  prod_enzymes <- .get_production_enzymes(db)
  prod_enzymes <- prod_enzymes[prod_enzymes$hmdb %in% common_mets, ]
  prod_enzymes <- prod_enzymes[prod_enzymes$gene_symbol %in% rownames(expr_data), ]

  # Get sensing proteins
  sensing_proteins <- .get_sensing_proteins(db)
  sensing_proteins <- sensing_proteins[sensing_proteins$hmdb %in% common_mets, ]
  sensing_proteins <- sensing_proteins[sensing_proteins$gene_symbol %in% rownames(expr_data), ]

  # For each metabolite, compute spot-level communication
  # Store as list of sparse matrices to save memory
  spot_comm <- list()

  if (verbose) pb <- txtProgressBar(min = 0, max = n_mets, style = 3)

  for (m_idx in seq_len(n_mets)) {
    met <- common_mets[m_idx]

    # Get production genes for this metabolite
    prod_genes <- prod_enzymes$gene_symbol[prod_enzymes$hmdb == met]
    # Get sensing genes for this metabolite
    sens_genes <- sensing_proteins$gene_symbol[sensing_proteins$hmdb == met]

    if (length(prod_genes) == 0 || length(sens_genes) == 0) {
      spot_comm[[met]] <- NULL
      if (verbose) setTxtProgressBar(pb, m_idx)
      next
    }

    # Compute spot-level production (mean of enzyme expression)
    prod_genes <- intersect(prod_genes, rownames(expr_data))
    if (length(prod_genes) == 0) {
      spot_comm[[met]] <- NULL
      if (verbose) setTxtProgressBar(pb, m_idx)
      next
    }
    spot_prod <- colMeans(as.matrix(expr_data[prod_genes, , drop = FALSE]))

    # Compute spot-level sensing
    sens_genes <- intersect(sens_genes, rownames(expr_data))
    if (length(sens_genes) == 0) {
      spot_comm[[met]] <- NULL
      if (verbose) setTxtProgressBar(pb, m_idx)
      next
    }
    spot_sens <- colMeans(as.matrix(expr_data[sens_genes, , drop = FALSE]))

    # Apply thresholds (FIXED: use absolute thresholds consistent with region-level)
    # Scale to 0-1 first for threshold comparison
    prod_max <- max(spot_prod)
    sens_max <- max(spot_sens)
    if (prod_max > 0) {
      spot_prod_scaled <- spot_prod / prod_max
      spot_prod[spot_prod_scaled < min_production] <- 0
    }
    if (sens_max > 0) {
      spot_sens_scaled <- spot_sens / sens_max
      spot_sens[spot_sens_scaled < min_sensing] <- 0
    }

    # Compute pairwise communication with spatial weighting
    if (comm_method == "geometric") {
      comm_mat <- sqrt(outer(spot_prod, spot_sens)) * weight_matrix
    } else if (comm_method == "product") {
      comm_mat <- outer(spot_prod, spot_sens) * weight_matrix
    } else {
      ps <- outer(spot_prod, spot_sens)
      sum_ps <- outer(spot_prod, rep(1, n_spots)) + outer(rep(1, n_spots), spot_sens)
      comm_mat <- 2 * ps / (sum_ps + 1e-10) * weight_matrix
    }

    # Convert to sparse matrix to save memory
    comm_mat[comm_mat < 1e-6] <- 0
    spot_comm[[met]] <- Matrix::Matrix(comm_mat, sparse = TRUE)

    if (verbose) setTxtProgressBar(pb, m_idx)
  }

  if (verbose) close(pb)

  spot_comm
}

#' Spatial Permutation Test - Consistent with Non-Spatial Version
#'
#' @description Permutation test for spatial communication significance.
#'   This version RECALCULATES production/sensing scores for each permutation,
#'   consistent with the non-spatial version in communication.R.
#'
#' @details
#' The null hypothesis is: "The observed communication pattern between cell
#' types is not associated with the actual cell type labels."
#'
#' This is tested by:
#' 1. Shuffling cell type labels
#' 2. Recalculating expression profiles for each cell type
#' 3. Recalculating production and sensing scores (with all adjustments)
#' 4. Computing spatial communication with the new scores
#' 5. Keeping spatial weight matrix UNCHANGED (preserving spatial structure)
#'
#' This approach is scientifically correct and consistent with the non-spatial
#' permutation test in communication.R.
#'
#' @keywords internal
.run_spatial_permutation_consistent <- function(object, comm_scores, weight_matrix,
                                                common_mets, n_permutations, n_cores,
                                                comm_method, min_production,
                                                min_sensing, verbose) {
  expr_data <- object@expression_data
  cell_meta <- object@cell_meta
  cell_type_col <- object@cell_type_column
  db <- object@database
  params <- object@parameters

  cell_types <- dimnames(comm_scores)[[1]]
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
  # PERMUTATION FUNCTION
  # ============================================================
  run_single_perm <- function(perm_id) {
    # Shuffle cell type labels (NOT spatial positions)
    shuffled_labels <- sample(cell_meta[[cell_type_col]])
    shuffled_meta <- cell_meta
    shuffled_meta[[cell_type_col]] <- shuffled_labels

    # Recalculate expression profiles with shuffled labels
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

    # ---- PRODUCTION SCORES (same logic as inferProduction) ----
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

    # ---- SENSING SCORES (same logic as inferSensing) ----
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

    # ---- SPATIAL COMMUNICATION ----
    # Use same weight_matrix (spatial structure preserved)
    # But with shuffled cell type assignments
    perm_result <- .compute_spatial_communication_region(
      prod_scores = perm_prod,
      sens_scores = perm_sens,
      weight_matrix = weight_matrix,  # UNCHANGED - preserves spatial structure
      cell_meta = shuffled_meta,      # SHUFFLED cell type labels
      cell_type_col = cell_type_col,
      cell_types = cell_types,
      common_mets = common_mets,
      comm_method = comm_method,
      min_production = min_production,
      min_sensing = min_sensing,
      verbose = FALSE
    )

    # Return comparison result
    perm_result$comm_scores >= comm_scores
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

  # Calculate p-values using standard permutation formula
  # Adding 1 to both numerator and denominator avoids p=0
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

# =============================================================================
# Utility Functions
# =============================================================================

#' @title Get Spatial Distance Statistics
#' @description Calculate statistics about spatial distances in the dataset
#'
#' @param object A scMetaLink object with spatial information
#' @return A list with distance statistics
#' @export
getSpatialDistanceStats <- function(object) {
  if (!isTRUE(object@parameters$is_spatial)) {
    stop("Object does not contain spatial information")
  }

  coords <- object@parameters$spatial_coords
  scale_factors <- object@parameters$scale_factors
  pixels_per_um <- scale_factors$pixels_per_um %||% 1

  # Convert to micrometers
  coords_um <- coords / pixels_per_um

  # Compute distances
  dist_vec <- as.vector(dist(coords_um))

  list(
    n_spots = nrow(coords),
    min_distance_um = min(dist_vec),
    max_distance_um = max(dist_vec),
    mean_distance_um = mean(dist_vec),
    median_distance_um = median(dist_vec),
    q25_distance_um = quantile(dist_vec, 0.25),
    q75_distance_um = quantile(dist_vec, 0.75),
    pixels_per_um = pixels_per_um,
    coord_unit = if (pixels_per_um == 1) "assumed micrometers" else "converted to micrometers"
  )
}

#' @title Identify Communication Hotspots
#' @description Find spatial regions with high metabolite communication activity.
#'   **Note**: This function requires running computeSpatialCommunication() with
#'   analysis_level='spot' first.
#'
#' @param object A scMetaLink object with spatial communication results
#' @param metabolite Character. Metabolite ID or name to analyze (NULL for aggregate)
#' @param type Character. "sender" or "receiver" to identify production or sensing hotspots
#' @param n_hotspots Integer. Number of hotspot regions to identify
#' @param method Character. Method for hotspot detection: "density" or "clustering"
#'
#' @return A data.frame with hotspot information
#' @export
identifyCommunicationHotspots <- function(object,
                                          metabolite = NULL,
                                          type = "sender",
                                          n_hotspots = 5,
                                          method = "density") {
  if (!isTRUE(object@parameters$is_spatial)) {
    stop("Object does not contain spatial information")
  }

  if (is.null(object@parameters$spot_communication)) {
    stop(paste0(
      "Spot-level communication not computed. ",
      "Run computeSpatialCommunication(obj, analysis_level='spot') first. ",
      "Note: The default analysis_level is 'region'. You need to explicitly set ",
      "analysis_level='spot' for hotspot identification."
    ))
  }

  if (!type %in% c("sender", "receiver")) {
    stop("type must be 'sender' or 'receiver'")
  }

  coords <- object@parameters$spatial_coords
  spot_comm <- object@parameters$spot_communication

  # Get communication scores for specified metabolite
  if (is.null(metabolite)) {
    # Aggregate across all metabolites
    all_scores <- lapply(spot_comm, function(m) {
      if (is.null(m)) return(NULL)
      if (type == "sender") Matrix::rowSums(m) else Matrix::colSums(m)
    })
    scores <- Reduce(`+`, Filter(Negate(is.null), all_scores))
  } else {
    if (!metabolite %in% names(spot_comm)) {
      stop("Metabolite not found in spot communication results")
    }
    m <- spot_comm[[metabolite]]
    if (is.null(m)) {
      stop("No communication data for this metabolite")
    }
    scores <- if (type == "sender") Matrix::rowSums(m) else Matrix::colSums(m)
  }

  # Identify hotspots
  if (method == "density") {
    # Simple approach: top scoring spots
    n_top <- min(n_hotspots * 10, length(scores))
    top_idx <- order(scores, decreasing = TRUE)[seq_len(n_top)]
    hotspot_coords <- coords[top_idx, ]
    hotspot_scores <- scores[top_idx]

    # Cluster nearby hotspots
    if (nrow(hotspot_coords) > n_hotspots) {
      km <- kmeans(hotspot_coords, centers = n_hotspots, nstart = 10)
      hotspots <- data.frame(
        hotspot_id = 1:n_hotspots,
        center_x = km$centers[, 1],
        center_y = km$centers[, 2],
        n_spots = as.vector(table(km$cluster)),
        mean_score = tapply(hotspot_scores, km$cluster, mean),
        stringsAsFactors = FALSE
      )
    } else {
      hotspots <- data.frame(
        hotspot_id = seq_along(top_idx),
        center_x = hotspot_coords$x,
        center_y = hotspot_coords$y,
        n_spots = 1,
        mean_score = hotspot_scores,
        stringsAsFactors = FALSE
      )
    }
  } else {
    stop("method must be 'density'. Other methods not yet implemented.")
  }

  hotspots[order(-hotspots$mean_score), ]
}

#' Check if object has spatial information
#' @keywords internal
.is_spatial <- function(object) {
  isTRUE(object@parameters$is_spatial)
}
