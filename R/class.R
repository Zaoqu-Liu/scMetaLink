#' @title scMetaLink Class Definition
#' @description S4 class for storing metabolite-mediated cell communication analysis results
#' @name scMetaLink-class
#' @rdname scMetaLink-class
#' @exportClass scMetaLink
#'
#' @slot expression_data Matrix. Normalized gene expression matrix (genes x cells)
#' @slot cell_meta data.frame. Cell metadata with cell type annotations
#' @slot cell_type_column Character. Column name for cell type in cell_meta
#' @slot production_scores Matrix. Metabolite production potential scores (metabolites x cell_types)
#' @slot sensing_scores Matrix. Metabolite sensing capability scores (metabolites x cell_types)
#' @slot communication_scores Array. Communication scores (sender x receiver x metabolite)
#' @slot communication_pvalues Array. P-values from permutation test
#' @slot significant_interactions data.frame. Filtered significant interactions
#' @slot pathway_aggregated data.frame. Pathway-level aggregated results
#' @slot parameters list. Analysis parameters
#' @slot database list. MetalinksDB reference data
#'
setClass(
  "scMetaLink",
  slots = list(
    expression_data = "ANY",
    cell_meta = "data.frame",
    cell_type_column = "character",
    production_scores = "ANY",
    sensing_scores = "ANY",
    communication_scores = "ANY",
    communication_pvalues = "ANY",
    significant_interactions = "data.frame",
    pathway_aggregated = "data.frame",
    parameters = "list",
    database = "list"
  ),
  prototype = list(
    expression_data = NULL,
    cell_meta = data.frame(),
    cell_type_column = "cell_type",
    production_scores = NULL,
    sensing_scores = NULL,
    communication_scores = NULL,
    communication_pvalues = NULL,
    significant_interactions = data.frame(),
    pathway_aggregated = data.frame(),
    parameters = list(),
    database = list()
  )
)

#' @title Create scMetaLink Object
#' @description Initialize a scMetaLink object from expression data and cell metadata
#'
#' @param expression_data A matrix or dgCMatrix of normalized expression values (genes x cells)
#' @param cell_meta A data.frame containing cell metadata
#' @param cell_type_column Character. Column name in cell_meta containing cell type annotations
#' @param min_cells Integer. Minimum number of cells per cell type (default: 10)
#' @param min_genes Integer. Minimum number of genes detected per cell (default: 200)
#'
#' @return A scMetaLink object
#' @export
#'
#' @examples
#' \donttest{
#' # Create from expression matrix
#' expr_mat <- matrix(rpois(1000, 5), nrow = 100, ncol = 10)
#' rownames(expr_mat) <- paste0("Gene", 1:100)
#' colnames(expr_mat) <- paste0("Cell", 1:10)
#' meta <- data.frame(cell_type = rep(c("TypeA", "TypeB"), each = 5))
#' rownames(meta) <- colnames(expr_mat)
#' obj <- createScMetaLink(expr_mat, meta, "cell_type")
#' }
createScMetaLink <- function(expression_data,
                             cell_meta,
                             cell_type_column = "cell_type",
                             min_cells = 10,
                             min_genes = 200) {
  # Input validation
  if (!is.matrix(expression_data) && !inherits(expression_data, "dgCMatrix")) {
    stop("expression_data must be a matrix or sparse dgCMatrix")
  }

  if (!is.data.frame(cell_meta)) {
    stop("cell_meta must be a data.frame")
  }

  if (!cell_type_column %in% colnames(cell_meta)) {
    stop(paste("Column", cell_type_column, "not found in cell_meta"))
  }

  if (is.null(rownames(expression_data))) {
    stop("expression_data must have row names (gene symbols)")
  }

  if (is.null(colnames(expression_data))) {
    stop("expression_data must have column names (cell IDs)")
  }

  if (is.null(rownames(cell_meta))) {
    stop("cell_meta must have row names (cell IDs)")
  }

  if (min_cells < 1) {
    stop("min_cells must be at least 1")
  }

  # Match cells
  common_cells <- intersect(colnames(expression_data), rownames(cell_meta))
  if (length(common_cells) == 0) {
    stop("No matching cells between expression_data and cell_meta")
  }

  # Filter cells
  expression_data <- expression_data[, common_cells, drop = FALSE]
  cell_meta <- cell_meta[common_cells, , drop = FALSE]

  # Filter cell types with minimum cells
  cell_type_counts <- table(cell_meta[[cell_type_column]])
  valid_types <- names(cell_type_counts)[cell_type_counts >= min_cells]

  if (length(valid_types) == 0) {
    stop(paste("No cell types with at least", min_cells, "cells"))
  }

  keep_cells <- cell_meta[[cell_type_column]] %in% valid_types
  expression_data <- expression_data[, keep_cells, drop = FALSE]
  cell_meta <- cell_meta[keep_cells, , drop = FALSE]

  message(sprintf(
    "Created scMetaLink object with %d genes, %d cells, %d cell types",
    nrow(expression_data), ncol(expression_data), length(valid_types)
  ))

  # Load database
  db <- .load_metalinksdb()

  # Create object
  new("scMetaLink",
    expression_data = expression_data,
    cell_meta = cell_meta,
    cell_type_column = cell_type_column,
    parameters = list(
      min_cells = min_cells,
      min_genes = min_genes,
      created_at = Sys.time()
    ),
    database = db
  )
}

#' @title Create scMetaLink from Seurat Object
#' @description Initialize a scMetaLink object from a Seurat object
#'
#' @param seurat_obj A Seurat object
#' @param cell_type_column Character. Column name in meta.data for cell type
#' @param assay Character. Assay to use (default: "RNA")
#' @param slot Character. Slot to use (default: "data" for normalized data)
#' @param min_cells Integer. Minimum cells per cell type
#'
#' @return A scMetaLink object
#' @export
createScMetaLinkFromSeurat <- function(seurat_obj,
                                       cell_type_column = "cell_type",
                                       assay = "RNA",
                                       slot = "data",
                                       min_cells = 10) {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' is required. Please install it.")
  }

  if (!slot %in% c("data", "counts")) {
    stop("slot must be 'data' or 'counts'")
  }

  # Extract expression data
  if (slot == "data") {
    expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = "data")
  } else {
    expr_data <- Seurat::GetAssayData(seurat_obj, assay = assay, slot = "counts")
  }

  # Extract metadata
  cell_meta <- seurat_obj@meta.data

  # Create scMetaLink object
  createScMetaLink(
    expression_data = expr_data,
    cell_meta = cell_meta,
    cell_type_column = cell_type_column,
    min_cells = min_cells
  )
}

#' @title Show Method for scMetaLink
#' @param object scMetaLink object
#' @rdname scMetaLink-class
#' @aliases show,scMetaLink-method
#' @exportMethod show
setMethod("show", "scMetaLink", function(object) {
  cat("scMetaLink Object\n")
  cat("=================\n")
  cat(sprintf("Genes: %d\n", nrow(object@expression_data)))
  cat(sprintf("Cells: %d\n", ncol(object@expression_data)))

  cell_types <- unique(object@cell_meta[[object@cell_type_column]])
  cat(sprintf(
    "Cell types: %d (%s%s)\n",
    length(cell_types),
    paste(head(cell_types, 3), collapse = ", "),
    if (length(cell_types) > 3) ", ..." else ""
  ))

  if (!is.null(object@production_scores)) {
    cat(sprintf("Metabolites with production scores: %d\n", nrow(object@production_scores)))
  }
  if (!is.null(object@sensing_scores)) {
    cat(sprintf("Metabolites with sensing scores: %d\n", nrow(object@sensing_scores)))
  }
  if (!is.null(object@communication_scores)) {
    cat(sprintf(
      "Communication array: %s\n",
      paste(dim(object@communication_scores), collapse = " x ")
    ))
  }
  if (nrow(object@significant_interactions) > 0) {
    cat(sprintf("Significant interactions: %d\n", nrow(object@significant_interactions)))
  }
  if (nrow(object@pathway_aggregated) > 0) {
    cat(sprintf("Pathway aggregations: %d\n", nrow(object@pathway_aggregated)))
  }
})

#' @title Accessor Functions for scMetaLink Objects
#' @description Functions to extract data from scMetaLink objects.
#' @name accessors
#' @rdname accessors
#' @aliases getProductionScores getSensingScores getCommunicationScores
#' @aliases getSignificantInteractions getPathwayAggregated getParameters
#' @param object A scMetaLink object
#' @return The requested data from the scMetaLink object:
#' \itemize{
#'   \item \code{getProductionScores}: Matrix of metabolite production scores
#'   \item \code{getSensingScores}: Matrix of metabolite sensing scores
#'   \item \code{getCommunicationScores}: 3D array of communication scores
#'   \item \code{getSignificantInteractions}: data.frame of significant interactions
#'   \item \code{getPathwayAggregated}: data.frame of pathway-aggregated results
#'   \item \code{getParameters}: list of analysis parameters
#' }
#' @examples
#' \donttest{
#' data(crc_example)
#' obj <- createScMetaLink(crc_expr, crc_meta, "cell_type")
#' obj <- inferProduction(obj)
#' prod_scores <- getProductionScores(obj)
#' }
#' @export
setGeneric("getProductionScores", function(object) standardGeneric("getProductionScores"))

#' @rdname accessors
#' @exportMethod getProductionScores
setMethod("getProductionScores", "scMetaLink", function(object) {
  object@production_scores
})

#' @rdname accessors
#' @export
setGeneric("getSensingScores", function(object) standardGeneric("getSensingScores"))

#' @rdname accessors
#' @exportMethod getSensingScores
setMethod("getSensingScores", "scMetaLink", function(object) {
  object@sensing_scores
})

#' @rdname accessors
#' @export
setGeneric("getCommunicationScores", function(object) standardGeneric("getCommunicationScores"))

#' @rdname accessors
#' @exportMethod getCommunicationScores
setMethod("getCommunicationScores", "scMetaLink", function(object) {
  object@communication_scores
})

#' @rdname accessors
#' @export
setGeneric("getSignificantInteractions", function(object) standardGeneric("getSignificantInteractions"))

#' @rdname accessors
#' @exportMethod getSignificantInteractions
setMethod("getSignificantInteractions", "scMetaLink", function(object) {
  object@significant_interactions
})

#' @rdname accessors
#' @export
setGeneric("getPathwayAggregated", function(object) standardGeneric("getPathwayAggregated"))

#' @rdname accessors
#' @exportMethod getPathwayAggregated
setMethod("getPathwayAggregated", "scMetaLink", function(object) {
  object@pathway_aggregated
})

#' @rdname accessors
#' @export
setGeneric("getParameters", function(object) standardGeneric("getParameters"))

#' @rdname accessors
#' @exportMethod getParameters
setMethod("getParameters", "scMetaLink", function(object) {
  object@parameters
})
