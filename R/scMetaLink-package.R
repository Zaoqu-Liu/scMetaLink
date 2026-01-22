#' @title scMetaLink: Single-Cell Metabolite-Mediated Cell Communication Analysis
#'
#' @description
#' A comprehensive framework for inferring metabolite-mediated cell-cell
#' communication from single-cell transcriptomic data. scMetaLink integrates
#' metabolite production potential via enzyme expression, metabolite sensing
#' capability via receptor and transporter expression, and secretion potential
#' to construct intercellular metabolic communication networks.
#'
#' @details
#' The package provides the following main functions:
#' \itemize{
#'   \item \code{\link{createScMetaLink}}: Create analysis object from expression data
#'   \item \code{\link{inferProduction}}: Infer metabolite production potential
#'   \item \code{\link{inferSensing}}: Infer metabolite sensing capability
#'   \item \code{\link{computeCommunication}}: Compute cell-cell communication scores
#'   \item \code{\link{filterSignificantInteractions}}: Filter significant interactions
#'   \item \code{\link{aggregateByPathway}}: Aggregate by metabolic pathways
#'   \item \code{\link{runScMetaLink}}: Run complete analysis pipeline
#' }
#'
#' For spatial transcriptomics data:
#' \itemize{
#'   \item \code{\link{createScMetaLinkFromSpatial}}: Create object from spatial data
#'   \item \code{\link{computeSpatialCommunication}}: Compute spatially-weighted communication
#' }
#'
#' @section Database:
#' scMetaLink utilizes MetalinksDB, containing:
#' \itemize{
#'   \item 41,894 metabolite-protein interactions
#'   \item 1,128 metabolites
#'   \item 4,374 proteins/genes
#'   \item 157,741 pathway associations
#' }
#'
#' @author Zaoqu Liu \email{liuzaoqu@@163.com}
#'
#' @seealso
#' Useful links:
#' \itemize{
#'   \item \url{https://github.com/Zaoqu-Liu/scMetaLink}
#'   \item \url{https://Zaoqu-Liu.github.io/scMetaLink/}
#'   \item Report bugs at \url{https://github.com/Zaoqu-Liu/scMetaLink/issues}
#' }
#'
#' @docType package
#' @name scMetaLink-package
#' @aliases scMetaLink-package
#'
#' @import methods
#' @import ggplot2
#' @importFrom Matrix Matrix sparseMatrix rowSums colSums rowMeans colMeans
#' @importFrom stats sd quantile p.adjust p.adjust.methods phyper median kmeans aggregate dist setNames
#' @importFrom utils head tail txtProgressBar setTxtProgressBar packageVersion write.csv
#' @importFrom grDevices colorRampPalette
#' @importFrom graphics title
#'
"_PACKAGE"

# Global variables declaration to avoid R CMD check NOTEs
# These are column names used in data.frame operations with non-standard evaluation
utils::globalVariables(c(

  # Data variables
  "metalinksdb",
  

  # Column names used in ggplot2 aes()
  "x", "y", "score", "cell_type", "name", "communication_score",
  "pvalue_adjusted", "fold_enrichment", "n_overlap", "change",
  "mean_score", "significant", "pathway", "n_metabolites",
  "strength_norm", "ct", "pair", "distance", "metabolite",
  
  # Column names in data manipulation
  "sender", "receiver", "metabolite_id", "metabolite_name",
  "hmdb", "uniprot", "gene_symbol", "protein_type", "type",
  "mor", "combined_score", "transport_direction",
  
  # ggraph aesthetics
  "node_size",
  
  # Spatial plot variables
  "center_x", "center_y", "hotspot_id", "label",
  "x_start", "y_start", "x_end", "y_end",
  "strength_norm",
  
  # Other variables
  ".data"
))
