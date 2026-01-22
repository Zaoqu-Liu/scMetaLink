# Helper functions and test data setup for scMetaLink tests

#' Create a minimal test expression matrix
#' @return A sparse matrix for testing
create_test_expression <- function(n_genes = 100, n_cells = 50, seed = 42) {
  set.seed(seed)
  
  # Create a simple expression matrix
  mat <- matrix(
    rpois(n_genes * n_cells, lambda = 2),
    nrow = n_genes,
    ncol = n_cells
  )
  
  # Add some structure - make some genes highly expressed in specific cells
  # Use safe indices based on actual dimensions
  half_cells <- n_cells %/% 2
  if (half_cells > 0 && n_genes >= 20) {
    mat[1:10, 1:half_cells] <- mat[1:10, 1:half_cells] + 5
    mat[11:20, (half_cells + 1):n_cells] <- mat[11:20, (half_cells + 1):n_cells] + 5
  }
  
  # Convert to sparse matrix
  mat <- Matrix::Matrix(mat, sparse = TRUE)
  
  # Add gene names (use some actual MetalinksDB genes if available)
  rownames(mat) <- paste0("Gene", seq_len(n_genes))
  colnames(mat) <- paste0("Cell", seq_len(n_cells))
  
  mat
}

#' Create test cell metadata
#' @return A data.frame with cell metadata
create_test_metadata <- function(n_cells = 50) {
  data.frame(
    cell_type = rep(c("TypeA", "TypeB"), each = n_cells / 2),
    condition = rep(c("Control", "Treatment"), n_cells / 2),
    row.names = paste0("Cell", seq_len(n_cells)),
    stringsAsFactors = FALSE
  )
}

#' Create test spatial coordinates
#' @return A data.frame with spatial coordinates
create_test_spatial_coords <- function(n_spots = 50, seed = 42) {
  set.seed(seed)
  data.frame(
    x = runif(n_spots, 0, 1000),
    y = runif(n_spots, 0, 1000),
    row.names = paste0("Cell", seq_len(n_spots))
  )
}
