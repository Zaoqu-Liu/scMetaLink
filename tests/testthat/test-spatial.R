# Tests for spatial analysis

test_that("createScMetaLinkFromSpatial works", {
  skip_on_cran()
  
  # Load spatial example data
  data("st_expr", package = "scMetaLink")
  data("st_meta", package = "scMetaLink")
  data("st_scalefactors", package = "scMetaLink")
  
  obj <- createScMetaLinkFromSpatial(
    expression_data = st_expr,
    spatial_coords = st_meta[, c("x", "y")],
    cell_meta = st_meta,
    cell_type_column = "cell_type",
    scale_factors = st_scalefactors,
    min_cells = 5
  )
  
  expect_s4_class(obj, "scMetaLink")
  expect_true(obj@parameters$is_spatial)
  expect_false(is.null(obj@parameters$spatial_coords))
})

test_that("createScMetaLinkFromSpatial validates input", {
  expr <- create_test_expression()
  meta <- create_test_metadata()
  coords <- create_test_spatial_coords()
  
  # Missing spatial coords
  expect_error(
    createScMetaLinkFromSpatial(
      expression_data = expr,
      spatial_coords = NULL,
      cell_meta = meta,
      cell_type_column = "cell_type"
    ),
    "spatial_coords must be provided"
  )
  
  # Invalid coords structure
  expect_error(
    createScMetaLinkFromSpatial(
      expression_data = expr,
      spatial_coords = data.frame(x = 1:10),  # Only 1 column
      cell_meta = meta,
      cell_type_column = "cell_type"
    ),
    "at least 2 columns"
  )
})

test_that("computeSpatialCommunication works with knn method", {
  skip_on_cran()
  
  # Load spatial example data
  data("st_expr", package = "scMetaLink")
  data("st_meta", package = "scMetaLink")
  data("st_scalefactors", package = "scMetaLink")
  
  obj <- createScMetaLinkFromSpatial(
    expression_data = st_expr,
    spatial_coords = st_meta[, c("x", "y")],
    cell_meta = st_meta,
    cell_type_column = "cell_type",
    scale_factors = st_scalefactors,
    min_cells = 5
  )
  
  obj <- inferProduction(obj, verbose = FALSE)
  obj <- inferSensing(obj, verbose = FALSE)
  
  obj <- computeSpatialCommunication(
    obj,
    method = "knn",
    k_neighbors = 6,
    n_permutations = 10,
    verbose = FALSE
  )
  
  expect_false(is.null(obj@communication_scores))
  expect_true("spatial_communication" %in% names(obj@parameters))
})

test_that("computeSpatialCommunication works with different methods", {
  skip_on_cran()
  
  # Load spatial example data
  data("st_expr", package = "scMetaLink")
  data("st_meta", package = "scMetaLink")
  data("st_scalefactors", package = "scMetaLink")
  
  obj <- createScMetaLinkFromSpatial(
    expression_data = st_expr,
    spatial_coords = st_meta[, c("x", "y")],
    cell_meta = st_meta,
    cell_type_column = "cell_type",
    scale_factors = st_scalefactors,
    min_cells = 5
  )
  
  obj <- inferProduction(obj, verbose = FALSE)
  obj <- inferSensing(obj, verbose = FALSE)
  
  methods <- c("gaussian", "exponential", "threshold")
  
  for (method in methods) {
    obj_test <- computeSpatialCommunication(
      obj,
      method = method,
      n_permutations = 0,
      verbose = FALSE
    )
    expect_false(is.null(obj_test@communication_scores), 
                 info = paste("Method:", method))
  }
})

test_that("getSpatialDistanceStats works", {
  skip_on_cran()
  
  # Load spatial example data
  data("st_expr", package = "scMetaLink")
  data("st_meta", package = "scMetaLink")
  data("st_scalefactors", package = "scMetaLink")
  
  obj <- createScMetaLinkFromSpatial(
    expression_data = st_expr,
    spatial_coords = st_meta[, c("x", "y")],
    cell_meta = st_meta,
    cell_type_column = "cell_type",
    scale_factors = st_scalefactors,
    min_cells = 5
  )
  
  stats <- getSpatialDistanceStats(obj)
  
  expect_type(stats, "list")
  expect_true("n_spots" %in% names(stats))
  expect_true("mean_distance_um" %in% names(stats))
})

test_that("getSpatialDistanceStats fails on non-spatial object", {
  expr <- create_test_expression()
  meta <- create_test_metadata()
  
  obj <- createScMetaLink(expr, meta, "cell_type", min_cells = 5)
  
  expect_error(
    getSpatialDistanceStats(obj),
    "does not contain spatial information"
  )
})
