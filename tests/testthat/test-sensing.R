# Tests for metabolite sensing inference

test_that("inferSensing runs without error", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(
    expression_data = crc_expr,
    cell_meta = crc_meta,
    cell_type_column = "cell_type",
    min_cells = 10
  )
  
  obj <- inferSensing(obj, verbose = FALSE)
  
  expect_false(is.null(obj@sensing_scores))
  expect_true(is.matrix(obj@sensing_scores))
  expect_true(nrow(obj@sensing_scores) > 0)
})

test_that("inferSensing validates input", {
  expr <- create_test_expression()
  meta <- create_test_metadata()
  
  obj <- createScMetaLink(expr, meta, "cell_type", min_cells = 5)
  
  # Invalid method
  expect_error(
    inferSensing(obj, method = "invalid"),
    "method must be"
  )
  
  # Invalid hill_n
  expect_error(
    inferSensing(obj, use_hill = TRUE, hill_n = -1),
    "hill_n must be positive"
  )
  
  # Invalid hill_Kh
  expect_error(
    inferSensing(obj, use_hill = TRUE, hill_Kh = 0),
    "hill_Kh must be between"
  )
})

test_that("Hill transformation works", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj1 <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj2 <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  
  obj1 <- inferSensing(obj1, use_hill = FALSE, verbose = FALSE)
  obj2 <- inferSensing(obj2, use_hill = TRUE, verbose = FALSE)
  
  # Results should differ when Hill is applied
  expect_false(identical(obj1@sensing_scores, obj2@sensing_scores))
})

test_that("getTopSensors works correctly", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj <- inferSensing(obj, verbose = FALSE)
  
  # Get a metabolite from the results
  met_id <- rownames(obj@sensing_scores)[1]
  
  top <- getTopSensors(obj, met_id, top_n = 3)
  
  expect_s3_class(top, "data.frame")
  expect_equal(nrow(top), 3)
  expect_true("cell_type" %in% names(top))
  expect_true("sensing_score" %in% names(top))
})

test_that("getMetaboliteReceptors works", {
  skip_on_cran()
  
  # Use a known metabolite
  receptors <- getMetaboliteReceptors("L-Glutamic acid")
  
  expect_s3_class(receptors, "data.frame")
  expect_true("gene_symbol" %in% names(receptors))
})
