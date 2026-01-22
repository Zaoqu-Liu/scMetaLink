# Tests for metabolite production inference

test_that("inferProduction runs without error", {
  skip_on_cran()
  
  # Load example data
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(
    expression_data = crc_expr,
    cell_meta = crc_meta,
    cell_type_column = "cell_type",
    min_cells = 10
  )
  
  obj <- inferProduction(obj, verbose = FALSE)
  
  expect_false(is.null(obj@production_scores))
  expect_true(is.matrix(obj@production_scores))
  expect_true(nrow(obj@production_scores) > 0)
})

test_that("inferProduction validates input", {
  expr <- create_test_expression()
  meta <- create_test_metadata()
  
  obj <- createScMetaLink(expr, meta, "cell_type", min_cells = 5)
  
  # Invalid method
  expect_error(
    inferProduction(obj, method = "invalid"),
    "method must be"
  )
  
  # Invalid mean_method
  expect_error(
    inferProduction(obj, mean_method = "invalid"),
    "mean_method must be"
  )
  
  # Invalid min_pct
  expect_error(
    inferProduction(obj, min_pct = 1.5),
    "min_pct must be"
  )
})

test_that("inferProduction methods produce different results", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj1 <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj2 <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj3 <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  
  obj1 <- inferProduction(obj1, method = "mean", verbose = FALSE)
  obj2 <- inferProduction(obj2, method = "proportion", verbose = FALSE)
  obj3 <- inferProduction(obj3, method = "combined", verbose = FALSE)
  
  # Results should differ
  expect_false(identical(obj1@production_scores, obj2@production_scores))
  expect_false(identical(obj2@production_scores, obj3@production_scores))
})

test_that("getTopProducers works correctly", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj <- inferProduction(obj, verbose = FALSE)
  
  # Get a metabolite from the results
  met_id <- rownames(obj@production_scores)[1]
  
  top <- getTopProducers(obj, met_id, top_n = 3)
  
  expect_s3_class(top, "data.frame")
  expect_equal(nrow(top), 3)
  expect_true("cell_type" %in% names(top))
  expect_true("production_score" %in% names(top))
})
