# Tests for workflow functions

test_that("runScMetaLink runs complete pipeline", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  # Run with no permutations and lenient threshold for testing
  obj <- runScMetaLink(
    expression_data = crc_expr,
    cell_meta = crc_meta,
    cell_type_column = "cell_type",
    n_permutations = 0,  # Skip permutation for faster testing
    pvalue_threshold = 1.0,  # Accept all interactions
    verbose = FALSE
  )
  
  expect_s4_class(obj, "scMetaLink")
  expect_false(is.null(obj@production_scores))
  expect_false(is.null(obj@sensing_scores))
  expect_false(is.null(obj@communication_scores))
})

test_that("runScMetaLink validates parameters", {
  data("crc_example", package = "scMetaLink")
  
  expect_error(
    runScMetaLink(crc_expr, crc_meta, "cell_type", method = "invalid"),
    "method must be"
  )
  
  expect_error(
    runScMetaLink(crc_expr, crc_meta, "cell_type", min_pct = 2),
    "min_pct must be"
  )
  
  expect_error(
    runScMetaLink(crc_expr, crc_meta, "cell_type", pvalue_threshold = 0),
    "pvalue_threshold must be"
  )
})

test_that("compareCommunication works", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  # Create two objects with guaranteed significant interactions
  obj1 <- runScMetaLink(
    expression_data = crc_expr,
    cell_meta = crc_meta,
    cell_type_column = "cell_type",
    n_permutations = 0,
    pvalue_threshold = 1.0,  # Accept all
    verbose = FALSE
  )
  
  obj2 <- obj1  # Same for testing
  
  if (nrow(obj1@significant_interactions) > 0) {
    diff <- compareCommunication(obj1, obj2, 
                                 condition_names = c("Cond1", "Cond2"),
                                 method = "log2fc")
    
    expect_s3_class(diff, "data.frame")
    expect_true("change" %in% names(diff))
  }
})

test_that("compareCommunication validates input", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  # Test with objects that have significant interactions
  obj <- runScMetaLink(crc_expr, crc_meta, "cell_type", 
                       n_permutations = 0, pvalue_threshold = 1.0, verbose = FALSE)
  
  if (nrow(obj@significant_interactions) > 0) {
    # Invalid condition_names length
    expect_error(
      compareCommunication(obj, obj, condition_names = c("A")),
      "length 2"
    )
    
    # Invalid method
    expect_error(
      compareCommunication(obj, obj, method = "invalid"),
      "method must be"
    )
  }
})

test_that("identifyCellTypeSpecificMetabolites works", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj <- inferProduction(obj, verbose = FALSE)
  
  specific <- identifyCellTypeSpecificMetabolites(obj, type = "production")
  
  expect_s3_class(specific, "data.frame")
  if (nrow(specific) > 0) {
    expect_true("cell_type" %in% names(specific))
    expect_true("z_score" %in% names(specific))
  }
})

test_that("getSummaryStats works", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- runScMetaLink(
    expression_data = crc_expr,
    cell_meta = crc_meta,
    cell_type_column = "cell_type",
    n_permutations = 0,
    pvalue_threshold = 1.0,
    verbose = FALSE
  )
  
  stats <- getSummaryStats(obj)
  
  expect_type(stats, "list")
  expect_true("n_genes" %in% names(stats))
  expect_true("n_cells" %in% names(stats))
  expect_true("n_cell_types" %in% names(stats))
})
