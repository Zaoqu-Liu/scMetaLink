# Tests for communication computation

test_that("computeCommunication runs without error", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj <- inferProduction(obj, verbose = FALSE)
  obj <- inferSensing(obj, verbose = FALSE)
  
  # Run with few permutations for speed
  obj <- computeCommunication(obj, n_permutations = 10, verbose = FALSE)
  
  expect_false(is.null(obj@communication_scores))
  expect_true(is.array(obj@communication_scores))
  expect_equal(length(dim(obj@communication_scores)), 3)
})

test_that("computeCommunication validates input", {
  skip_on_cran()
  
  # Use real data for validation tests
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  
  # No production scores
  expect_error(
    computeCommunication(obj),
    "Run inferProduction"
  )
  
  obj <- inferProduction(obj, verbose = FALSE)
  
  # No sensing scores
  expect_error(
    computeCommunication(obj),
    "Run inferSensing"
  )
})

test_that("computeCommunication methods produce different results", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj1 <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj1 <- inferProduction(obj1, verbose = FALSE)
  obj1 <- inferSensing(obj1, verbose = FALSE)
  
  obj2 <- obj1
  obj3 <- obj1
  
  obj1 <- computeCommunication(obj1, method = "geometric", n_permutations = 0, verbose = FALSE)
  obj2 <- computeCommunication(obj2, method = "product", n_permutations = 0, verbose = FALSE)
  obj3 <- computeCommunication(obj3, method = "harmonic", n_permutations = 0, verbose = FALSE)
  
  # Results should differ
  expect_false(identical(obj1@communication_scores, obj2@communication_scores))
  expect_false(identical(obj2@communication_scores, obj3@communication_scores))
})

test_that("filterSignificantInteractions works", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj <- inferProduction(obj, verbose = FALSE)
  obj <- inferSensing(obj, verbose = FALSE)
  obj <- computeCommunication(obj, n_permutations = 10, verbose = FALSE)
  
  # Filter with lenient threshold
  obj <- filterSignificantInteractions(obj, pvalue_threshold = 0.5, adjust_method = "none")
  
  sig <- getSignificantInteractions(obj)
  expect_s3_class(sig, "data.frame")
  
  if (nrow(sig) > 0) {
    expect_true("sender" %in% names(sig))
    expect_true("receiver" %in% names(sig))
    expect_true("metabolite_id" %in% names(sig))
    expect_true("communication_score" %in% names(sig))
  }
})

test_that("summarizeCommunicationPairs works", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj <- inferProduction(obj, verbose = FALSE)
  obj <- inferSensing(obj, verbose = FALSE)
  obj <- computeCommunication(obj, n_permutations = 10, verbose = FALSE)
  obj <- filterSignificantInteractions(obj, pvalue_threshold = 0.5, adjust_method = "none")
  
  if (nrow(obj@significant_interactions) > 0) {
    summary <- summarizeCommunicationPairs(obj, aggregate_method = "sum")
    expect_s3_class(summary, "data.frame")
  }
})

test_that("getCommunicationMatrix works", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj <- inferProduction(obj, verbose = FALSE)
  obj <- inferSensing(obj, verbose = FALSE)
  obj <- computeCommunication(obj, n_permutations = 10, verbose = FALSE)
  obj <- filterSignificantInteractions(obj, pvalue_threshold = 0.5, adjust_method = "none")
  
  if (nrow(obj@significant_interactions) > 0) {
    mat <- getCommunicationMatrix(obj)
    expect_true(is.matrix(mat))
  }
})
