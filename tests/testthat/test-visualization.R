# Tests for visualization functions

test_that("plotCommunicationHeatmap works", {
  skip_on_cran()
  skip_if_not_installed("ComplexHeatmap")
  
  data("crc_example", package = "scMetaLink")
  
  obj <- runScMetaLink(
    expression_data = crc_expr,
    cell_meta = crc_meta,
    cell_type_column = "cell_type",
    n_permutations = 10,
    pvalue_threshold = 0.5,
    verbose = FALSE
  )
  
  if (nrow(obj@significant_interactions) > 0) {
    ht <- plotCommunicationHeatmap(obj)
    expect_s4_class(ht, "Heatmap")
  }
})

test_that("plotCommunicationCircle runs without error", {
  skip_on_cran()
  skip_if_not_installed("circlize")
  
  data("crc_example", package = "scMetaLink")
  
  obj <- runScMetaLink(
    expression_data = crc_expr,
    cell_meta = crc_meta,
    cell_type_column = "cell_type",
    n_permutations = 10,
    pvalue_threshold = 0.5,
    verbose = FALSE
  )
  
  if (nrow(obj@significant_interactions) > 0) {
    # Should run without error
    expect_silent(plotCommunicationCircle(obj, top_n = 20))
  }
})

test_that("plotCommunicationNetwork works", {
  skip_on_cran()
  skip_if_not_installed("ggraph")
  skip_if_not_installed("igraph")
  
  data("crc_example", package = "scMetaLink")
  
  obj <- runScMetaLink(
    expression_data = crc_expr,
    cell_meta = crc_meta,
    cell_type_column = "cell_type",
    n_permutations = 10,
    pvalue_threshold = 0.5,
    verbose = FALSE
  )
  
  if (nrow(obj@significant_interactions) > 0) {
    p <- plotCommunicationNetwork(obj, min_score = 0)
    expect_s3_class(p, "ggplot")
  }
})

test_that("plotMetaboliteProfile works", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj <- inferProduction(obj, verbose = FALSE)
  obj <- inferSensing(obj, verbose = FALSE)
  
  # Get a metabolite from the scores
  met_id <- rownames(obj@production_scores)[1]
  
  p <- plotMetaboliteProfile(obj, met_id)
  expect_s3_class(p, "ggplot")
})

test_that("plotTopInteractions works", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- runScMetaLink(
    expression_data = crc_expr,
    cell_meta = crc_meta,
    cell_type_column = "cell_type",
    n_permutations = 10,
    pvalue_threshold = 0.5,
    verbose = FALSE
  )
  
  if (nrow(obj@significant_interactions) > 0) {
    p <- plotTopInteractions(obj, top_n = 10)
    expect_s3_class(p, "ggplot")
  }
})

test_that("plotDifferentialCommunication works", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj1 <- runScMetaLink(
    expression_data = crc_expr,
    cell_meta = crc_meta,
    cell_type_column = "cell_type",
    n_permutations = 0,
    pvalue_threshold = 1.0,
    verbose = FALSE
  )
  
  if (nrow(obj1@significant_interactions) > 0) {
    diff <- compareCommunication(obj1, obj1)
    
    if (nrow(diff) > 0) {
      p <- plotDifferentialCommunication(diff, top_n = 10, type = "bar")
      expect_s3_class(p, "ggplot")
    }
  }
})

test_that("visualization functions validate input", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  
  # No communication scores - should error
  expect_error(
    plotCommunicationHeatmap(obj),
    "Communication scores not calculated|significant_interactions|No significant"
  )
})
