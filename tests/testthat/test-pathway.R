# Tests for pathway analysis

test_that("aggregateByPathway runs without error", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj <- inferProduction(obj, verbose = FALSE)
  obj <- inferSensing(obj, verbose = FALSE)
  obj <- computeCommunication(obj, n_permutations = 10, verbose = FALSE)
  obj <- filterSignificantInteractions(obj, pvalue_threshold = 0.5, adjust_method = "none")
  
  if (nrow(obj@significant_interactions) > 0) {
    obj <- aggregateByPathway(obj)
    
    pw_agg <- getPathwayAggregated(obj)
    expect_s3_class(pw_agg, "data.frame")
  }
})

test_that("enrichPathways works", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj <- inferProduction(obj, verbose = FALSE)
  obj <- inferSensing(obj, verbose = FALSE)
  obj <- computeCommunication(obj, n_permutations = 10, verbose = FALSE)
  obj <- filterSignificantInteractions(obj, pvalue_threshold = 0.5, adjust_method = "none")
  
  if (nrow(obj@significant_interactions) > 0) {
    enriched <- enrichPathways(obj, pvalue_threshold = 1, min_overlap = 1)
    expect_s3_class(enriched, "data.frame")
  }
})

test_that("getPathwayMetabolites works", {
  skip_on_cran()
  
  mets <- getPathwayMetabolites("glycolysis", only_signaling = FALSE)
  expect_s3_class(mets, "data.frame")
})

test_that("listTopPathways validates input", {
  expr <- create_test_expression()
  meta <- create_test_metadata()
  
  obj <- createScMetaLink(expr, meta, "cell_type", min_cells = 5)
  
  # No pathway aggregation done
  expect_error(
    listTopPathways(obj),
    "Pathway aggregation not done"
  )
})
