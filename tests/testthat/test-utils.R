# Tests for utility functions

test_that("getDatabaseInfo returns valid info", {
  info <- getDatabaseInfo()
  
  expect_s3_class(info, "data.frame")
  expect_true("Component" %in% names(info))
  expect_true("Count" %in% names(info))
  expect_true(all(info$Count > 0))
})

test_that("searchMetabolite works", {
  skip_on_cran()
  
  # Search for a common metabolite
  results <- searchMetabolite("glucose")
  
  expect_s3_class(results, "data.frame")
  expect_true(nrow(results) > 0)
})

test_that("searchMetabolite validates input", {
  expect_error(searchMetabolite(), "must be provided")
  expect_error(searchMetabolite(""), "must be provided")
})

test_that("searchGene works", {
  skip_on_cran()
  
  # Search for a known gene
  results <- searchGene("SLC")
  
  expect_s3_class(results, "data.frame")
})

test_that("searchGene validates input", {
  expect_error(searchGene(), "must be provided")
  expect_error(searchGene(""), "must be provided")
})

test_that("listMetabolites works", {
  skip_on_cran()
  
  all_mets <- listMetabolites(type = "all")
  signaling_mets <- listMetabolites(type = "signaling")
  metabolic_mets <- listMetabolites(type = "metabolic")
  
  expect_s3_class(all_mets, "data.frame")
  expect_s3_class(signaling_mets, "data.frame")
  expect_s3_class(metabolic_mets, "data.frame")
  
  # Signaling should be a subset
  expect_true(nrow(signaling_mets) <= nrow(all_mets))
})

test_that("listMetabolites validates input", {
  expect_error(listMetabolites(type = "invalid"), "must be")
})

test_that("listGenes works", {
  skip_on_cran()
  
  all_genes <- listGenes(role = "all")
  receptors <- listGenes(role = "receptor")
  transporters <- listGenes(role = "transporter")
  
  expect_s3_class(all_genes, "data.frame")
  expect_s3_class(receptors, "data.frame")
  expect_s3_class(transporters, "data.frame")
})

test_that("listGenes validates input", {
  expect_error(listGenes(role = "invalid"), "must be")
})

test_that("exportResults creates files", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  obj <- inferProduction(obj, verbose = FALSE)
  
  # Create temp directory
  temp_dir <- tempdir()
  files <- exportResults(obj, output_dir = temp_dir, prefix = "test")
  
  expect_true(length(files) > 0)
  expect_true(all(file.exists(files)))
  
  # Cleanup
  unlink(files)
})

test_that("saveScMetaLink and loadScMetaLink work", {
  skip_on_cran()
  
  data("crc_example", package = "scMetaLink")
  
  obj <- createScMetaLink(crc_expr, crc_meta, "cell_type", min_cells = 10)
  
  # Save to temp file
  temp_file <- tempfile(fileext = ".rds")
  saveScMetaLink(obj, temp_file)
  
  expect_true(file.exists(temp_file))
  
  # Load back
  loaded <- loadScMetaLink(temp_file)
  
  expect_s4_class(loaded, "scMetaLink")
  expect_equal(ncol(loaded@expression_data), ncol(obj@expression_data))
  
  # Cleanup
  unlink(temp_file)
})
