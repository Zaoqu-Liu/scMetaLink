# Tests for scMetaLink class and object creation

test_that("createScMetaLink creates valid object", {
  expr <- create_test_expression()
  meta <- create_test_metadata()
  
  obj <- createScMetaLink(
    expression_data = expr,
    cell_meta = meta,
    cell_type_column = "cell_type",
    min_cells = 5
  )
  
  expect_s4_class(obj, "scMetaLink")
  expect_equal(ncol(obj@expression_data), 50)
  expect_equal(nrow(obj@cell_meta), 50)
  expect_equal(obj@cell_type_column, "cell_type")
})

test_that("createScMetaLink validates input correctly", {
  expr <- create_test_expression()
  meta <- create_test_metadata()
  

  # Test invalid expression data type
  expect_error(
    createScMetaLink(
      expression_data = as.data.frame(as.matrix(expr)),
      cell_meta = meta,
      cell_type_column = "cell_type"
    ),
    "must be a matrix"
  )
  
  # Test missing cell type column
  expect_error(
    createScMetaLink(
      expression_data = expr,
      cell_meta = meta,
      cell_type_column = "nonexistent"
    ),
    "not found"
  )
  
  # Test missing row names (results in no matching cells)
  meta_no_rownames <- meta
  rownames(meta_no_rownames) <- NULL
  expect_error(
    createScMetaLink(
      expression_data = expr,
      cell_meta = meta_no_rownames,
      cell_type_column = "cell_type"
    ),
    "row names|No matching"
  )
})

test_that("createScMetaLink filters cell types correctly", {
  expr <- create_test_expression(n_cells = 30)
  meta <- data.frame(
    cell_type = c(rep("TypeA", 20), rep("TypeB", 5), rep("TypeC", 5)),
    row.names = paste0("Cell", 1:30),
    stringsAsFactors = FALSE
  )
  
  # With min_cells = 10, TypeB and TypeC should be filtered out
  obj <- createScMetaLink(
    expression_data = expr,
    cell_meta = meta,
    cell_type_column = "cell_type",
    min_cells = 10
  )
  
  expect_equal(ncol(obj@expression_data), 20)
  expect_equal(unique(obj@cell_meta$cell_type), "TypeA")
})

test_that("show method works", {
  expr <- create_test_expression()
  meta <- create_test_metadata()
  
  obj <- createScMetaLink(expr, meta, "cell_type", min_cells = 5)
  
  # show method should not error
  expect_output(show(obj), "scMetaLink Object")
})

test_that("accessor methods work correctly", {
  expr <- create_test_expression()
  meta <- create_test_metadata()
  
  obj <- createScMetaLink(expr, meta, "cell_type", min_cells = 5)
  
  # Initially should be NULL
  expect_null(getProductionScores(obj))
  expect_null(getSensingScores(obj))
  expect_null(getCommunicationScores(obj))
  
  # Parameters should be a list
  params <- getParameters(obj)
  expect_type(params, "list")
  expect_true("min_cells" %in% names(params))
})
