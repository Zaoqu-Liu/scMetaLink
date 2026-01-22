# Tests for lactate signaling analysis module

# ============================================================
# Test 1: Gene set completeness and accuracy
# ============================================================
test_that("lactate gene sets are scientifically accurate", {
  genes <- getLactateGenes()

  # Check structure
  expect_type(genes, "list")
  expect_true("production" %in% names(genes))
  expect_true("degradation" %in% names(genes))
  expect_true("direct_sensing" %in% names(genes))
  expect_true("indirect_sensing" %in% names(genes))
  expect_true("uptake" %in% names(genes))

  # Production - synthesis
  expect_true("LDHA" %in% genes$production$synthesis)
  expect_false("LDHB" %in% genes$production$synthesis)

  # Production - export
  expect_true("SLC16A3" %in% genes$production$export)
  expect_true("BSG" %in% genes$production$export)

  # Degradation
  expect_true("LDHB" %in% genes$degradation$enzymes)
  expect_true("LDHD" %in% genes$degradation$enzymes)

  # Direct sensing
  expect_equal(genes$direct_sensing$receptor, "HCAR1")

  # Indirect sensing - proton GPCRs only
  expect_true(all(c("GPR4", "GPR65", "GPR68", "GPR132") %in%
                  genes$indirect_sensing$proton_receptors))
  expect_false("HTR2B" %in% genes$indirect_sensing$proton_receptors)
  expect_false("TLR9" %in% genes$indirect_sensing$proton_receptors)

  # Uptake
  expect_true("SLC16A1" %in% genes$uptake$import)
  expect_true("BSG" %in% genes$uptake$import)
})


# ============================================================
# Test 2: ALDH exclusion verification
# ============================================================
test_that("ALDH family is correctly excluded", {
  genes <- getLactateGenes()
  all_genes <- unlist(genes, use.names = FALSE)

  aldh_genes <- c("ALDH1A1", "ALDH1A2", "ALDH1A3", "ALDH2", "ALDH3A1")
  expect_false(any(aldh_genes %in% all_genes))
})


# ============================================================
# Test 3: getLactateGenes category filtering
# ============================================================
test_that("getLactateGenes returns correct categories", {
  # All categories
  all_genes <- getLactateGenes("all")
  expect_type(all_genes, "list")
  expect_equal(length(all_genes), 5)

  # Production only
  prod <- getLactateGenes("production")
  expect_true("synthesis" %in% names(prod))
  expect_true("export" %in% names(prod))

  # Direct sensing only
  direct <- getLactateGenes("direct_sensing")
  expect_true("receptor" %in% names(direct))

  # Invalid category
  expect_error(getLactateGenes("invalid"))
})


# ============================================================
# Test 4: Create test data with lactate genes
# ============================================================
create_lactate_test_data <- function(n_cells = 100, seed = 42) {
  set.seed(seed)

  # Get lactate genes
  gene_sets <- getLactateGenes()
  lactate_genes <- unique(c(
    gene_sets$production$synthesis,
    gene_sets$production$export,
    gene_sets$degradation$enzymes,
    gene_sets$direct_sensing$receptor,
    gene_sets$indirect_sensing$proton_receptors,
    gene_sets$uptake$import
  ))

  n_genes <- length(lactate_genes) + 50

  # Create expression matrix
  mat <- matrix(
    rpois(n_genes * n_cells, lambda = 2),
    nrow = n_genes,
    ncol = n_cells
  )

  # Name genes
  rownames(mat) <- c(lactate_genes, paste0("OtherGene", 1:50))
  colnames(mat) <- paste0("Cell", seq_len(n_cells))

  # Add cell type-specific expression patterns
  # TypeA (tumor-like): high LDHA, high SLC16A3 (MCT4)
  mat["LDHA", 1:50] <- mat["LDHA", 1:50] + 10
  mat["SLC16A3", 1:50] <- mat["SLC16A3", 1:50] + 8

  # TypeB (immune-like): high GPR65, high GPR132
  mat["GPR65", 51:100] <- mat["GPR65", 51:100] + 8
  mat["GPR132", 51:100] <- mat["GPR132", 51:100] + 6

  # Convert to sparse
  mat <- Matrix::Matrix(mat, sparse = TRUE)

  # Create metadata
  meta <- data.frame(
    cell_type = rep(c("TumorLike", "ImmuneLike"), each = 50),
    row.names = colnames(mat),
    stringsAsFactors = FALSE
  )

  list(expr = mat, meta = meta)
}


# ============================================================
# Test 5: inferLactateSignaling basic execution
# ============================================================
test_that("inferLactateSignaling runs without error", {
  skip_if_not_installed("Matrix")

  test_data <- create_lactate_test_data()
  obj <- createScMetaLink(test_data$expr, test_data$meta, "cell_type", min_cells = 5)

  # Run with no permutations for speed
  result <- inferLactateSignaling(obj, n_permutations = 0, verbose = FALSE)

  expect_s4_class(result, "scMetaLink")
  expect_true(!is.null(result@parameters$lactate_signaling))
})


# ============================================================
# Test 6: Lactate signaling results structure
# ============================================================
test_that("inferLactateSignaling returns correct structure", {
  skip_if_not_installed("Matrix")

  test_data <- create_lactate_test_data()
  obj <- createScMetaLink(test_data$expr, test_data$meta, "cell_type", min_cells = 5)
  obj <- inferLactateSignaling(obj, n_permutations = 0, verbose = FALSE)

  lactate_res <- obj@parameters$lactate_signaling

  # Check all expected components
  expect_true("production" %in% names(lactate_res))
  expect_true("direct_sensing" %in% names(lactate_res))
  expect_true("indirect_sensing" %in% names(lactate_res))
  expect_true("direct_communication" %in% names(lactate_res))
  expect_true("indirect_communication" %in% names(lactate_res))
  expect_true("combined_communication" %in% names(lactate_res))
  expect_true("gene_contributions" %in% names(lactate_res))
  expect_true("parameters" %in% names(lactate_res))

  # Check dimensions
  n_types <- 2
  expect_equal(length(lactate_res$production), n_types)
  expect_equal(nrow(lactate_res$direct_communication), n_types)
  expect_equal(ncol(lactate_res$direct_communication), n_types)
})


# ============================================================
# Test 7: Direct pathway only
# ============================================================
test_that("inferLactateSignaling works with direct only", {
  skip_if_not_installed("Matrix")

  test_data <- create_lactate_test_data()
  obj <- createScMetaLink(test_data$expr, test_data$meta, "cell_type", min_cells = 5)
  obj <- inferLactateSignaling(obj,
                                include_direct = TRUE,
                                include_indirect = FALSE,
                                n_permutations = 0,
                                verbose = FALSE)

  lactate_res <- obj@parameters$lactate_signaling

  expect_true(!is.null(lactate_res$direct_communication))
  expect_null(lactate_res$indirect_communication)
})


# ============================================================
# Test 8: Indirect pathway only
# ============================================================
test_that("inferLactateSignaling works with indirect only", {
  skip_if_not_installed("Matrix")

  test_data <- create_lactate_test_data()
  obj <- createScMetaLink(test_data$expr, test_data$meta, "cell_type", min_cells = 5)
  obj <- inferLactateSignaling(obj,
                                include_direct = FALSE,
                                include_indirect = TRUE,
                                n_permutations = 0,
                                verbose = FALSE)

  lactate_res <- obj@parameters$lactate_signaling

  expect_null(lactate_res$direct_communication)
  expect_true(!is.null(lactate_res$indirect_communication))
})


# ============================================================
# Test 9: MCT direction logic validation
# ============================================================
test_that("MCT transporters are correctly categorized by direction", {
  genes <- getLactateGenes()

  # MCT4 (SLC16A3) should be in export (production)
  expect_true("SLC16A3" %in% genes$production$export)

  # MCT1 (SLC16A1) should be in import (uptake)
  expect_true("SLC16A1" %in% genes$uptake$import)
})


# ============================================================
# Test 10: Helper functions work correctly
# ============================================================
test_that("getTopLactateProducers works correctly", {
  skip_if_not_installed("Matrix")

  test_data <- create_lactate_test_data()
  obj <- createScMetaLink(test_data$expr, test_data$meta, "cell_type", min_cells = 5)
  obj <- inferLactateSignaling(obj, n_permutations = 0, verbose = FALSE)

  top_producers <- getTopLactateProducers(obj, top_n = 2)

  expect_s3_class(top_producers, "data.frame")
  expect_equal(nrow(top_producers), 2)
  expect_true("cell_type" %in% names(top_producers))
  expect_true("production_score" %in% names(top_producers))
  expect_true("rank" %in% names(top_producers))

  # TumorLike should have higher production due to high LDHA
  expect_equal(top_producers$cell_type[1], "TumorLike")
})


# ============================================================
# Test 11: getTopLactateSensors works correctly
# ============================================================
test_that("getTopLactateSensors works correctly", {
  skip_if_not_installed("Matrix")

  test_data <- create_lactate_test_data()
  obj <- createScMetaLink(test_data$expr, test_data$meta, "cell_type", min_cells = 5)
  obj <- inferLactateSignaling(obj, n_permutations = 0, verbose = FALSE)

  # Test indirect pathway
  top_sensors <- getTopLactateSensors(obj, pathway = "indirect", top_n = 2)

  expect_s3_class(top_sensors, "data.frame")
  expect_equal(nrow(top_sensors), 2)

  # ImmuneLike should have higher indirect sensing due to high GPR65
  expect_equal(top_sensors$cell_type[1], "ImmuneLike")
})


# ============================================================
# Test 12: getLactateSignalingSummary works correctly
# ============================================================
test_that("getLactateSignalingSummary works correctly", {
  skip_if_not_installed("Matrix")

  test_data <- create_lactate_test_data()
  obj <- createScMetaLink(test_data$expr, test_data$meta, "cell_type", min_cells = 5)
  obj <- inferLactateSignaling(obj, n_permutations = 0, verbose = FALSE)

  summary <- getLactateSignalingSummary(obj, pathway = "combined", top_n = 5)

  expect_s3_class(summary, "data.frame")
  expect_true("sender" %in% names(summary))
  expect_true("receiver" %in% names(summary))
  expect_true("communication_score" %in% names(summary))
})


# ============================================================
# Test 13: checkLactateGenes works correctly
# ============================================================
test_that("checkLactateGenes works correctly", {
  skip_if_not_installed("Matrix")

  test_data <- create_lactate_test_data()
  obj <- createScMetaLink(test_data$expr, test_data$meta, "cell_type", min_cells = 5)

  # Should run without error and return invisibly
  result <- checkLactateGenes(obj)

  expect_s3_class(result, "data.frame")
  expect_true("category" %in% names(result))
  expect_true("gene" %in% names(result))
  expect_true("available" %in% names(result))
})


# ============================================================
# Test 14: Input validation
# ============================================================
test_that("inferLactateSignaling validates inputs correctly", {
  skip_if_not_installed("Matrix")

  test_data <- create_lactate_test_data()
  obj <- createScMetaLink(test_data$expr, test_data$meta, "cell_type", min_cells = 5)

  # Both pathways FALSE should error
  expect_error(
    inferLactateSignaling(obj, include_direct = FALSE, include_indirect = FALSE),
    "At least one"
  )

  # Invalid method should error
  expect_error(
    inferLactateSignaling(obj, method = "invalid"),
    "method must be one of"
  )

  # Invalid comm_method should error
  expect_error(
    inferLactateSignaling(obj, comm_method = "invalid"),
    "comm_method must be one of"
  )
})


# ============================================================
# Test 15: Permutation test (basic)
# ============================================================
test_that("permutation test produces valid p-values", {
  skip_if_not_installed("Matrix")
  skip_on_cran()

  test_data <- create_lactate_test_data()
  obj <- createScMetaLink(test_data$expr, test_data$meta, "cell_type", min_cells = 5)
  obj <- inferLactateSignaling(obj, n_permutations = 10, verbose = FALSE)

  lactate_res <- obj@parameters$lactate_signaling

  # Check p-values exist and are valid
  expect_true("pvalues" %in% names(lactate_res))
  expect_true(!is.null(lactate_res$pvalues$direct))
  expect_true(!is.null(lactate_res$pvalues$indirect))

  # P-values should be between 0 and 1
  expect_true(all(lactate_res$pvalues$direct >= 0 & lactate_res$pvalues$direct <= 1))
  expect_true(all(lactate_res$pvalues$indirect >= 0 & lactate_res$pvalues$indirect <= 1))
})


# ============================================================
# Test 16: Biological validation - tumor-immune axis
# ============================================================
test_that("indirect pathway captures tumor-immune communication", {
  skip_if_not_installed("Matrix")

  test_data <- create_lactate_test_data()
  obj <- createScMetaLink(test_data$expr, test_data$meta, "cell_type", min_cells = 5)
  obj <- inferLactateSignaling(obj, n_permutations = 0, verbose = FALSE)

  lactate_res <- obj@parameters$lactate_signaling

  # TumorLike should have high production
  expect_gt(lactate_res$production["TumorLike"], lactate_res$production["ImmuneLike"])

  # ImmuneLike should have high indirect sensing
  expect_gt(lactate_res$indirect_sensing["ImmuneLike"],
            lactate_res$indirect_sensing["TumorLike"])

  # Tumor -> Immune communication via indirect pathway should be strong
  indirect_comm <- lactate_res$indirect_communication
  expect_true(is.matrix(indirect_comm))

  # This specific communication should be positive
  expect_gt(indirect_comm["TumorLike", "ImmuneLike"], 0)
})


# ============================================================
# Test 17: Proton receptor weights
# ============================================================
test_that("proton receptor weights are correctly defined", {
  genes <- getLactateGenes()

  # Check weights exist
  expect_true("weights" %in% names(genes$indirect_sensing))

  weights <- genes$indirect_sensing$weights

  # GPR4, GPR65, GPR68 should have weight 1.0
  expect_equal(unname(weights["GPR4"]), 1.0)
  expect_equal(unname(weights["GPR65"]), 1.0)
  expect_equal(unname(weights["GPR68"]), 1.0)

  # GPR132 should have lower weight (0.5) due to weak pH sensitivity
  expect_equal(unname(weights["GPR132"]), 0.5)
})
