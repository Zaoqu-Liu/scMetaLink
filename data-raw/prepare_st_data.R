#!/usr/bin/env Rscript
# =============================================================================
# Prepare Spatial Transcriptomics Example Data for scMetaLink
# =============================================================================
# Author: Zaoqu Liu
# Date: 2026-01-22
# Purpose: Create a compact ST dataset for testing spatial communication analysis
# =============================================================================

# Load required libraries
suppressPackageStartupMessages({
  library(Matrix)
  library(scMetaLink)
})

cat("=================================================================\n")
cat("Preparing ST Example Data for scMetaLink\n")
cat("=================================================================\n\n")

# =============================================================================
# 1. Load MetalinksDB to get relevant genes
# =============================================================================
cat("1. Loading MetalinksDB...\n")
db <- scMetaLink:::.load_metalinksdb()
metalinks_genes <- unique(db$proteins$gene_symbol)
cat(sprintf("   MetalinksDB genes: %d\n", length(metalinks_genes)))

# =============================================================================
# 2. Load ST data
# =============================================================================
cat("\n2. Loading ST data...\n")
st_path <- "/Users/liuzaoqu/Desktop/tmp/ST/ST-colon1"

# Read expression matrix (Matrix Market format)
expr_mtx <- Matrix::readMM(file.path(st_path, "matrix.mtx"))
barcodes <- read.table(file.path(st_path, "barcodes.tsv"), stringsAsFactors = FALSE)$V1
genes <- read.table(file.path(st_path, "genes.tsv"), stringsAsFactors = FALSE)

# Use gene symbols (column 2)
gene_symbols <- genes$V2

# Set dimnames
rownames(expr_mtx) <- gene_symbols
colnames(expr_mtx) <- barcodes

cat(sprintf("   Original data: %d genes x %d spots\n", nrow(expr_mtx), ncol(expr_mtx)))

# =============================================================================
# 3. Load spatial coordinates
# =============================================================================
cat("\n3. Loading spatial coordinates...\n")

# Read coordinates
coords <- read.table(
  file.path(st_path, "spatial/coordinates.tsv"),
  header = TRUE,
  row.names = 1,
  stringsAsFactors = FALSE
)

cat(sprintf("   Coordinates loaded: %d spots\n", nrow(coords)))
cat(sprintf("   Columns: %s\n", paste(colnames(coords), collapse = ", ")))

# Read scale factors
scale_factors <- jsonlite::fromJSON(
  file.path(st_path, "spatial/scalefactors_json.json")
)
cat(sprintf("   Spot diameter: %.2f pixels\n", scale_factors$spot_diameter_fullres))

# =============================================================================
# 4. Filter to MetalinksDB genes
# =============================================================================
cat("\n4. Filtering genes...\n")

# Find overlap
overlap_genes <- intersect(gene_symbols, metalinks_genes)
cat(sprintf("   Genes in MetalinksDB: %d (%.1f%%)\n", 
            length(overlap_genes), 100 * length(overlap_genes) / length(gene_symbols)))

# Filter expression matrix
expr_filtered <- expr_mtx[overlap_genes, , drop = FALSE]
cat(sprintf("   Filtered matrix: %d genes x %d spots\n", nrow(expr_filtered), ncol(expr_filtered)))

# =============================================================================
# 5. Subsample spots (for package size)
# =============================================================================
cat("\n5. Subsampling spots...\n")

# Target: ~1000 spots for reasonable package size
set.seed(42)
n_target <- 1000

# Match spots between expression and coordinates
common_spots <- intersect(colnames(expr_filtered), rownames(coords))
cat(sprintf("   Common spots: %d\n", length(common_spots)))

# Subsample if needed
if (length(common_spots) > n_target) {
  selected_spots <- sample(common_spots, n_target)
} else {
  selected_spots <- common_spots
}

cat(sprintf("   Selected spots: %d\n", length(selected_spots)))

# Subset data
st_expr <- expr_filtered[, selected_spots, drop = FALSE]
st_coords <- coords[selected_spots, ]

cat(sprintf("   Final expression matrix: %d genes x %d spots\n", 
            nrow(st_expr), ncol(st_expr)))

# =============================================================================
# 6. Create spatial metadata
# =============================================================================
cat("\n6. Creating spatial metadata...\n")

st_meta <- data.frame(
  row.names = selected_spots,
  # Spatial coordinates (use image coordinates for visualization)
  x = st_coords$imagecol,
  y = st_coords$imagerow,
  # Array coordinates
  array_row = st_coords$row,
  array_col = st_coords$col,
  # Placeholder for cell type (to be filled by user via deconvolution)
  cell_type = "Unknown",
  stringsAsFactors = FALSE
)

# Add some mock cell types based on spatial location for testing
# In real analysis, users would use deconvolution methods
cat("   Adding mock cell types for testing...\n")

# Simple spatial clustering based on coordinates
set.seed(123)
kmeans_result <- kmeans(st_meta[, c("x", "y")], centers = 6)
cell_type_names <- c("Epithelial", "Immune", "Stromal", "Tumor", "Fibroblast", "Endothelial")
st_meta$cell_type <- cell_type_names[kmeans_result$cluster]

cat(sprintf("   Cell type distribution:\n"))
print(table(st_meta$cell_type))

# =============================================================================
# 7. Create scale factors object
# =============================================================================
cat("\n7. Creating scale factors...\n")

st_scalefactors <- list(
  spot_diameter_fullres = scale_factors$spot_diameter_fullres,
  tissue_hires_scalef = scale_factors$tissue_hires_scalef,
  tissue_lowres_scalef = scale_factors$tissue_lowres_scalef,
  # Calculated values
  spot_diameter_um = 55,  # Standard Visium spot diameter
  pixels_per_um = scale_factors$spot_diameter_fullres / 55
)

cat(sprintf("   Spot diameter: 55 μm (%.2f pixels)\n", st_scalefactors$spot_diameter_fullres))
cat(sprintf("   Resolution: %.2f pixels/μm\n", st_scalefactors$pixels_per_um))

# =============================================================================
# 8. Check data size
# =============================================================================
cat("\n8. Checking data sizes...\n")

expr_size <- object.size(st_expr)
meta_size <- object.size(st_meta)
scale_size <- object.size(st_scalefactors)

cat(sprintf("   st_expr: %.2f MB\n", expr_size / 1024^2))
cat(sprintf("   st_meta: %.2f KB\n", meta_size / 1024))
cat(sprintf("   st_scalefactors: %.2f KB\n", scale_size / 1024))
cat(sprintf("   Total: %.2f MB\n", (expr_size + meta_size + scale_size) / 1024^2))

# =============================================================================
# 9. Save data
# =============================================================================
cat("\n9. Saving data...\n")

# Create data directory if needed
data_dir <- "/Users/liuzaoqu/Desktop/develop/scMetaLink/data"
if (!dir.exists(data_dir)) {
  dir.create(data_dir, recursive = TRUE)
}

# Convert to dgCMatrix for better compatibility
st_expr <- as(st_expr, "CsparseMatrix")

# Save as single .rda file
save(st_expr, st_meta, st_scalefactors,
     file = file.path(data_dir, "st_colon.rda"),
     compress = "xz")

cat(sprintf("   Saved to: %s\n", file.path(data_dir, "st_colon.rda")))

# Check final file size
file_size <- file.size(file.path(data_dir, "st_colon.rda"))
cat(sprintf("   File size: %.2f MB\n", file_size / 1024^2))

# =============================================================================
# 10. Summary
# =============================================================================
cat("\n")
cat("=================================================================\n")
cat("Summary\n")
cat("=================================================================\n")
cat(sprintf("Dataset: ST-colon (Colon Spatial Transcriptomics)\n"))
cat(sprintf("Genes: %d (MetalinksDB-filtered)\n", nrow(st_expr)))
cat(sprintf("Spots: %d\n", ncol(st_expr)))
cat(sprintf("Cell types: %d (mock clusters for testing)\n", length(unique(st_meta$cell_type))))
cat(sprintf("Spatial resolution: Visium (55 μm spots)\n"))
cat(sprintf("File size: %.2f MB\n", file_size / 1024^2))
cat("=================================================================\n")

cat("\nST example data prepared successfully!\n")
