# Create scMetaLink Object

Initialize a scMetaLink object from expression data and cell metadata

## Usage

``` r
createScMetaLink(
  expression_data,
  cell_meta,
  cell_type_column = "cell_type",
  min_cells = 10,
  min_genes = 200
)
```

## Arguments

- expression_data:

  A matrix or dgCMatrix of normalized expression values (genes x cells)

- cell_meta:

  A data.frame containing cell metadata

- cell_type_column:

  Character. Column name in cell_meta containing cell type annotations

- min_cells:

  Integer. Minimum number of cells per cell type (default: 10)

- min_genes:

  Integer. Minimum number of genes detected per cell (default: 200)

## Value

A scMetaLink object

## Examples

``` r
# \donttest{
# Create from expression matrix
expr_mat <- matrix(rpois(1000, 5), nrow = 100, ncol = 10)
rownames(expr_mat) <- paste0("Gene", 1:100)
colnames(expr_mat) <- paste0("Cell", 1:10)
meta <- data.frame(cell_type = rep(c("TypeA", "TypeB"), each = 5))
rownames(meta) <- colnames(expr_mat)
obj <- createScMetaLink(expr_mat, meta, "cell_type")
#> Error in createScMetaLink(expr_mat, meta, "cell_type"): No cell types with at least 10 cells
# }
```
