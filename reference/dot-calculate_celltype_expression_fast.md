# Calculate Cell Type Expression Profiles

Efficiently calculate mean expression and percentage of expressing cells
for each gene across cell types using vectorized operations.

## Usage

``` r
.calculate_celltype_expression_fast(
  expr_data,
  cell_meta,
  cell_type_col,
  min_expression = 0,
  mean_method = "arithmetic"
)
```

## Arguments

- expr_data:

  Expression matrix (genes x cells)

- cell_meta:

  Cell metadata data.frame

- cell_type_col:

  Column name for cell type annotation

- min_expression:

  Minimum expression threshold

- mean_method:

  Method for calculating mean: "arithmetic" or "trimean"

## Value

List with mean_expr and pct_expr matrices

## Details

The trimean is defined as: (Q1 + 2\*Q2 + Q3) / 4, where Q1, Q2, Q3 are
the 25th, 50th, and 75th percentiles respectively. This provides a
robust estimate of central tendency that is less sensitive to outliers
and the high dropout rates common in single-cell RNA-seq data.
