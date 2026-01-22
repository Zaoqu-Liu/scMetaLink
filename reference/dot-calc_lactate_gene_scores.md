# Calculate Cell Type Expression Profiles for Lactate Genes

Calculates mean expression and percentage of expressing cells for
lactate-related genes across cell types.

## Usage

``` r
.calc_lactate_gene_scores(
  expr_data,
  cell_meta,
  cell_type_col,
  genes,
  min_expression = 0,
  method = "combined"
)
```

## Arguments

- expr_data:

  Expression matrix (genes x cells)

- cell_meta:

  Cell metadata data.frame

- cell_type_col:

  Column name for cell type annotation

- genes:

  Character vector of genes to analyze

- min_expression:

  Minimum expression threshold

- method:

  Scoring method: "combined", "mean", or "proportion"

## Value

Matrix of gene scores (genes x cell_types)
