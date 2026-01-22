# Calculate Lactate Production Score

Calculates lactate production potential for each cell type based on:
Production = mean(Synthesis) \* mean(Export) - alpha \*
mean(Degradation)

## Usage

``` r
.calc_lactate_production(
  gene_scores,
  gene_sets,
  consider_export = TRUE,
  consider_degradation = TRUE,
  degradation_weight = 0.3
)
```

## Arguments

- gene_scores:

  Matrix of gene scores (genes x cell_types)

- gene_sets:

  Lactate gene sets from .get_lactate_gene_sets()

- consider_export:

  Logical. Weight by export potential

- consider_degradation:

  Logical. Subtract degradation

- degradation_weight:

  Numeric. Weight for degradation subtraction (default 0.3)

## Value

Named numeric vector of production scores per cell type
