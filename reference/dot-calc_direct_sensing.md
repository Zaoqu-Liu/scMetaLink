# Calculate Direct Lactate Sensing Score (HCAR1)

Calculates direct lactate sensing capability based on HCAR1 expression,
optionally weighted by uptake transporter expression.

DirectSensing = HCAR1_expression + beta \* mean(Uptake_transporters)

## Usage

``` r
.calc_direct_sensing(
  gene_scores,
  gene_sets,
  include_uptake = TRUE,
  uptake_weight = 0.5
)
```

## Arguments

- gene_scores:

  Matrix of gene scores (genes x cell_types)

- gene_sets:

  Lactate gene sets from .get_lactate_gene_sets()

- include_uptake:

  Logical. Include MCT uptake transporters

- uptake_weight:

  Numeric. Weight for uptake contribution (default 0.5)

## Value

Named numeric vector of direct sensing scores per cell type
