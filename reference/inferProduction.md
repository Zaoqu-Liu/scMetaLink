# Infer Metabolite Production Potential

Infer metabolite production potential for each cell type based on enzyme
expression patterns. Production potential reflects a cell type's
capacity to synthesize and secrete metabolites for intercellular
communication.

## Usage

``` r
inferProduction(
  object,
  method = "combined",
  mean_method = "arithmetic",
  min_expression = 0,
  min_pct = 0.1,
  consider_degradation = TRUE,
  consider_secretion = TRUE,
  normalize = TRUE,
  verbose = TRUE
)
```

## Arguments

- object:

  A scMetaLink object

- method:

  Character. Scoring method: "mean", "proportion", or "combined"

- mean_method:

  Character. Method for calculating mean expression: "arithmetic"
  (standard mean) or "trimean" (more robust to outliers and dropout).
  Trimean is recommended for single-cell data with high dropout rates.

- min_expression:

  Numeric. Minimum expression threshold

- min_pct:

  Numeric. Minimum percentage of expressing cells (0-1)

- consider_degradation:

  Logical. Subtract degradation enzyme expression

- consider_secretion:

  Logical. Weight by secretion potential

- normalize:

  Logical. Normalize scores across cell types

- verbose:

  Logical. Print progress messages

## Value

Updated scMetaLink object with production_scores slot filled
