# Identify Cell Type Specific Metabolites

Find metabolites specifically produced or sensed by cell types

## Usage

``` r
identifyCellTypeSpecificMetabolites(
  object,
  type = "production",
  specificity_threshold = 1.5
)
```

## Arguments

- object:

  scMetaLink object

- type:

  Character. "production" or "sensing"

- specificity_threshold:

  Numeric. Z-score threshold for specificity

## Value

data.frame with cell type specific metabolites
