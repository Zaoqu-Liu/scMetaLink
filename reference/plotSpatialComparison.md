# Plot Spatial Communication Comparison

Compare communication patterns for multiple metabolites

## Usage

``` r
plotSpatialComparison(
  object,
  metabolites,
  type = "production",
  ncol = 2,
  point_size = 1
)
```

## Arguments

- object:

  A scMetaLink object with spatial communication results

- metabolites:

  Character vector. Metabolites to compare

- type:

  Character. "production" or "sensing"

- ncol:

  Integer. Number of columns in facet grid

- point_size:

  Numeric. Size of spot points

## Value

A ggplot2 object
