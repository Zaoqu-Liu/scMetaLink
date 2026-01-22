# Plot Spatial Distance Distribution

Visualize the distribution of spatial distances between spots

## Usage

``` r
plotSpatialDistanceDistribution(
  object,
  by_celltype = FALSE,
  max_distance = NULL,
  bins = 50
)
```

## Arguments

- object:

  A scMetaLink object with spatial information

- by_celltype:

  Logical. Show distances stratified by cell type pairs

- max_distance:

  Numeric. Maximum distance to show (NULL for all)

- bins:

  Integer. Number of histogram bins

## Value

A ggplot2 object
