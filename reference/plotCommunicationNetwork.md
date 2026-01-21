# Plot Communication Network

Create a network visualization of cell-cell communication

## Usage

``` r
plotCommunicationNetwork(
  object,
  metabolite = NULL,
  min_score = 0,
  layout = "fr",
  node_size_by = "degree",
  edge_width_scale = 2,
  colors = NULL
)
```

## Arguments

- object:

  scMetaLink object

- metabolite:

  Character. Specific metabolite (NULL for aggregated)

- min_score:

  Numeric. Minimum score threshold

- layout:

  Character. Network layout algorithm

- node_size_by:

  Character. Size nodes by "degree" or "centrality"

- edge_width_scale:

  Numeric. Scale factor for edge widths

- colors:

  Named vector. Colors for cell types

## Value

A ggplot object
