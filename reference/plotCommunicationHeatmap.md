# Plot Communication Heatmap

Create a heatmap of cell-cell communication

## Usage

``` r
plotCommunicationHeatmap(
  object,
  metabolite = NULL,
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_values = FALSE,
  colors = NULL,
  title = NULL
)
```

## Arguments

- object:

  scMetaLink object

- metabolite:

  Character. Specific metabolite to show (NULL for aggregated)

- cluster_rows:

  Logical. Cluster rows

- cluster_cols:

  Logical. Cluster columns

- show_values:

  Logical. Show values in cells

- colors:

  Character vector. Color palette

- title:

  Character. Plot title

## Value

A ComplexHeatmap object

## Examples

``` r
# \donttest{
# Plot aggregated communication heatmap
plotCommunicationHeatmap(obj)
#> Error: object 'obj' not found

# Plot specific metabolite
plotCommunicationHeatmap(obj, metabolite = "HMDB0000148")
#> Error: object 'obj' not found
# }
```
