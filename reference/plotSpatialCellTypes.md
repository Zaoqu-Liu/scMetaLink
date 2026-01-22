# Plot Spatial Cell Types

Visualize cell type distribution on spatial coordinates

## Usage

``` r
plotSpatialCellTypes(
  object,
  point_size = 1.5,
  alpha = 0.8,
  colors = NULL,
  show_legend = TRUE
)
```

## Arguments

- object:

  A scMetaLink object with spatial information

- point_size:

  Numeric. Size of spot points

- alpha:

  Numeric. Transparency (0-1)

- colors:

  Character vector. Colors for cell types (NULL for default)

- show_legend:

  Logical. Show legend

## Value

A ggplot2 object
