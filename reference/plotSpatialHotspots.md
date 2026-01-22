# Plot Spatial Hotspots

Visualize communication hotspots on spatial coordinates

## Usage

``` r
plotSpatialHotspots(
  object,
  metabolite = NULL,
  type = "sender",
  point_size = 1,
  hotspot_size = 5,
  n_hotspots = 5
)
```

## Arguments

- object:

  A scMetaLink object with spatial communication results

- metabolite:

  Character. Metabolite to analyze (NULL for aggregate)

- type:

  Character. "sender" or "receiver"

- point_size:

  Numeric. Size of spot points

- hotspot_size:

  Numeric. Size of hotspot markers

- n_hotspots:

  Integer. Number of hotspots to highlight

## Value

A ggplot2 object
