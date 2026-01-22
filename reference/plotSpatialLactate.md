# Plot Spatial Lactate Signaling

Visualizes spatial distribution of lactate production and sensing
scores.

## Usage

``` r
plotSpatialLactate(object, type = "production", point_size = 2, title = NULL)
```

## Arguments

- object:

  scMetaLink object with spatial lactate signaling results

- type:

  Character. What to plot: "production", "direct_sensing",
  "indirect_sensing", or "all". Default "production".

- point_size:

  Numeric. Size of spots. Default 2.

- title:

  Character. Plot title. Default auto-generated.

## Value

A ggplot2 object or base R plot
