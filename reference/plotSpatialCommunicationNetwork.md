# Plot Spatial Communication Network

Visualize cell-cell communication with spatial context

## Usage

``` r
plotSpatialCommunicationNetwork(
  object,
  metabolite = NULL,
  sender_type = NULL,
  receiver_type = NULL,
  top_n = 20,
  arrow_scale = 1,
  point_size = 5,
  show_labels = TRUE
)
```

## Arguments

- object:

  A scMetaLink object with spatial communication results

- metabolite:

  Character. Metabolite to visualize (NULL for aggregate)

- sender_type:

  Character. Sender cell type (NULL for all)

- receiver_type:

  Character. Receiver cell type (NULL for all)

- top_n:

  Integer. Number of top interactions to show arrows for

- arrow_scale:

  Numeric. Scale factor for arrow thickness

- point_size:

  Numeric. Size of cell type center points

- show_labels:

  Logical. Show cell type labels

## Value

A ggplot2 object
