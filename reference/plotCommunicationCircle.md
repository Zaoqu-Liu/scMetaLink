# Plot Communication Circle

Create a chord diagram of cell-cell communication

## Usage

``` r
plotCommunicationCircle(
  object,
  top_n = 50,
  metabolite = NULL,
  colors = NULL,
  transparency = 0.5,
  title = NULL
)
```

## Arguments

- object:

  scMetaLink object

- top_n:

  Integer. Number of top interactions to show

- metabolite:

  Character. Specific metabolite (NULL for aggregated)

- colors:

  Named vector. Colors for cell types

- transparency:

  Numeric. Link transparency (0-1)

- title:

  Character. Plot title

## Value

Invisibly returns NULL, plots chord diagram
