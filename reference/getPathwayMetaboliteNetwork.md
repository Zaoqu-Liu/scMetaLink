# Get Pathway-Metabolite Network

Extract pathway-metabolite relationships for network visualization

## Usage

``` r
getPathwayMetaboliteNetwork(object, pathways = NULL, top_n = 10)
```

## Arguments

- object:

  scMetaLink object

- pathways:

  Character vector. Specific pathways to include (NULL for top pathways)

- top_n:

  Integer. Number of top pathways if pathways is NULL

## Value

data.frame with pathway-metabolite edges
