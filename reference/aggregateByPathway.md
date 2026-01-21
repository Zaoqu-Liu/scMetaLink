# Aggregate Communication by Pathway

Aggregate metabolite-mediated communication at the pathway level. Due to
the large number of pathway associations, this function uses a
simplified approach focusing on top pathways.

## Usage

``` r
aggregateByPathway(object, top_pathways = 50, min_metabolites = 3)
```

## Arguments

- object:

  A scMetaLink object with significant interactions

- top_pathways:

  Integer. Number of top pathways to analyze (default: 50)

- min_metabolites:

  Integer. Minimum metabolites per pathway (default: 3)

## Value

Updated scMetaLink object with pathway_aggregated slot filled
