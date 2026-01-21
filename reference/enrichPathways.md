# Pathway Enrichment Analysis

Perform hypergeometric test to identify enriched pathways among
significant metabolite-mediated interactions

## Usage

``` r
enrichPathways(
  object,
  pvalue_threshold = 0.05,
  min_overlap = 2,
  adjust_method = "BH"
)
```

## Arguments

- object:

  scMetaLink object with significant interactions

- pvalue_threshold:

  Numeric. P-value cutoff for enrichment (default: 0.05)

- min_overlap:

  Integer. Minimum overlap between significant metabolites and pathway
  (default: 2)

- adjust_method:

  Character. P-value adjustment method (default: "BH")

## Value

data.frame with enriched pathways and statistics

## Examples

``` r
# \donttest{
# Perform pathway enrichment
enriched <- enrichPathways(result)
#> Error: object 'result' not found
head(enriched)
#> Error: object 'enriched' not found
# }
```
