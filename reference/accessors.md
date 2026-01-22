# Accessor Functions for scMetaLink Objects

Functions to extract data from scMetaLink objects.

## Usage

``` r
getProductionScores(object)

# S4 method for class 'scMetaLink'
getProductionScores(object)

getSensingScores(object)

# S4 method for class 'scMetaLink'
getSensingScores(object)

getCommunicationScores(object)

# S4 method for class 'scMetaLink'
getCommunicationScores(object)

getSignificantInteractions(object)

# S4 method for class 'scMetaLink'
getSignificantInteractions(object)

getPathwayAggregated(object)

# S4 method for class 'scMetaLink'
getPathwayAggregated(object)

getParameters(object)

# S4 method for class 'scMetaLink'
getParameters(object)
```

## Arguments

- object:

  A scMetaLink object

## Value

The requested data from the scMetaLink object:

- `getProductionScores`: Matrix of metabolite production scores

- `getSensingScores`: Matrix of metabolite sensing scores

- `getCommunicationScores`: 3D array of communication scores

- `getSignificantInteractions`: data.frame of significant interactions

- `getPathwayAggregated`: data.frame of pathway-aggregated results

- `getParameters`: list of analysis parameters

## Examples

``` r
# \donttest{
data(crc_example)
obj <- createScMetaLink(crc_expr, crc_meta, "cell_type")
#> Created scMetaLink object with 4210 genes, 2850 cells, 15 cell types
obj <- inferProduction(obj)
#> Inferring production potential for 15 cell types...
#>   Calculating cell type expression profiles...
#>   Building metabolite-gene mapping...
#>   Found 15606 production enzyme-metabolite pairs
#>   Computing production scores (matrix multiplication)...
#>   Adjusting for degradation...
#>   Applying secretion potential weights...
#>   Normalizing scores...
#>   Computed production scores for 790 metabolites
#> Done!
prod_scores <- getProductionScores(obj)
# }
```
