# Check Lactate Gene Availability

Checks which lactate signaling genes are available in the expression
data. Useful for quality control before running analysis.

## Usage

``` r
checkLactateGenes(object)
```

## Arguments

- object:

  scMetaLink object

## Value

data.frame showing gene availability by category

## Examples

``` r
# \donttest{
data(crc_example)
obj <- createScMetaLink(crc_expr, crc_meta, "cell_type")
#> Created scMetaLink object with 4210 genes, 2850 cells, 15 cell types
checkLactateGenes(obj)
#> Lactate Gene Availability Summary:
#> ==================================
#>          category    available
#>       degradation   2/2 (100%)
#>    direct_sensing   1/1 (100%)
#>  indirect_sensing   4/4 (100%)
#>        production 10/10 (100%)
#>            uptake   3/3 (100%)
#> 
# }
```
