# List Available Metabolites

Get a list of all metabolites in the database with their properties

## Usage

``` r
listMetabolites(type = "all")
```

## Arguments

- type:

  Character. Filter by interaction type: "all", "signaling" (lr only),
  or "metabolic" (pd only)

## Value

data.frame with metabolite information

## Examples

``` r
# \donttest{
# Get all metabolites
mets <- listMetabolites()
head(mets)
#> # A tibble: 6 × 4
#>   hmdb        metabolite        pubchem metabolite_subclass                     
#>   <chr>       <chr>             <chr>   <chr>                                   
#> 1 HMDB0001448 Sulfate           1117    Non-metal sulfates                      
#> 2 HMDB0000208 Oxoglutaric acid  51      Gamma-keto acids and derivatives        
#> 3 HMDB0000464 Calcium           271     NA                                      
#> 4 HMDB0000131 Glycerol          753     Carbohydrates and carbohydrate conjugat…
#> 5 HMDB0000429 17alpha-Estradiol 68570   Estrane steroids                        
#> 6 HMDB0000990 Acetaldehyde      177     Carbonyl compounds                      

# Get only signaling metabolites (those with receptors)
signaling_mets <- listMetabolites(type = "signaling")
# }
```
