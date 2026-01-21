# List Available Genes

Get a list of all genes in the database with their roles

## Usage

``` r
listGenes(role = "all")
```

## Arguments

- role:

  Character. Filter by role: "all", "enzyme", "receptor", "transporter"

## Value

data.frame with gene information

## Examples

``` r
# \donttest{
# Get all genes
genes <- listGenes()

# Get only receptor genes
receptors <- listGenes(role = "receptor")
# }
```
