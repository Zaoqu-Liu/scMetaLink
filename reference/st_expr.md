# Colon Spatial Transcriptomics Expression Data

Example 10x Visium spatial transcriptomics expression matrix from colon
tissue

## Usage

``` r
st_expr
```

## Format

A sparse dgCMatrix with 4,284 genes (rows) x 1,000 spots (columns). Only
genes present in MetalinksDB are included to reduce file size.

## Source

Colon spatial transcriptomics dataset

## Details

This is a subset of colon Visium data containing 1,000 randomly sampled
spots. The data is suitable for testing spatial metabolite communication
analysis.

Technical specifications:

- Platform: 10x Genomics Visium

- Spot diameter: 55 um

- Resolution: ~2.37 pixels/um

- Genes: 4,284 (filtered to MetalinksDB)

- Spots: 1,000 (subsampled from 3,313)

## See also

[`st_meta`](https://Zaoqu-Liu.github.io/scMetaLink/reference/st_meta.md),
[`st_scalefactors`](https://Zaoqu-Liu.github.io/scMetaLink/reference/st_scalefactors.md)

## Examples

``` r
# \donttest{
data(st_expr)
data(st_meta)
data(st_scalefactors)

# Check dimensions
dim(st_expr)
#> [1] 4284 1000

# View spatial coordinates
head(st_meta)
#>                        x     y array_row array_col   cell_type
#> TAGTCCCGGAGACCAC-1  6642 13195        66        48 Endothelial
#> GTGCGTGTATATGAGC-1 10097  8796        41        83     Stromal
#> CCGGTATCTGGCGACT-1 10371 15069        77        85  Fibroblast
#> CCCAAGAATGCACGGT-1 10519  1993         2        88      Immune
#> CGAAACATAGATGGCA-1  9435  3574        11        77      Immune
#> TTCTTATCCGCTGGGT-1 11917 10169        49       101  Fibroblast

# Run spatial scMetaLink analysis
obj <- createScMetaLinkFromSpatial(
  expression_data = st_expr,
  spatial_coords = st_meta[, c("x", "y")],
  cell_meta = st_meta,
  cell_type_column = "cell_type"
)
#> Created spatial scMetaLink object with 4284 genes, 1000 spots, 6 cell types
#> Warning: Coordinates appear to be in pixels (range: 13511). For Visium data, please provide scale_factors with 'pixels_per_um'. Without correct scaling, distance parameters (sigma, threshold) will be interpreted as pixels instead of micrometers, leading to incorrect results. Assuming coordinates are in micrometers for now (pixels_per_um = 1).
# }
```
