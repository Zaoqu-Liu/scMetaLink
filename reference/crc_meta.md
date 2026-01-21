# CRC Example Cell Metadata

Cell metadata for the CRC example dataset

## Usage

``` r
crc_meta
```

## Format

A data.frame with 2,850 rows (cells) and 3 columns:

- cell_type:

  Cell type annotation (15 types)

- tumor_normal:

  Tumor or Normal tissue origin

- tissue_region:

  Tissue region

## Examples

``` r
# \donttest{
data(crc_meta)
#> Warning: data set ‘crc_meta’ not found
table(crc_meta$cell_type)
#> 
#>                 B               CAF       Endothelial          Gliacyte 
#>               150               200               100                20 
#>              Mast          Monocyte Normal Epithelial Normal Fibroblast 
#>                30               120               300               200 
#> Normal Macrophage          Pericyte            Plasma               SMC 
#>               150                50               250                30 
#>                 T               TAM  Tumor Epithelial 
#>               500               150               600 
# }
```
