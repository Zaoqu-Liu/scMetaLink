# CRC Example Expression Data

Example single-cell RNA-seq expression matrix from colorectal cancer

## Usage

``` r
crc_expr
```

## Format

A sparse dgCMatrix with 4,210 genes (rows) x 2,850 cells (columns). Only
genes present in MetalinksDB are included to reduce file size.

## Source

CellScope package example data

## Details

This is a subset of colorectal cancer single-cell data containing:

- Tumor cells: Tumor Epithelial (600 cells)

- Immune cells: T (500), Plasma (250), B (150), TAM (150), Monocyte
  (120), Normal Macrophage (150), Mast (30)

- Stromal cells: CAF (200), Normal Fibroblast (200), Endothelial (100),
  Pericyte (50), SMC (30)

- Epithelial: Normal Epithelial (300)

- Other: Gliacyte (20)

## Examples

``` r
# \donttest{
data(crc_expr)
#> Warning: data set ‘crc_expr’ not found
data(crc_meta)
#> Warning: data set ‘crc_meta’ not found

# Run scMetaLink analysis
obj <- createScMetaLink(crc_expr, crc_meta, "cell_type")
#> Created scMetaLink object with 4210 genes, 2850 cells, 15 cell types
# }
```
