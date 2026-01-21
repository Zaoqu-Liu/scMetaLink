# Run scMetaLink from Seurat Object

Convenience function to run analysis from Seurat object

## Usage

``` r
runScMetaLinkSeurat(
  seurat_obj,
  cell_type_column = "cell_type",
  assay = "RNA",
  slot = "data",
  ...
)
```

## Arguments

- seurat_obj:

  A Seurat object

- cell_type_column:

  Character. Column in meta.data for cell type

- assay:

  Character. Assay to use

- slot:

  Character. Slot to use ("data" or "counts")

- ...:

  Additional arguments passed to runScMetaLink

## Value

A scMetaLink object
