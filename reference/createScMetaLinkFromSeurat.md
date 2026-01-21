# Create scMetaLink from Seurat Object

Initialize a scMetaLink object from a Seurat object

## Usage

``` r
createScMetaLinkFromSeurat(
  seurat_obj,
  cell_type_column = "cell_type",
  assay = "RNA",
  slot = "data",
  min_cells = 10
)
```

## Arguments

- seurat_obj:

  A Seurat object

- cell_type_column:

  Character. Column name in meta.data for cell type

- assay:

  Character. Assay to use (default: "RNA")

- slot:

  Character. Slot to use (default: "data" for normalized data)

- min_cells:

  Integer. Minimum cells per cell type

## Value

A scMetaLink object
