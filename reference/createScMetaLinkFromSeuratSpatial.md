# Create scMetaLink from Seurat Spatial Object

Initialize a scMetaLink object from a Seurat object with spatial data

## Usage

``` r
createScMetaLinkFromSeuratSpatial(
  seurat_obj,
  cell_type_column = "cell_type",
  assay = "Spatial",
  slot = "data",
  image = NULL,
  min_cells = 5
)
```

## Arguments

- seurat_obj:

  A Seurat object with spatial assay and coordinates

- cell_type_column:

  Character. Column name in meta.data for cell type

- assay:

  Character. Assay to use (default: "Spatial")

- slot:

  Character. Slot to use (default: "data" for normalized data)

- image:

  Character. Name of the spatial image to use (default: first available)

- min_cells:

  Integer. Minimum spots per cell type

## Value

A scMetaLink object with spatial information
