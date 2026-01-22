# Create scMetaLink Object from Spatial Data

Initialize a scMetaLink object from spatial transcriptomics data with
spatial coordinate information.

## Usage

``` r
createScMetaLinkFromSpatial(
  expression_data,
  spatial_coords,
  cell_meta,
  cell_type_column = "cell_type",
  scale_factors = NULL,
  min_cells = 5
)
```

## Arguments

- expression_data:

  A matrix or dgCMatrix of normalized expression values (genes x spots)

- spatial_coords:

  A data.frame or matrix with spatial coordinates (spots x 2). Must have
  row names matching column names of expression_data.

- cell_meta:

  A data.frame containing spot metadata (e.g., cell type from
  deconvolution)

- cell_type_column:

  Character. Column name in cell_meta containing cell type annotations

- scale_factors:

  List. Optional scale factors for coordinate conversion. Should contain
  'pixels_per_um' for distance calculations. \*\*IMPORTANT\*\*: For
  Visium data, coordinates are in pixels, not micrometers. Without
  correct scale_factors, distance-based parameters will be wrong.

- min_cells:

  Integer. Minimum number of spots per cell type (default: 5)

## Value

A scMetaLink object with spatial information

## Details

This function creates a scMetaLink object with spatial information
stored in additional slots. The spatial coordinates are used for
distance-weighted communication analysis.

\*\*Important: Cell Type Annotation for Visium Data\*\*

For 10x Visium data, each spot (55 micrometer diameter) typically
contains 1-10 cells of potentially different types. Therefore, cell type
annotations should ideally come from deconvolution methods such as:

- RCTD (spacexr package)

- cell2location

- SPOTlight

- CARD

If using dominant cell type assignment per spot, be aware that this is a
simplification and may miss important cell type heterogeneity within
spots.

\*\*Important: Coordinate Units\*\*

For 10x Visium data, coordinates are typically in pixels. You MUST
provide scale_factors with 'pixels_per_um' to convert to micrometers for
biologically meaningful distance calculations. Without this, sigma=50
will be interpreted as 50 pixels instead of 50 micrometers.

## Examples

``` r
# \donttest{
data(st_expr)
data(st_meta)
data(st_scalefactors)

obj <- createScMetaLinkFromSpatial(
  expression_data = st_expr,
  spatial_coords = st_meta[, c("x", "y")],
  cell_meta = st_meta,
  cell_type_column = "cell_type",
  scale_factors = st_scalefactors
)
#> Created spatial scMetaLink object with 4284 genes, 1000 spots, 6 cell types
# }
```
