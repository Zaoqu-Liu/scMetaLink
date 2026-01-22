# Colon Spatial Transcriptomics Scale Factors

Scale factors for the colon Visium example dataset

## Usage

``` r
st_scalefactors
```

## Format

A list containing:

- spot_diameter_fullres:

  Spot diameter in full resolution pixels (~130)

- tissue_hires_scalef:

  Scale factor for high-res image (~0.12)

- tissue_lowres_scalef:

  Scale factor for low-res image (~0.037)

- spot_diameter_um:

  Spot diameter in micrometers (55 um)

- pixels_per_um:

  Pixels per micrometer (~2.37)

## Details

These scale factors are essential for converting between pixel
coordinates and physical distances (micrometers). The standard Visium
spot diameter is 55 um, which is used to calculate the pixels_per_um
conversion factor.

For spatial communication analysis, distances should be converted to
micrometers to enable biologically meaningful distance thresholds.

## See also

[`st_expr`](https://Zaoqu-Liu.github.io/scMetaLink/reference/st_expr.md),
[`st_meta`](https://Zaoqu-Liu.github.io/scMetaLink/reference/st_meta.md)

## Examples

``` r
# \donttest{
data(st_scalefactors)

# Convert pixel distance to micrometers
pixel_distance <- 500
um_distance <- pixel_distance / st_scalefactors$pixels_per_um
cat("Distance:", um_distance, "um\n")
#> Distance: 211.1615 um
# }
```
