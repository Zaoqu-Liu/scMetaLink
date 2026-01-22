# Compute Spot-Level Spatial Communication (Simplified)

Simplified spot-level communication for exploratory analysis. This does
NOT include: degradation adjustment, secretion weighting, Hill
transformation, or full normalization.

## Usage

``` r
.compute_spatial_communication_spot(
  object,
  weight_matrix,
  common_mets,
  comm_method,
  min_production,
  min_sensing,
  verbose
)
```
