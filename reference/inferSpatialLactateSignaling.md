# Infer Spatial Lactate Signaling

Infers lactate signaling in spatial transcriptomics data with
distance-weighted communication scores. Supports both direct (HCAR1) and
indirect (proton-sensing GPCRs) pathways with spatial context.

## Usage

``` r
inferSpatialLactateSignaling(
  object,
  max_distance = 200,
  distance_decay = "gaussian",
  sigma = 50,
  include_direct = TRUE,
  include_indirect = TRUE,
  method = "combined",
  comm_method = "geometric",
  aggregate_by = "cell_type",
  min_production = 0,
  min_sensing = 0,
  normalize = TRUE,
  verbose = TRUE
)
```

## Arguments

- object:

  A spatial scMetaLink object (created with createScMetaLinkFromSpatial)

- max_distance:

  Numeric. Maximum communication distance in micrometers. Spot pairs
  beyond this distance are considered non-interacting. Default 200 um.

- distance_decay:

  Character. Distance decay function:

  - "gaussian": Gaussian decay exp(-d^2/(2\*sigma^2)) (default)

  - "exponential": Exponential decay exp(-d/sigma)

  - "linear": Linear decay max(0, 1 - d/max_distance)

  - "none": No distance weighting, use cell type level aggregation

- sigma:

  Numeric. Decay parameter for gaussian/exponential (in um). Default 50
  um. Literature suggests lactate has medium-range diffusion (~50-80
  um).

- include_direct:

  Logical. Include direct lactate-HCAR1 signaling. Default TRUE.

- include_indirect:

  Logical. Include indirect lactate-H+-GPCR signaling. Default TRUE.

- method:

  Character. Scoring method: "combined", "mean", or "proportion".
  Default "combined".

- comm_method:

  Character. Communication score method: "geometric", "product",
  "harmonic". Default "geometric".

- aggregate_by:

  Character. Aggregation level:

  - "cell_type": Aggregate by cell type (default)

  - "spot": Return spot-level scores (memory intensive)

  - "both": Return both levels

- min_production:

  Numeric. Minimum production score threshold. Default 0.

- min_sensing:

  Numeric. Minimum sensing score threshold. Default 0.

- normalize:

  Logical. Normalize scores. Default TRUE.

- verbose:

  Logical. Print progress messages. Default TRUE.

## Value

Updated scMetaLink object with spatial_lactate_signaling results stored
in the parameters slot, containing:

- spot_production:

  Spot-level production scores

- spot_direct_sensing:

  Spot-level direct sensing scores

- spot_indirect_sensing:

  Spot-level indirect sensing scores

- celltype_communication:

  Cell type level communication (if aggregate_by includes "cell_type")

- spot_communication:

  Spot level communication (if aggregate_by includes "spot")

- parameters:

  Analysis parameters

## Examples

``` r
# \donttest{
data(st_expr)
data(st_meta)
data(st_scalefactors)

# Create spatial object
obj <- createScMetaLinkFromSpatial(
  expression_data = st_expr,
  spatial_coords = st_meta[, c("x", "y")],
  cell_meta = st_meta,
  cell_type_column = "cell_type",
  scale_factors = st_scalefactors
)
#> Created spatial scMetaLink object with 4284 genes, 1000 spots, 6 cell types

# Run spatial lactate analysis
obj <- inferSpatialLactateSignaling(obj)
#> Inferring spatial lactate signaling...
#>   Spots: 1000
#>   Cell types: 6
#>   Distance decay: gaussian (sigma=50 um, max=200 um)
#>   Available lactate genes: 17/17
#>   Calculating spot-level expression scores...
#>   Calculating spot-level production scores...
#>   Calculating spot-level direct sensing scores...
#>   Calculating spot-level indirect sensing scores...
#>   Computing spatial distance weights...
#>   Aggregating by cell type...
#> Done!
#>   Spot production scores: 1000 spots
#>   Cell type communication: 6 x 6

# Access results
spatial_results <- obj@parameters$spatial_lactate_signaling
# }
```
