# Compute Communication Scores (Vectorized)

Compute Communication Scores (Vectorized)

## Usage

``` r
.compute_communication_scores(
  prod_scores,
  sens_scores,
  method,
  min_production,
  min_sensing,
  pop_weights = NULL
)
```

## Arguments

- prod_scores:

  Production score matrix (metabolites x cell_types)

- sens_scores:

  Sensing score matrix (metabolites x cell_types)

- method:

  Communication method

- min_production:

  Minimum production threshold

- min_sensing:

  Minimum sensing threshold

- pop_weights:

  Optional population weight matrix (cell_types x cell_types)

## Value

3D array of communication scores (sender x receiver x metabolite)
