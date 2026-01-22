# Run Permutation Test for Lactate Communication

Performs permutation testing by shuffling cell type labels.

## Usage

``` r
.run_lactate_permutation(
  object,
  comm_scores,
  pathway = "direct",
  n_permutations = 100,
  method = "combined",
  comm_method = "geometric",
  min_production = 0,
  min_sensing = 0,
  include_uptake = TRUE,
  use_weights = TRUE,
  verbose = TRUE
)
```

## Arguments

- object:

  scMetaLink object

- comm_scores:

  Communication score matrix

- pathway:

  Character. "direct" or "indirect"

- n_permutations:

  Number of permutations

- method:

  Scoring method

- min_production:

  Minimum production threshold

- min_sensing:

  Minimum sensing threshold

- include_uptake:

  Include uptake in direct sensing

- use_weights:

  Use weights in indirect sensing

- verbose:

  Print progress

## Value

Matrix of p-values
