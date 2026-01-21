# Permutation Test - Consistent with inferProduction/inferSensing

This function replicates EXACTLY the same computation as inferProduction
and inferSensing to ensure statistical validity

## Usage

``` r
.run_permutation_test_consistent(
  object,
  comm_scores,
  common_mets,
  n_permutations,
  n_cores,
  method,
  min_production,
  min_sensing,
  population.size = FALSE,
  verbose
)
```
