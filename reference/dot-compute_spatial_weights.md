# Compute Spatial Weight Matrix

Compute Spatial Weight Matrix

## Usage

``` r
.compute_spatial_weights(
  dist_matrix,
  method,
  threshold,
  sigma,
  lambda,
  k = 6,
  symmetric = TRUE
)
```

## Arguments

- dist_matrix:

  Distance matrix (spots x spots)

- method:

  Spatial weighting method

- threshold:

  Distance threshold in micrometers

- sigma:

  Sigma parameter for Gaussian decay

- lambda:

  Lambda parameter for exponential decay

- k:

  Number of nearest neighbors for knn method

- symmetric:

  Whether to symmetrize the KNN weight matrix
