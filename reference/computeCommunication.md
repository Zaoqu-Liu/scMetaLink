# Compute Metabolite-Mediated Cell Communication

Calculate communication scores between all cell type pairs mediated by
metabolites. Communication strength represents the potential for signal
transmission from sender to receiver cells via metabolites.

## Usage

``` r
computeCommunication(
  object,
  method = "geometric",
  min_production = 0.1,
  min_sensing = 0.1,
  population.size = FALSE,
  n_permutations = 100,
  n_cores = 1,
  seed = 42,
  verbose = TRUE
)
```

## Arguments

- object:

  A scMetaLink object with production and sensing scores

- method:

  Character. Communication score method: "geometric" (default),
  "product", "harmonic"

- min_production:

  Numeric. Minimum production score threshold (0-1)

- min_sensing:

  Numeric. Minimum sensing score threshold (0-1)

- population.size:

  Logical. Whether to weight communication by cell type population
  sizes. When TRUE, communication strength is scaled by the relative
  abundance of sender and receiver cell types, reflecting the biological
  reality that larger populations contribute more to overall tissue
  signaling.

- n_permutations:

  Integer. Number of permutations for significance testing (0 to skip)

- n_cores:

  Integer. Number of cores for parallel computing

- seed:

  Integer. Random seed for reproducibility

- verbose:

  Logical. Print progress messages

## Value

Updated scMetaLink object with communication_scores and p-values

## Details

The communication score combines production and sensing capabilities:

- geometric: sqrt(production \* sensing) - balanced measure

- product: production \* sensing - emphasizes strong bilateral signals

- harmonic: 2\*production\*sensing/(production+sensing) - penalizes
  imbalance

When `population.size = TRUE`, scores are multiplied by:
sqrt(n_sender/n_total \* n_receiver/n_total), which accounts for the
relative contribution of each cell type to tissue-level communication.
