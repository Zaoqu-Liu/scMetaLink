# Run Complete scMetaLink Analysis

One-step function to run the complete metabolite-mediated cell
communication analysis pipeline

## Usage

``` r
runScMetaLink(
  expression_data,
  cell_meta,
  cell_type_column = "cell_type",
  method = "combined",
  min_cells = 10,
  min_pct = 0.1,
  n_permutations = 1000,
  pvalue_threshold = 0.05,
  n_cores = 1,
  verbose = TRUE
)
```

## Arguments

- expression_data:

  A matrix or sparse matrix of normalized expression

- cell_meta:

  A data.frame with cell metadata

- cell_type_column:

  Character. Column name for cell type annotations

- method:

  Character. Expression scoring method: "mean", "proportion", "combined"

- min_cells:

  Integer. Minimum cells per cell type

- min_pct:

  Numeric. Minimum percentage of expressing cells (0-1)

- n_permutations:

  Integer. Number of permutations for significance testing

- pvalue_threshold:

  Numeric. P-value threshold for significance

- n_cores:

  Integer. Number of cores for parallel computing

- verbose:

  Logical. Print progress messages

## Value

A scMetaLink object with all analysis completed

## Details

This function runs the complete scMetaLink pipeline:

1.  Create scMetaLink object from expression data

2.  Infer metabolite production potential (MPP)

3.  Infer metabolite sensing capability (MSC)

4.  Compute cell-cell communication scores

5.  Perform permutation-based significance testing

6.  Filter significant interactions

7.  Aggregate by pathway

## Examples

``` r
# \donttest{
# Run complete analysis
result <- runScMetaLink(
  expression_data = expr_matrix,
  cell_meta = cell_metadata,
  cell_type_column = "cell_type",
  n_permutations = 1000
)
#> === Step 1/6: Creating scMetaLink object ===
#> Error: object 'expr_matrix' not found

# View significant interactions
sig <- getSignificantInteractions(result)
#> Error in h(simpleError(msg, call)): error in evaluating the argument 'object' in selecting a method for function 'getSignificantInteractions': object 'result' not found
head(sig)
#> Error: object 'sig' not found
# }
```
