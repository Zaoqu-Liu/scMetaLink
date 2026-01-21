# Filter Significant Interactions

Filter Significant Interactions

## Usage

``` r
filterSignificantInteractions(
  object,
  pvalue_threshold = 0.05,
  adjust_method = "metabolite_stratified",
  min_score = 0
)
```

## Arguments

- object:

  scMetaLink object

- pvalue_threshold:

  Numeric. P-value threshold (default: 0.05)

- adjust_method:

  Character. Multiple testing correction method: "BH", "bonferroni",
  "holm", "none", or "metabolite_stratified" (recommended, performs BH
  within each metabolite)

- min_score:

  Numeric. Minimum communication score (default: 0)

## Value

Updated scMetaLink object

## Details

The `adjust_method` parameter controls how p-values are corrected for
multiple testing:

- "metabolite_stratified" (recommended): Performs BH correction within
  each metabolite, then combines results. This is less conservative than
  global correction and more biologically meaningful since metabolites
  are independent biological signals.

- "BH", "bonferroni", "holm": Standard global correction methods. Can be
  very conservative when testing many interactions.

- "none": No correction, uses raw p-values. Use with caution.
