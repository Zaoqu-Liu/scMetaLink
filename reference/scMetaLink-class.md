# scMetaLink Class Definition

S4 class for storing metabolite-mediated cell communication analysis
results

## Usage

``` r
# S4 method for class 'scMetaLink'
show(object)
```

## Arguments

- object:

  scMetaLink object

## Slots

- `expression_data`:

  Matrix. Normalized gene expression matrix (genes x cells)

- `cell_meta`:

  data.frame. Cell metadata with cell type annotations

- `cell_type_column`:

  Character. Column name for cell type in cell_meta

- `production_scores`:

  Matrix. Metabolite production potential scores (metabolites x
  cell_types)

- `sensing_scores`:

  Matrix. Metabolite sensing capability scores (metabolites x
  cell_types)

- `communication_scores`:

  Array. Communication scores (sender x receiver x metabolite)

- `communication_pvalues`:

  Array. P-values from permutation test

- `significant_interactions`:

  data.frame. Filtered significant interactions

- `pathway_aggregated`:

  data.frame. Pathway-level aggregated results

- `parameters`:

  list. Analysis parameters

- `database`:

  list. MetalinksDB reference data
