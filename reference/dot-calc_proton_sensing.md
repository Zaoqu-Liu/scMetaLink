# Calculate Indirect Proton Sensing Score (GPR4/65/68/132)

Calculates indirect lactate sensing capability through proton-sensing
GPCRs. Uses weighted mean based on pH sensitivity strength of each
receptor.

## Usage

``` r
.calc_proton_sensing(gene_scores, gene_sets, use_weights = TRUE)
```

## Arguments

- gene_scores:

  Matrix of gene scores (genes x cell_types)

- gene_sets:

  Lactate gene sets from .get_lactate_gene_sets()

- use_weights:

  Logical. Use receptor-specific weights

## Value

Named numeric vector of indirect sensing scores per cell type
