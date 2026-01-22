# Get Lactate Signaling Gene Sets

Returns curated, literature-validated gene sets involved in lactate
signaling pathways. Useful for pathway analysis, visualization, and
validation.

**Gene selection criteria**:

- Synthesis: LDHA family (excludes LDHB which prefers reverse reaction)

- Export: MCT4 (SLC16A3) as primary exporter, plus other MCTs and BSG
  chaperone

- Direct sensing: HCAR1 only (the sole confirmed lactate GPCR)

- Indirect sensing: Classic proton-sensing GPCRs (GPR4, GPR65, GPR68,
  GPR132)

- Uptake: MCT1 (SLC16A1) as primary importer, plus BSG chaperone

## Usage

``` r
getLactateGenes(category = "all")
```

## Arguments

- category:

  Character. Gene category to return:

  - "all" (default): All gene sets

  - "production": Synthesis and export genes

  - "degradation": Lactate consumption enzymes

  - "direct_sensing": HCAR1 receptor

  - "indirect_sensing": Proton-sensing GPCRs

  - "uptake": Import transporters

## Value

A named list of gene vectors. If category is "all", returns the complete
nested list structure. Otherwise, returns the specific category.

## Examples

``` r
# Get all gene sets
genes <- getLactateGenes()
names(genes)
#> [1] "production"       "degradation"      "direct_sensing"   "indirect_sensing"
#> [5] "uptake"          

# Get synthesis enzymes
genes$production$synthesis
#> [1] "LDHA"    "LDHC"    "LDHAL6A" "LDHAL6B"
# [1] "LDHA" "LDHC" "LDHAL6A" "LDHAL6B"

# Get proton-sensing receptors
genes$indirect_sensing$proton_receptors
#> [1] "GPR4"   "GPR65"  "GPR68"  "GPR132"
# [1] "GPR4" "GPR65" "GPR68" "GPR132"

# Get only production genes
prod_genes <- getLactateGenes("production")

# Get direct sensing receptor
getLactateGenes("direct_sensing")$receptor
#> [1] "HCAR1"
# [1] "HCAR1"
```
