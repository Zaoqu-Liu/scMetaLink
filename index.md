# scMetaLink

**scMetaLink** is an R package for inferring metabolite-mediated
cell-cell communication from single-cell RNA sequencing data.

📖 **Documentation**: <https://Zaoqu-Liu.github.io/scMetaLink/>

## Background

Metabolites serve as critical signaling molecules in the cellular
microenvironment. While tools like CellChat and CellPhoneDB have
advanced our understanding of ligand-receptor interactions, the
metabolic dimension of intercellular communication remains
underexplored. scMetaLink addresses this gap by integrating enzyme
expression, transporter activity, and receptor-metabolite interactions
to reconstruct metabolite-mediated communication networks.

## Features

- **Metabolite Production Inference**: Estimate metabolite production
  potential based on biosynthetic enzyme expression
- **Metabolite Sensing Analysis**: Quantify sensing capability via
  receptors (GPCRs, nuclear receptors) and transporters
- **Statistical Framework**: Permutation-based significance testing with
  multiple hypothesis correction
- **Pathway Integration**: Aggregate communication patterns at the
  pathway level
- **Visualization**: Heatmaps, chord diagrams, and network plots

## Database

scMetaLink utilizes MetalinksDB, containing:

- 41,894 metabolite-protein interactions
- 1,128 metabolites
- 4,374 proteins/genes
- 157,741 pathway associations

## Installation

``` r
# Install from GitHub
devtools::install_github("Zaoqu-Liu/scMetaLink")
```

## Quick Start

``` r
library(scMetaLink)

# Create scMetaLink object
obj <- createScMetaLink(
  expression_data = expr_matrix,
  cell_meta = cell_metadata,
  cell_type_column = "cell_type"
)

# Infer metabolite production and sensing
obj <- inferProduction(obj)
obj <- inferSensing(obj)

# Compute communication with significance testing
obj <- computeCommunication(obj, n_permutations = 1000)
obj <- filterSignificantInteractions(obj, pvalue_threshold = 0.05)

# Visualize results
plotCommunicationHeatmap(obj)
plotCommunicationCircle(obj)
```

For Seurat users:

``` r
result <- runScMetaLinkSeurat(seurat_obj, cell_type_column = "cell_type")
```

## Methodology

### Metabolite Production Potential (MPP)

$$MPP(m,c) = \frac{\sum\limits_{e \in E_{m}^{+}}w_{e} \cdot \text{expr}(e,c)}{\left| E_{m}^{+} \right|} \cdot S(m)$$

where $E_{m}^{+}$ denotes biosynthetic enzymes for metabolite $m$,
$w_{e}$ is the evidence weight, and $S(m)$ represents secretion
potential.

### Metabolite Sensing Capability (MSC)

$$MSC(m,c) = \sum\limits_{r \in R_{m}}w_{r} \cdot \tau(r) \cdot \text{expr}(r,c)$$

where $R_{m}$ includes receptors and transporters for metabolite $m$,
and $\tau(r)$ is the receptor type weight.

### Communication Score

$$C\left( s\rightarrow r,m \right) = \sqrt{MPP(m,s) \times MSC(m,r)}$$

Statistical significance is assessed via permutation testing with FDR
correction.

## Tutorials

| Tutorial                                                                                           | Description            |
|:---------------------------------------------------------------------------------------------------|:-----------------------|
| [Quick Start](https://Zaoqu-Liu.github.io/scMetaLink/articles/01-quick-start.html)                 | Basic workflow         |
| [Theory](https://Zaoqu-Liu.github.io/scMetaLink/articles/02-theory.html)                           | Mathematical framework |
| [Production & Sensing](https://Zaoqu-Liu.github.io/scMetaLink/articles/03-production-sensing.html) | Inference details      |
| [Communication](https://Zaoqu-Liu.github.io/scMetaLink/articles/04-communication.html)             | Statistical analysis   |
| [Visualization](https://Zaoqu-Liu.github.io/scMetaLink/articles/05-visualization.html)             | Plotting functions     |
| [Applications](https://Zaoqu-Liu.github.io/scMetaLink/articles/06-applications.html)               | Case studies           |

## Contact

- **Author**: Zaoqu Liu
- **Email**: <liuzaoqu@163.com>
- **GitHub**: <https://github.com/Zaoqu-Liu/scMetaLink>
