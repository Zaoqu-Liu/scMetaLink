# scMetaLink

<!-- badges: start -->
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![R version](https://img.shields.io/badge/R-%3E%3D4.0-blue.svg)](https://www.r-project.org/)
[![Platform](https://img.shields.io/badge/platform-linux%20%7C%20macOS%20%7C%20windows-lightgrey.svg)](https://github.com/Zaoqu-Liu/scMetaLink)
[![GitHub stars](https://img.shields.io/github/stars/Zaoqu-Liu/scMetaLink?style=social)](https://github.com/Zaoqu-Liu/scMetaLink/stargazers)
<!-- badges: end -->

**scMetaLink** is an R package for inferring metabolite-mediated cell-cell communication from single-cell RNA sequencing data.

📖 **Documentation**: https://Zaoqu-Liu.github.io/scMetaLink/

## Background

Metabolites serve as critical signaling molecules in the cellular microenvironment. While tools like CellChat and CellPhoneDB have advanced our understanding of ligand-receptor interactions, the metabolic dimension of intercellular communication remains underexplored. scMetaLink addresses this gap by integrating enzyme expression, transporter activity, and receptor-metabolite interactions to reconstruct metabolite-mediated communication networks.

## Workflow

<p align="center">
  <img src="man/figures/flowchart.png" width="90%" />
</p>

## Features

- **Metabolite Production Inference**: Estimate metabolite production potential based on biosynthetic enzyme expression
- **Metabolite Sensing Analysis**: Quantify sensing capability via receptors (GPCRs, nuclear receptors) and transporters
- **Statistical Framework**: Permutation-based significance testing with multiple hypothesis correction
- **Pathway Integration**: Aggregate communication patterns at the pathway level
- **Visualization**: Heatmaps, chord diagrams, and network plots

## Database

scMetaLink utilizes MetalinksDB, containing:

- 41,894 metabolite-protein interactions
- 1,128 metabolites
- 4,374 proteins/genes
- 157,741 pathway associations

## Installation

```r
# Install from GitHub
devtools::install_github("Zaoqu-Liu/scMetaLink")
```

## Quick Start

```r
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

```r
result <- runScMetaLinkSeurat(seurat_obj, cell_type_column = "cell_type")
```

## Methodology

### Metabolite Production Potential (MPP)

$$MPP(m, c) = \frac{\sum_{e \in E_m^+} w_e \cdot \text{expr}(e, c)}{|E_m^+|} \cdot S(m)$$

where $E_m^+$ denotes biosynthetic enzymes for metabolite $m$, $w_e$ is the evidence weight, and $S(m)$ represents secretion potential.

### Metabolite Sensing Capability (MSC)

$$MSC(m, c) = \sum_{r \in R_m} w_r \cdot \tau(r) \cdot \text{expr}(r, c)$$

where $R_m$ includes receptors and transporters for metabolite $m$, and $\tau(r)$ is the receptor type weight.

### Communication Score

$$C(s \rightarrow r, m) = \sqrt{MPP(m, s) \times MSC(m, r)}$$

Statistical significance is assessed via permutation testing with FDR correction.

## Tutorials

| Tutorial | Description |
|:--|:--|
| [Quick Start](https://Zaoqu-Liu.github.io/scMetaLink/articles/quick-start.html) | Basic workflow |
| [Theory](https://Zaoqu-Liu.github.io/scMetaLink/articles/theory.html) | Mathematical framework |
| [Production & Sensing](https://Zaoqu-Liu.github.io/scMetaLink/articles/production-sensing.html) | Inference details |
| [Communication](https://Zaoqu-Liu.github.io/scMetaLink/articles/communication.html) | Statistical analysis |
| [Spatial Analysis](https://Zaoqu-Liu.github.io/scMetaLink/articles/spatial-analysis.html) | Spatial transcriptomics |
| [Visualization](https://Zaoqu-Liu.github.io/scMetaLink/articles/visualization.html) | Plotting functions |
| [Applications](https://Zaoqu-Liu.github.io/scMetaLink/articles/applications.html) | Case studies |

## Contact

- **Author**: Zaoqu Liu
- **Email**: liuzaoqu@163.com
- **GitHub**: https://github.com/Zaoqu-Liu/scMetaLink
