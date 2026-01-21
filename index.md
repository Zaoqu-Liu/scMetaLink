# scMetaLink ![](reference/figures/logo.png)

## Single-Cell Metabolite-Mediated Cell Communication Analysis

**scMetaLink** is a comprehensive R package for inferring
metabolite-mediated cell-cell communication from single-cell
transcriptomic data. While existing tools like CellChat and CellPhoneDB
focus on ligand-receptor interactions, scMetaLink fills a critical gap
by enabling the study of metabolite signaling in the cellular
microenvironment.

### Key Features

🧬 **Metabolite Production Inference**: Infer metabolite production
potential from enzyme expression patterns

🎯 **Metabolite Sensing Analysis**: Quantify metabolite sensing
capability via receptor/transporter expression

📊 **Statistical Framework**: Permutation-based significance testing
with multiple hypothesis correction

🔬 **Pathway Integration**: Aggregate results at the pathway level for
biological interpretation

📈 **Publication-Ready Visualization**: Heatmaps, chord diagrams, and
network plots

### Database

scMetaLink leverages **MetalinksDB**, a comprehensive database
containing: - **41,894** metabolite-protein interactions - **1,128**
metabolites - **4,374** proteins/genes - **157,741** pathway
associations

## Installation

``` r
# Install from GitHub
devtools::install_github("liuzaoqu/scMetaLink")

# Or install from local
install.packages("scMetaLink", repos = NULL, type = "source")
```

## Quick Start

### One-Step Analysis

``` r
library(scMetaLink)

# Run complete analysis pipeline
result <- runScMetaLink(
  expression_data = expr_matrix,     # genes x cells matrix
  cell_meta = cell_metadata,         # data.frame with cell annotations
  cell_type_column = "cell_type",    # column name for cell type
  n_permutations = 1000              # for significance testing
)

# View significant interactions
sig <- getSignificantInteractions(result)
head(sig)
```

### From Seurat Object

``` r
# Direct analysis from Seurat
result <- runScMetaLinkSeurat(
  seurat_obj,
  cell_type_column = "seurat_clusters"
)
```

### Step-by-Step Analysis

``` r
# Create object
obj <- createScMetaLink(expr_matrix, cell_meta, "cell_type")

# Infer metabolite production
obj <- inferProduction(obj)

# Infer metabolite sensing
obj <- inferSensing(obj)

# Compute communication
obj <- computeCommunication(obj, n_permutations = 1000)

# Filter significant
obj <- filterSignificantInteractions(obj, pvalue_threshold = 0.05)

# Pathway analysis
obj <- aggregateByPathway(obj)
```

## Visualization

``` r
# Communication heatmap
plotCommunicationHeatmap(result)

# Chord diagram
plotCommunicationCircle(result)

# Network plot
plotCommunicationNetwork(result)

# Top interactions
plotTopInteractions(result)

# Pathway analysis
plotPathwayCommunication(result)
```

## Methodology

### Metabolite Production Potential (MPP)

$$MPP(m,c) = \frac{\sum\limits_{e \in E_{m}^{+}}w_{e} \cdot expr(e,c)}{\left| E_{m}^{+} \right|} \cdot S(m)$$

Where: - $E_{m}^{+}$: Enzymes producing metabolite $m$ - $w_{e}$:
Confidence weight for enzyme - $S(m)$: Secretion potential

### Metabolite Sensing Capability (MSC)

$$MSC(m,c) = \sum\limits_{r \in R_{m}}w_{r} \cdot \tau(r) \cdot expr(r,c)$$

Where: - $R_{m}$: Receptors for metabolite $m$ - $\tau(r)$: Receptor
type weight

### Communication Score

$$C\left( s\rightarrow r,m \right) = \sqrt{MPP(m,s) \times MSC(m,r)}$$

## Comparison with Existing Tools

| Feature                   | CellChat | CellPhoneDB | **scMetaLink** |
|---------------------------|----------|-------------|----------------|
| Ligand-Receptor           | ✓        | ✓           | ✓              |
| **Metabolite Ligands**    | ✗        | ✗           | **✓**          |
| **Metabolite Production** | ✗        | ✗           | **✓**          |
| **Secretion Potential**   | ✗        | ✗           | **✓**          |
| Pathway Integration       | Limited  | Limited     | **Extensive**  |

## Applications

- **Tumor Microenvironment**: Study metabolic crosstalk between cancer
  and immune cells
- **Immune Metabolism**: Understand metabolite-mediated immune
  regulation
- **Development**: Track metabolic communication during differentiation
- **Disease**: Identify dysregulated metabolic signaling

## Citation

    Liu Z, et al. (2026). scMetaLink: Inferring metabolite-mediated cell-cell
    communication from single-cell transcriptomic data. [Journal Name].

## License

MIT License

## Contact

- **Author**: Zaoqu Liu
- **Email**: <liuzaoqu@163.com>
- **GitHub**: <https://github.com/liuzaoqu/scMetaLink>
