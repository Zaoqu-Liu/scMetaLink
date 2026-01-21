# Theory & Methods

## Overview

scMetaLink models metabolite-mediated cell communication as a two-step
process:

![\*\*Figure 1: scMetaLink Workflow.\*\* Schematic overview of the
analysis pipeline for inferring metabolite-mediated cell-cell
communication.](flowchart.png)

**Figure 1: scMetaLink Workflow.** Schematic overview of the analysis
pipeline for inferring metabolite-mediated cell-cell communication.

## The MetalinksDB Knowledge Base

scMetaLink uses **MetalinksDB**, a comprehensive database of
metabolite-protein interactions.

``` r
library(scMetaLink)

# Load and explore the database
db <- scMetaLink:::.load_metalinksdb()

cat("=== MetalinksDB Statistics ===\n")
#> === MetalinksDB Statistics ===
cat("Metabolites:", nrow(db$metabolites), "\n")
#> Metabolites: 1128
cat("Proteins:", nrow(db$proteins), "\n")
#> Proteins: 4374
cat("Interactions:", nrow(db$edges), "\n")
#> Interactions: 41894
cat("Pathways:", length(unique(db$pathway$pathway)), "\n")
#> Pathways: 36806
```

### Interaction Types

``` r
# Two main interaction types
cat("=== Interaction Types ===\n")
#> === Interaction Types ===
print(table(db$edges$type))
#> 
#>    lr    pd 
#> 11869 30025

cat("\n- lr (ligand-receptor): Direct metabolite-receptor binding\n")
#> 
#> - lr (ligand-receptor): Direct metabolite-receptor binding
cat("- pd (produce-degrade): Enzymatic production/degradation\n")
#> - pd (produce-degrade): Enzymatic production/degradation
```

### Mode of Regulation (MOR)

For `pd` type interactions, MOR indicates the direction:

``` r
pd_edges <- db$edges[db$edges$type == "pd", ]
cat("=== Mode of Regulation for pd interactions ===\n")
#> === Mode of Regulation for pd interactions ===
print(table(pd_edges$mor))
#> 
#>    -1     1 
#> 13144 16881

cat("\n- mor = 1: Enzyme produces/secretes the metabolite\n")
#> 
#> - mor = 1: Enzyme produces/secretes the metabolite
cat("- mor = -1: Enzyme degrades/consumes the metabolite\n")
#> - mor = -1: Enzyme degrades/consumes the metabolite
cat("- mor = 0: Enzyme binds without direction change\n")
#> - mor = 0: Enzyme binds without direction change
```

### Protein Types

``` r
cat("=== Protein Classifications ===\n")
#> === Protein Classifications ===
print(table(db$proteins$protein_type, useNA = "ifany"))
#> 
#> catalytic_receptor             enzyme               gpcr               lgic 
#>                179                918                306                 71 
#>                nhr           other_ic      other_protein        transporter 
#>                 43                 39                 28                393 
#>               vgic               <NA> 
#>                 73               2324

cat("\n")
cat("Receptors: gpcr, lgic, nhr, vgic, catalytic_receptor, other_ic\n")
#> Receptors: gpcr, lgic, nhr, vgic, catalytic_receptor, other_ic
cat("Enzymes: enzyme (or NA in protein_type)\n")
#> Enzymes: enzyme (or NA in protein_type)
cat("Transporters: transporter\n")
#> Transporters: transporter
```

## Mathematical Framework

### 1. Metabolite Production Potential (MPP)

The production potential for metabolite $m$ in cell type $c$ is:

$$\text{MPP}_{m,c} = \frac{1}{\left| E_{m} \right|}\sum\limits_{e \in E_{m}}\text{Score}(e,c) - \alpha \cdot \frac{1}{\left| D_{m} \right|}\sum\limits_{d \in D_{m}}\text{Score}(d,c)$$

Where: - $E_{m}$: Set of enzymes that produce metabolite $m$ - $D_{m}$:
Set of enzymes that degrade metabolite $m$ - $\alpha = 0.5$: Degradation
weight factor - $\text{Score}(g,c)$: Gene expression score for gene $g$
in cell type $c$

#### Gene Expression Scoring

Three methods are available:

| Method       | Formula                                          | Use Case                       |
|--------------|--------------------------------------------------|--------------------------------|
| `mean`       | ${\bar{x}}_{g}$                                  | Simple average expression      |
| `proportion` | $P\left( x_{g} > 0 \right)$                      | Proportion of expressing cells |
| `combined`   | ${\bar{x}}_{g} \times P\left( x_{g} > 0 \right)$ | Balanced (recommended)         |

| Method     | Formula           | Value |
|:-----------|:------------------|------:|
| mean       | mean(expr)        | 0.610 |
| proportion | mean(expr \> 0)   | 0.300 |
| combined   | mean x proportion | 0.183 |

#### Trimean Option

For single-cell data with high dropout, **trimean** provides a robust
alternative:

$$\text{Trimean} = \frac{Q_{1} + 2Q_{2} + Q_{3}}{4}$$

``` r
# Example: Trimean vs Arithmetic Mean
expr_with_outlier <- c(0, 0, 0, 1, 2, 3, 100) # Outlier = 100

cat("Data:", paste(expr_with_outlier, collapse = ", "), "\n")
#> Data: 0, 0, 0, 1, 2, 3, 100
cat("Arithmetic mean:", round(mean(expr_with_outlier), 2), "\n")
#> Arithmetic mean: 15.14
cat("Median:", median(expr_with_outlier), "\n")
#> Median: 1

# Trimean calculation
q <- quantile(expr_with_outlier, c(0.25, 0.5, 0.75))
trimean <- (q[1] + 2 * q[2] + q[3]) / 4
cat("Trimean:", round(trimean, 2), "\n")
#> Trimean: 1.12
cat("\nTrimean is robust to the outlier!\n")
#> 
#> Trimean is robust to the outlier!
```

### 2. Metabolite Sensing Capability (MSC)

The sensing capability for metabolite $m$ in cell type $c$ is:

$$\text{MSC}_{m,c} = \frac{1}{\left| R_{m} \right|}\sum\limits_{r \in R_{m}}w_{r} \cdot \text{Score}(r,c)$$

Where: - $R_{m}$: Set of receptors/transporters for metabolite $m$ -
$w_{r}$: Affinity weight based on interaction confidence score

#### Affinity Weighting

Receptor-metabolite binding affinities from MetalinksDB are used as
weights:

``` r
lr_edges <- db$edges[db$edges$type == "lr", ]
cat("Affinity score distribution (0-1000):\n")
#> Affinity score distribution (0-1000):
summary(lr_edges$combined_score)
#>    Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#>   202.0   287.0   413.0   553.5   921.0   999.0     470
```

#### Hill Function (Optional)

To model receptor saturation kinetics, an optional Hill transformation
can be applied:

$$P = \frac{E^{n}}{K_{h}^{n} + E^{n}}$$

Where: - $E$: Expression level (normalized) - $n$: Hill coefficient
(cooperativity) - $K_{h}$: Half-maximal threshold

``` r
# Hill function visualization
x <- seq(0, 1, 0.01)
hill <- function(x, n, Kh) x^n / (Kh^n + x^n)

plot(x, x,
  type = "l", lty = 2, col = "gray",
  xlab = "Expression (normalized)", ylab = "Response",
  main = "Hill Function Transformation"
)
lines(x, hill(x, n = 1, Kh = 0.5), col = "blue", lwd = 2)
lines(x, hill(x, n = 2, Kh = 0.5), col = "red", lwd = 2)
lines(x, hill(x, n = 0.5, Kh = 0.5), col = "green", lwd = 2)
legend("bottomright",
  legend = c("Linear", "n=1 (Michaelis-Menten)", "n=2 (Cooperative)", "n=0.5"),
  col = c("gray", "blue", "red", "green"),
  lty = c(2, 1, 1, 1), lwd = 2
)
abline(h = 0.5, v = 0.5, lty = 3, col = "gray")
```

![\*\*Figure 1: Hill Function Transformation.\*\* Different Hill
coefficients (n) produce different response curves. n=1 gives standard
Michaelis-Menten kinetics. n\>1 indicates positive cooperativity
(steeper curve). The half-maximal response occurs at
x=Kh.](02-theory_files/figure-html/hill_demo-1.png)

**Figure 1: Hill Function Transformation.** Different Hill coefficients
(n) produce different response curves. n=1 gives standard
Michaelis-Menten kinetics. n\>1 indicates positive cooperativity
(steeper curve). The half-maximal response occurs at x=Kh.

### 3. Communication Score

The communication score from sender $s$ to receiver $r$ via metabolite
$m$:

$$\text{Comm}_{s\rightarrow r}^{m} = f\left( \text{MPP}_{m,s},\text{MSC}_{m,r} \right)$$

Three aggregation methods are available:

| Method      | Formula                                                                 | Properties                          |
|-------------|-------------------------------------------------------------------------|-------------------------------------|
| `geometric` | $\sqrt{\text{MPP} \times \text{MSC}}$                                   | Balanced, penalizes imbalance       |
| `product`   | $\text{MPP} \times \text{MSC}$                                          | Emphasizes strong bilateral signals |
| `harmonic`  | $\frac{2 \times \text{MPP} \times \text{MSC}}{\text{MPP} + \text{MSC}}$ | Strongly penalizes imbalance        |

![\*\*Figure 2: Communication Score Methods.\*\* Heat maps showing how
different methods combine production (x-axis) and sensing (y-axis)
scores. Geometric mean (left) provides balanced weighting. Product
(center) emphasizes strong bilateral signals. Harmonic mean (right)
penalizes imbalanced
communication.](02-theory_files/figure-html/comm_methods-1.png)

**Figure 2: Communication Score Methods.** Heat maps showing how
different methods combine production (x-axis) and sensing (y-axis)
scores. Geometric mean (left) provides balanced weighting. Product
(center) emphasizes strong bilateral signals. Harmonic mean (right)
penalizes imbalanced communication.

### 4. Population Size Correction (Optional)

When `population.size = TRUE`, communication is weighted by cell type
abundance:

$$\text{Comm}_{s\rightarrow r}^{\text{weighted}} = \text{Comm}_{s\rightarrow r} \times \sqrt{\frac{n_{s}}{N} \times \frac{n_{r}}{N}}$$

Where $n_{s},n_{r}$ are cell counts and $N$ is total cells.

**Rationale**: Larger cell populations contribute more to tissue-level
signaling.

## Statistical Framework

### Permutation Test

To assess significance, cell type labels are randomly shuffled:

1.  Shuffle cell type labels $B$ times (default: 100-1000)
2.  Recalculate production, sensing, and communication for each
    permutation
3.  Compute empirical p-value:

$$p = \frac{1 + \sum\limits_{b = 1}^{B}\mathbb{1}\left\lbrack \text{Comm}_{b} \geq \text{Comm}_{\text{obs}} \right\rbrack}{B + 1}$$

### Multiple Testing Correction

Three strategies for handling multiple comparisons:

| Method                  | Description                                        |
|-------------------------|----------------------------------------------------|
| `metabolite_stratified` | BH correction within each metabolite (recommended) |
| `BH`                    | Global Benjamini-Hochberg                          |
| `bonferroni`            | Global Bonferroni (very conservative)              |

**Why `metabolite_stratified`?**

- Metabolites are independent biological signals
- Global correction is overly conservative
- Per-metabolite correction preserves biological interpretability

## Data Processing Pipeline

![](02-theory_files/figure-html/pipeline-1.png)

## Key Assumptions

1.  **Gene expression correlates with protein activity**:
    Enzyme/receptor expression levels reflect functional activity
2.  **Metabolites can diffuse between cells**: Produced metabolites are
    available for sensing by neighboring cells
3.  **Cell type labels are accurate**: Incorrect annotations will affect
    results
4.  **Adequate cell sampling**: Each cell type needs sufficient cells
    for robust estimation

## Comparison with Ligand-Receptor Methods

| Aspect                | Ligand-Receptor (e.g., CellChat) | Metabolite-Mediated (scMetaLink)    |
|-----------------------|----------------------------------|-------------------------------------|
| **Signal type**       | Proteins                         | Small molecules                     |
| **Diffusion range**   | Short (juxtacrine/paracrine)     | Variable (paracrine/endocrine)      |
| **Database**          | CellChatDB, CellPhoneDB          | MetalinksDB                         |
| **Key genes**         | Ligands, receptors               | Enzymes, transporters, receptors    |
| **Signal complexity** | Direct binding                   | Synthesis -\> Secretion -\> Sensing |

## References

1.  Schafer, S., et al. (2023). MetalinksDB: a knowledgebase of
    metabolite-centric signaling. *Nature Communications*.

2.  Xiao, Z., et al. (2022). Metabolite-mediated intercellular
    communication in tumor microenvironment. *Frontiers in Cell and
    Developmental Biology*.

## Next

- **[Production & Sensing
  Analysis](https://Zaoqu-Liu.github.io/scMetaLink/articles/03-production-sensing.md)**:
  Detailed parameter tuning
- **[Communication
  Analysis](https://Zaoqu-Liu.github.io/scMetaLink/articles/04-communication.md)**:
  Advanced communication analysis
