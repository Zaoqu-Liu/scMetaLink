# Changelog

## scMetaLink 0.99.0

### New Features

- Initial Bioconductor submission
- [`createScMetaLink()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/createScMetaLink.md) -
  Create scMetaLink object from expression data
- [`createScMetaLinkFromSeurat()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/createScMetaLinkFromSeurat.md) -
  Create from Seurat object
- `createScMetaLinkFromSCE()` - Create from SingleCellExperiment
- [`createScMetaLinkFromSpatial()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/createScMetaLinkFromSpatial.md) -
  Create from spatial transcriptomics data
- [`inferProduction()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/inferProduction.md) -
  Infer metabolite production potential
- [`inferSensing()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/inferSensing.md) -
  Infer metabolite sensing capability
- [`computeCommunication()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/computeCommunication.md) -
  Compute cell-cell communication scores
- [`computeSpatialCommunication()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/computeSpatialCommunication.md) -
  Spatial communication analysis
- [`filterSignificantInteractions()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/filterSignificantInteractions.md) -
  Filter by statistical significance
- [`aggregateByPathway()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/aggregateByPathway.md) -
  Pathway-level aggregation

### Lactate Signaling Module

- [`inferLactateSignaling()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/inferLactateSignaling.md) -
  Specialized lactate communication analysis
- [`checkLactateGenes()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/checkLactateGenes.md) -
  Check lactate-related gene availability
- [`getLactateGenes()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/getLactateGenes.md) -
  Get curated lactate gene sets

### Visualization

- [`plotCommunicationHeatmap()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotCommunicationHeatmap.md) -
  Communication heatmap
- [`plotCommunicationCircle()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotCommunicationCircle.md) -
  Chord diagram
- [`plotCommunicationNetwork()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotCommunicationNetwork.md) -
  Network visualization
- `plotMetaboliteFlow()` - Metabolite flow diagram
- [`plotSpatialCellTypes()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotSpatialCellTypes.md) -
  Spatial cell type distribution
- [`plotSpatialFeature()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotSpatialFeature.md) -
  Spatial feature maps
- [`plotSpatialCommunicationNetwork()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotSpatialCommunicationNetwork.md) -
  Spatial network overlay
