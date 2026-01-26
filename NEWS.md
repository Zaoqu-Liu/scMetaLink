# scMetaLink 0.99.1

## Bug Fixes and Improvements

* Added `createScMetaLinkFromSCE()` function for SingleCellExperiment support
* Cleaned repository structure for Bioconductor submission
* Reduced package size by optimizing vignettes

# scMetaLink 0.99.0

## New Features

* Initial Bioconductor submission
* `createScMetaLink()` - Create scMetaLink object from expression data
* `createScMetaLinkFromSeurat()` - Create from Seurat object
* `createScMetaLinkFromSCE()` - Create from SingleCellExperiment
* `createScMetaLinkFromSpatial()` - Create from spatial transcriptomics data
* `inferProduction()` - Infer metabolite production potential
* `inferSensing()` - Infer metabolite sensing capability
* `computeCommunication()` - Compute cell-cell communication scores
* `computeSpatialCommunication()` - Spatial communication analysis
* `filterSignificantInteractions()` - Filter by statistical significance
* `aggregateByPathway()` - Pathway-level aggregation

## Lactate Signaling Module

* `inferLactateSignaling()` - Specialized lactate communication analysis
* `checkLactateGenes()` - Check lactate-related gene availability
* `getLactateGenes()` - Get curated lactate gene sets

## Visualization

* `plotCommunicationHeatmap()` - Communication heatmap
* `plotCommunicationCircle()` - Chord diagram
* `plotCommunicationNetwork()` - Network visualization
* `plotMetaboliteFlow()` - Metabolite flow diagram
* `plotSpatialCellTypes()` - Spatial cell type distribution
* `plotSpatialFeature()` - Spatial feature maps
* `plotSpatialCommunicationNetwork()` - Spatial network overlay
