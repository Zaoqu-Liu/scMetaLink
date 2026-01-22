# Package index

## Main Workflow

High-level functions for running scMetaLink analysis

- [`createScMetaLink()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/createScMetaLink.md)
  : Create scMetaLink Object
- [`createScMetaLinkFromSeurat()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/createScMetaLinkFromSeurat.md)
  : Create scMetaLink from Seurat Object
- [`runScMetaLink()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/runScMetaLink.md)
  : Run Complete scMetaLink Analysis
- [`runScMetaLinkSeurat()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/runScMetaLinkSeurat.md)
  : Run scMetaLink from Seurat Object

## Production Analysis

Infer metabolite production potential

- [`inferProduction()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/inferProduction.md)
  : Infer Metabolite Production Potential
- [`getTopProducers()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/getTopProducers.md)
  : Get Top Producing Cell Types for a Metabolite

## Sensing Analysis

Infer metabolite sensing capability

- [`inferSensing()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/inferSensing.md)
  : Infer Metabolite Sensing Capability
- [`getTopSensors()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/getTopSensors.md)
  : Get Top Sensing Cell Types for a Metabolite
- [`getMetaboliteReceptors()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/getMetaboliteReceptors.md)
  : Get Receptors for a Metabolite

## Communication Analysis

Compute cell-cell communication scores

- [`computeCommunication()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/computeCommunication.md)
  : Compute Metabolite-Mediated Cell Communication
- [`filterSignificantInteractions()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/filterSignificantInteractions.md)
  : Filter Significant Interactions
- [`summarizeCommunicationPairs()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/summarizeCommunicationPairs.md)
  : Summarize Communication by Cell Type Pairs
- [`getCommunicationMatrix()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/getCommunicationMatrix.md)
  : Get Communication Summary Matrix
- [`compareCommunication()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/compareCommunication.md)
  : Compare Two Conditions
- [`identifyCellTypeSpecificMetabolites()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/identifyCellTypeSpecificMetabolites.md)
  : Identify Cell Type Specific Metabolites

## Spatial Analysis

Spatial transcriptomics support

- [`createScMetaLinkFromSpatial()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/createScMetaLinkFromSpatial.md)
  : Create scMetaLink Object from Spatial Data
- [`createScMetaLinkFromSeuratSpatial()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/createScMetaLinkFromSeuratSpatial.md)
  : Create scMetaLink from Seurat Spatial Object
- [`computeSpatialCommunication()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/computeSpatialCommunication.md)
  : Compute Spatial Communication
- [`getSpatialDistanceStats()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/getSpatialDistanceStats.md)
  : Get Spatial Distance Statistics
- [`identifyCommunicationHotspots()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/identifyCommunicationHotspots.md)
  : Identify Communication Hotspots

## Spatial Visualization

Spatial plotting functions

- [`plotSpatialCellTypes()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotSpatialCellTypes.md)
  : Plot Spatial Cell Types
- [`plotSpatialFeature()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotSpatialFeature.md)
  : Plot Spatial Communication Heatmap
- [`plotSpatialCommunicationNetwork()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotSpatialCommunicationNetwork.md)
  : Plot Spatial Communication Network
- [`plotSpatialDistanceDistribution()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotSpatialDistanceDistribution.md)
  : Plot Spatial Distance Distribution
- [`plotSpatialHotspots()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotSpatialHotspots.md)
  : Plot Spatial Hotspots
- [`plotSpatialComparison()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotSpatialComparison.md)
  : Plot Spatial Communication Comparison

## Pathway Analysis

Pathway-level aggregation and enrichment

- [`aggregateByPathway()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/aggregateByPathway.md)
  : Aggregate Communication by Pathway
- [`enrichPathways()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/enrichPathways.md)
  : Pathway Enrichment Analysis
- [`getPathwayCommunicationMatrix()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/getPathwayCommunicationMatrix.md)
  : Get Pathway Communication Matrix
- [`getPathwayMetabolites()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/getPathwayMetabolites.md)
  : Get Metabolites in Pathway
- [`getPathwayMetaboliteNetwork()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/getPathwayMetaboliteNetwork.md)
  : Get Pathway-Metabolite Network
- [`summarizePathwayActivity()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/summarizePathwayActivity.md)
  : Summarize Pathway Activity
- [`listTopPathways()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/listTopPathways.md)
  : List Top Pathways

## Visualization

General plotting functions

- [`plotCommunicationHeatmap()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotCommunicationHeatmap.md)
  : Plot Communication Heatmap
- [`plotCommunicationCircle()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotCommunicationCircle.md)
  : Plot Communication Circle
- [`plotCommunicationNetwork()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotCommunicationNetwork.md)
  : Plot Communication Network
- [`plotMetaboliteProfile()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotMetaboliteProfile.md)
  : Plot Metabolite Profile
- [`plotTopInteractions()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotTopInteractions.md)
  : Plot Top Interactions
- [`plotPathwayCommunication()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotPathwayCommunication.md)
  : Plot Pathway Communication
- [`plotDifferentialCommunication()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotDifferentialCommunication.md)
  : Plot Differential Communication
- [`plotEnrichedPathways()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/plotEnrichedPathways.md)
  : Plot Enriched Pathways

## Database Utilities

Functions to query MetalinksDB

- [`listMetabolites()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/listMetabolites.md)
  : List Available Metabolites
- [`listGenes()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/listGenes.md)
  : List Available Genes
- [`searchMetabolite()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/searchMetabolite.md)
  : Search Metabolite in Database
- [`searchGene()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/searchGene.md)
  : Search Gene in Database
- [`getDatabaseInfo()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/getDatabaseInfo.md)
  : Get Database Information

## Export & Utilities

Export results and helper functions

- [`exportResults()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/exportResults.md)
  : Export Results to CSV
- [`saveScMetaLink()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/saveScMetaLink.md)
  : Save scMetaLink Object
- [`loadScMetaLink()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/loadScMetaLink.md)
  : Load scMetaLink Object
- [`getSummaryStats()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/getSummaryStats.md)
  : Get Summary Statistics
- [`citationScMetaLink()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/citationScMetaLink.md)
  : Show Package Citation

## Data Accessors

S4 methods and class definition

- [`show(`*`<scMetaLink>`*`)`](https://Zaoqu-Liu.github.io/scMetaLink/reference/scMetaLink-class.md)
  : scMetaLink Class Definition
- [`getProductionScores()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/accessors.md)
  [`getSensingScores()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/accessors.md)
  [`getCommunicationScores()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/accessors.md)
  [`getSignificantInteractions()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/accessors.md)
  [`getPathwayAggregated()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/accessors.md)
  [`getParameters()`](https://Zaoqu-Liu.github.io/scMetaLink/reference/accessors.md)
  : Accessor Functions for scMetaLink Objects

## Example Data

Built-in datasets for tutorials

- [`crc_expr`](https://Zaoqu-Liu.github.io/scMetaLink/reference/crc_expr.md)
  : CRC Example Expression Data
- [`crc_meta`](https://Zaoqu-Liu.github.io/scMetaLink/reference/crc_meta.md)
  : CRC Example Cell Metadata
- [`metalinksdb`](https://Zaoqu-Liu.github.io/scMetaLink/reference/metalinksdb.md)
  : MetalinksDB Database
- [`st_expr`](https://Zaoqu-Liu.github.io/scMetaLink/reference/st_expr.md)
  : Colon Spatial Transcriptomics Expression Data
- [`st_meta`](https://Zaoqu-Liu.github.io/scMetaLink/reference/st_meta.md)
  : Colon Spatial Transcriptomics Metadata
- [`st_scalefactors`](https://Zaoqu-Liu.github.io/scMetaLink/reference/st_scalefactors.md)
  : Colon Spatial Transcriptomics Scale Factors

## Package

Package documentation

- [`scMetaLink`](https://Zaoqu-Liu.github.io/scMetaLink/reference/scMetaLink-package.md)
  [`scMetaLink-package`](https://Zaoqu-Liu.github.io/scMetaLink/reference/scMetaLink-package.md)
  : scMetaLink: Single-Cell Metabolite-Mediated Cell Communication
  Analysis
