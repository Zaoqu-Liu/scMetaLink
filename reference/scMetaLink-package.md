# scMetaLink: Single-Cell Metabolite-Mediated Cell Communication Analysis

A comprehensive framework for inferring metabolite-mediated cell-cell
communication from single-cell transcriptomic data. scMetaLink
integrates metabolite production potential via enzyme expression,
metabolite sensing capability via receptor and transporter expression,
and secretion potential to construct intercellular metabolic
communication networks.

## Details

The package provides the following main functions:

- [`createScMetaLink`](https://Zaoqu-Liu.github.io/scMetaLink/reference/createScMetaLink.md):
  Create analysis object from expression data

- [`inferProduction`](https://Zaoqu-Liu.github.io/scMetaLink/reference/inferProduction.md):
  Infer metabolite production potential

- [`inferSensing`](https://Zaoqu-Liu.github.io/scMetaLink/reference/inferSensing.md):
  Infer metabolite sensing capability

- [`computeCommunication`](https://Zaoqu-Liu.github.io/scMetaLink/reference/computeCommunication.md):
  Compute cell-cell communication scores

- [`filterSignificantInteractions`](https://Zaoqu-Liu.github.io/scMetaLink/reference/filterSignificantInteractions.md):
  Filter significant interactions

- [`aggregateByPathway`](https://Zaoqu-Liu.github.io/scMetaLink/reference/aggregateByPathway.md):
  Aggregate by metabolic pathways

- [`runScMetaLink`](https://Zaoqu-Liu.github.io/scMetaLink/reference/runScMetaLink.md):
  Run complete analysis pipeline

For spatial transcriptomics data:

- [`createScMetaLinkFromSpatial`](https://Zaoqu-Liu.github.io/scMetaLink/reference/createScMetaLinkFromSpatial.md):
  Create object from spatial data

- [`computeSpatialCommunication`](https://Zaoqu-Liu.github.io/scMetaLink/reference/computeSpatialCommunication.md):
  Compute spatially-weighted communication

## Database

scMetaLink utilizes MetalinksDB, containing:

- 41,894 metabolite-protein interactions

- 1,128 metabolites

- 4,374 proteins/genes

- 157,741 pathway associations

## See also

Useful links:

- <https://github.com/Zaoqu-Liu/scMetaLink>

- <https://Zaoqu-Liu.github.io/scMetaLink/>

- Report bugs at <https://github.com/Zaoqu-Liu/scMetaLink/issues>

## Author

Zaoqu Liu <liuzaoqu@163.com>
