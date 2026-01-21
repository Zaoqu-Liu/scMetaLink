# Get Production Enzymes for Metabolites

Extract enzymes and secretion transporters involved in metabolite
production and release. This function filters for proteins that
contribute to making metabolites available for intercellular
communication.

## Usage

``` r
.get_production_enzymes(db)
```

## Arguments

- db:

  MetalinksDB database list

## Value

data.frame of production relationships

## Details

The production potential consists of two components:

- Synthesis enzymes: Proteins that catalyze metabolite biosynthesis

- Secretion transporters: Membrane transporters that release metabolites
  to the extracellular space (transport_direction = "out")

Uptake transporters (transport_direction = "in") are excluded from
production as they represent metabolite sensing/uptake rather than
production/release.
