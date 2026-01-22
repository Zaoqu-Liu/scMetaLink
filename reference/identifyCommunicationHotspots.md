# Identify Communication Hotspots

Find spatial regions with high metabolite communication activity.
\*\*Note\*\*: This function requires running
computeSpatialCommunication() with analysis_level='spot' first.

## Usage

``` r
identifyCommunicationHotspots(
  object,
  metabolite = NULL,
  type = "sender",
  n_hotspots = 5,
  method = "density"
)
```

## Arguments

- object:

  A scMetaLink object with spatial communication results

- metabolite:

  Character. Metabolite ID or name to analyze (NULL for aggregate)

- type:

  Character. "sender" or "receiver" to identify production or sensing
  hotspots

- n_hotspots:

  Integer. Number of hotspot regions to identify

- method:

  Character. Method for hotspot detection: "density" or "clustering"

## Value

A data.frame with hotspot information
