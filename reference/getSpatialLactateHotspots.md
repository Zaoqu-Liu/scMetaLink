# Get Spatial Lactate Hotspots

Identifies spatial hotspots of lactate production and sensing.

## Usage

``` r
getSpatialLactateHotspots(object, type = "all", top_n = 20)
```

## Arguments

- object:

  scMetaLink object with spatial lactate signaling results

- type:

  Character. Type of hotspot: "production", "direct_sensing",
  "indirect_sensing", or "all". Default "all".

- top_n:

  Integer. Number of top spots to return. Default 20.

## Value

data.frame with spot IDs, coordinates, scores, and cell types
