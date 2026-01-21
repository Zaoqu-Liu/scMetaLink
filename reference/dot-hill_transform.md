# Hill Function Transformation

Apply Hill function to model receptor-ligand binding saturation. The
Hill function captures the sigmoidal relationship between expression and
biological response, reflecting receptor saturation at high
concentrations.

## Usage

``` r
.hill_transform(x, n = 1, Kh = 0.5)
```

## Arguments

- x:

  Numeric matrix or vector of expression values (should be 0-1 scaled)

- n:

  Hill coefficient (cooperativity parameter). n=1 gives standard
  Michaelis-Menten kinetics; n\>1 indicates positive cooperativity.

- Kh:

  Half-maximal response threshold. The value of x at which response =
  0.5.

## Value

Transformed values in the same structure as input
