# Calculate Lactate Communication Scores

Combines production and sensing scores into cell-cell communication
scores.

## Usage

``` r
.calc_lactate_communication(
  production_scores,
  sensing_scores,
  method = "geometric",
  min_production = 0,
  min_sensing = 0
)
```

## Arguments

- production_scores:

  Named numeric vector of production scores

- sensing_scores:

  Named numeric vector of sensing scores

- method:

  Character. Communication method: "geometric", "product", "harmonic"

- min_production:

  Numeric. Minimum production threshold

- min_sensing:

  Numeric. Minimum sensing threshold

## Value

Matrix of communication scores (sender x receiver)
