# Calculate Row-wise Trimean

Calculate the trimean for each row of a matrix. Trimean = (Q1 + 2\*Q2 +
Q3) / 4, providing a robust central tendency estimate.

## Usage

``` r
.row_trimean(mat)
```

## Arguments

- mat:

  Matrix (genes x cells), can be sparse or dense

## Value

Numeric vector of trimean values for each row
