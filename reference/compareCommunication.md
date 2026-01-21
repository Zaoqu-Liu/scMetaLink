# Compare Two Conditions

Compare metabolite-mediated communication between two conditions

## Usage

``` r
compareCommunication(
  object1,
  object2,
  condition_names = c("Condition1", "Condition2"),
  method = "log2fc"
)
```

## Arguments

- object1:

  scMetaLink object for condition 1

- object2:

  scMetaLink object for condition 2

- condition_names:

  Character vector of length 2 for condition names

- method:

  Character. Comparison method: "difference", "ratio", or "log2fc"

## Value

data.frame with differential communication results

## Examples

``` r
# \donttest{
# Compare tumor vs normal
diff <- compareCommunication(tumor_obj, normal_obj,
                             condition_names = c("Tumor", "Normal"))
#> Error: object 'tumor_obj' not found
# }
```
