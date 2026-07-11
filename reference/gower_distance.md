# Weighted Gower's dissimilarity measure (Trials' comparability for transitivity evaluation)

`gower_distance` calculate the weighted Gower's dissimilarity
coefficient for all pairs of trials included in a network of
interventions, considering several characteristics measured at trial
level (Spineli et al., 2025). It takes values from 0 to 1, with 0
implying complete similarity and 1 complete dissimilarity.

## Usage

``` r
gower_distance(input, weight)
```

## Arguments

- input:

  A data-frame in the long arm-based format. Two-arm trials occupy one
  row in the data-frame. Multi-arm trials occupy as many rows as the
  number of possible comparisons among the interventions. The first two
  columns refer to the trial name, and the pairwise comparison,
  respectively. The remaining columns refer to summary characteristics.
  See 'Details' for the specification of the columns.

- weight:

  A vector of non-negative numbers to define the weight contribution of
  each characteristic. The default is a vector of 1s for all
  characteristics.

## Value

`gower_distance` returns the following list of elements:

- Dissimilarity_table:

  A lower off-diagonal matrix of 'dist' class with the dissimilarities
  of all pairs of trials.

- Types_used:

  A data-frame with type mode (i.e., double or integer) of each
  characteristic.

- Total_missing:

  The percentage of missing cases in the comparison, calculated as the
  ratio of total missing cases to the product of the number of studies
  with the number of characteristics.

## Details

The correct type mode of columns in `input` must be ensured to use the
function `gower_distance`. The first two columns referring to the trial
name, and pairwise comparison, respectively, must be **character**. The
remaining columns referring to the characteristics must be **double** or
**integer** depending on whether the corresponding characteristic refers
to a quantitative or qualitative variable. The type mode of each column
is assessed by `gower_distance` using the base function `typeof`. Note
that `gower_distance` invites unordered and ordered variables; for the
latter, add the argument `ordered = TRUE` in the base function
**factor()**.

`gower_distance` is integrated in the function
[`comp_clustering`](https://loukiaspin.github.io/rnmamod/reference/comp_clustering.md).

## References

Gower J. General Coefficient of Similarity and Some of Its Properties.
*Biometrics* 1971;**27**(4):857–71. doi: 10.2307/2528823

Spineli LM, Papadimitropoulou K, Kalyvas C. Exploring the Transitivity
Assumption in Network Meta-Analysis: A Novel Approach and Its
Implications. *Stat Med* 2025;**44**(7):e70068. doi: 10.1002/sim.70068.

## See also

[`comp_clustering`](https://loukiaspin.github.io/rnmamod/reference/comp_clustering.md)

## Author

Loukia M. Spineli

## Examples

``` r
# Fictional dataset
data_set <- data.frame(Trial_name = as.character(1:7),
                      comp = c("1vs2", "1vs2", "1vs2", "1vs3", "1vs3", "2vs3", "2vs3"),
                      sample = c(140, 145, 150, 40, 45, 75, 80),
                      age = c(18, 18, 18, 48, 48, 35, 35),
                      blinding = as.integer(c("yes", "yes", "yes", "no", "no", "no", "no")))
#> Warning: NAs introduced by coercion

# Calculate the weighted Gower dissimilarity of all study pairs in the network
gower_distance(input = data_set)
#> $Dissimilarity_table
#>            [,1]       [,2]      [,3]       [,4]      [,5]       [,6] [,7]
#> [1,] 0.00000000         NA        NA         NA        NA         NA   NA
#> [2,] 0.02272727 0.00000000        NA         NA        NA         NA   NA
#> [3,] 0.04545455 0.02272727 0.0000000         NA        NA         NA   NA
#> [4,] 0.95454545 0.97727273 1.0000000 0.00000000        NA         NA   NA
#> [5,] 0.93181818 0.95454545 0.9772727 0.02272727 0.0000000         NA   NA
#> [6,] 0.57878788 0.60151515 0.6242424 0.37575758 0.3530303 0.00000000   NA
#> [7,] 0.55606061 0.57878788 0.6015152 0.39848485 0.3757576 0.02272727    0
#> 
#> $Types_used
#>   characteristic    type
#> 1         sample  double
#> 2            age  double
#> 3       blinding integer
#> 
#> $Total_missing
#> [1] "33.33%"
#> 
#> $Variable_dissimilarities
#>           [,1]      [,2] [,3]
#> 1-2 0.04545455 0.0000000   NA
#> 1-3 0.09090909 0.0000000   NA
#> 1-4 0.90909091 1.0000000   NA
#> 1-5 0.86363636 1.0000000   NA
#> 1-6 0.59090909 0.5666667   NA
#> 1-7 0.54545455 0.5666667   NA
#> 2-3 0.04545455 0.0000000   NA
#> 2-4 0.95454545 1.0000000   NA
#> 2-5 0.90909091 1.0000000   NA
#> 2-6 0.63636364 0.5666667   NA
#> 2-7 0.59090909 0.5666667   NA
#> 3-4 1.00000000 1.0000000   NA
#> 3-5 0.95454545 1.0000000   NA
#> 3-6 0.68181818 0.5666667   NA
#> 3-7 0.63636364 0.5666667   NA
#> 4-5 0.04545455 0.0000000   NA
#> 4-6 0.31818182 0.4333333   NA
#> 4-7 0.36363636 0.4333333   NA
#> 5-6 0.27272727 0.4333333   NA
#> 5-7 0.31818182 0.4333333   NA
#> 6-7 0.04545455 0.0000000   NA
#> 
```
