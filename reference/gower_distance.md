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
