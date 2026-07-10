# Pattern-mixture model with Taylor series for continuous outcome

Applies the pattern-mixture model under a specific assumption about the
informative missingness parameter in trial-arms with **continuous**
missing participant outcome data and uses the Taylor series to obtain
the effect size and standard error for each trial (Mavridis et al.,
2015).

## Usage

``` r
taylor_continuous(data, measure, mean_value, var_value, rho)
```

## Format

The columns of the data-frame in the argument `data` refer to the
following ordered elements for a continuous outcome:

|  |  |
|----|----|
| **id** | A unique identifier for each trial. |
| **y1** | The observed mean outcome in the first arm of the comparison. |
| **y2** | The observed mean outcome in the second arm of the comparison. |
| **sd1** | The observed standard deviation of the outcome in the first arm of the comparison. |
| **sd2** | The observed standard deviation of the outcome in the second arm of the comparison. |
| **m1** | The number of missing participants in the first arm of the comparison. |
| **m2** | The number of missing participants in the second arm of the comparison. |
| **n1** | The number randomised in the first arm of the comparison. |
| **n2** | The number randomised in the second arm of the comparison. |
| **t1** | An identifier for the intervention in the first arm of the comparison. |
| **t2** | An identifier for the intervention in the second arm of the comparison. |

## Arguments

- data:

  A data-frame in the long arm-based format. Two-arm trials occupy one
  row in the data-frame. Multi-arm trials occupy as many rows as the
  number of possible comparisons among the interventions. See 'Format'
  for the specification of the columns.

- measure:

  Character string indicating the effect measure with values `"MD"`,
  `"SMD"`, or `"ROM"` for the mean difference, standardised mean
  difference, and ratio of means, respectively.

- mean_value:

  A numeric value for the mean of the normal distribution of the
  informative missingness parameter. The same value is considered for
  all trial-arms of the dataset. The default argument is 0 and
  corresponds to the missing-at-random assumption. For the informative
  missingness ratio of means, the mean value is defined in the
  logarithmic scale.

- var_value:

  A positive non-zero number for the variance of the normal distribution
  of the informative missingness parameter. When the `measure` is
  `"MD"`, or `"SMD"` the default argument is 1; when the `measure` is
  `"ROM"` the default argument is 0.04. The same value is considered for
  all trial-arms of the dataset.

- rho:

  A numeric value in the interval \[-1, 1\] that indicates the
  correlation coefficient between two informative missingness parameters
  in a trial. The same value is considered across all trials of the
  dataset. The default argument is 0 and corresponds to uncorrelated
  missingness parameters.

## Value

A data-frame that additionally includes the following elements:

- EM:

  The effect size adjusted for the missing participants and obtained
  using the Taylor series.

- se.EM:

  The standard error of the effect size adjusted for the missing
  participants and obtained using the Taylor series.

## Details

The `taylor_continuous` function is integrated in the
[`unrelated_effects_plot`](https://loukiaspin.github.io/rnmamod/reference/unrelated_effects_plot.md)
function.

## References

Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
uncertainty due to missing continuous outcome data in pairwise and
network meta-analysis. *Stat Med* 2015;**34**(5):721–41. doi:
10.1002/sim.6365

## See also

[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`unrelated_effects_plot`](https://loukiaspin.github.io/rnmamod/reference/unrelated_effects_plot.md)

## Author

Loukia M. Spineli
