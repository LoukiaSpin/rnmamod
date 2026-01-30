# Pattern-mixture model with Taylor series for a binary outcome

Applies the pattern-mixture model under a specific assumption about the
informative missingness odds ratio in trial-arms with **binary** missing
participant outcome data and uses the Taylor series to obtain the odds
ratio (in the logarithmic scale) and standard error for each trial
(White et al., 2008).

## Usage

``` r
taylor_imor(data, mean_value, var_value, rho)
```

## Format

The columns of the data-frame in the argument `data` refer to the
following ordered elements for a binary outcome:

|        |                                                                            |
|--------|----------------------------------------------------------------------------|
| **id** | A unique identifier for each trial.                                        |
| **r1** | The observed number of events in the first arm of the comparison.          |
| **r2** | The observed number of events in the second arm of the comparison.         |
| **m1** | The number of missing participants in the first arm of the comparison.     |
| **m2** | The number of missing participants in the second arm of the comparison.    |
| **n1** | The number of participants randomised in the first arm of the comparison.  |
| **n2** | The number of participants randomised in the second arm of the comparison. |
| **t1** | An identifier for the intervention in the first arm of the comparison.     |
| **t2** | An identifier for the intervention in the second arm of the comparison.    |

## Arguments

- data:

  A data-frame in the long arm-based format. Two-arm trials occupy one
  row in the data-frame. Multi-arm trials occupy as many rows as the
  number of possible comparisons among the interventions. See 'Format'
  for the specification of the columns.

- mean_value:

  A numeric value for the mean of the normal distribution of the
  informative missingness odds ratio in the logarithmic scale. The same
  value is considered for all trial-arms of the dataset. The default
  argument is 0 and corresponds to the missing-at-random assumption.

- var_value:

  A positive non-zero number for the variance of the normal distribution
  of the informative missingness odds ratio in the logarithmic scale.
  The default argument is 1.

- rho:

  A numeric value in the interval \[-1, 1\] that indicates the
  correlation coefficient between two missingness parameters in a trial.
  The same value is considered across all trials of the dataset. The
  default argument is 0 and corresponds to uncorrelated missingness
  parameters.

## Value

A data-frame that additionally includes the following elements:

- EM:

  The odds ratio in the logarithmic scale (log OR) adjusted for missing
  participants and obtained using the Taylor series.

- se.EM:

  The standard error of the log OR adjusted for missing participants and
  obtained using the Taylor series.

## Details

The `taylor_imor` function is integrated in the
[`unrelated_effects_plot`](https://loukiaspin.github.io/rnmamod/reference/unrelated_effects_plot.md)
function.

## References

White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing
data in meta-analysis–part 1: two-stage methods. *Stat Med*
2008;**27**(5):711–27. doi: 10.1002/sim.3008

## See also

[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`unrelated_effects_plot`](https://loukiaspin.github.io/rnmamod/reference/unrelated_effects_plot.md)

## Author

Loukia M. Spineli
