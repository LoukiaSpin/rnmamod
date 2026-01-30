# End-user-ready results for unrelated trial effects model

Performs the unrelated trial effects model (also known as fixed effects
model) and illustrates the results of each trial and corresponding
pairwise comparison.

## Usage

``` r
unrelated_effects_plot(
  data,
  measure,
  char,
  drug_names,
  trial_names,
  mean_misspar,
  var_misspar,
  rho,
  save_xls
)
```

## Arguments

- data:

  A data-frame of a one-trial-per-row format containing arm-level data
  of each trial. See 'Format' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

- measure:

  Character string indicating the effect measure with values `"OR"`,
  `"MD"`, `"SMD"`, or `"ROM"` for the odds ratio, mean difference,
  standardised mean difference and ratio of means, respectively.

- char:

  A data-frame of three columns and number of rows equal to the number
  of trials in `data`. Each column refers to a trial-characteristic with
  **nominal** elements.

- drug_names:

  A vector of labels with the name of the interventions in the order
  they appear in the argument `data`. If `drug_names` is not defined,
  the order of the interventions as they appear in `data` is used,
  instead.

- trial_names:

  A vector of labels with the name of the trials in the order they
  appear in the argument `data`. If `trial_names` is not defined, the
  order of the trials as they appear in `data` is used, instead.

- mean_misspar:

  A numeric value for the mean of the normal distribution of the
  informative missingness parameter (see 'Details'). The default
  argument is 0 and corresponds to the missing-at-random assumption. The
  same value is considered across all trials of the dataset.

- var_misspar:

  A positive non-zero number for the variance of the normal distribution
  of the informative missingness parameter. When the `measure` is
  `"OR"`, `"MD"`, or `"SMD"` the default argument is 1. When the
  `measure` is `"ROM"` the default argument is 0.04. The same value is
  considered across all trials of the dataset.

- rho:

  A numeric value in the interval \[-1, 1\] that indicates the
  correlation coefficient between two informative missingness parameters
  in a trial. The same value is considered across all trials of the
  dataset. The default argument is 0 and corresponds to uncorrelated
  missingness parameters.

- save_xls:

  Logical to indicate whether to export the tabulated results to an
  'xlsx' file (via the
  [`write_xlsx`](https://docs.ropensci.org/writexl//reference/write_xlsx.html)
  function of the R-package
  [writexl](https://CRAN.R-project.org/package=writexl)) to the working
  directory of the user. The default is `FALSE` (do not export).

## Value

A panel of interval plots for each observed comparison in the network,
when there are up to 15 trials in the `data`. Otherwise,
`unrelated_effects_plot` exports a data-frame to an 'xlsx' file at the
working directory of the user. This data-frame includes the `data` in
the long format, the within-trial effect measure and 95% confidence
interval of the corresponding comparisons, the interventions compared,
and the three characteristics (as defined in `char`). For datasets with
more than 15 trials, the plot becomes cluttered and it is difficult to
identify the trial-names. Hence, exporting the results in an Excel file
is a viable alternative.

## Details

The unrelated trial effects model may be an alternative to network
meta-analysis, when the latter is not deemed appropriate (e.g., there is
considerable statistical heterogeneity, or substantial intransitivity).
In the presence of missing participant outcome data, the effect size and
standard error are adjusted by applying the pattern-mixture model with
Taylor series in trial-arms with reported missing participants (Mavridis
et al., 2015; White et al., 2008). The `unrelated_effects_plot` function
calls the
[`taylor_imor`](https://loukiaspin.github.io/rnmamod/reference/taylor_imor.md)
and
[`taylor_continuous`](https://loukiaspin.github.io/rnmamod/reference/taylor_continuous.md)
functions (for a binary and continuous outcome, respectively) to employ
pattern-mixture model with Taylor series. The `unrelated_effects_plot`
function considers the informative missingness odds ratio in the
logarithmic scale for binary outcome data (White et al., 2008), the
informative missingness difference of means when `measure` is `"MD"` or
`"SMD"`, and the informative missingness ratio of means in the
logarithmic scale when `measure` is `"ROM"` (Mavridis et al., 2015).

The number of interval plots equals the number of observed comparisons
in the network. In each interval plot, the y-axis refers to all trials
of the network and x-axis refers to the selected effect measure. The
odds ratio and ratio of means are calculated in the logarithmic scale
but they are reported in their original scale after exponentiation.

`unrelated_effects_plot` depicts all three characteristics for each
trial using different colours, line-types and point-shapes for the
corresponding 95% confidence interval and point estimate. Ideally, each
characteristic should have no more than three categories; otherwise, the
plot becomes cluttered. For now, the `unrelated_effects_plot` function
uses the default colour palette, line-types and point-shapes.

## References

Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
uncertainty due to missing continuous outcome data in pairwise and
network meta-analysis. *Stat Med* 2015;**34**(5):721–41. doi:
10.1002/sim.6365

White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing
data in meta-analysis–part 1: two-stage methods. *Stat Med*
2008;**27**(5):711–27. doi: 10.1002/sim.3008

## See also

[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`taylor_continuous`](https://loukiaspin.github.io/rnmamod/reference/taylor_continuous.md),
[`taylor_imor`](https://loukiaspin.github.io/rnmamod/reference/taylor_imor.md),
[`write_xlsx`](https://docs.ropensci.org/writexl//reference/write_xlsx.html)

## Author

Loukia M. Spineli
