# Predictive distributions for the between-study variance in a future meta-analysis on odds ratio or standardised mean difference

A table with the hyperparameters of the predictive distributions for the
between-study variance developed by Turner et al. (2015) and Rhodes et
al. (2015): log-normal distribution and t-distribution (with 5 degrees
of freedom) when the outcome data are analysed in the odds ratio or
standardised mean difference scale, respectively.

## Usage

``` r
table_tau2_prior(measure, area)
```

## Arguments

- measure:

  Character string indicating the effect measure with possible values
  `"OR"` for odds ratio and `"SMD"` for standardised mean difference.

- area:

  Character string indicating the medical area relating to the
  predictive distributions for standardised mean difference with
  possible values `"cancer"` for medical areas of cancer,
  `"respiratory"` for medical areas of respiratory diseases, and
  `"other"` for medical areas other than cancer or respiratory diseases.
  The argument is not relevant for odds ratio.

## Value

A cross-sectional table as a heatmap showing the hyperparameters (mean
and standard deviation) of the corresponding predictive distribution for
all combinations between the outcome types and treatment-comparison
types and according to the selected medical area (only relevant with
standardised mean difference) as defined by Turner et al. (2015) and
Rhodes et al. (2015). The tiles are coloured with different shades
according to the corresponding median value: the larger the median, the
darker the colour.

## Details

This table aids in selecting the hyperparameters for the function
[`heterogeneity_param_prior`](https://loukiaspin.github.io/rnmamod/reference/heterogeneity_param_prior.md)
when considering an informative prior distribution for the between-study
variance parameter based on the two publications mentioned above
(relevant for the function
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
to conduct random-effects network meta-analysis).

## References

Rhodes KM, Turner RM, Higgins JP. Predictive distributions were
developed for the extent of heterogeneity in meta-analyses of continuous
outcome data. *J Clin Epidemiol* 2015;**68**(1):52–60. doi:
10.1016/j.jclinepi.2014.08.012

Turner RM, Jackson D, Wei Y, Thompson SG, Higgins JP. Predictive
distributions for between-study heterogeneity and simple methods for
their application in Bayesian meta-analysis. *Stat Med*
2015;**34**(6):984–98. doi: 10.1002/sim.6381

## See also

[`heterogeneity_param_prior`](https://loukiaspin.github.io/rnmamod/reference/heterogeneity_param_prior.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)

## Author

Loukia M. Spineli
