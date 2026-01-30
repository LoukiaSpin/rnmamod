# Function for the hyper-parameters of the prior distribution of the inconsistency variance (network meta-analysis with random inconsistency effects)

Calculates the mean and standard deviation of the log-normal
distribution and location-scale t-distribution of the inconsistency
variance in the log-odds ratio and standardised mean difference scales,
respectively, based on corresponding empirical distributions for the
between-study variance proposed by Turner et al. (2015) and Rhodes et
al. (2015). It also return the median value of the inconsistency
standard deviation.

## Usage

``` r
inconsistency_variance_prior(mean_tau2, sd_tau2, mean_scale, measure)
```

## Arguments

- mean_tau2:

  Mean value from the empirical prior distribution for the between-study
  variance.

- sd_tau2:

  Standard deviation value from the empirical prior distribution for the
  between-study variance.

- mean_scale:

  Positive (non-zero value) as a scaling factor of `mean_tau2`. See Law
  et al. (2016).

- measure:

  Character string indicating the effect measure. For a binary outcome,
  use only `"OR"` for the odds ratio. For a continuous outcome, use only
  `"SMD"` for standardised mean difference.

## Value

A list of three elements: the mean and standard deviation for the prior
distribution for the inconsistency variance, and the median
inconsistency standard deviation according to the selected empirical
prior distribution for the between-study variance.

## Details

Law et al. (2016) suggested using the proposed empirical prior
distributions for between-study variance to construct a prior
distribution for the inconsistency variance. The authors provided the
formulas for the hyper-parameters of the inconsistency variance for a
binary outcome measured in the log odds ratio scale. We extended the
idea for a continuous outcome measured in the standardised mean
difference scale. Currently, the empirical prior distributions for the
between-study variance have been proposed for these effect measures only
(Turner et al. (2015), Rhodes et al. (2015)).

## References

Law M, Jackson D, Turner R, Rhodes K, Viechtbauer W. Two new methods to
fit models for network meta-analysis with random inconsistency effects.
*BMC Med Res Methodol* 2016;**16**:87. doi: 10.1186/s12874-016-0184-5

Rhodes KM, Turner RM, Higgins JP. Predictive distributions were
developed for the extent of heterogeneity in meta-analyses of continuous
outcome data. *J Clin Epidemiol* 2015;**68**(1):52–60. doi:
10.1016/j.jclinepi.2014.08.012

Turner RM, Jackson D, Wei Y, Thompson SG, Higgins JP. Predictive
distributions for between-study heterogeneity and simple methods for
their application in Bayesian meta-analysis. *Stat Med*
2015;**34**(6):984–98. doi: 10.1002/sim.6381

## Author

Loukia M. Spineli
