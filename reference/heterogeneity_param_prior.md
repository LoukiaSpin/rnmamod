# Determine the prior distribution for the heterogeneity parameter

Generates the prior distribution (weakly informative or
empirically-based) for the heterogeneity parameter.
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
inherits `heterogeneity_param_prior` via the argument `heter_prior`.

## Usage

``` r
heterogeneity_param_prior(measure, model, heter_prior)
```

## Arguments

- measure:

  Character string indicating the effect measure. For a binary outcome,
  the following can be considered: `"OR"`, `"RR"` or `"RD"` for the odds
  ratio, relative risk, and risk difference, respectively. For a
  continuous outcome, the following can be considered: `"MD"`, `"SMD"`,
  or `"ROM"` for mean difference, standardised mean difference and ratio
  of means, respectively.

- model:

  Character string indicating the analysis model with values `"RE"`, or
  `"FE"` for the random-effects and fixed-effect model, respectively.
  The default argument is `"RE"`.

- heter_prior:

  A list of three elements with the following order: 1) a character
  string indicating the distribution with (currently available) values
  `"halfnormal"`, `"uniform"`, `"lognormal"`, or `"logt"`; 2) two
  numeric values that refer to the parameters of the selected
  distribution. For `"lognormal"`, and `"logt"` these numbers refer to
  the mean and precision, respectively. For `"halfnormal"`, these
  numbers refer to zero and the scale parameter (equal to 4 or 1 being
  the corresponding precision of the scale parameter 4 or 1). For
  `"uniform"`, these numbers refer to the minimum and maximum value of
  the distribution.

## Value

A value to be passed to
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

## Details

The names of the (current) prior distributions follow the JAGS syntax.
The mean and precision of `"lognormal"` and `"logt"` should align with
the values proposed by Turner et al. (2015) and Rhodes et al. (2015) for
the corresponding empirically-based prior distributions when `measure`
is `"OR"` or `"SMD"`, respectively. The users may refer to Dias et al.
(2013) to determine the minimum and maximum value of the uniform
distribution, and to Friede et al. (2017) to determine the mean and
precision of the half-normal distribution. When `model` is `"FE"`,
`heterogeneity_param_prior` is ignored in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

## References

Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
making 2: a generalized linear modeling framework for pairwise and
network meta-analysis of randomized controlled trials. *Med Decis
Making* 2013;**33**(5):607–17. doi: 10.1177/0272989X12458724

Friede T, Roever C, Wandel S, Neuenschwander B. Meta-analysis of two
studies in the presence of heterogeneity with applications in rare
diseases. *Biom J* 2017;**59**(4):658–71. doi: 10.1002/bimj.201500236

Rhodes KM, Turner RM, Higgins JP. Predictive distributions were
developed for the extent of heterogeneity in meta-analyses of continuous
outcome data. *J Clin Epidemiol* 2015;**68**(1):52–60. doi:
10.1016/j.jclinepi.2014.08.012

Turner RM, Jackson D, Wei Y, Thompson SG, Higgins JP. Predictive
distributions for between-study heterogeneity and simple methods for
their application in Bayesian meta-analysis. *Stat Med*
2015;**34**(6):984–98. doi: 10.1002/sim.6381

## See also

[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)

## Author

Loukia M. Spineli
