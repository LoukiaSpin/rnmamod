# Define the mean value of the normal distribution of the missingness parameter

Generates the mean value of the normal distribution of the missingness
parameter in the proper format depending on the assumed structure of the
missingness parameter.
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
inherits `missingness_param_prior` through the argument `mean_misspar`
(see 'Argument' in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).

## Usage

``` r
missingness_param_prior(assumption, mean_misspar)
```

## Arguments

- assumption:

  Character string indicating the structure of the informative
  missingness parameter. Set `assumption` equal to one of the following:
  `"HIE-COMMON"`, `"HIE-TRIAL"`, `"HIE-ARM"`, `"IDE-COMMON"`,
  `"IDE-TRIAL"`, `"IDE-ARM"`, `"IND-CORR"`, or `"IND-UNCORR"`. The
  default argument is `"IDE-ARM"`. The abbreviations `"IDE"`, `"HIE"`,
  and `"IND"` stand for identical, hierarchical and independent,
  respectively. `"CORR"` and `"UNCORR"` stand for correlated and
  uncorrelated, respectively.

- mean_misspar:

  A numeric value or a vector of two numeric values for the mean of the
  normal distribution of the informative missingness parameter (see
  'Details'). The default argument is 0 and corresponds to the
  missing-at-random assumption for `assumption = "IDE-ARM"`.

## Value

A scalar or numeric vector to be passed to
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

## Details

[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
considers the informative missingness odds ratio in the logarithmic
scale for binary outcome data (Spineli, 2019a; Turner et al., 2015;
White et al., 2008), the informative missingness difference of means
when `measure` is `"MD"` or `"SMD"`, and the informative missingness
ratio of means in the logarithmic scale when `measure` is `"ROM"`
(Spineli et al., 2021; Mavridis et al., 2015).

When `assumption` is trial-specific (i.e., `"IDE-TRIAL"` or
`"HIE-TRIAL"`), or independent (i.e., `"IND-CORR"` or `"IND-UNCORR"`),
only one numeric value can be assigned to `mean_misspar` because the
same missingness scenario is applied to all trials and trial-arms of the
dataset, respectively. When `assumption` is `"IDE-ARM"` or `"HIE-ARM"`,
a maximum of two *different or identical* numeric values can be assigned
as a vector to `mean_misspars`: the first value refers to the
experimental arm, and the second value refers to the control arm of a
trial. In the case of a network, the first value is considered for all
non-reference interventions and the second value is considered for the
reference intervention of the network (see 'Argument' `ref` in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).
This is necessary to ensure transitivity in the assumptions for the
missingness parameter across the comparisons in the network (Spineli,
2019b).

Currently, there are no empirically-based prior distributions for the
informative missingness parameters. The users may refer to Mavridis et
al. (2015) and Spineli (2019) to determine `mean_misspar` for an
informative missingness parameter.

## References

Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
uncertainty due to missing continuous outcome data in pairwise and
network meta-analysis. *Stat Med* 2015;**34**(5):721–41. doi:
10.1002/sim.6365

Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing
outcome data in network meta-analysis: a one-stage pattern-mixture model
approach. *Stat Methods Med Res* 2021;**30**(4):958–75. doi:
10.1177/0962280220983544

Spineli LM. An empirical comparison of Bayesian modelling strategies for
missing binary outcome data in network meta-analysis. *BMC Med Res
Methodol* 2019a;**19**(1):86. doi: 10.1186/s12874-019-0731-y

Spineli LM. Modeling missing binary outcome data while preserving
transitivity assumption yielded more credible network meta-analysis
results. *J Clin Epidemiol* 2019b;**105**:19–26. doi:
10.1016/j.jclinepi.2018.09.002

Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account
for uncertainty due to missing binary outcome data in pairwise
meta-analysis. *Stat Med* 2015;**34**(12):2062–80. doi: 10.1002/sim.6475

White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing
data in meta-analysis–part 1: two-stage methods. *Stat Med*
2008;**27**(5):711–27. doi: 10.1002/sim.3008

## See also

[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)

## Author

Loukia M. Spineli
