# WinBUGS code for Bayesian pairwise or network meta-analysis and meta-regression

The WinBUGS code, as written by Dias et al. (2013) to run a one-stage
Bayesian network meta-analysis, extended to incorporate the
pattern-mixture model for binary or continuous missing participant
outcome data (Spineli et al., 2021; Spineli, 2019). The model has been
also extended to incorporate a trial-level covariate to apply
meta-regression (Cooper et al., 2009). In the case of two interventions,
the code boils down to a one-stage Bayesian pairwise meta-analysis with
pattern-mixture model (Turner et al., 2015; Spineli et al, 2021).

## Usage

``` r
prepare_model(measure, model, covar_assumption, assumption)
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

- covar_assumption:

  Character string indicating the structure of the
  intervention-by-covariate interaction, as described in Cooper et al.,
  (2009). Set `covar_assumption` equal to one of the following, when
  meta-regression is performed: `"exchangeable"`, `"independent"`, and
  `"common"`. Assign `"NO"` to perform pairwise or network
  meta-analysis.

- assumption:

  Character string indicating the structure of the informative
  missingness parameter. Set `assumption` equal to one of the following:
  `"HIE-COMMON"`, `"HIE-TRIAL"`, `"HIE-ARM"`, `"IDE-COMMON"`,
  `"IDE-TRIAL"`, `"IDE-ARM"`, `"IND-CORR"`, or `"IND-UNCORR"`. The
  default argument is `"IDE-ARM"`. The abbreviations `"IDE"`, `"HIE"`,
  and `"IND"` stand for identical, hierarchical and independent,
  respectively. `"CORR"` and `"UNCORR"` stand for correlated and
  uncorrelated, respectively.

## Value

An R character vector object to be passed to
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
and
[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md)
through the
[`textConnection`](https://rdrr.io/r/base/textconnections.html) function
as the argument `object`.

## Details

`prepare_model` creates the model in the JAGS dialect of the BUGS
language. The output of this function constitutes the argument
`model.file` of the [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html)
function (in the R-package
[R2jags](https://CRAN.R-project.org/package=R2jags)) via the
[`textConnection`](https://rdrr.io/r/base/textconnections.html)
function.

## References

Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing
between-study heterogeneity and inconsistency in mixed treatment
comparisons: Application to stroke prevention treatments in individuals
with non-rheumatic atrial fibrillation. *Stat Med*
2009;**28**(14):1861–81. doi: 10.1002/sim.3594

Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
making 2: a generalized linear modeling framework for pairwise and
network meta-analysis of randomized controlled trials. *Med Decis
Making* 2013;**33**(5):607–17. doi: 10.1177/0272989X12458724

Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing
outcome data in network meta-analysis: a one-stage pattern-mixture model
approach. *Stat Methods Med Res* 2021;**30**(4):958–75. doi:
10.1177/0962280220983544

Spineli LM. An empirical comparison of Bayesian modelling strategies for
missing binary outcome data in network meta-analysis. *BMC Med Res
Methodol* 2019;**19**(1):86. doi: 10.1186/s12874-019-0731-y

Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account
for uncertainty due to missing binary outcome data in pairwise
meta-analysis. *Stat Med* 2015;**34**(12):2062–80. doi: 10.1002/sim.6475

## See also

[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`jags`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`textConnection`](https://rdrr.io/r/base/textconnections.html)

## Author

Loukia M. Spineli
