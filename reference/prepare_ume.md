# WinBUGS code for the unrelated mean effects model

The WinBUGS code, as proposed by Dias et al. (2013) to run a one-stage
Bayesian unrelated mean effects model, refined (Spineli, 2021), and
extended to incorporate the pattern-mixture model for binary or
continuous missing participant outcome data (Spineli et al., 2021;
Spineli, 2019).

## Usage

``` r
prepare_ume(measure, model, assumption, connected)
```

## Arguments

- measure:

  Character string indicating the effect measure with values `"OR"`,
  `"MD"`, `"SMD"`, or `"ROM"` for the odds ratio, mean difference,
  standardised mean difference and ratio of means, respectively.

- model:

  Character string indicating the analysis model with values `"RE"`, or
  `"FE"` for the random-effects and fixed-effect model, respectively.
  The default argument is `"RE"`.

- assumption:

  Character string indicating the structure of the informative
  missingness parameter. Set `assumption` equal to one of the following:
  `"HIE-COMMON"`, `"HIE-TRIAL"`, `"HIE-ARM"`, `"IDE-COMMON"`,
  `"IDE-TRIAL"`, `"IDE-ARM"`, `"IND-CORR"`, or `"IND-UNCORR"`. The
  default argument is `"IDE-ARM"`. The abbreviations `"IDE"`, `"HIE"`,
  and `"IND"` stand for identical, hierarchical and independent,
  respectively. `"CORR"` and `"UNCORR"` stand for correlated and
  uncorrelated, respectively.

- connected:

  An integer equal to one or larger that indicates the number of
  subnetworks.

## Value

An R character vector object to be passed to
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md)
through the
[`textConnection`](https://rdrr.io/r/base/textconnections.html) function
as the argument `object`.

## Details

This functions creates the model in the JAGS dialect of the BUGS
language. The output of this function constitutes the argument
`model.file` of [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) (in
the R-package [R2jags](https://CRAN.R-project.org/package=R2jags)) via
the [`textConnection`](https://rdrr.io/r/base/textconnections.html)
function.

`prepare_ume` inherits `measure`, `model`, and `assumption` from the
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
function. For a binary outcome, when `measure` is "RR" (relative risk)
or "RD" (risk difference) in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
`prepare_ume` currently considers the WinBUGS code for the odds ratio.

## References

Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence
synthesis for decision making 4: inconsistency in networks of evidence
based on randomized controlled trials. *Med Decis Making*
2013;**33**(5):641–56. doi: 10.1177/0272989X12455847

Spineli LM. A revised framework to evaluate the consistency assumption
globally in a network of interventions. *Med Decis Making* 2021. doi:
10.1177/0272989X211068005

Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing
outcome data in network meta-analysis: a one-stage pattern-mixture model
approach. *Stat Methods Med Res* 2021;**30**(4):958–75. doi:
10.1177/0962280220983544

Spineli LM. An empirical comparison of Bayesian modelling strategies for
missing binary outcome data in network meta-analysis. *BMC Med Res
Methodol* 2019;**19**(1):86. doi: 10.1186/s12874-019-0731-y

## See also

[`jags`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md),
[`textConnection`](https://rdrr.io/r/base/textconnections.html)

## Author

Loukia M. Spineli
