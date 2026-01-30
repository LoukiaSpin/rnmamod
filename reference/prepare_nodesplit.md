# WinBUGS code for the node-splitting approach

The WinBUGS code, as written by Dias et al. (2010) to run a one-stage
Bayesian node-splitting model, extended to incorporate the
pattern-mixture model for binary or continuous missing participant
outcome data (Spineli et al., 2021; Spineli, 2019).

## Usage

``` r
prepare_nodesplit(measure, model, assumption)
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
[`run_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/run_nodesplit.md)
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

`prepare_nodesplit` inherits `measure`, `model`, and `assumption` from
the
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
function. For a binary outcome, when `measure` is "RR" (relative risk)
or "RD" (risk difference) in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
`prepare_nodesplit` currently considers the WinBUGS code for the odds
ratio.

The split nodes have been automatically selected via the
[`mtc.nodesplit.comparisons`](https://rdrr.io/pkg/gemtc/man/mtc.nodesplit.html)
function of the R-package
[gemtc](https://CRAN.R-project.org/package=gemtc). See 'Details' in
[`run_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/run_nodesplit.md).

## References

Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
treatment comparison meta-analysis. *Stat Med* 2010;**29**(7-8):932–44.
doi: 10.1002/sim.3767

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
[`mtc.nodesplit.comparisons`](https://rdrr.io/pkg/gemtc/man/mtc.nodesplit.html),
[`run_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/run_nodesplit.md),
[`textConnection`](https://rdrr.io/r/base/textconnections.html)

## Author

Loukia M. Spineli
