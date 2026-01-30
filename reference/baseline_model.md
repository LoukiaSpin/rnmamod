# The baseline model for binary outcome

To process the elements in the argument `base_risk` of the
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
function. It also runs the hierarchical baseline model, separately from
the relative effects model as described in Dias et al. (2018) and Dias
et al. (2013b). The output is to be passed to
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
and
[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md)
to obtain the (unadjusted and adjusted, respectively) absolute risks for
each intervention in the dataset.

## Usage

``` r
baseline_model(base_risk, n_chains, n_iter, n_burnin, n_thin)
```

## Arguments

- base_risk:

  A scalar, a vector of length three with elements sorted in ascending
  order, or a matrix with two columns and number of rows equal to the
  number of relevant trials. In the case of a scalar or vector, the
  elements should be in the interval (0, 1). For the matrix, the first
  column refers to the number of events and the second column to the
  sample size of the trials comprising the dataset for the baseline
  model. See 'Details' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  This argument is only relevant for a binary outcome.

- n_chains:

  Positive integer specifying the number of chains for the MCMC
  sampling; an argument of the
  [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) function of the
  R-package [R2jags](https://CRAN.R-project.org/package=R2jags). The
  default argument is 2.

- n_iter:

  Positive integer specifying the number of Markov chains for the MCMC
  sampling; an argument of the
  [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) function of the
  R-package [R2jags](https://CRAN.R-project.org/package=R2jags). The
  default argument is 10000.

- n_burnin:

  Positive integer specifying the number of iterations to discard at the
  beginning of the MCMC sampling; an argument of the
  [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) function of the
  R-package [R2jags](https://CRAN.R-project.org/package=R2jags). The
  default argument is 1000.

- n_thin:

  Positive integer specifying the thinning rate for the MCMC sampling;
  an argument of the [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html)
  function of the R-package
  [R2jags](https://CRAN.R-project.org/package=R2jags). The default
  argument is 1.

## Value

When `base_risk` is scalar (fixed baseline), the function returns the
user-defined baseline for the selected reference intervention in the
logit scale. When `base_risk` is a vector (random baseline), the
function returns a vector with the calculated logit of an event for the
selected reference intervention and its precision. Finally, when
`base_risk` is a matrix (predicted baseline), the function returns the
following elements:

- ref_base:

  A vector with the posterior mean and precision of the **predicted**
  logit of an event for the selected reference intervention. This vector
  is be passed to
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
  and
  [`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md).

- figure:

  A forest plot on the trial-specific observed and estimated baseline
  risk. See 'Details'.

- table_baseline:

  A table with the posterior and predictive distribution of the summary
  baseline mean and the posterior distribution of the between-trial
  standard deviation in baseline. All results are in the logit scale.

## Details

If `base_risk` is a matrix, `baseline_model` creates the hierarchical
baseline model in the JAGS dialect of the BUGS language. The output of
this function (see 'Value') constitutes the posterior mean and precision
of the predicted logit of an event for the selected reference
intervention and it is plugged in the WinBUGS code for the relative
effects model (Dias et al., 2013a) via the
[`prepare_model`](https://loukiaspin.github.io/rnmamod/reference/prepare_model.md)
function. Following (Dias et al., 2013a), a uniform prior distribution
is assigned on the between-trial standard deviation with upper and lower
limit equal to 0 and 5, respectively.

When `base_risk` is a matrix, the function also returns a forest plot
with the estimated trial-specific probability of an event and 95%
credible intervals (the random effects) alongside the corresponding
observed probability of an event for the selected reference
intervention. A grey rectangular illustrates the summary mean and 95%
credible interval of the random effects.

When `base_risk` is a matrix (predicted baseline), the model is updated
until convergence using the
[`autojags`](https://rdrr.io/pkg/R2jags/man/autojags.html) function of
the R-package [R2jags](https://CRAN.R-project.org/package=R2jags) with 2
updates and number of iterations and thinning equal to `n_iter` and
`n_thin`, respectively.

## References

Dias S, Ades AE, Welton NJ, Jansen JP, Sutton AJ. Network Meta-Analysis
for Decision Making. Chichester (UK): Wiley; 2018.

Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
making 2: a generalized linear modeling framework for pairwise and
network meta-analysis of randomized controlled trials. *Med Decis
Making* 2013a;**33**(5):607–17. doi: 10.1177/0272989X12458724

Dias S, Welton NJ, Sutton AJ, Ades AE. Evidence synthesis for decision
making 5: the baseline natural history model. *Med Decis Making*
2013b;**33**(5):657–70. doi: 10.1177/0272989X13485155

## See also

[`prepare_model`](https://loukiaspin.github.io/rnmamod/reference/prepare_model.md),
[`autojags`](https://rdrr.io/pkg/R2jags/man/autojags.html),
[`jags`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)

## Author

Loukia M. Spineli
