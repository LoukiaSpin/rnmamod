# Perform a series of Bayesian pairwise meta-analyses

Performs a Bayesian pairwise meta-analysis for each pairwise comparison
with at least two trials in the network.

## Usage

``` r
run_series_meta(full, n_chains, n_iter, n_burnin, n_thin, inits = NULL)
```

## Arguments

- full:

  An object of S3 class
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  See 'Value' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

- n_chains:

  Integer specifying the number of chains for the MCMC sampling; an
  argument of the [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html)
  function of the R-package
  [R2jags](https://CRAN.R-project.org/package=R2jags). The default
  argument is 2.

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

- inits:

  A list with the initial values for the parameters; an argument of the
  [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) function of the
  R-package [R2jags](https://CRAN.R-project.org/package=R2jags). The
  default argument is `NULL`, and JAGS generates the initial values.

## Value

An R2jags output on the summaries of the posterior distribution, and the
Gelman-Rubin convergence diagnostic (Gelman et al., 1992) of the
following monitored parameters:

- EM:

  The summary effect estimate (according to the argument `measure`
  defined in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md))
  of each observed pairwise comparison with at least two trials in the
  network.

- tau:

  The between-trial standard deviation for pairwise comparisons with at
  least two trials, when the random-effects model has been specified.

- single:

  A binary vector that indicates the comparisons in `EM` with one trial.

## Details

`run_series_meta` inherits the arguments `data`, `measure`, `model`,
`assumption`, `heter_prior`, `mean_misspar`, and `var_misspar` from
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
(now contained in the argument `full`). This prevents specifying a
different Bayesian model from that considered in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
Therefore, the user needs first to apply
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
and then use `run_series_meta` (see 'Examples').

For a binary outcome, when `measure` is "RR" (relative risk) or "RD"
(risk difference) in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
`run_series_meta` currently performs a series of pairwise meta-analysis
using the odds ratio as effect measure for being the **base-case**
effect measure in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
for a binary outcome (see also 'Details' in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).

`run_series_meta` runs a series of Bayesian pairwise meta-analyses in
`JAGS`. The progress of the simulation appears on the R console. The
number of times the function is used is also printed on the console (in
red) and is equal to the number of observed pairwise comparisons in the
network (see 'Examples'). The model is updated until convergence using
the [`autojags`](https://rdrr.io/pkg/R2jags/man/autojags.html) function
of the R-package [R2jags](https://CRAN.R-project.org/package=R2jags)
with 2 updates and number of iterations and thinning equal to `n_iter`
and `n_thin`, respectively.

The output of `run_series_meta` is not end-user-ready. The
[`series_meta_plot`](https://loukiaspin.github.io/rnmamod/reference/series_meta_plot.md)
function inherits the output of `run_series_meta` as an S3 object and
processes it further to provide an end-user-ready output.

`run_series_meta` can be used only for a network of interventions. In
the case of two interventions, the execution of the function will be
stopped and an error message will be printed on the R console.

## References

Gelman A, Rubin DB. Inference from iterative simulation using multiple
sequences. *Stat Sci* 1992;**7**(4):457–72. doi: 10.1214/ss/1177011136

## See also

[`jags`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`series_meta_plot`](https://loukiaspin.github.io/rnmamod/reference/series_meta_plot.md)

## Author

Loukia M. Spineli

## Examples

``` r
data("nma.dogliotti2014")

# Show the first six trials of the dataset (one-trial-per-row format)
head(nma.dogliotti2014)
#>              study t1 t2 t3   r1   r2 r3  m1  m2 m3   n1   n2 n3
#> 1     BAATAF, 1990  1  7 NA  195  188 NA   0  21 NA  208  212 NA
#> 2     SPINAF, 1992  1  7 NA  186  172 NA  56  81 NA  265  260 NA
#> 3    SPAF-II, 1994  2  7 NA  480  478 NA  23  38 NA  545  555 NA
#> 4      PATAF, 1999  2  7 NA  243   86 NA  54  42 NA  319  131 NA
#> 5 ACTIVE (W), 2006  3  7 NA 2825 3089 NA 410 223 NA 3335 3371 NA
#> 6       JAST, 2006  1  2 NA  338  313 NA  89  96 NA  445  426 NA

# Read results from 'run_model' (using the default arguments)
res <- readRDS(system.file('extdata/res_dogliotti.rds', package = 'rnmamod'))

# \donttest{
# Run separate random-effects pairwise meta-analyses
# Note: Ideally, set 'n_iter' to 10000 and 'n_burnin' to 1000
run_series_meta(full = res,
                n_chains = 3,
                n_iter = 1000,
                n_burnin = 100,
                n_thin = 1)
#> JAGS generates initial values for the parameters.
#> 1 out of 11 observed comparisons (and model updating until convergence)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 20
#>    Unobserved stochastic nodes: 29
#>    Total graph size: 529
#> 
#> Initializing model
#> 
#> 2 out of 11 observed comparisons (and model updating until convergence)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 28
#>    Unobserved stochastic nodes: 37
#>    Total graph size: 696
#> 
#> Initializing model
#> 
#> 3 out of 11 observed comparisons (and model updating until convergence)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 13
#>    Total graph size: 192
#> 
#> Initializing model
#> 
#> 4 out of 11 observed comparisons (and model updating until convergence)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 25
#>    Total graph size: 444
#> 
#> Initializing model
#> 
#> 5 out of 11 observed comparisons (and model updating until convergence)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 13
#>    Total graph size: 192
#> 
#> Initializing model
#> 
#> 6 out of 11 observed comparisons (and model updating until convergence)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 13
#>    Total graph size: 192
#> 
#> Initializing model
#> 
#> 7 out of 11 observed comparisons (and model updating until convergence)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 13
#>    Total graph size: 192
#> 
#> Initializing model
#> 
#> 8 out of 11 observed comparisons (and model updating until convergence)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 13
#>    Total graph size: 192
#> 
#> Initializing model
#> 
#> 9 out of 11 observed comparisons (and model updating until convergence)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 13
#>    Total graph size: 192
#> 
#> Initializing model
#> 
#> 10 out of 11 observed comparisons (and model updating until convergence)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 13
#>    Total graph size: 192
#> 
#> Initializing model
#> 
#> 11 out of 11 observed comparisons (and model updating until convergence)
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 4
#>    Unobserved stochastic nodes: 13
#>    Total graph size: 192
#> 
#> Initializing model
#> 
#> $EM
#>    t1 t2        mean        sd       2.5%          25%         50%       75%
#> 1   1  7  1.14720837 0.4414257  0.2662862  0.875972006  1.15272998 1.4180296
#> 2   2  7  0.50962570 0.4106020 -0.2759662  0.243483964  0.50282502 0.7620356
#> 3   3  7  0.65912448 0.9836207 -1.4072967  0.222874012  0.63554987 1.0856227
#> 4   1  2  0.35124473 0.3874048 -0.5403764  0.147158627  0.37156562 0.5892799
#> 5   2  3  0.39831713 1.0179364 -1.8903159 -0.001420263  0.40416854 0.8366353
#> 6   6  7 -0.15961173 1.0008721 -2.1353603 -0.679548047 -0.18712012 0.2952978
#> 7   7  8  0.22496861 1.0061226 -1.9882075 -0.152483727  0.22671009 0.5827302
#> 8   2  8  0.79353330 1.0304747 -1.3781104  0.328469165  0.80759058 1.2601346
#> 9   4  5  0.34846684 1.0658458 -2.0286782 -0.138776975  0.37767385 0.8323036
#> 10  4  7 -0.05029033 1.0941014 -2.4185521 -0.564956653 -0.04419221 0.4696194
#> 11  5  7 -0.27253577 1.0643378 -2.6973038 -0.707568242 -0.29205356 0.1887347
#>       97.5%     Rhat n.eff
#> 1  2.018798 1.002606   950
#> 2  1.339162 1.060231    39
#> 3  2.902930 1.000962  3000
#> 4  1.012881 1.025843   390
#> 5  2.591209 1.004379   920
#> 6  2.058547 1.011551   230
#> 7  2.413475 1.002548   980
#> 8  3.052401 1.000660  3000
#> 9  2.592051 1.007375   420
#> 10 2.237218 1.001422  2200
#> 11 1.984132 1.006760   320
#> 
#> $tau
#>   t1 t2      mean        sd       2.5%       25%       50%       75%    97.5%
#> 1  1  7 0.5283772 0.3559817 0.03293219 0.2781101 0.4665814 0.7033896 1.385347
#> 2  2  7 0.5849285 0.3203662 0.15112306 0.3706364 0.5289242 0.7219508 1.352580
#> 4  1  2 0.3288302 0.3060678 0.01951776 0.1178705 0.2355982 0.4437842 1.113341
#>       Rhat n.eff
#> 1 1.076793    86
#> 2 1.036147    66
#> 4 1.047846    61
#> 
#> $measure
#> [1] "OR"
#> 
#> $model
#> [1] "RE"
#> 
#> $single
#>  [1] 0 0 1 0 1 1 1 1 1 1 1
#> 
#> $n_chains
#> [1] 3
#> 
#> $n_iter
#> [1] 1000
#> 
#> $n_burnin
#> [1] 100
#> 
#> $n_thin
#> [1] 1
#> 
#> attr(,"class")
#> [1] "run_series_meta"
# }
```
