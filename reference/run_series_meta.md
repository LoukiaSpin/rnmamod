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
#>    t1 t2       mean        sd       2.5%        25%        50%        75%
#> 1   1  7  1.0206227 0.4174857  0.1485436  0.7647818  1.0339332 1.26676073
#> 2   2  7  0.4814008 0.3851172 -0.2699264  0.2374215  0.4710247 0.70925941
#> 3   3  7  0.6128957 1.0883626 -1.6687571  0.1743553  0.6213914 1.11748338
#> 4   1  2  0.3853745 0.3315957 -0.2324264  0.1509196  0.3667991 0.62759298
#> 5   2  3  0.3410806 1.0819657 -1.9431615 -0.1241617  0.3276144 0.77489569
#> 6   6  7 -0.1311163 1.0542260 -2.4232290 -0.6605076 -0.1463187 0.38065654
#> 7   7  8  0.2262757 0.9709205 -1.9097594 -0.1205395  0.2354991 0.58009804
#> 8   2  8  0.7370633 1.0978805 -1.6007735  0.2563489  0.7654580 1.22779705
#> 9   4  5  0.3915420 1.0225743 -1.6851867 -0.1013331  0.3917442 0.85568626
#> 10  4  7 -0.1938794 1.0407208 -2.4089293 -0.6658701 -0.1849046 0.30340363
#> 11  5  7 -0.4310873 1.0480695 -2.6692603 -0.8630550 -0.4013567 0.02603885
#>       97.5%     Rhat n.eff
#> 1  1.862921 1.004616   490
#> 2  1.278263 1.024333   410
#> 3  2.783911 1.006684  1100
#> 4  1.022280 1.131041    20
#> 5  2.666976 1.002026  1300
#> 6  2.095051 1.015413   200
#> 7  2.326279 1.004925  1900
#> 8  3.004225 1.002160  1300
#> 9  2.644165 1.000831  3000
#> 10 2.069731 1.002093  2800
#> 11 1.663752 1.002505  1700
#> 
#> $tau
#>   t1 t2      mean        sd       2.5%        25%       50%       75%     97.5%
#> 1  1  7 0.5219022 0.3391181 0.04704861 0.27811246 0.4640576 0.6998448 1.3803677
#> 2  2  7 0.5574011 0.2821244 0.14494251 0.36097865 0.5078711 0.6988280 1.2400000
#> 4  1  2 0.2677859 0.2201964 0.01622410 0.09905223 0.2203420 0.3730381 0.8353221
#>       Rhat n.eff
#> 1 1.031029   110
#> 2 1.015841   170
#> 4 1.019281   110
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
