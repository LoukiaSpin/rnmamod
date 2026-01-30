# Perform the node-splitting approach

Performs the Bayesian node-splitting approach of Dias et al. (2010)
extended to address aggregate binary and continuous missing participant
outcome data via the pattern-mixture model (Spineli et al., 2021;
Spineli, 2019). This model offers a local evaluation of the plausibility
of the consistency assumption in the network (Dias et al., 2010).

## Usage

``` r
run_nodesplit(full, n_chains, n_iter, n_burnin, n_thin, inits = NULL)
```

## Arguments

- full:

  An object of S3 class
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  See 'Value' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

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

- inits:

  A list with the initial values for the parameters; an argument of the
  [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) function of the
  R-package [R2jags](https://CRAN.R-project.org/package=R2jags). The
  default argument is `NULL`, and JAGS generates the initial values.

## Value

An R2jags output on the summaries of the posterior distribution, and the
Gelman-Rubin convergence diagnostic of the following monitored
parameters:

- direct:

  The summary effect measure (according to the argument `measure`
  defined in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md))
  of each split node based on the corresponding trials.

- indirect:

  The indirect summary effect measure (according to the argument
  `measure` defined in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md))
  of each split node based on the remaining network after removing
  (splitting) the corresponding node.

- diff:

  The inconsistency parameter for each split node defined as the
  difference between the direct and indirect effect of the corresponding
  split node.

- p_value:

  The two-sided Bayesian p-value (based on the posterior mean on
  step(diff)) for each split node defined.

- tau:

  The between-trial standard deviation after each split node, when the
  random-effects model has been specified.

Furthermore, the output includes the following element:

- model_assessment:

  A data-frame on the measures of model assessment after each split
  node: deviance information criterion, total residual deviance, and
  number of effective parameters.

## Details

`run_nodesplit` inherits the arguments `data`, `measure`, `model`,
`assumption`, `heter_prior`, `mean_misspar`, `var_misspar`, `ref`, and
`indic` from
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
(now contained in the argument `full`). This prevents specifying a
different Bayesian model from that considered in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
Therefore, the user needs first to apply
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
and then use `run_nodesplit` (see 'Examples').

For a binary outcome, when `measure` is "RR" (relative risk) or "RD"
(risk difference) in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
`run_nodesplit` currently performs node-splitting using the odds ratio
as effect measure for being the **base-case** effect measure in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
for a binary outcome (see also 'Details' in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).

To perform the Bayesian node-splitting approach, the
[`prepare_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/prepare_nodesplit.md)
function is called which contains the WinBUGS code as written by Dias et
al. (2010) for binomial and normal likelihood to analyse binary and
continuous outcome data, respectively.
[`prepare_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/prepare_nodesplit.md)
has been extended to incorporate the pattern-mixture model with
informative missingness parameters for binary and continuous outcome
data (see 'Details' in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).

`run_nodesplit` runs the Bayesian node-splitting approach in `JAGS`. The
progress of the simulation appears on the R console. The number of times
`run_nodesplit` is used appears on the R console as a text in red and it
equals the number of split nodes (see 'Examples'). If there are no split
nodes in the network, the execution of the function will be stopped and
an error message will be printed on the R console. The model is updated
until convergence using the
[`autojags`](https://rdrr.io/pkg/R2jags/man/autojags.html) function of
the R-package [R2jags](https://CRAN.R-project.org/package=R2jags) with 2
updates and number of iterations and thinning equal to `n_iter` and
`n_thin`, respectively.

`run_nodesplit` uses the
[`mtc.nodesplit.comparisons`](https://rdrr.io/pkg/gemtc/man/mtc.nodesplit.html)
function of the R-package
[gemtc](https://CRAN.R-project.org/package=gemtc) to obtain
automatically the nodes to split based on the decision rule of van
Valkenhoef et al. (2016). `run_nodesplit` uses the option (1) in van
Valkenhoef et al. (2016) to parameterise multi-arm trials that contain
the node-to-split. In contrast,
[`mtc.nodesplit.comparisons`](https://rdrr.io/pkg/gemtc/man/mtc.nodesplit.html)
uses the option (3) in van Valkenhoef et al. (2016). Option (1) keeps
the baseline arm of the node-to-split in the corresponding multi-arms.
Option (3) excludes both arms of the node-to-split from the
corresponding multi-arm trials.

The output of `run_nodesplit` is not end-user-ready. The
[`nodesplit_plot`](https://loukiaspin.github.io/rnmamod/reference/nodesplit_plot.md)
function inherits the output of `run_nodesplit` as an S3 object and
processes it further to provide an end-user-ready output.

`run_nodesplit` can be used only for a network of interventions. In the
case of two interventions, the execution of the function will be stopped
and an error message will be printed on the R console.

## References

Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
treatment comparison meta-analysis. *Stat Med* 2010;**29**(7-8):932–44.
doi: 10.1002/sim.3767

Gelman A, Rubin DB. Inference from iterative simulation using multiple
sequences. *Stat Sci* 1992;**7**(4):457–72. doi: 10.1214/ss/1177011136

Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing
outcome data in network meta-analysis: a one-stage pattern-mixture model
approach. *Stat Methods Med Res* 2021;**30**(4):958–75. doi:
10.1177/0962280220983544

Spineli LM. An empirical comparison of Bayesian modelling strategies for
missing binary outcome data in network meta-analysis. *BMC Med Res
Methodol* 2019;**19**(1):86. doi: 10.1186/s12874-019-0731-y

van Valkenhoef G, Dias S, Ades AE, Welton NJ. Automated generation of
node-splitting models for assessment of inconsistency in network
meta-analysis. *Res Synth Methods* 2016;**7**(1):80–93. doi:
10.1002/jrsm.1167

## See also

[`autojags`](https://rdrr.io/pkg/R2jags/man/autojags.html),
[`jags`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`mtc.nodesplit.comparisons`](https://rdrr.io/pkg/gemtc/man/mtc.nodesplit.html),
[`nodesplit_plot`](https://loukiaspin.github.io/rnmamod/reference/nodesplit_plot.md),
[`prepare_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/prepare_nodesplit.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)

## Author

Loukia M. Spineli

## Examples

``` r
data("nma.baker2009")

# Read results from 'run_model' (using the default arguments)
res <- readRDS(system.file('extdata/res_baker.rds', package = 'rnmamod'))

# \donttest{
# Run random-effects node-splitting approach
# Note: Ideally, set 'n_iter' to 10000 and 'n_burnin' to 1000
run_nodesplit(full = res,
              n_chains = 3,
              n_iter = 1000,
              n_burnin = 100,
              n_thin = 1)
#> JAGS generates initial values for the parameters.
#> 1 out of 8 split nodes
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 79
#>    Unobserved stochastic nodes: 141
#>    Total graph size: 2526
#> 
#> Initializing model
#> 
#> Updating model for split node 1 until convergence
#> 2 out of 8 split nodes
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 79
#>    Unobserved stochastic nodes: 141
#>    Total graph size: 2526
#> 
#> Initializing model
#> 
#> Updating model for split node 2 until convergence
#> 3 out of 8 split nodes
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 79
#>    Unobserved stochastic nodes: 141
#>    Total graph size: 2526
#> 
#> Initializing model
#> 
#> Updating model for split node 3 until convergence
#> 4 out of 8 split nodes
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 79
#>    Unobserved stochastic nodes: 141
#>    Total graph size: 2535
#> 
#> Initializing model
#> 
#> Updating model for split node 4 until convergence
#> 5 out of 8 split nodes
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 79
#>    Unobserved stochastic nodes: 141
#>    Total graph size: 2535
#> 
#> Initializing model
#> 
#> Updating model for split node 5 until convergence
#> 6 out of 8 split nodes
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 79
#>    Unobserved stochastic nodes: 141
#>    Total graph size: 2535
#> 
#> Initializing model
#> 
#> Updating model for split node 6 until convergence
#> 7 out of 8 split nodes
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 79
#>    Unobserved stochastic nodes: 141
#>    Total graph size: 2526
#> 
#> Initializing model
#> 
#> Updating model for split node 7 until convergence
#> 8 out of 8 split nodes
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 79
#>    Unobserved stochastic nodes: 141
#>    Total graph size: 2531
#> 
#> Initializing model
#> 
#> Updating model for split node 8 until convergence
#> $direct
#>   treat1 treat2       mean        sd       2.5%       97.5%     Rhat n.eff
#> 1      6      1 -0.4764911 0.2381483 -0.9264469 -0.01406681 1.415765     8
#> 2      7      1 -0.3745740 0.2696855 -0.9424315  0.09418342 1.079514    38
#> 3      8      1 -0.4664331 0.2582532 -0.9715420 -0.04043404 2.061279     5
#> 4      5      4 -0.2280042 0.2406259 -0.6935732  0.22504783 1.042505    52
#> 5      7      4 -0.5827774 0.2283905 -1.0534659 -0.13836142 1.055178   280
#> 6      7      5 -0.3770219 0.2278810 -0.8344568  0.04813429 1.028349    84
#> 7      8      6 -1.2160887 0.6983302 -2.6804946  0.07677250 1.167290    19
#> 8      8      7 -0.8195170 0.3265053 -1.5663669 -0.19284957 1.133646    21
#> 
#> $indirect
#>   treat1 treat2         mean        sd       2.5%       97.5%     Rhat n.eff
#> 1      6      1  0.942626515 0.7152297 -0.4308032  2.42841815 1.031156   110
#> 2      7      1  0.007416371 0.7308743 -1.4404710  1.50464085 1.172967    18
#> 3      8      1 -1.310061728 0.5889717 -2.4076024 -0.31119339 1.165150    18
#> 4      5      4 -0.782616933 0.4322692 -1.6951215 -0.01684944 1.037721    69
#> 5      7      4 -0.224408912 0.2038462 -0.6462462  0.12099295 1.156643    18
#> 6      7      5 -0.356059158 0.2981125 -0.8934331  0.24034101 1.067769    40
#> 7      8      6 -0.103020577 0.2223173 -0.5891935  0.35522987 1.043316   100
#> 8      8      7 -0.036842314 0.1979648 -0.4407604  0.35640351 1.021216   100
#> 
#> $diff
#>   treat1 treat2       mean        sd        2.5%        97.5%     Rhat n.eff
#> 1      6      1 -1.4191176 0.7432236 -2.91328816  0.057931461 1.100386    27
#> 2      7      1 -0.3819904 0.7380319 -1.90395303  0.885808392 1.222356    13
#> 3      8      1  0.8436286 0.5124235 -0.08009813  1.897718957 1.016184   190
#> 4      5      4  0.5546128 0.4498314 -0.30779607  1.475054163 1.017742   490
#> 5      7      4 -0.3583685 0.3012351 -1.01576870  0.221783417 1.073042    51
#> 6      7      5 -0.0209627 0.4207426 -0.85822525  0.717934625 1.057364    44
#> 7      8      6 -1.1130681 0.7220204 -2.66363363  0.247829733 1.118242    26
#> 8      8      7 -0.7826746 0.3952756 -1.65657253 -0.005988471 1.074462    35
#> 
#> $p_value
#>   treat1 treat2    p_value
#> 1      6      1 0.06866667
#> 2      7      1 0.62933333
#> 3      8      1 0.08066667
#> 4      5      4 0.20466667
#> 5      7      4 0.24666667
#> 6      7      5 0.97066667
#> 7      8      6 0.15333333
#> 8      8      7 0.04866667
#> 
#> $tau
#>   treat1 treat2       50%         sd         2.5%     97.5%     Rhat n.eff
#> 1      6      1 0.1372632 0.08624545 0.0182242886 0.3353547 1.381202    11
#> 2      7      1 0.1467188 0.08777046 0.0232173297 0.3534738 1.131333    24
#> 3      8      1 0.1352372 0.09991365 0.0154244468 0.3653359 1.221620    15
#> 4      5      4 0.1311072 0.08410075 0.0291202896 0.3372812 1.101564    25
#> 5      7      4 0.1174756 0.08255381 0.0063959691 0.2976880 1.212370    16
#> 6      7      5 0.1749594 0.09848069 0.0367016911 0.4113188 1.178860    16
#> 7      8      6 0.1326717 0.10307751 0.0008973746 0.3754917 1.367326    12
#> 8      8      7 0.1647816 0.09657632 0.0296912124 0.3814366 1.129143    25
#> 
#> $model
#> [1] "RE"
#> 
#> $model_assessment
#>   treat1 treat2 deviance      DIC       pD
#> 1      6      1 53.06219 89.05718 35.99499
#> 2      7      1 53.37384 91.06124 37.68740
#> 3      8      1 53.41883 90.81611 37.39729
#> 4      5      4 54.21706 92.14342 37.92635
#> 5      7      4 52.39749 88.96004 36.56255
#> 6      7      5 53.34828 91.44039 38.09211
#> 7      8      6 53.57174 89.77703 36.20529
#> 8      8      7 53.30739 91.30597 37.99857
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
#> [1] "run_nodesplit"
# }
```
