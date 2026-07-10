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
#> 1      6      1 -0.3475195 0.2602996 -0.8250056  0.14883633 1.127628    22
#> 2      7      1 -0.6430262 0.2319460 -1.0959832 -0.21210808 1.037831    71
#> 3      8      1 -0.4050210 0.1678153 -0.7639951 -0.08004837 1.021618   140
#> 4      5      4 -0.1899555 0.2080949 -0.6273254  0.19779258 1.015704   190
#> 5      7      4 -0.5709905 0.2112652 -0.9607394 -0.15696131 1.136719    25
#> 6      7      5 -0.3608537 0.2313994 -0.8171038  0.08420671 1.031627    73
#> 7      8      6 -1.5415404 0.6420013 -2.7872876 -0.04214797 1.039615   330
#> 8      8      7 -0.8964306 0.3862680 -1.7145885 -0.15782588 1.247387    13
#> 
#> $indirect
#>   treat1 treat2        mean        sd       2.5%       97.5%     Rhat n.eff
#> 1      6      1  1.18185777 1.0172878 -0.4561300  3.06036788 2.138137     5
#> 2      7      1 -0.44282277 0.8149694 -2.0758064  0.86855275 1.081094    31
#> 3      8      1 -1.29301628 0.5283317 -2.2837218 -0.32314054 1.228704    13
#> 4      5      4 -0.67293295 0.3679087 -1.4270310  0.02705463 1.015950   390
#> 5      7      4 -0.18342210 0.1997755 -0.5461401  0.21901942 1.124704    21
#> 6      7      5 -0.39620448 0.2758297 -0.9447325  0.18625857 1.053474    57
#> 7      8      6 -0.06677049 0.2194732 -0.5204410  0.31541549 1.075757    34
#> 8      8      7 -0.02523183 0.1966435 -0.4200616  0.34142787 1.002945   820
#> 
#> $diff
#>   treat1 treat2        mean        sd        2.5%        97.5%     Rhat n.eff
#> 1      6      1 -1.52937728 0.9812317 -3.47376423  0.077611225 1.982397     5
#> 2      7      1 -0.20020344 0.7192733 -1.51415053  1.227400009 1.095414    29
#> 3      8      1  0.88799524 0.5463173 -0.04556213  1.899368146 1.213016    14
#> 4      5      4  0.48297748 0.3760610 -0.28218680  1.209727041 1.012588   280
#> 5      7      4 -0.38756835 0.2897400 -0.94496509  0.142615551 1.023809   840
#> 6      7      5  0.03535077 0.3912771 -0.82790747  0.789128152 1.038185    96
#> 7      8      6 -1.47476990 0.6700360 -2.63068712 -0.000297547 1.048805    83
#> 8      8      7 -0.87119878 0.4358669 -1.69174050  0.075655267 1.169309    17
#> 
#> $p_value
#>   treat1 treat2    p_value
#> 1      6      1 0.06600000
#> 2      7      1 0.73200000
#> 3      8      1 0.06266667
#> 4      5      4 0.21133333
#> 5      7      4 0.14066667
#> 6      7      5 0.84800000
#> 7      8      6 0.05000000
#> 8      8      7 0.05800000
#> 
#> $tau
#>   treat1 treat2       50%         sd         2.5%     97.5%     Rhat n.eff
#> 1      6      1 0.1215437 0.08246478 0.0082049142 0.3260568 1.103753    49
#> 2      7      1 0.1405358 0.08750373 0.0147942634 0.3361953 1.086028    34
#> 3      8      1 0.1217684 0.09416532 0.0148473826 0.3615808 1.314259    11
#> 4      5      4 0.1027161 0.08058699 0.0292031041 0.2880754 1.005836   770
#> 5      7      4 0.0870589 0.08095747 0.0006200626 0.2728706 1.389290    11
#> 6      7      5 0.1558986 0.10610953 0.0070656223 0.3937198 1.467487     9
#> 7      8      6 0.1240954 0.07847724 0.0141844211 0.3071928 1.161875    34
#> 8      8      7 0.1575697 0.09560836 0.0068516860 0.3738418 1.155318    35
#> 
#> $model
#> [1] "RE"
#> 
#> $model_assessment
#>   treat1 treat2 deviance      DIC       pD
#> 1      6      1 53.92127 90.42423 36.50296
#> 2      7      1 54.23491 91.29573 37.06082
#> 3      8      1 53.34470 88.85501 35.51031
#> 4      5      4 53.71467 91.13127 37.41660
#> 5      7      4 54.19413 91.75235 37.55822
#> 6      7      5 54.46978 92.80598 38.33620
#> 7      8      6 52.78928 88.71401 35.92473
#> 8      8      7 54.19800 92.15368 37.95568
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
