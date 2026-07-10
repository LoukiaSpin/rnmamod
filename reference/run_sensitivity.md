# Perform sensitivity analysis for missing participant outcome data

Performs a sensitivity analysis by applying pairwise meta-analysis or
network meta-analysis for a series of different scenarios about the
informative missingness parameter.

## Usage

``` r
run_sensitivity(
  full,
  assumption,
  mean_scenarios,
  var_misspar,
  n_chains,
  n_iter,
  n_burnin,
  n_thin,
  inits = NULL
)
```

## Arguments

- full:

  An object of S3 class
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  See 'Value' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

- assumption:

  Character string indicating the structure of the informative
  missingness parameter. Set `assumption` equal to one of the following
  two: `"HIE-ARM"`, or `"IDE-ARM"` (see 'Details'). The default argument
  is `"IDE-ARM"`. The abbreviations `"IDE"`, and `"HIE"` stand for
  identical, and hierarchical, respectively.

- mean_scenarios:

  A vector with numeric values for the mean of the normal distribution
  of the informative missingness parameter (see 'Details'). The vector
  should have a length equal to 5 or larger. The missing-at-random (MAR)
  assumption should be the median of the vector, so that the same number
  of informative scenarios appear before and after the MAR. The default
  scenarios are c(-log(3), -log(2), log(0.9999), log(2), log(3)) and
  c(-2, -1, 0, 1, 2) for binary and continuous outcome data,
  respectively.

- var_misspar:

  A positive non-zero number for the variance of the normal distribution
  of the informative missingness parameter. When the `measure` (defined
  in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md))
  is `"OR"`, `"MD"`, or `"SMD"` the default argument is 1. When the
  `measure` is `"ROM"`, the default argument is 0.04.

- n_chains:

  Integer specifying the number of chains for the MCMC sampling; an
  argument of the [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html)
  function of the R-package
  [R2jags](https://CRAN.R-project.org/package=R2jags). The default
  argument is 2.

- n_iter:

  Integer specifying the number of Markov chains for the MCMC sampling;
  an argument of the [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html)
  function of the R-package
  [R2jags](https://CRAN.R-project.org/package=R2jags). The default
  argument is 10000.

- n_burnin:

  Integer specifying the number of iterations to discard at the
  beginning of the MCMC sampling; an argument of the
  [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) function of the
  R-package [R2jags](https://CRAN.R-project.org/package=R2jags). The
  default argument is 1000.

- n_thin:

  Integer specifying the thinning rate for the MCMC sampling; an
  argument of the [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html)
  function of the R-package
  [R2jags](https://CRAN.R-project.org/package=R2jags). The default
  argument is 1.

- inits:

  A list with the initial values for the parameters; an argument of the
  [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) function of the
  R-package [R2jags](https://CRAN.R-project.org/package=R2jags). The
  default argument is `NULL`, and JAGS generates the initial values.

## Value

A list of R2jags outputs on the summaries of the posterior distribution,
and the Gelman-Rubin convergence diagnostic (Gelman et al., 1992) of the
following monitored parameters for a random-effects pairwise
meta-analysis:

- EM:

  The estimated summary effect measure (according to the argument
  `measure` defined in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).

- EM_pred:

  The predicted summary effect measure (according to the argument
  `measure` defined in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).
  This element does not appear in the case of a fixed-effect
  meta-analysis.

- EM_LOR:

  The estimated summary odd ratio in the logarithmic scale when
  `measure = "RR"` or `measure = "RD"`.

- EM_LOR_pred:

  The predicted summary odd ratio in the logarithmic scale when
  `measure = "RR"` or `measure = "RD"`. This element does not appear in
  the case of a fixed-effect meta-analysis.

- tau:

  The between-trial standard deviation. This element does not appear in
  the case of a fixed-effect meta-analysis.

In a random-effects network meta-analysis, `EM` refer to all possible
pairwise comparisons of interventions in the network. Furthermore, `tau`
is typically assumed to be common for all observed comparisons in the
network.

## Details

The model runs in `JAGS` and the progress of the simulation appears on
the R console. The number of times `run_sensitivity` is used appears on
the R console as a text in red and it equals the **number of scenarios**
defined as *the square of the length of the vector specified in
`mean_scenarios`* (see 'Examples'). The output of `run_sensitivity` is
used as an S3 object by other functions of the package to be processed
further and provide an end-user-ready output.

In the case of pairwise meta-analysis, `EM` and `tau` are estimated as
many times as the number of scenarios considered. In the case of network
meta-analysis, each possible pairwise comparison is estimated as many
times as the number of scenarios considered.

The informative missingness parameter is assumed to differ only across
the interventions of the dataset. Therefore, the user can specify the
informative missingness parameter to be arm-specific *and* identical
(`assumption = "IDE-ARM"`), or arm-specific *and* hierarchical
(`assumption = "HIE-ARM"`) (Spineli et al., 2021).

The length of the vector specified in argument `mean_scenarios` should
be equal to or more than 5 (a positive odd integer) to allow for an
adequate number of scenarios. It is important that the number
corresponding to the MAR assumption is the middle of the numbers in the
vector specified in argument `mean_scenarios`. The **MAR assumption**
constitutes the **primary analysis**. Under the informative missingness
difference of means parameter (relevant for the raw and standardised
mean diffenre), the MAR assumption equals 0. Under the informative
missingness odds ratio parameter (IMOR; relevant for the odds ratio) and
the informative missingness ratio of means (IMRoM; relevant for the
ratio of means) parameter, the MAR assumption equals 1; however, both
parameters are analysed in the logarithmic scale. We advise using the
value `0.999` rather than `1` in `mean_scenarios` for the IMOR and IMRoM
parameters; otherwise, the execution of the function will be stopped and
the error 'Invalid parent values' will be printed on the R console.

Currently, there are no empirically-based prior distributions for the
informative missingness parameters. The users may refer to Spineli
(2019), Mavridis et al. (2015), Turner et al. (2015), and White et al.
(2008) to determine `mean_scenarios` for an informative missingness
mechanism and select a proper value for `var_misspar`.

`run_sensitivity` inherits the arguments `data`, `measure`, `model`,
`heter_prior`, `D`, `indic`, `base_risk`, and `ref` from
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
(now contained in the argument `full`). This prevents specifying a
different Bayesian model from that considered in the primary analysis
(via
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md))–an
exception in the `assumption` argument as it is restricted to only two
character strings. Therefore, the user needs first to apply
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
and then use `run_sensitivity` (see 'Examples').

The `run_sensitivity` function also returns the arguments `measure`,
`scenarios`, `D`, `heter`, `n_chains`, `n_iter`, `n_burnin`, and
`n_thin` as specified by the user to be inherited by other relevant
functions of the package.

The model is updated until convergence using the
[`autojags`](https://rdrr.io/pkg/R2jags/man/autojags.html) function of
the R-package [R2jags](https://CRAN.R-project.org/package=R2jags) with 2
updates and number of iterations and thinning equal to `n_iter` and
`n_thin`, respectively.

`run_sensitivity` can be used only when missing participant outcome data
have been extracted for at least one trial. Otherwise, the execution of
the function will be stopped and an error message will be printed on the
R console.

## References

Gelman, A, Rubin, DB. Inference from iterative simulation using multiple
sequences. *Stat Sci* 1992;**7**(4):457–72. doi: 10.1214/ss/1177011136

Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
uncertainty due to missing continuous outcome data in pairwise and
network meta-analysis. *Stat Med* 2015;**34**(5):721–41. doi:
10.1002/sim.6365

Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness
of primary analysis results: A case study on missing outcome data in
pairwise and network meta-analysis. *Res Synth Methods*
2021;**12**(4):475–90. doi: 10.1002/jrsm.1478

Spineli LM. An empirical comparison of Bayesian modelling strategies for
missing binary outcome data in network meta-analysis. *BMC Med Res
Methodol* 2019;**19**(1):86. doi: 10.1186/s12874-019-0731-y

Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account
for uncertainty due to missing binary outcome data in pairwise
meta-analysis. *Stat Med* 2015;**34**(12):2062–80. doi: 10.1002/sim.6475

White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing
data in meta-analysis–part 1: two-stage methods. *Stat Med*
2008;**27**(5):711–27. doi: 10.1002/sim.3008

## See also

[`autojags`](https://rdrr.io/pkg/R2jags/man/autojags.html),
[`jags`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)

## Author

Loukia M. Spineli

## Examples

``` r
data("pma.taylor2004")

# Read results from 'run_model' (using the default arguments)
res <- readRDS(system.file('extdata/res_taylor.rds', package = 'rnmamod'))

# \donttest{
# Perform the sensitivity analysis (default arguments)
# Note: Ideally, set 'n_iter' to 10000 and 'n_burnin' to 1000
run_sensitivity(full = res,
                assumption = "IDE-ARM",
                var_misspar = 1,
                n_chains = 3,
                n_iter = 1000,
                n_burnin = 100,
                n_thin = 5)
#> The default scenarios were considered: c(-2, -1, 0, 1, 2).
#> JAGS generates initial values for the parameters.
#> 1 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 1 until convergence
#> 2 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 2 until convergence
#> 3 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 3 until convergence
#> 4 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 4 until convergence
#> 5 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 5 until convergence
#> 6 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 6 until convergence
#> 7 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 7 until convergence
#> 8 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 8 until convergence
#> 9 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 9 until convergence
#> 10 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 10 until convergence
#> 11 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 11 until convergence
#> 12 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 12 until convergence
#> 13 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 13 until convergence
#> 14 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 14 until convergence
#> 15 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 15 until convergence
#> 16 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 16 until convergence
#> 17 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 17 until convergence
#> 18 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 18 until convergence
#> 19 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 19 until convergence
#> 20 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 20 until convergence
#> 21 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 21 until convergence
#> 22 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 22 until convergence
#> 23 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 23 until convergence
#> 24 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 24 until convergence
#> 25 out of 25 total scenarios
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 16
#>    Unobserved stochastic nodes: 28
#>    Total graph size: 343
#> 
#> Initializing model
#> 
#> Updating model for scenario 25 until convergence
#> $EM
#>               mean        sd       2.5%        25%          50%         75%
#>  [1,] -0.098585142 0.2477695 -0.5448514 -0.2511135 -0.105667108  0.06666307
#>  [2,] -0.122173224 0.2834465 -0.5960651 -0.3081796 -0.125425163  0.04643160
#>  [3,] -0.140355224 0.2753617 -0.6427526 -0.3098226 -0.117604138  0.01232324
#>  [4,] -0.160552623 0.2527240 -0.6560794 -0.3039149 -0.158758356 -0.01294888
#>  [5,] -0.173947065 0.2617689 -0.6751999 -0.3449300 -0.175092438 -0.02170691
#>  [6,] -0.061175482 0.2622071 -0.5837536 -0.2162821 -0.070993556  0.11008269
#>  [7,] -0.075977635 0.2867725 -0.6024684 -0.2318131 -0.074291056  0.08331781
#>  [8,] -0.132448490 0.2380859 -0.6125521 -0.2772930 -0.135351510  0.01833532
#>  [9,] -0.117621307 0.2584321 -0.6900217 -0.2552826 -0.101700843  0.02961161
#> [10,] -0.147559029 0.2568479 -0.6868824 -0.2989574 -0.153037347  0.01048411
#> [11,] -0.028235912 0.2476901 -0.5247728 -0.1734267 -0.005740572  0.12745670
#> [12,] -0.100897704 0.2449501 -0.5713168 -0.2578322 -0.106142106  0.05178246
#> [13,] -0.097751434 0.2561717 -0.5511514 -0.2511608 -0.105182015  0.06272879
#> [14,] -0.127487597 0.2620980 -0.6508298 -0.2782384 -0.113433582  0.02056332
#> [15,] -0.135539424 0.2556825 -0.6620708 -0.2740940 -0.124402058  0.01608418
#> [16,] -0.023232345 0.2483834 -0.5197527 -0.1735213 -0.026653561  0.14197808
#> [17,] -0.039584071 0.2605575 -0.5735668 -0.1914476 -0.037859733  0.13168181
#> [18,] -0.057639112 0.2691562 -0.5760693 -0.2132433 -0.065759158  0.10168698
#> [19,] -0.067984470 0.2583847 -0.5710132 -0.2318370 -0.064363375  0.08142817
#> [20,] -0.115503093 0.2708791 -0.6448837 -0.2807997 -0.117078278  0.04842236
#> [21,]  0.009629076 0.2520274 -0.4976840 -0.1511609  0.025190790  0.17293649
#> [22,] -0.013778359 0.2374024 -0.4978035 -0.1587481 -0.015786390  0.12802566
#> [23,] -0.066611019 0.2390253 -0.5366574 -0.2237952 -0.057299230  0.09362090
#> [24,] -0.069179779 0.2721143 -0.6200773 -0.2396625 -0.074214109  0.09444577
#> [25,] -0.084587507 0.2521505 -0.5687310 -0.2374860 -0.088594252  0.08138683
#>           97.5%     Rhat n.eff
#>  [1,] 0.3615704 1.009198   210
#>  [2,] 0.4586432 1.006852   600
#>  [3,] 0.3884170 1.022334   170
#>  [4,] 0.3332049 1.010474   370
#>  [5,] 0.2832294 1.001417   600
#>  [6,] 0.4009676 1.015431   390
#>  [7,] 0.4105424 1.041476    80
#>  [8,] 0.3259568 1.040007    52
#>  [9,] 0.3328339 1.024181   150
#> [10,] 0.3224797 1.049614    61
#> [11,] 0.4140850 1.001252   600
#> [12,] 0.3614750 1.000299   600
#> [13,] 0.3823749 1.007284   470
#> [14,] 0.3693723 1.027623    72
#> [15,] 0.3580304 1.021584   180
#> [16,] 0.4405524 1.016846   350
#> [17,] 0.4990078 1.001661   600
#> [18,] 0.4708527 1.017649   130
#> [19,] 0.4441300 1.002221   600
#> [20,] 0.3673389 1.020872   190
#> [21,] 0.4436620 1.003259   600
#> [22,] 0.4507411 1.014672   240
#> [23,] 0.3599848 1.008879   280
#> [24,] 0.4662736 1.000479   600
#> [25,] 0.3776312 1.003403   410
#> 
#> $EM_pred
#>              mean        sd       2.5%        25%          50%        75%
#>  [1,] -0.10401285 0.3733000 -0.9046155 -0.3110445 -0.089701923 0.10895335
#>  [2,] -0.13303577 0.5198664 -1.0831250 -0.3516214 -0.124796139 0.08093764
#>  [3,] -0.17094184 0.4832055 -1.2268919 -0.3758341 -0.129666701 0.05072319
#>  [4,] -0.13530303 0.4338663 -0.9818621 -0.3169703 -0.140727304 0.06083527
#>  [5,] -0.19676043 0.4316295 -1.0877829 -0.3780687 -0.165875361 0.01247477
#>  [6,] -0.05583495 0.4254428 -0.8742355 -0.2553714 -0.067687740 0.17015312
#>  [7,] -0.08955995 0.4531946 -0.9661055 -0.3037336 -0.066995178 0.13543428
#>  [8,] -0.13937064 0.3799901 -0.9123243 -0.3145656 -0.120426161 0.05667444
#>  [9,] -0.14277894 0.4555477 -1.0257072 -0.3100469 -0.106326741 0.06157884
#> [10,] -0.16027558 0.4042255 -1.0410098 -0.3719204 -0.169094865 0.05827288
#> [11,] -0.02860951 0.4069202 -0.7678804 -0.2464405 -0.018143562 0.18649954
#> [12,] -0.11499760 0.4195468 -1.0313512 -0.3261744 -0.114661164 0.09495771
#> [13,] -0.06587474 0.4625878 -0.8137655 -0.3332975 -0.080286998 0.14827991
#> [14,] -0.13556110 0.4686706 -1.0501913 -0.3442697 -0.127565169 0.07376001
#> [15,] -0.13711391 0.4170586 -1.1031621 -0.3201304 -0.124164067 0.08255056
#> [16,] -0.02279512 0.3880477 -0.8594341 -0.2052677 -0.014935061 0.17230248
#> [17,] -0.04549223 0.4302355 -0.8834158 -0.2528167 -0.038816872 0.16970570
#> [18,] -0.06913113 0.4726580 -1.0248836 -0.2667871 -0.066444602 0.12836887
#> [19,] -0.05833253 0.4018164 -0.8814351 -0.2654280 -0.057963203 0.13712049
#> [20,] -0.12354569 0.4231540 -0.9865541 -0.3210973 -0.120102151 0.10604299
#> [21,]  0.02031559 0.3752779 -0.7422987 -0.1849218  0.009033804 0.24059589
#> [22,] -0.02977679 0.4210429 -0.7725389 -0.2072706 -0.019288948 0.17923785
#> [23,] -0.05899007 0.3962546 -0.8439239 -0.2742125 -0.050574915 0.12939704
#> [24,] -0.04897328 0.4098437 -0.8580556 -0.2694113 -0.055732817 0.15507599
#> [25,] -0.09585334 0.4251563 -0.9065438 -0.3076704 -0.094446836 0.11851809
#>           97.5%     Rhat n.eff
#>  [1,] 0.5989984 1.004212   380
#>  [2,] 0.6364095 1.084590   450
#>  [3,] 0.7091741 1.022615   420
#>  [4,] 0.8477200 1.012556   560
#>  [5,] 0.5214141 1.022888   600
#>  [6,] 0.7061002 1.017046   310
#>  [7,] 0.7932607 1.029567   160
#>  [8,] 0.6681611 1.030157    75
#>  [9,] 0.6009639 1.015619   600
#> [10,] 0.6287712 1.017966   110
#> [11,] 0.6801038 1.007697   460
#> [12,] 0.6976434 1.002784   470
#> [13,] 0.8613740 1.024311   600
#> [14,] 0.7054536 1.012191   150
#> [15,] 0.6112558 1.015406   290
#> [16,] 0.7004359 1.007324   340
#> [17,] 0.8318758 1.007208   580
#> [18,] 0.7098425 1.016090   490
#> [19,] 0.7581055 1.002322   600
#> [20,] 0.5656467 1.014618   160
#> [21,] 0.7659916 1.002846   600
#> [22,] 0.7351277 1.005583   600
#> [23,] 0.8375734 1.004431   420
#> [24,] 0.8340633 1.007240   570
#> [25,] 0.6522140 1.009818   500
#> 
#> $tau
#>            mean        sd        2.5%        25%       50%       75%     97.5%
#>  [1,] 0.2380954 0.2153363 0.013500658 0.08330434 0.1757870 0.3197826 0.7692328
#>  [2,] 0.2611240 0.2871277 0.020409472 0.09426574 0.1858401 0.3377040 0.8606076
#>  [3,] 0.2799059 0.2677395 0.019163284 0.10497913 0.2017856 0.3597663 1.0328539
#>  [4,] 0.2568161 0.2133799 0.011851635 0.11082470 0.2034187 0.3483924 0.7334966
#>  [5,] 0.2278381 0.2124609 0.016737509 0.08367280 0.1686403 0.3102389 0.7309463
#>  [6,] 0.2424450 0.2159332 0.016689730 0.10513529 0.1840408 0.3121328 0.8594661
#>  [7,] 0.2470998 0.2403946 0.007621982 0.08910924 0.1781861 0.3327522 0.8027401
#>  [8,] 0.2187730 0.1953566 0.007684777 0.08889238 0.1645052 0.2958507 0.7393782
#>  [9,] 0.2545337 0.2461588 0.010088629 0.08756754 0.1917536 0.3450815 0.8824341
#> [10,] 0.2470150 0.2160814 0.019031335 0.10732125 0.1940555 0.3104467 0.8702732
#> [11,] 0.2537160 0.2427514 0.012560056 0.09030012 0.1899642 0.3378698 0.9154208
#> [12,] 0.2436171 0.2209712 0.022775443 0.10123559 0.1868818 0.3226801 0.7983610
#> [13,] 0.2535471 0.2195511 0.021794880 0.10648953 0.2023536 0.3372929 0.7729696
#> [14,] 0.2516386 0.2353316 0.022690202 0.09882553 0.1896556 0.3361160 0.8278557
#> [15,] 0.2427511 0.2053660 0.020906521 0.10009628 0.1817961 0.3257088 0.8080272
#> [16,] 0.2370021 0.2114464 0.014949896 0.09158281 0.1774504 0.3188572 0.8118558
#> [17,] 0.2493036 0.2398025 0.021958679 0.09546764 0.1861069 0.3217919 0.8072718
#> [18,] 0.2494615 0.2556767 0.007675567 0.08725075 0.1764193 0.3293491 0.8900274
#> [19,] 0.2389104 0.2159729 0.013128587 0.08630301 0.1801127 0.3244172 0.7958209
#> [20,] 0.2426813 0.2546914 0.005535893 0.08563427 0.1827922 0.3165899 0.7863569
#> [21,] 0.2364303 0.2243880 0.011131865 0.09055316 0.1790937 0.3116436 0.7670977
#> [22,] 0.2188739 0.2490662 0.010269351 0.07643701 0.1693234 0.2901998 0.6910289
#> [23,] 0.2368482 0.2057956 0.022285154 0.09028180 0.1793596 0.3152875 0.7674719
#> [24,] 0.2414290 0.2059855 0.013362208 0.09341286 0.1917661 0.3357936 0.7509146
#> [25,] 0.2248440 0.2127396 0.006184910 0.07549843 0.1664587 0.2945672 0.9150064
#>           Rhat n.eff
#>  [1,] 1.023887   600
#>  [2,] 1.094453    31
#>  [3,] 1.008864   600
#>  [4,] 1.002374   520
#>  [5,] 1.008517   600
#>  [6,] 1.012482   160
#>  [7,] 1.050396   110
#>  [8,] 1.076901    47
#>  [9,] 1.028398   260
#> [10,] 1.019464   130
#> [11,] 1.012460   410
#> [12,] 1.012384   570
#> [13,] 1.000525   600
#> [14,] 1.007930   280
#> [15,] 1.023112   110
#> [16,] 1.015944   120
#> [17,] 1.005138   300
#> [18,] 1.032400    93
#> [19,] 1.016889   600
#> [20,] 1.006836   240
#> [21,] 1.021354   140
#> [22,] 1.009349   600
#> [23,] 1.011602   600
#> [24,] 1.012395   600
#> [25,] 1.007432   530
#> 
#> $measure
#> [1] "SMD"
#> 
#> $model
#> [1] "RE"
#> 
#> $scenarios
#> [1] -2 -1  0  1  2
#> 
#> $D
#> [1] 0
#> 
#> $heter
#> [1] -2.9900000  0.2143347  4.0000000
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
#> [1] 5
#> 
#> attr(,"class")
#> [1] "run_sensitivity"
# }
```
