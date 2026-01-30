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
#>                mean        sd       2.5%        25%          50%          75%
#>  [1,] -0.0973579904 0.2546467 -0.5939338 -0.2460538 -0.087936317  0.056160733
#>  [2,] -0.1052182942 0.2709751 -0.6070855 -0.2699012 -0.104952534  0.079010065
#>  [3,] -0.1343641942 0.2614908 -0.6600825 -0.2915702 -0.136110575  0.034959003
#>  [4,] -0.1657205011 0.2561071 -0.6187473 -0.3262107 -0.169518086 -0.021256345
#>  [5,] -0.2204974571 0.2691333 -0.7324628 -0.3609978 -0.225233487 -0.071885432
#>  [6,] -0.0800973678 0.2343642 -0.5571773 -0.2176328 -0.075405563  0.069969553
#>  [7,] -0.1279648161 0.2485966 -0.6154290 -0.3061032 -0.123339395  0.037209110
#>  [8,] -0.1021070784 0.2696719 -0.6311652 -0.2637169 -0.106838169  0.070584562
#>  [9,] -0.1161461393 0.2577142 -0.6236957 -0.2821637 -0.120158927  0.048999432
#> [10,] -0.1723745071 0.2646703 -0.7168304 -0.3205019 -0.161367927 -0.007390435
#> [11,] -0.0637500264 0.2746212 -0.6316931 -0.2350423 -0.035312019  0.110106111
#> [12,] -0.0876210328 0.2540161 -0.6431850 -0.2451833 -0.079866769  0.076607172
#> [13,] -0.0687171719 0.2691161 -0.5742566 -0.2279099 -0.063535429  0.096921278
#> [14,] -0.1228397590 0.2424198 -0.5763775 -0.2921612 -0.120403535  0.026316535
#> [15,] -0.1350450194 0.2616705 -0.5760420 -0.3032858 -0.142338865  0.023951428
#> [16,] -0.0169691743 0.2415629 -0.4698544 -0.1519638 -0.011075738  0.131749002
#> [17,] -0.0178099804 0.2356225 -0.4712291 -0.1698828 -0.015790089  0.131907241
#> [18,] -0.0759187476 0.2485094 -0.5460830 -0.2370044 -0.073450765  0.087781443
#> [19,] -0.0921415141 0.2449487 -0.5276367 -0.2501207 -0.097311545  0.071044373
#> [20,] -0.1052836203 0.2634321 -0.6603374 -0.2677533 -0.096732394  0.061267966
#> [21,] -0.0006551035 0.2437898 -0.4763413 -0.1529490 -0.006242237  0.147177937
#> [22,] -0.0250827329 0.2622401 -0.5046998 -0.1897100 -0.021104401  0.136804048
#> [23,] -0.0704940502 0.2547580 -0.5292606 -0.2204118 -0.065884568  0.073562989
#> [24,] -0.0550279038 0.2572022 -0.5439474 -0.2122852 -0.037758890  0.084398225
#> [25,] -0.0757789612 0.2505398 -0.5704946 -0.2326497 -0.066314481  0.089536812
#>           97.5%      Rhat n.eff
#>  [1,] 0.3868343 1.0024649   550
#>  [2,] 0.3394279 0.9994329   600
#>  [3,] 0.3559661 1.0039836   480
#>  [4,] 0.3761893 1.0027441   520
#>  [5,] 0.3366788 1.0168359   390
#>  [6,] 0.3404073 1.0015798   600
#>  [7,] 0.3766154 1.0045460   440
#>  [8,] 0.4510190 1.0078434   600
#>  [9,] 0.3479595 1.0108278   220
#> [10,] 0.3241280 1.0059342   300
#> [11,] 0.3953833 1.0022962   600
#> [12,] 0.3461018 1.0071103   310
#> [13,] 0.4629967 1.0085406   240
#> [14,] 0.3501632 1.0200873   220
#> [15,] 0.3797150 1.0036896   430
#> [16,] 0.4591470 1.0138180   190
#> [17,] 0.4184134 1.0353663    64
#> [18,] 0.3924223 1.0390267    55
#> [19,] 0.3718940 1.0032038   550
#> [20,] 0.4164982 1.0114614   600
#> [21,] 0.4629188 1.0010881   600
#> [22,] 0.5061437 1.0161961   600
#> [23,] 0.4197884 1.0027136   480
#> [24,] 0.4176639 1.0031185   600
#> [25,] 0.4124910 1.0198472   100
#> 
#> $EM_pred
#>               mean        sd       2.5%        25%          50%         75%
#>  [1,] -0.101626836 0.4213674 -0.8996424 -0.3106107 -0.100668746  0.11100649
#>  [2,] -0.112106494 0.4517032 -1.2101535 -0.3180577 -0.086350329  0.13844520
#>  [3,] -0.145785455 0.3803537 -0.9351691 -0.3580299 -0.116166807  0.08256562
#>  [4,] -0.131937323 0.4270649 -0.8692487 -0.3614889 -0.168881763  0.05324081
#>  [5,] -0.218241072 0.4305724 -1.1841876 -0.3919359 -0.224621524 -0.03633561
#>  [6,] -0.100872732 0.3688023 -0.8319359 -0.2654218 -0.084333169  0.09428305
#>  [7,] -0.121765540 0.4032927 -0.9312353 -0.3268987 -0.107782031  0.08772561
#>  [8,] -0.099087169 0.4250696 -0.9808517 -0.3163975 -0.122574177  0.12488952
#>  [9,] -0.101162883 0.4128599 -0.9003753 -0.3193588 -0.127145431  0.10713890
#> [10,] -0.170918784 0.4276518 -1.0129692 -0.3667567 -0.165598728  0.03939704
#> [11,] -0.074258432 0.4133049 -1.0372877 -0.2725687 -0.054941427  0.14624424
#> [12,] -0.115680027 0.3879262 -0.8894224 -0.3447570 -0.081479708  0.10626107
#> [13,] -0.074447775 0.4036275 -0.8917092 -0.2863055 -0.071520035  0.13384528
#> [14,] -0.130973243 0.3968200 -0.9008704 -0.3524987 -0.125904431  0.06917630
#> [15,] -0.130893097 0.4491334 -0.9393694 -0.3357972 -0.140732646  0.05272187
#> [16,] -0.002222782 0.4031206 -0.8856529 -0.1758811 -0.008929197  0.16371005
#> [17,] -0.027441936 0.3463463 -0.7585571 -0.2243715 -0.022440599  0.16817426
#> [18,] -0.056336235 0.4099155 -0.8465791 -0.2735734 -0.074446055  0.15498026
#> [19,] -0.083342391 0.3751447 -0.8492176 -0.2888712 -0.087184078  0.11555434
#> [20,] -0.117444129 0.4119419 -0.9546386 -0.3120880 -0.103609266  0.08830900
#> [21,] -0.026798196 0.3941460 -0.7866214 -0.2001650 -0.026027361  0.16089360
#> [22,] -0.030395677 0.4509904 -0.8332859 -0.2424463 -0.028277433  0.18152012
#> [23,] -0.043883551 0.4408230 -0.8244007 -0.2732453 -0.064048357  0.14199114
#> [24,] -0.083007638 0.4451453 -1.0416071 -0.2734333 -0.051052212  0.15565662
#> [25,] -0.085364666 0.3840182 -0.9114035 -0.2862948 -0.081993788  0.13792788
#>           97.5%      Rhat n.eff
#>  [1,] 0.7637538 1.0070666   600
#>  [2,] 0.7174865 1.0073224   600
#>  [3,] 0.4781399 1.0085464   240
#>  [4,] 0.9145718 1.0032607   550
#>  [5,] 0.7346174 1.0025522   490
#>  [6,] 0.5472201 1.0073898   600
#>  [7,] 0.6335772 1.0023639   600
#>  [8,] 0.7984146 1.0015773   600
#>  [9,] 0.7352443 1.0190806   170
#> [10,] 0.5719676 1.0158937   360
#> [11,] 0.8341144 1.0087916   600
#> [12,] 0.5907775 1.0060247   600
#> [13,] 0.8153435 1.0071597   230
#> [14,] 0.6438688 1.0157898   210
#> [15,] 0.8133746 1.0110221   600
#> [16,] 0.8725114 1.0077834   600
#> [17,] 0.6378544 1.0214202    94
#> [18,] 0.7610488 1.0181256   180
#> [19,] 0.6628234 0.9986588   600
#> [20,] 0.7424318 1.0083366   600
#> [21,] 0.6501878 1.0058007   560
#> [22,] 0.7490570 1.0145364   600
#> [23,] 0.8813103 1.0151707   600
#> [24,] 0.6223160 1.0096009   600
#> [25,] 0.6273587 1.0200288   190
#> 
#> $tau
#>            mean        sd        2.5%        25%       50%       75%     97.5%
#>  [1,] 0.2470356 0.2345611 0.022094470 0.09778269 0.1781678 0.3218848 0.8776256
#>  [2,] 0.2337189 0.2030537 0.019106951 0.08830734 0.1765247 0.3119168 0.7762553
#>  [3,] 0.2434104 0.2202254 0.015735372 0.09151982 0.1884932 0.3125342 0.8524733
#>  [4,] 0.2505638 0.2234751 0.018382233 0.09490715 0.1812520 0.3391095 0.8224866
#>  [5,] 0.2591985 0.2670077 0.013982614 0.08958559 0.1874765 0.3353690 1.0195637
#>  [6,] 0.2339714 0.2081216 0.009531287 0.09126980 0.1773621 0.3184116 0.7869249
#>  [7,] 0.2415613 0.2067599 0.018630893 0.09702733 0.1787244 0.3247771 0.7845238
#>  [8,] 0.2489288 0.2181019 0.019704795 0.09159432 0.1854637 0.3356410 0.8014655
#>  [9,] 0.2463036 0.2111778 0.015475924 0.09766240 0.1950898 0.3434188 0.7289854
#> [10,] 0.2496112 0.2293513 0.024920175 0.10274427 0.1946577 0.3227955 0.9288785
#> [11,] 0.2366712 0.2299573 0.005635760 0.09671367 0.1867084 0.3080924 0.8070852
#> [12,] 0.2470876 0.2001686 0.014251883 0.09953360 0.1986750 0.3331938 0.7749875
#> [13,] 0.2414738 0.2496836 0.016396450 0.09884915 0.1720365 0.3121035 0.8228996
#> [14,] 0.2515945 0.2122072 0.032006074 0.11923016 0.1912771 0.3254277 0.7692744
#> [15,] 0.2478150 0.2521699 0.009495180 0.08559089 0.1715397 0.3335054 0.9497043
#> [16,] 0.2279020 0.2227109 0.011076895 0.07833710 0.1703105 0.2968139 0.7736991
#> [17,] 0.2186878 0.1798216 0.015722574 0.08858736 0.1715917 0.2946677 0.6831957
#> [18,] 0.2462919 0.2179819 0.022018152 0.09783550 0.1889458 0.3266151 0.7929396
#> [19,] 0.2323055 0.1772659 0.016357267 0.10706460 0.1960828 0.3096833 0.6972008
#> [20,] 0.2394037 0.2141484 0.013867392 0.08309020 0.1735556 0.3362784 0.7932552
#> [21,] 0.2266067 0.1911928 0.025645312 0.08977633 0.1770006 0.2991644 0.7055894
#> [22,] 0.2433070 0.2145409 0.010484699 0.10576161 0.1818798 0.3141801 0.8453104
#> [23,] 0.2407115 0.2307713 0.006050818 0.08484236 0.1753980 0.3203851 0.8335994
#> [24,] 0.2212248 0.2123137 0.009085150 0.08080785 0.1715390 0.2908365 0.7535461
#> [25,] 0.2227302 0.2021759 0.015155842 0.08481658 0.1738648 0.2908239 0.7161720
#>           Rhat n.eff
#>  [1,] 1.001556   600
#>  [2,] 1.019169   210
#>  [3,] 1.009175   600
#>  [4,] 1.008148   540
#>  [5,] 1.040063   600
#>  [6,] 1.019445   180
#>  [7,] 1.008807   420
#>  [8,] 1.006000   600
#>  [9,] 1.009632   180
#> [10,] 1.019628   110
#> [11,] 1.019962   130
#> [12,] 1.006182   600
#> [13,] 1.006822   600
#> [14,] 1.007817   220
#> [15,] 1.031809   190
#> [16,] 1.031903   180
#> [17,] 1.017271   130
#> [18,] 1.012093   170
#> [19,] 1.058366    53
#> [20,] 1.013292   230
#> [21,] 1.010476   260
#> [22,] 1.086491    53
#> [23,] 1.048939    94
#> [24,] 1.065189    70
#> [25,] 1.004650   480
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
