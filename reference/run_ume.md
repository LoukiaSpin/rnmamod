# Perform the unrelated mean effects model

Performs the unrelated mean effects model of Dias et al. (2013) that has
been refined (Spineli, 2021) and extended to address aggregate binary
and continuous missing participant outcome data via the pattern-mixture
model (Spineli et al. 2021; Spineli, 2019). This model offers a global
evaluation of the plausibility of the consistency assumption in the
network.

## Usage

``` r
run_ume(full, n_iter, n_burnin, n_chains, n_thin, inits = NULL)
```

## Arguments

- full:

  An object of S3 class
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  See 'Value' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

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

- n_chains:

  Positive integer specifying the number of chains for the MCMC
  sampling; an argument of the
  [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) function of the
  R-package [R2jags](https://CRAN.R-project.org/package=R2jags). The
  default argument is 2.

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
  for each pairwise comparison observed in the network.

- dev_o:

  The deviance contribution of each trial-arm based on the observed
  outcome.

- hat_par:

  The fitted outcome at each trial-arm.

- tau:

  The between-trial standard deviation (assumed common across the
  observed pairwise comparisons) for the whole network, when a
  random-effects model has been specified.

- m_tau:

  The between-trial standard deviation (assumed common across the
  observed pairwise comparisons) for the subset of multi-arm trials,
  when a random-effects model has been specified.

The output also includes the following elements:

- leverage_o:

  The leverage for the observed outcome at each trial-arm.

- sign_dev_o:

  The sign of the difference between observed and fitted outcome at each
  trial-arm.

- model_assessment:

  A data-frame on the measures of model assessment: deviance information
  criterion, number of effective parameters, and total residual
  deviance.

- jagsfit:

  An object of S3 class
  [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) with the posterior
  results on all monitored parameters to be used in the
  [`mcmc_diagnostics`](https://loukiaspin.github.io/rnmamod/reference/mcmc_diagnostics.md)
  function.

Furthermore, `run_ume` returns a character vector with the pairwise
comparisons observed in the network, `obs_comp`, and a character vector
with comparisons between the non-baseline interventions observed in
multi-arm trials only, `frail_comp`. Both vectors are used in
[`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md)
function.

## Details

`run_ume` inherits the arguments `data`, `measure`, `model`,
`assumption`, `heter_prior`, `mean_misspar`, `var_misspar`, and `ref`
from
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
This prevents specifying a different Bayesian model from that considered
in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).Therefore,
the user needs first to apply
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
and then use `run_ume` (see 'Examples').

The `run_ume` function also returns the arguments `data`, `model`,
`measure`, `assumption`, `n_chains`, `n_iter`, `n_burnin`, and `n_thin`
as specified by the user to be inherited by other relevant functions of
the package.

Initially, `run_ume` calls the
[`improved_ume`](https://loukiaspin.github.io/rnmamod/reference/improved_ume.md)
function to identify the *frail comparisons*, that is, comparisons
between non-baseline interventions in multi-arm trials not investigated
in any two-arm or multi-arm trial of the network (Spineli, 2021). The
'original' model of Dias et al. (2013) omits the frail comparisons from
the estimation process. Consequently, the number of estimated summary
effects is less than those obtained by performing separate pairwise
meta-analyses (see
[`run_series_meta`](https://loukiaspin.github.io/rnmamod/reference/run_series_meta.md)).

For a binary outcome, when `measure` is "RR" (relative risk) or "RD"
(risk difference) in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
`run_ume` currently considers the odds ratio as effect measure for being
the **base-case** effect measure in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
for a binary outcome (see also 'Details' in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).

`run_ume` calls the
[`prepare_ume`](https://loukiaspin.github.io/rnmamod/reference/prepare_ume.md)
function which contains the WinBUGS code as written by Dias et al.
(2013) for binomial and normal likelihood to analyse binary and
continuous outcome data, respectively.
[`prepare_ume`](https://loukiaspin.github.io/rnmamod/reference/prepare_ume.md)
has been extended to incorporate the pattern-mixture model with
informative missingness parameters for binary and continuous outcome
data (see 'Details' in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).
[`prepare_ume`](https://loukiaspin.github.io/rnmamod/reference/prepare_ume.md)
has also been refined to account for the multi-arm trials by assigning
conditional univariate normal distributions on the underlying
trial-specific effect size of comparisons with the baseline arm of the
multi-arm trial (Spineli, 2021).

`run_ume` runs Bayesian unrelated mean effects model in `JAGS`. The
progress of the simulation appears on the R console. The model is
updated until convergence using the
[`autojags`](https://rdrr.io/pkg/R2jags/man/autojags.html) function of
the R-package [R2jags](https://CRAN.R-project.org/package=R2jags) with 2
updates and number of iterations and thinning equal to `n_iter` and
`n_thin`, respectively.

The output of `run_ume` is not end-user-ready. The
[`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md)
function uses the output of `run_ume` as an S3 object and processes it
further to provide an end-user-ready output.

`run_ume` can be used only for a network of interventions. In the case
of two interventions, the execution of the function will be stopped and
an error message will be printed on the R console.

## References

Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence
synthesis for decision making 4: inconsistency in networks of evidence
based on randomized controlled trials. *Med Decis Making*
2013;**33**(5):641–56. doi: 10.1177/0272989X12455847

Gelman A, Rubin DB. Inference from iterative simulation using multiple
sequences. *Stat Sci* 1992;**7**(4):457–72. doi: 10.1214/ss/1177011136

Spineli LM. A Revised Framework to Evaluate the Consistency Assumption
Globally in a Network of Interventions. *Med Decis Making* 2021. doi:
10.1177/0272989X211068005

Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing
outcome data in network meta-analysis: a one-stage pattern-mixture model
approach. *Stat Methods Med Res* 2021;**30**(4):958–75. doi:
10.1177/0962280220983544

Spineli LM. An empirical comparison of Bayesian modelling strategies for
missing binary outcome data in network meta-analysis. *BMC Med Res
Methodol* 2019;**19**(1):86. doi: 10.1186/s12874-019-0731-y

## See also

[`autojags`](https://rdrr.io/pkg/R2jags/man/autojags.html),
[`jags`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`prepare_ume`](https://loukiaspin.github.io/rnmamod/reference/prepare_ume.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`run_series_meta`](https://loukiaspin.github.io/rnmamod/reference/run_series_meta.md),
[`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md)

## Author

Loukia M. Spineli

## Examples

``` r
data("nma.liu2013")

# Read results from 'run_model' (using the default arguments)
res <- readRDS(system.file('extdata/res_liu.rds', package = 'rnmamod'))

# \donttest{
# Run random-effects unrelated mean effects model
# Note: Ideally, set 'n_iter' to 10000 and 'n_burnin' to 1000
run_ume(full = res,
        n_chains = 3,
        n_iter = 1000,
        n_burnin = 100,
        n_thin = 1)
#> JAGS generates initial values for the parameters.
#> Running the model ...
#> Compiling model graph
#>    Resolving undeclared variables
#>    Allocating nodes
#> Graph information:
#>    Observed stochastic nodes: 59
#>    Unobserved stochastic nodes: 80
#>    Total graph size: 1291
#> 
#> Initializing model
#> 
#> ... Updating the model until convergence
#> $EM
#>               mean        sd       2.5%        25%        50%        75%
#> EM[2,1]  0.6222651 0.8692068 -1.2826497  0.1592340  0.6701062  1.0774445
#> EM[3,1]  0.1785900 0.7302484 -1.2391343 -0.2945278  0.1923722  0.6535423
#> EM[4,1]  0.7048289 0.6239851 -0.6803186  0.3323913  0.7588838  1.1252183
#> EM[5,1]  1.7565404 0.8503496  0.2491505  1.1558146  1.7789509  2.3079050
#> EM[4,2] -0.4242400 0.9904817 -2.4814668 -0.9998214 -0.3994952  0.1833293
#> EM[6,2] -1.4649326 1.2264658 -3.9236798 -2.2194803 -1.4138484 -0.6909996
#> EM[4,3]  1.1965780 0.9856488 -0.8578727  0.6203434  1.2428284  1.7883028
#> EM[5,4] -1.1398820 1.3441896 -3.4785908 -2.0416571 -1.2738581 -0.3243877
#>             97.5%     Rhat n.eff
#> EM[2,1] 2.3419981 1.004551  1700
#> EM[3,1] 1.6456126 1.024294   130
#> EM[4,1] 1.8593559 1.107807    25
#> EM[5,1] 3.4184985 1.015384   230
#> EM[4,2] 1.4938821 1.003308  1900
#> EM[6,2] 0.9271331 1.008605  1500
#> EM[4,3] 3.0165227 1.009325   340
#> EM[5,4] 1.7984066 1.086650    38
#> 
#> $dev_o
#>                    mean         sd         2.5%        25%       50%
#> dev.o[1,1]  0.976003854 1.39607393 0.0011460044 0.10198725 0.4306728
#> dev.o[2,1]  0.871112577 1.28244294 0.0005982561 0.08359298 0.3722219
#> dev.o[3,1]  1.012839717 1.28769021 0.0010607370 0.11491937 0.5149106
#> dev.o[4,1]  1.006625466 1.39261857 0.0010794329 0.09353413 0.4673917
#> dev.o[5,1]  1.084174616 1.55459406 0.0010117277 0.10448676 0.4899437
#> dev.o[6,1]  1.067250743 1.41538097 0.0012672402 0.11100549 0.5209381
#> dev.o[7,1]  0.982751040 1.37240725 0.0014792830 0.10016528 0.4581158
#> dev.o[8,1]  0.004754167 0.05801256 0.0000000000 0.00000000 0.0000000
#> dev.o[9,1]  0.884455648 1.21515789 0.0012114430 0.08883597 0.4165631
#> dev.o[10,1] 0.926157225 1.27169874 0.0012500972 0.09447921 0.4457951
#> dev.o[11,1] 0.875404437 1.23549165 0.0009995200 0.09264021 0.4100028
#> dev.o[1,2]  0.976626458 1.39151757 0.0008428999 0.10140718 0.4577621
#> dev.o[2,2]  1.121945740 1.50437490 0.0011486483 0.12069271 0.5090218
#> dev.o[3,2]  1.035483482 1.53764933 0.0009423255 0.10341816 0.4734266
#> dev.o[4,2]  0.986717742 1.32939164 0.0009795902 0.11093454 0.4675963
#> dev.o[5,2]  1.239220480 1.73317596 0.0011687460 0.12158472 0.5352552
#> dev.o[6,2]  1.184696892 1.51161206 0.0013559866 0.14204325 0.6234204
#> dev.o[7,2]  0.840302115 1.17276054 0.0006249377 0.08772113 0.3860450
#> dev.o[8,2]  0.008773052 0.09801368 0.0000000000 0.00000000 0.0000000
#> dev.o[9,2]  1.022415547 1.43820852 0.0013067480 0.11263849 0.4595847
#> dev.o[10,2] 1.346960286 1.67058401 0.0015720712 0.15724162 0.7103110
#> dev.o[11,2] 1.003767779 1.36913040 0.0009935364 0.09850465 0.4688047
#> dev.o[9,3]  1.239462278 1.60772213 0.0012578733 0.12642030 0.6210999
#> dev.o[10,3] 1.024712072 1.38329319 0.0009337266 0.09562252 0.4951711
#> dev.o[11,3] 0.925211777 1.21735756 0.0004953771 0.09541532 0.4615575
#>                      75%      97.5%     Rhat n.eff
#> dev.o[1,1]  1.262698e+00 5.04196552 1.002145  1600
#> dev.o[2,1]  1.131356e+00 4.69754548 1.001454  2100
#> dev.o[3,1]  1.457124e+00 4.61884275 1.004786   900
#> dev.o[4,1]  1.340516e+00 5.16874191 1.000808  3000
#> dev.o[5,1]  1.411297e+00 5.34265663 1.002291  1100
#> dev.o[6,1]  1.444541e+00 5.27409078 1.002632   940
#> dev.o[7,1]  1.315281e+00 4.98523872 1.002468  1000
#> dev.o[8,1]  1.065814e-14 0.01508521 1.135482  1800
#> dev.o[9,1]  1.181916e+00 4.31036097 1.004013  1400
#> dev.o[10,1] 1.287612e+00 4.36733772 1.005622   390
#> dev.o[11,1] 1.144427e+00 4.41962804 1.006728   320
#> dev.o[1,2]  1.267749e+00 5.04021427 1.000523  3000
#> dev.o[2,2]  1.567881e+00 5.48292482 1.001445  2100
#> dev.o[3,2]  1.352937e+00 5.30072320 1.007514   290
#> dev.o[4,2]  1.329288e+00 4.73407766 1.001334  2400
#> dev.o[5,2]  1.658457e+00 6.29615277 1.003127   760
#> dev.o[6,2]  1.697070e+00 5.42660073 1.005507   410
#> dev.o[7,2]  1.107023e+00 4.12548514 1.000854  3000
#> dev.o[8,2]  2.664535e-14 0.04779659 1.153474  3000
#> dev.o[9,2]  1.340815e+00 5.34468338 1.000770  3000
#> dev.o[10,2] 1.929184e+00 5.87402768 1.014574   200
#> dev.o[11,2] 1.332556e+00 5.03517099 1.002733   900
#> dev.o[9,3]  1.718315e+00 5.71985767 1.000656  3000
#> dev.o[10,3] 1.398576e+00 4.78993128 1.003662   630
#> dev.o[11,3] 1.282196e+00 4.30264923 1.016013   130
#> 
#> $hat_par
#>                    mean         sd       2.5%       25%       50%       75%
#> hat.par[1,1]  26.662334 4.49353591 18.2276579 23.554113 26.466278 29.557409
#> hat.par[2,1]   2.743415 1.34821374  0.7710154  1.755521  2.544492  3.509565
#> hat.par[3,1]  10.139600 1.17655407  7.2667402  9.498921 10.386991 11.036401
#> hat.par[4,1]  22.865137 2.51344298 17.5015652 21.199076 23.051587 24.662264
#> hat.par[5,1]   8.145112 2.15501404  4.1355217  6.579734  8.110985  9.584765
#> hat.par[6,1]   3.294189 0.97510995  1.3617792  2.616085  3.306796  3.996286
#> hat.par[7,1]   2.619953 1.30884848  0.6675262  1.676552  2.430515  3.360208
#> hat.par[8,1]  11.997657 0.02804988 11.9924598 12.000000 12.000000 12.000000
#> hat.par[9,1]  15.647277 2.57537212 10.4975275 13.887042 15.655014 17.465822
#> hat.par[10,1]  3.477230 1.27419879  1.3563752  2.525444  3.397467  4.295578
#> hat.par[11,1]  4.768414 1.53313256  2.0334519  3.672171  4.698947  5.737874
#> hat.par[1,2]  38.247122 5.05649690 28.7869946 34.705492 38.134221 41.616114
#> hat.par[2,2]   4.251885 1.65638865  1.5439555  3.033635  4.088644  5.286088
#> hat.par[3,2]   7.889164 1.39733776  4.9964309  6.927980  7.964043  8.941345
#> hat.par[4,2]  16.222297 2.45068976 11.2254307 14.544803 16.322644 17.934953
#> hat.par[5,2]   2.803788 1.53859782  0.5239361  1.599391  2.608436  3.729748
#> hat.par[6,2]   3.734348 0.93602908  1.7213532  3.117686  3.813339  4.449655
#> hat.par[7,2]   2.340276 1.20160731  0.4621918  1.430560  2.196674  3.101845
#> hat.par[8,2]  14.995692 0.04690426 14.9761207 15.000000 15.000000 15.000000
#> hat.par[9,2]  16.909992 2.75074591 11.7351359 14.917093 16.897924 18.852387
#> hat.par[10,2]  3.126981 1.34051348  0.8040481  2.140357  3.002424  4.021114
#> hat.par[11,2]  7.033723 1.47500240  4.0627918  6.031211  7.051608  8.086801
#> hat.par[9,3]  21.354232 2.26556712 16.7861116 19.836315 21.472542 22.994460
#> hat.par[10,3]  8.390388 1.47518089  5.3610433  7.397323  8.558933  9.476101
#> hat.par[11,3] 11.244381 1.65432490  7.8045619 10.161343 11.391136 12.489974
#>                   97.5%     Rhat n.eff
#> hat.par[1,1]  35.947908 1.004134  1100
#> hat.par[2,1]   6.048701 1.001430  2100
#> hat.par[3,1]  11.656457 1.065176    58
#> hat.par[4,1]  27.245950 1.002279  1100
#> hat.par[5,1]  12.567582 1.002344  1600
#> hat.par[6,1]   5.074586 1.006880   360
#> hat.par[7,1]   5.605217 1.044541    50
#> hat.par[8,1]  12.000000 1.135482  1800
#> hat.par[9,1]  20.616886 1.011655   270
#> hat.par[10,1]  6.174334 1.035384    64
#> hat.par[11,1]  8.011958 1.018654   490
#> hat.par[1,2]  48.513713 1.000854  3000
#> hat.par[2,2]   7.945940 1.000543  3000
#> hat.par[3,2]  10.307465 1.012666   190
#> hat.par[4,2]  20.784389 1.002847   850
#> hat.par[5,2]   6.294247 1.002597  2700
#> hat.par[6,2]   5.306243 1.030664   120
#> hat.par[7,2]   5.022927 1.080634    36
#> hat.par[8,2]  15.000000 1.153474  3000
#> hat.par[9,2]  22.463739 1.004174   600
#> hat.par[10,2]  5.892114 1.119005    27
#> hat.par[11,2]  9.782071 1.000729  3000
#> hat.par[9,3]  25.407572 1.004801   470
#> hat.par[10,3] 10.876516 1.017648   180
#> hat.par[11,3] 14.012652 1.012897   480
#> 
#> $leverage_o
#>  [1]  0.970672037  0.618638263  1.000686208  1.003830639  1.079642517
#>  [6]  0.723193503  0.919975696 -0.008962488  0.833453335  0.814030242
#> [11]  0.692919832  0.974314842  0.965433848  1.029936803  0.979103745
#> [16]  1.223074410  0.813964813  0.777488574  0.003313173  0.924926052
#> [21]  0.726549865  0.711013615  0.780462644  0.872159425  0.907513141
#> 
#> $sign_dev_o
#>  [1]  1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1
#> 
#> $tau
#>        mean          sd        2.5%         25%         50%         75% 
#>  0.68984565  0.41890551  0.04736705  0.36601168  0.62758949  1.00237176 
#>       97.5%        Rhat       n.eff 
#>  1.50714572  1.03188264 90.00000000 
#> 
#> $m_tau
#>        mean          sd        2.5%         25%         50%         75% 
#>  0.72462102  0.49998405  0.04744372  0.30360757  0.66055672  1.03020817 
#>       97.5%        Rhat       n.eff 
#>  1.94814690  1.09200846 36.00000000 
#> 
#> $model_assessment
#>        DIC       pD      dev
#> 1 43.98516 20.33733 23.64783
#> 
#> $obs_comp
#> [1] "2vs1" "3vs1" "4vs1" "5vs1" "4vs2" "6vs2" "4vs3" "5vs4"
#> 
#> $jagsfit
#> Inference for Bugs model at "14", fit using jags,
#>  3 chains, each with 1000 iterations (first 0 discarded)
#>  n.sims = 3000 iterations saved. Running time = secs
#>               mu.vect sd.vect   2.5%    25%    50%    75%  97.5%  Rhat n.eff
#> EM[2,1]         0.622   0.869 -1.283  0.159  0.670  1.077  2.342 1.005  1700
#> EM[3,1]         0.179   0.730 -1.239 -0.295  0.192  0.654  1.646 1.024   130
#> EM[4,1]         0.705   0.624 -0.680  0.332  0.759  1.125  1.859 1.108    25
#> EM[5,1]         1.757   0.850  0.249  1.156  1.779  2.308  3.418 1.015   230
#> EM[4,2]        -0.424   0.990 -2.481 -1.000 -0.399  0.183  1.494 1.003  1900
#> EM[6,2]        -1.465   1.226 -3.924 -2.219 -1.414 -0.691  0.927 1.009  1500
#> EM[4,3]         1.197   0.986 -0.858  0.620  1.243  1.788  3.017 1.009   340
#> EM[5,4]        -1.140   1.344 -3.479 -2.042 -1.274 -0.324  1.798 1.087    38
#> dev.o[1,1]      0.976   1.396  0.001  0.102  0.431  1.263  5.042 1.002  1600
#> dev.o[2,1]      0.871   1.282  0.001  0.084  0.372  1.131  4.698 1.001  2100
#> dev.o[3,1]      1.013   1.288  0.001  0.115  0.515  1.457  4.619 1.005   900
#> dev.o[4,1]      1.007   1.393  0.001  0.094  0.467  1.341  5.169 1.001  3000
#> dev.o[5,1]      1.084   1.555  0.001  0.104  0.490  1.411  5.343 1.002  1100
#> dev.o[6,1]      1.067   1.415  0.001  0.111  0.521  1.445  5.274 1.003   940
#> dev.o[7,1]      0.983   1.372  0.001  0.100  0.458  1.315  4.985 1.002  1000
#> dev.o[8,1]      0.005   0.058  0.000  0.000  0.000  0.000  0.015 1.135  1800
#> dev.o[9,1]      0.884   1.215  0.001  0.089  0.417  1.182  4.310 1.004  1400
#> dev.o[10,1]     0.926   1.272  0.001  0.094  0.446  1.288  4.367 1.006   390
#> dev.o[11,1]     0.875   1.235  0.001  0.093  0.410  1.144  4.420 1.007   320
#> dev.o[1,2]      0.977   1.392  0.001  0.101  0.458  1.268  5.040 1.001  3000
#> dev.o[2,2]      1.122   1.504  0.001  0.121  0.509  1.568  5.483 1.001  2100
#> dev.o[3,2]      1.035   1.538  0.001  0.103  0.473  1.353  5.301 1.008   290
#> dev.o[4,2]      0.987   1.329  0.001  0.111  0.468  1.329  4.734 1.001  2400
#> dev.o[5,2]      1.239   1.733  0.001  0.122  0.535  1.658  6.296 1.003   760
#> dev.o[6,2]      1.185   1.512  0.001  0.142  0.623  1.697  5.427 1.006   410
#> dev.o[7,2]      0.840   1.173  0.001  0.088  0.386  1.107  4.125 1.001  3000
#> dev.o[8,2]      0.009   0.098  0.000  0.000  0.000  0.000  0.048 1.153  3000
#> dev.o[9,2]      1.022   1.438  0.001  0.113  0.460  1.341  5.345 1.001  3000
#> dev.o[10,2]     1.347   1.671  0.002  0.157  0.710  1.929  5.874 1.015   200
#> dev.o[11,2]     1.004   1.369  0.001  0.099  0.469  1.333  5.035 1.003   900
#> dev.o[9,3]      1.239   1.608  0.001  0.126  0.621  1.718  5.720 1.001  3000
#> dev.o[10,3]     1.025   1.383  0.001  0.096  0.495  1.399  4.790 1.004   630
#> dev.o[11,3]     0.925   1.217  0.000  0.095  0.462  1.282  4.303 1.016   130
#> hat.par[1,1]   26.662   4.494 18.228 23.554 26.466 29.557 35.948 1.004  1100
#> hat.par[2,1]    2.743   1.348  0.771  1.756  2.544  3.510  6.049 1.001  2100
#> hat.par[3,1]   10.140   1.177  7.267  9.499 10.387 11.036 11.656 1.065    58
#> hat.par[4,1]   22.865   2.513 17.502 21.199 23.052 24.662 27.246 1.002  1100
#> hat.par[5,1]    8.145   2.155  4.136  6.580  8.111  9.585 12.568 1.002  1600
#> hat.par[6,1]    3.294   0.975  1.362  2.616  3.307  3.996  5.075 1.007   360
#> hat.par[7,1]    2.620   1.309  0.668  1.677  2.431  3.360  5.605 1.045    50
#> hat.par[8,1]   11.998   0.028 11.992 12.000 12.000 12.000 12.000 1.135  1800
#> hat.par[9,1]   15.647   2.575 10.498 13.887 15.655 17.466 20.617 1.012   270
#> hat.par[10,1]   3.477   1.274  1.356  2.525  3.397  4.296  6.174 1.035    64
#> hat.par[11,1]   4.768   1.533  2.033  3.672  4.699  5.738  8.012 1.019   490
#> hat.par[1,2]   38.247   5.056 28.787 34.705 38.134 41.616 48.514 1.001  3000
#> hat.par[2,2]    4.252   1.656  1.544  3.034  4.089  5.286  7.946 1.001  3000
#> hat.par[3,2]    7.889   1.397  4.996  6.928  7.964  8.941 10.307 1.013   190
#> hat.par[4,2]   16.222   2.451 11.225 14.545 16.323 17.935 20.784 1.003   850
#> hat.par[5,2]    2.804   1.539  0.524  1.599  2.608  3.730  6.294 1.003  2700
#> hat.par[6,2]    3.734   0.936  1.721  3.118  3.813  4.450  5.306 1.031   120
#> hat.par[7,2]    2.340   1.202  0.462  1.431  2.197  3.102  5.023 1.081    36
#> hat.par[8,2]   14.996   0.047 14.976 15.000 15.000 15.000 15.000 1.153  3000
#> hat.par[9,2]   16.910   2.751 11.735 14.917 16.898 18.852 22.464 1.004   600
#> hat.par[10,2]   3.127   1.341  0.804  2.140  3.002  4.021  5.892 1.119    27
#> hat.par[11,2]   7.034   1.475  4.063  6.031  7.052  8.087  9.782 1.001  3000
#> hat.par[9,3]   21.354   2.266 16.786 19.836 21.473 22.994 25.408 1.005   470
#> hat.par[10,3]   8.390   1.475  5.361  7.397  8.559  9.476 10.877 1.018   180
#> hat.par[11,3]  11.244   1.654  7.805 10.161 11.391 12.490 14.013 1.013   480
#> m.tau           0.725   0.500  0.047  0.304  0.661  1.030  1.948 1.092    36
#> resdev.o[1]     1.953   1.999  0.055  0.547  1.292  2.684  7.337 1.002  1400
#> resdev.o[2]     1.993   1.890  0.052  0.594  1.442  2.742  6.976 1.001  3000
#> resdev.o[3]     2.048   1.977  0.061  0.633  1.487  2.881  7.090 1.014   150
#> resdev.o[4]     1.993   1.894  0.053  0.603  1.433  2.831  6.971 1.002  2200
#> resdev.o[5]     2.323   2.296  0.057  0.672  1.533  3.232  8.386 1.009   240
#> resdev.o[6]     2.252   1.878  0.111  0.884  1.837  3.074  7.081 1.014   190
#> resdev.o[7]     1.823   1.862  0.047  0.565  1.269  2.484  7.163 1.002  2800
#> resdev.o[8]     0.014   0.134  0.000  0.000  0.000  0.000  0.078 1.035  3000
#> resdev.o[9]     3.146   2.482  0.237  1.296  2.481  4.371  9.702 1.002  1100
#> resdev.o[10]    3.298   2.403  0.235  1.480  2.813  4.509  9.122 1.033    71
#> resdev.o[11]    2.804   2.097  0.256  1.256  2.351  3.819  8.098 1.017   130
#> tau             0.690   0.419  0.047  0.366  0.628  1.002  1.507 1.032    90
#> totresdev.o    23.648   6.595 12.122 18.963 23.260 27.719 37.834 1.030    70
#> 
#> For each parameter, n.eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
#> 
#> $data
#>                study t1 t2 t3 r1 r2 r3 m1 m2 m3  n1  n2 n3
#> 356    Richard, 2012  1  3  4 15 16 23  6  8  4  39  42 34
#> 357     Barone, 2010  1  2 NA 27 38 NA 19 20 NA 152 144 NA
#> 358 Weinbtraub, 2010  1  3 NA  2  5 NA  6  6 NA  27  28 NA
#> 359      Menza, 2009  1  4  5  4  2  9  6  7  5  17  18 17
#> 360      Devos, 2008  1  4  5  4  8 11  0  2  1  16  15 17
#> 361   Antonini, 2006  4  5 NA 10  8 NA  4  4 NA  16  15 NA
#> 362     Barone, 2006  2  4 NA 23 16 NA  1  7 NA  33  34 NA
#> 363  Rektorova, 2003  2  6 NA  8  3 NA  3  2 NA  22  19 NA
#> 364  Leentjens, 2003  1  4 NA  4  3 NA  0  0 NA   6   6 NA
#> 365    Wermuth, 1998  1  4 NA  3  2 NA  2  5 NA  19  18 NA
#> 366      Rabey, 1996  4  5 NA 12 15 NA  8 12 NA  20  27 NA
#> 
#> $model
#> [1] "RE"
#> 
#> $measure
#> [1] "OR"
#> 
#> $assumption
#> [1] "IDE-ARM"
#> 
#> $phi
#> NULL
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
#> $m_tau
#>        mean          sd        2.5%         25%         50%         75% 
#>  0.72462102  0.49998405  0.04744372  0.30360757  0.66055672  1.03020817 
#>       97.5%        Rhat       n.eff 
#>  1.94814690  1.09200846 36.00000000 
#> 
#> $frail_comp
#> [1] "4vs3"
#> 
#> attr(,"class")
#> [1] "run_ume"
# }
```
