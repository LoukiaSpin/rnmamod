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
#> EM[2,1]  0.6411756 0.8640783 -1.1278313  0.1485561  0.6549098  1.1130284
#> EM[3,1]  0.2138138 0.7277686 -1.2075876 -0.2617058  0.1969826  0.6857658
#> EM[4,1]  0.7857477 0.5874192 -0.5002016  0.4266032  0.8165053  1.1751339
#> EM[5,1]  1.8692829 0.8058200  0.2071046  1.3404155  1.8799027  2.3984321
#> EM[4,2] -0.4371592 1.0188729 -2.5550948 -1.0588785 -0.4172995  0.2310008
#> EM[6,2] -1.3695623 1.0670491 -3.5249186 -2.0700454 -1.3415512 -0.7024801
#> EM[4,3]  1.1275144 0.9990354 -0.8842449  0.5492763  1.1473832  1.7029166
#> EM[5,4] -1.2543557 1.5643457 -4.3514011 -2.4355531 -1.1403605 -0.1933647
#>             97.5%     Rhat n.eff
#> EM[2,1] 2.4165400 1.004742   590
#> EM[3,1] 1.6699797 1.009861   240
#> EM[4,1] 1.8133567 1.014485   290
#> EM[5,1] 3.4448129 1.013247  3000
#> EM[4,2] 1.4597717 1.009626   220
#> EM[6,2] 0.8093906 1.002493  1000
#> EM[4,3] 3.1397025 1.006827   650
#> EM[5,4] 1.7877831 1.035660    63
#> 
#> $dev_o
#>                    mean        sd         2.5%        25%       50%
#> dev.o[1,1]  1.087981797 1.4700033 0.0010246068 0.11492381 0.4905191
#> dev.o[2,1]  0.860650266 1.2173806 0.0008871066 0.08551090 0.3854807
#> dev.o[3,1]  1.401854355 1.9816491 0.0016828616 0.13838053 0.6074578
#> dev.o[4,1]  1.129676362 1.6053435 0.0014533556 0.11526583 0.5176434
#> dev.o[5,1]  1.019390392 1.4876829 0.0006432110 0.10105329 0.4296978
#> dev.o[6,1]  1.044869715 1.3608280 0.0008829905 0.11693039 0.5004370
#> dev.o[7,1]  1.034846204 1.3391678 0.0009033066 0.11024861 0.4984616
#> dev.o[8,1]  0.006036393 0.1063059 0.0000000000 0.00000000 0.0000000
#> dev.o[9,1]  0.891666859 1.2696417 0.0012679947 0.09723516 0.3965096
#> dev.o[10,1] 1.027408857 1.3053784 0.0013303410 0.11677908 0.5226921
#> dev.o[11,1] 0.870973830 1.2461908 0.0005384276 0.08471039 0.3877897
#> dev.o[1,2]  1.014184360 1.3843809 0.0010611204 0.10503386 0.4654616
#> dev.o[2,2]  1.166098235 1.5486182 0.0010256899 0.11866349 0.5588803
#> dev.o[3,2]  1.112410037 1.5901120 0.0007468840 0.10984866 0.4814913
#> dev.o[4,2]  1.081605776 1.4606211 0.0009950555 0.11670956 0.5146342
#> dev.o[5,2]  0.980369712 1.3315661 0.0009441174 0.10774468 0.4372487
#> dev.o[6,2]  1.148051934 1.4221610 0.0016532668 0.14568572 0.6016700
#> dev.o[7,2]  0.792879718 1.1097164 0.0007792466 0.08648384 0.3557735
#> dev.o[8,2]  0.013262129 0.1322449 0.0000000000 0.00000000 0.0000000
#> dev.o[9,2]  1.123605738 1.6248036 0.0011408239 0.09687636 0.4878892
#> dev.o[10,2] 1.406088275 1.7728873 0.0016861601 0.16705791 0.7491742
#> dev.o[11,2] 1.082492023 1.4734886 0.0010586106 0.10686120 0.5062203
#> dev.o[9,3]  1.094853460 1.5215918 0.0009251520 0.11277908 0.5299804
#> dev.o[10,3] 0.896401953 1.2136581 0.0012880443 0.10852684 0.4200497
#> dev.o[11,3] 0.991619333 1.4498019 0.0010066324 0.08708065 0.4508166
#>                      75%      97.5%     Rhat n.eff
#> dev.o[1,1]  1.512174e+00 5.48050170 1.001108  3000
#> dev.o[2,1]  1.143684e+00 4.26620120 1.001338  2400
#> dev.o[3,1]  1.852095e+00 7.36380260 1.012485   230
#> dev.o[4,1]  1.521927e+00 5.38805392 1.001300  3000
#> dev.o[5,1]  1.299893e+00 5.17458052 1.003216   770
#> dev.o[6,1]  1.492598e+00 4.82141773 1.001325  2400
#> dev.o[7,1]  1.413698e+00 4.75168842 1.001463  2100
#> dev.o[8,1]  3.197442e-14 0.01245678 1.200005  3000
#> dev.o[9,1]  1.186485e+00 4.45169332 1.000669  3000
#> dev.o[10,1] 1.397107e+00 4.81133473 1.002756   890
#> dev.o[11,1] 1.133162e+00 4.41842136 1.006116   460
#> dev.o[1,2]  1.402976e+00 5.11420444 1.000508  3000
#> dev.o[2,2]  1.591358e+00 5.71243576 1.001438  2100
#> dev.o[3,2]  1.494369e+00 5.45854660 1.002572  3000
#> dev.o[4,2]  1.472286e+00 5.29032044 1.002509  3000
#> dev.o[5,2]  1.369540e+00 4.61228484 1.000710  3000
#> dev.o[6,2]  1.670931e+00 4.99957367 1.004529   510
#> dev.o[7,2]  1.047903e+00 3.92595284 1.001446  2100
#> dev.o[8,2]  6.994405e-14 0.03395512 1.147125   660
#> dev.o[9,2]  1.503558e+00 5.71786291 1.001738  1700
#> dev.o[10,2] 2.027033e+00 6.30535070 1.001382  2300
#> dev.o[11,2] 1.486976e+00 5.36049366 1.001164  3000
#> dev.o[9,3]  1.438287e+00 5.32127557 1.006299   480
#> dev.o[10,3] 1.228097e+00 4.32514351 1.002036  1500
#> dev.o[11,3] 1.329610e+00 4.93126278 1.001104  3000
#> 
#> $hat_par
#>                    mean         sd       2.5%       25%       50%       75%
#> hat.par[1,1]  26.845878 4.78837744 18.0684546 23.510929 26.537944 29.979851
#> hat.par[2,1]   2.721465 1.33628820  0.7648237  1.747671  2.511186  3.523310
#> hat.par[3,1]  10.226134 1.23386829  7.3489904  9.528879 10.466885 11.182142
#> hat.par[4,1]  22.953304 2.64548645 17.4164863 21.207513 23.154867 24.820044
#> hat.par[5,1]   8.052856 2.09552480  4.1279610  6.573381  7.970784  9.405726
#> hat.par[6,1]   3.288968 0.95995239  1.4522109  2.589933  3.315218  3.995785
#> hat.par[7,1]   2.488775 1.26179284  0.6749404  1.560291  2.288980  3.202030
#> hat.par[8,1]  11.997093 0.04880353 11.9937732 12.000000 12.000000 12.000000
#> hat.par[9,1]  15.161336 2.64781492 10.0436460 13.328541 15.179412 16.933131
#> hat.par[10,1]  3.325602 1.27060929  1.2289497  2.391501  3.184565  4.140561
#> hat.par[11,1]  4.724338 1.55100410  2.0156863  3.636192  4.623136  5.694326
#> hat.par[1,2]  38.078742 5.13882373 28.3468268 34.485307 38.103953 41.523992
#> hat.par[2,2]   4.193637 1.66658573  1.5451085  2.973127  3.983172  5.200298
#> hat.par[3,2]   7.803955 1.47824818  4.5555184  6.843399  7.975981  8.886889
#> hat.par[4,2]  16.181571 2.56264961 11.0913833 14.392261 16.268928 18.014792
#> hat.par[5,2]   2.930561 1.45790136  0.7393941  1.835342  2.703530  3.828519
#> hat.par[6,2]   3.728065 0.92510664  1.7726629  3.095352  3.791322  4.425271
#> hat.par[7,2]   2.497733 1.19808864  0.6662640  1.601848  2.354538  3.195539
#> hat.par[8,2]  14.993513 0.06390693 14.9830320 15.000000 15.000000 15.000000
#> hat.par[9,2]  17.017037 2.85415906 11.5449462 15.113568 16.926601 18.981203
#> hat.par[10,2]  3.240577 1.31611663  1.0104264  2.256722  3.147384  4.123880
#> hat.par[11,2]  7.006276 1.53290123  3.9524878  5.947098  7.049734  8.076096
#> hat.par[9,3]  21.748889 2.28366106 16.9516925 20.269483 21.867515 23.389315
#> hat.par[10,3]  8.402404 1.36984637  5.5470297  7.447513  8.480455  9.467578
#> hat.par[11,3] 11.437128 1.62427132  7.9150357 10.448154 11.540144 12.608576
#>                   97.5%     Rhat n.eff
#> hat.par[1,1]  37.089281 1.001602  1800
#> hat.par[2,1]   5.800153 1.017302   130
#> hat.par[3,1]  11.854543 1.025668    93
#> hat.par[4,1]  27.559029 1.010784   340
#> hat.par[5,1]  12.400661 1.005619   750
#> hat.par[6,1]   5.141731 1.005754   600
#> hat.par[7,1]   5.523761 1.001077  3000
#> hat.par[8,1]  12.000000 1.200005  3000
#> hat.par[9,1]  20.386412 1.020162   110
#> hat.par[10,1]  6.082728 1.006205  1300
#> hat.par[11,1]  8.064451 1.002125  1200
#> hat.par[1,2]  48.572661 1.001683  1700
#> hat.par[2,2]   7.937909 1.005902   600
#> hat.par[3,2]  10.229482 1.011413   190
#> hat.par[4,2]  20.890703 1.010646   220
#> hat.par[5,2]   6.288531 1.003751   730
#> hat.par[6,2]   5.254138 1.004907   540
#> hat.par[7,2]   5.228095 1.003493   690
#> hat.par[8,2]  15.000000 1.147125   660
#> hat.par[9,2]  22.550706 1.001303  2500
#> hat.par[10,2]  6.040552 1.002452  1000
#> hat.par[11,2]  9.914263 1.004036   570
#> hat.par[9,3]  25.680903 1.003461   680
#> hat.par[10,3] 10.712773 1.002731  3000
#> hat.par[11,3] 14.174022 1.006412   720
#> 
#> $leverage_o
#>  [1]  1.086874797  0.621703407  1.369141051  1.129339977  1.018787904
#>  [6]  0.695834321  0.918140593 -0.004492432  0.888489137  0.837975523
#> [11]  0.708197291  1.013949276  0.983048373  1.095246780  1.076530228
#> [16]  0.978393958  0.783831521  0.663133915  0.011631843  1.001847886
#> [21]  0.666806462  0.773078597  0.824130650  0.749602194  0.934094743
#> 
#> $sign_dev_o
#>  [1]  1 -1 -1  1 -1  1  1 -1 -1  1 -1 -1  1  1 -1  1 -1 -1 -1 -1 -1  1  1  1 -1
#> 
#> $tau
#>        mean          sd        2.5%         25%         50%         75% 
#>  0.68600986  0.36465586  0.08759403  0.42396051  0.64759612  0.90332490 
#>       97.5%        Rhat       n.eff 
#>  1.43352676  1.11418436 42.00000000 
#> 
#> $m_tau
#>        mean          sd        2.5%         25%         50%         75% 
#>  0.76712939  0.51778577  0.07013626  0.36532302  0.69364887  1.07378061 
#>       97.5%        Rhat       n.eff 
#>  1.98854231  1.02921672 92.00000000 
#> 
#> $model_assessment
#>       DIC       pD      dev
#> 1 45.1046 20.82532 24.27928
#> 
#> $obs_comp
#> [1] "2vs1" "3vs1" "4vs1" "5vs1" "4vs2" "6vs2" "4vs3" "5vs4"
#> 
#> $jagsfit
#> Inference for Bugs model at "40", fit using jags,
#>  3 chains, each with 1000 iterations (first 0 discarded)
#>  n.sims = 3000 iterations saved. Running time = secs
#>               mu.vect sd.vect   2.5%    25%    50%    75%  97.5%  Rhat n.eff
#> EM[2,1]         0.641   0.864 -1.128  0.149  0.655  1.113  2.417 1.005   590
#> EM[3,1]         0.214   0.728 -1.208 -0.262  0.197  0.686  1.670 1.010   240
#> EM[4,1]         0.786   0.587 -0.500  0.427  0.817  1.175  1.813 1.014   290
#> EM[5,1]         1.869   0.806  0.207  1.340  1.880  2.398  3.445 1.013  3000
#> EM[4,2]        -0.437   1.019 -2.555 -1.059 -0.417  0.231  1.460 1.010   220
#> EM[6,2]        -1.370   1.067 -3.525 -2.070 -1.342 -0.702  0.809 1.002  1000
#> EM[4,3]         1.128   0.999 -0.884  0.549  1.147  1.703  3.140 1.007   650
#> EM[5,4]        -1.254   1.564 -4.351 -2.436 -1.140 -0.193  1.788 1.036    63
#> dev.o[1,1]      1.088   1.470  0.001  0.115  0.491  1.512  5.481 1.001  3000
#> dev.o[2,1]      0.861   1.217  0.001  0.086  0.385  1.144  4.266 1.001  2400
#> dev.o[3,1]      1.402   1.982  0.002  0.138  0.607  1.852  7.364 1.012   230
#> dev.o[4,1]      1.130   1.605  0.001  0.115  0.518  1.522  5.388 1.001  3000
#> dev.o[5,1]      1.019   1.488  0.001  0.101  0.430  1.300  5.175 1.003   770
#> dev.o[6,1]      1.045   1.361  0.001  0.117  0.500  1.493  4.821 1.001  2400
#> dev.o[7,1]      1.035   1.339  0.001  0.110  0.498  1.414  4.752 1.001  2100
#> dev.o[8,1]      0.006   0.106  0.000  0.000  0.000  0.000  0.012 1.200  3000
#> dev.o[9,1]      0.892   1.270  0.001  0.097  0.397  1.186  4.452 1.001  3000
#> dev.o[10,1]     1.027   1.305  0.001  0.117  0.523  1.397  4.811 1.003   890
#> dev.o[11,1]     0.871   1.246  0.001  0.085  0.388  1.133  4.418 1.006   460
#> dev.o[1,2]      1.014   1.384  0.001  0.105  0.465  1.403  5.114 1.001  3000
#> dev.o[2,2]      1.166   1.549  0.001  0.119  0.559  1.591  5.712 1.001  2100
#> dev.o[3,2]      1.112   1.590  0.001  0.110  0.481  1.494  5.459 1.003  3000
#> dev.o[4,2]      1.082   1.461  0.001  0.117  0.515  1.472  5.290 1.003  3000
#> dev.o[5,2]      0.980   1.332  0.001  0.108  0.437  1.370  4.612 1.001  3000
#> dev.o[6,2]      1.148   1.422  0.002  0.146  0.602  1.671  5.000 1.005   510
#> dev.o[7,2]      0.793   1.110  0.001  0.086  0.356  1.048  3.926 1.001  2100
#> dev.o[8,2]      0.013   0.132  0.000  0.000  0.000  0.000  0.034 1.147   660
#> dev.o[9,2]      1.124   1.625  0.001  0.097  0.488  1.504  5.718 1.002  1700
#> dev.o[10,2]     1.406   1.773  0.002  0.167  0.749  2.027  6.305 1.001  2300
#> dev.o[11,2]     1.082   1.473  0.001  0.107  0.506  1.487  5.360 1.001  3000
#> dev.o[9,3]      1.095   1.522  0.001  0.113  0.530  1.438  5.321 1.006   480
#> dev.o[10,3]     0.896   1.214  0.001  0.109  0.420  1.228  4.325 1.002  1500
#> dev.o[11,3]     0.992   1.450  0.001  0.087  0.451  1.330  4.931 1.001  3000
#> hat.par[1,1]   26.846   4.788 18.068 23.511 26.538 29.980 37.089 1.002  1800
#> hat.par[2,1]    2.721   1.336  0.765  1.748  2.511  3.523  5.800 1.017   130
#> hat.par[3,1]   10.226   1.234  7.349  9.529 10.467 11.182 11.855 1.026    93
#> hat.par[4,1]   22.953   2.645 17.416 21.208 23.155 24.820 27.559 1.011   340
#> hat.par[5,1]    8.053   2.096  4.128  6.573  7.971  9.406 12.401 1.006   750
#> hat.par[6,1]    3.289   0.960  1.452  2.590  3.315  3.996  5.142 1.006   600
#> hat.par[7,1]    2.489   1.262  0.675  1.560  2.289  3.202  5.524 1.001  3000
#> hat.par[8,1]   11.997   0.049 11.994 12.000 12.000 12.000 12.000 1.200  3000
#> hat.par[9,1]   15.161   2.648 10.044 13.329 15.179 16.933 20.386 1.020   110
#> hat.par[10,1]   3.326   1.271  1.229  2.392  3.185  4.141  6.083 1.006  1300
#> hat.par[11,1]   4.724   1.551  2.016  3.636  4.623  5.694  8.064 1.002  1200
#> hat.par[1,2]   38.079   5.139 28.347 34.485 38.104 41.524 48.573 1.002  1700
#> hat.par[2,2]    4.194   1.667  1.545  2.973  3.983  5.200  7.938 1.006   600
#> hat.par[3,2]    7.804   1.478  4.556  6.843  7.976  8.887 10.229 1.011   190
#> hat.par[4,2]   16.182   2.563 11.091 14.392 16.269 18.015 20.891 1.011   220
#> hat.par[5,2]    2.931   1.458  0.739  1.835  2.704  3.829  6.289 1.004   730
#> hat.par[6,2]    3.728   0.925  1.773  3.095  3.791  4.425  5.254 1.005   540
#> hat.par[7,2]    2.498   1.198  0.666  1.602  2.355  3.196  5.228 1.003   690
#> hat.par[8,2]   14.994   0.064 14.983 15.000 15.000 15.000 15.000 1.147   660
#> hat.par[9,2]   17.017   2.854 11.545 15.114 16.927 18.981 22.551 1.001  2500
#> hat.par[10,2]   3.241   1.316  1.010  2.257  3.147  4.124  6.041 1.002  1000
#> hat.par[11,2]   7.006   1.533  3.952  5.947  7.050  8.076  9.914 1.004   570
#> hat.par[9,3]   21.749   2.284 16.952 20.269 21.868 23.389 25.681 1.003   680
#> hat.par[10,3]   8.402   1.370  5.547  7.448  8.480  9.468 10.713 1.003  3000
#> hat.par[11,3]  11.437   1.624  7.915 10.448 11.540 12.609 14.174 1.006   720
#> m.tau           0.767   0.518  0.070  0.365  0.694  1.074  1.989 1.029    92
#> resdev.o[1]     2.102   2.010  0.051  0.598  1.479  3.053  7.222 1.001  3000
#> resdev.o[2]     2.027   1.894  0.056  0.642  1.506  2.837  7.127 1.006   390
#> resdev.o[3]     2.514   2.531  0.063  0.658  1.740  3.572  9.318 1.007   400
#> resdev.o[4]     2.211   2.160  0.065  0.663  1.593  3.054  7.893 1.002  3000
#> resdev.o[5]     2.000   1.981  0.053  0.569  1.372  2.838  7.146 1.003   730
#> resdev.o[6]     2.193   1.745  0.088  0.889  1.832  3.050  6.514 1.006   350
#> resdev.o[7]     1.828   1.712  0.058  0.575  1.324  2.537  6.418 1.003   830
#> resdev.o[8]     0.019   0.209  0.000  0.000  0.000  0.000  0.065 1.175   910
#> resdev.o[9]     3.110   2.566  0.241  1.239  2.489  4.210  9.850 1.001  2900
#> resdev.o[10]    3.330   2.346  0.315  1.589  2.838  4.494  9.076 1.002  1100
#> resdev.o[11]    2.945   2.403  0.211  1.218  2.314  3.993  9.320 1.011   400
#> tau             0.686   0.365  0.088  0.424  0.648  0.903  1.434 1.114    42
#> totresdev.o    24.279   6.906 12.751 19.240 23.655 28.503 39.696 1.007   320
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
#>  0.76712939  0.51778577  0.07013626  0.36532302  0.69364887  1.07378061 
#>       97.5%        Rhat       n.eff 
#>  1.98854231  1.02921672 92.00000000 
#> 
#> $frail_comp
#> [1] "4vs3"
#> 
#> attr(,"class")
#> [1] "run_ume"
# }
```
