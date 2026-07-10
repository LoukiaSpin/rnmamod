# Perform Bayesian pairwise or network meta-analysis

Performs a one-stage pairwise or network meta-analysis while addressing
aggregate binary or continuous missing participant outcome data via the
pattern-mixture model.

## Usage

``` r
run_model(
  data,
  measure,
  model,
  assumption,
  heter_prior,
  mean_misspar,
  var_misspar,
  D,
  ref,
  base_risk,
  n_chains,
  n_iter,
  n_burnin,
  n_thin,
  inits = NULL
)
```

## Format

The columns of the data-frame in the argument `data` refer to the
following elements for a continuous outcome:

|        |                                                             |
|--------|-------------------------------------------------------------|
| **t**  | An intervention identifier in each arm.                     |
|        |                                                             |
| **y**  | The observed mean value of the outcome in each arm.         |
|        |                                                             |
| **sd** | The observed standard deviation of the outcome in each arm. |
|        |                                                             |
| **m**  | The number of missing participant outcome data in each arm. |
|        |                                                             |
| **n**  | The number of randomised participants in each arm.          |

For a binary outcome, the columns of the data-frame in the argument
`data` refer to the following elements:

|       |                                                             |
|-------|-------------------------------------------------------------|
| **t** | An intervention identifier in each arm.                     |
|       |                                                             |
| **r** | The observed number of events of the outcome in each arm.   |
|       |                                                             |
| **m** | The number of missing participant outcome data in each arm. |
|       |                                                             |
| **n** | The number of randomised participants in each arm.          |

The number of rows in `data` equals the number of collected trials. Each
element appears in `data` as many times as the maximum number of
interventions compared in a trial of the dataset. In pairwise
meta-analysis, the maximum number of arms is inherently two. The same
holds for a network meta-analysis without multi-arm trials. In the case
of network meta-analysis with multi-arm trials, the maximum number of
arms exceeds two. See 'Examples' that illustrates the structure of
`data` for a network with a maximum number of four arms. It is not a
prerequisite of `run_model` that the multi-arm trials appear at the
bottom of the dataset.

## Arguments

- data:

  A data-frame of the one-trial-per-row format with arm-level data. See
  'Format' for the specification of the columns.

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

- heter_prior:

  A list of three elements with the following order: 1) a character
  string indicating the distribution with (currently available) values
  `"halfnormal"`, `"uniform"`, `"lognormal"`, or `"logt"`; 2) two
  numeric values that refer to the parameters of the selected
  distribution. For `"lognormal"`, and `"logt"` these numbers refer to
  the mean and precision, respectively. For `"halfnormal"`, these
  numbers refer to zero and the scale parameter (equal to 4 or 1 being
  the corresponding precision of the scale parameter 0.5 or 1). For
  `"uniform"`, these numbers refer to the minimum and maximum value of
  the distribution. See 'Details' in
  [`heterogeneity_param_prior`](https://loukiaspin.github.io/rnmamod/reference/heterogeneity_param_prior.md).

- mean_misspar:

  A scalar or numeric vector of two numeric values for the mean of the
  normal distribution of the informative missingness parameter (see
  'Details'). The default argument is 0 and corresponds to the
  missing-at-random assumption. See also 'Details' in
  [`missingness_param_prior`](https://loukiaspin.github.io/rnmamod/reference/missingness_param_prior.md).

- var_misspar:

  A positive non-zero number for the variance of the normal distribution
  of the informative missingness parameter. When the `measure` is
  `"OR"`, `"MD"`, or `"SMD"` the default argument is 1. When the
  `measure` is `"ROM"` the default argument is 0.04.

- D:

  A binary number for the direction of the outcome. Set `D = 1` for
  beneficial outcome and `D = 0` for harmful outcome.

- ref:

  An integer specifying the reference intervention. The number should
  match the intervention identifier under element **t** in `data` (See
  'Format').

- base_risk:

  A scalar, a vector of length three with elements sorted in ascending
  order, or a matrix with two columns and number of rows equal to the
  number of relevant trials. In the case of a scalar or vector, the
  elements should be in the interval (0, 1) (see 'Details'). If
  `base_risk` has not been defined, the function uses the median event
  risk for the reference intervention from the corresponding trials in
  `data`. This argument is only relevant for a binary outcome.

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

A list of R2jags output on the summaries of the posterior distribution,
and the Gelman-Rubin convergence diagnostic (Gelman et al., 1992) of the
following monitored parameters for a fixed-effect pairwise
meta-analysis:

- EM:

  The estimated summary effect measure (according to the argument
  `measure`).

- EM_LOR:

  The estimated summary odd ratio in the logarithmic scale when
  `measure = "RR"` or `measure = "RD"`.

- dev_o:

  The deviance contribution of each trial-arm based on the observed
  outcome.

- hat_par:

  The fitted outcome at each trial-arm.

- phi:

  The informative missingness parameter.

For a fixed-effect network meta-analysis, the output additionally
includes:

- SUCRA:

  The surface under the cumulative ranking curve for each intervention.

- SUCRA_LOR:

  The surface under the cumulative ranking curve for each intervention
  under the odds ratio effect measure when `measure = "RR"` or
  `measure = "RD"`.

- effectiveneness:

  The ranking probability of each intervention for every rank.

For a random-effects pairwise meta-analysis, the output additionally
includes the following elements:

- EM_pred:

  The predicted summary effect measure (according to the argument
  `measure`).

- EM_pred_LOR:

  The predicted summary odds ratio in the logarithmic scale when
  `measure = "RR"` or `measure = "RD"`.

- delta:

  The estimated trial-specific effect measure (according to the argument
  `measure`).

- tau:

  The between-trial standard deviation.

In network meta-analysis, `EM` and `EM_pred` refer to all possible
pairwise comparisons of interventions in the network. Furthermore, `tau`
is typically assumed to be common for all observed comparisons in the
network. For a multi-arm trial, we estimate a total of *T-1* `delta` for
comparisons with the baseline intervention of the trial (found in the
first column of the element **t**), with *T* being the number of
interventions in the trial.

Furthermore, the output includes the following elements:

- leverage_o:

  The leverage for the observed outcome at each trial-arm.

- sign_dev_o:

  The sign of the difference between observed and fitted outcome at each
  trial-arm.

- model_assessment:

  A data-frame on the measures of model assessment: deviance information
  criterion, number of effective parameters, and total residual
  deviance.

- indic:

  The sign of basic parameters in relation to the reference intervention
  as specified in argument `reg`

- jagsfit:

  An object of S3 class
  [`jags`](https://rdrr.io/pkg/R2jags/man/jags.html) with the posterior
  results on all monitored parameters to be used in the
  [`mcmc_diagnostics`](https://loukiaspin.github.io/rnmamod/reference/mcmc_diagnostics.md)
  function.

The `run_model` function also returns the arguments `data`, `measure`,
`model`, `assumption`, `heter_prior`, `mean_misspar`, `var_misspar`,
`D`, `ref`, `base_risk`, `n_chains`, `n_iter`, `n_burnin`, and `n_thin`
as specified by the user to be inherited by other functions of the
package.

## Details

The model runs in `JAGS` and the progress of the simulation appears on
the R console. The output of `run_model` is used as an S3 object by
other functions of the package to be processed further and provide an
end-user-ready output.

The
[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md)
function is called to prepare the data for the Bayesian analysis.
[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md)
creates the pseudo-data-frames `m_new`, and `I`, that have the same
dimensions with the element `N`. `m_new` takes the zero value for the
observed trial-arms with unreported missing participant outcome data
(i.e., `m` equals `NA` for the corresponding trial-arms), the same value
with `m` for the observed trial-arms with reported missing participant
outcome data, and `NA` for the unobserved trial-arms. `I` is a dummy
pseudo-data-frame and takes the value one for the observed trial-arms
with reported missing participant outcome data, the zero value for the
observed trial-arms with unreported missing participant outcome data
(i.e., `m_new` equals zero for the corresponding trial-arms), and `NA`
for the unobserved trial-arms. Thus, `I` indicates whether missing
participant outcome data have been collected for the observed
trial-arms. If the user has not defined the element **m** in `data`,
`m_new` and `I` take the zero value for all observed trial-arms to
indicate that no missing participant outcome data have been collected
for the analysed outcome. See 'Details' in
[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md).

Furthermore,
[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md)
sorts the interventions across the arms of each trial in an ascending
order and correspondingly the remaining elements in `data` (see
'Format').
[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md)
considers the first column in **t** as being the control arm for every
trial. Thus, this sorting ensures that interventions with a lower
identifier are consistently treated as the control arm in each trial.
This case is relevant in non-star-shaped networks.

The model is updated until convergence using the
[`autojags`](https://rdrr.io/pkg/R2jags/man/autojags.html) function of
the R-package [R2jags](https://CRAN.R-project.org/package=R2jags) with 2
updates and number of iterations and thinning equal to `n_iter` and
`n_thin`, respectively.

To perform a Bayesian pairwise or network meta-analysis, the
[`prepare_model`](https://loukiaspin.github.io/rnmamod/reference/prepare_model.md)
function is called which contains the WinBUGS code as written by Dias et
al. (2013a) for binomial and normal likelihood to analyse aggregate
binary and continuous outcome data, respectively.
[`prepare_model`](https://loukiaspin.github.io/rnmamod/reference/prepare_model.md)
uses the consistency model (as described in Lu and Ades (2006)) to
estimate all possible comparisons in the network. It also accounts for
the multi-arm trials by assigning conditional univariate normal
distributions on the underlying trial-specific effect size of
comparisons with the baseline arm of the multi-arm trial (Dias et al.,
2013a).

The code of Dias et al. (2013a) has been extended to incorporate the
pattern-mixture model to adjust the underlying outcome in each arm of
every trial for missing participant outcome data (Spineli et al., 2021;
Spineli, 2019a; Turner et al., 2015). The assumptions about the
missingness parameter are specified using the arguments `mean_misspar`
and `var_misspar`. Specifically, `run_model` considers the informative
missingness odds ratio in the logarithmic scale for binary outcome data
(Spineli, 2019a; Turner et al., 2015; White et al., 2008), the
informative missingness difference of means when `measure` is `"MD"` or
`"SMD"`, and the informative missingness ratio of means in the
logarithmic scale when `measure` is `"ROM"` (Spineli et al., 2021;
Mavridis et al., 2015).

When `assumption` is trial-specific (i.e., `"IDE-TRIAL"` or
`"HIE-TRIAL"`), or independent (i.e., `"IND-CORR"` or `"IND-UNCORR"`),
only one numeric value can be assigned to `mean_misspar` because the
same missingness scenario is applied to all trials and trial-arms of the
dataset, respectively. When `assumption` is `"IDE-ARM"` or `"HIE-ARM"`,
a maximum of two *different* or *identical* numeric values can be
assigned as a vector to `mean_misspars`: the first value refers to the
experimental arm, and the second value refers to the control arm of a
trial. In the case of a network, the first value is considered for all
non-reference interventions and the second value is considered for the
reference intervention of the network (i.e., the intervention with
identifier equal to `ref`). This is necessary to ensure transitivity in
the assumptions for the missingness parameter across the network
(Spineli, 2019b).

When there is at least one trial-arm with unreported missing participant
outcome data (i.e., `m` equals `NA` for the corresponding trial-arms) or
when missing participant outcome data have not been collected for the
analysed outcome (i.e., `m` is missing in `data`), `run_model` assigns
the assumption `"IND-UNCORR"` to `assumption`.

Currently, there are no empirically-based prior distributions for the
informative missingness parameters. The user may refer to Spineli
(2019), Turner et al. (2015), Mavridis et al. (2015), and White et al.
(2008) to determine `mean_misspar` and select a proper value for
`var_misspar`.

The scalar `base_risk` refers to a fixed baseline risk for the selected
reference intervention (as specified with `ref`). When `base_risk` is a
three-element vector, it refers to a random baseline risk and the
elements should be sorted in ascending order as they refer to the lower
bound, mean value, and upper bound of the 95% confidence interval for
the baseline risk for the selected reference intervention. The
[`baseline_model`](https://loukiaspin.github.io/rnmamod/reference/baseline_model.md)
function is called to calculate the mean and variance of the
approximately normal distribution of the logit of an event for `ref`
using these three elements (Dias et al., 2018). When `base_risk` is a
matrix, it refers to the predicted baseline risk with first column being
the number of events, and second column being the sample size of the
corresponding trials on the selected reference intervention. Then the
[`baseline_model`](https://loukiaspin.github.io/rnmamod/reference/baseline_model.md)
function is called that contains the WinBUGS code as written by Dias et
al. (2013b) for the hierarchical baseline model. The posterior mean and
precision of the predictive distribution of the logit of an event for
the selected reference intervention are plugged in the WinBUGS code for
the relative effects model (via the
[`prepare_model`](https://loukiaspin.github.io/rnmamod/reference/prepare_model.md)
function). The matrix `base_risk` should not comprise the trials in
`data` that include the `ref`, unless justified (Dias et al., 2018).

To obtain unique absolute risks for each intervention, the network
meta-analysis model has been extended to incorporate the transitive
risks framework, namely, an intervention has the same absolute risk
regardless of the comparator intervention(s) in a trial (Spineli et al.,
2017). The absolute risks are a function of the odds ratio (the
**base-case** effect measure for a binary outcome) and the selected
baseline risk for the reference intervention (`ref`) (Appendix in Dias
et al., 2013a). We advocate using the odds ratio as an effect measure
for its desired mathematical properties. Then, the relative risk and
risk difference can be obtained as a function of the absolute risks of
the corresponding interventions in the comparison of interest. Hence,
regardless of the selected `measure` for a binary outcome, `run_model`
performs pairwise or network meta-analysis based on the odds ratio.

## References

Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing
between-study heterogeneity and inconsistency in mixed treatment
comparisons: Application to stroke prevention treatments in individuals
with non-rheumatic atrial fibrillation. *Stat Med*
2009;**28**(14):1861–81. doi: 10.1002/sim.3594

Dias S, Ades AE, Welton NJ, Jansen JP, Sutton AJ. Network Meta-Analysis
for Decision Making. Chichester (UK): Wiley; 2018.

Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
making 2: a generalized linear modeling framework for pairwise and
network meta-analysis of randomized controlled trials. *Med Decis
Making* 2013a;**33**(5):607–17. doi: 10.1177/0272989X12458724

Dias S, Welton NJ, Sutton AJ, Ades AE. Evidence synthesis for decision
making 5: the baseline natural history model. *Med Decis Making*
2013b;**33**(5):657–70. doi: 10.1177/0272989X13485155

Gelman A, Rubin DB. Inference from iterative simulation using multiple
sequences. *Stat Sci* 1992;**7**(4):457–72. doi: 10.1214/ss/1177011136

Lu G, Ades AE. Assessing evidence inconsistency in mixed treatment
comparisons. *J Am Stat Assoc* 2006;**101**:447–59. doi:
10.1198/016214505000001302

Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for
uncertainty due to missing continuous outcome data in pairwise and
network meta-analysis. *Stat Med* 2015;**34**(5):721–41. doi:
10.1002/sim.6365

Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing
outcome data in network meta-analysis: a one-stage pattern-mixture model
approach. *Stat Methods Med Res* 2021;**30**(4):958–75. doi:
10.1177/0962280220983544

Spineli LM. An empirical comparison of Bayesian modelling strategies for
missing binary outcome data in network meta-analysis. *BMC Med Res
Methodol* 2019a;**19**(1):86. doi: 10.1186/s12874-019-0731-y

Spineli LM. Modeling missing binary outcome data while preserving
transitivity assumption yielded more credible network meta-analysis
results. *J Clin Epidemiol* 2019b;**105**:19–26. doi:
10.1016/j.jclinepi.2018.09.002

Spineli LM, Brignardello-Petersen R, Heen AF, Achille F, Brandt L,
Guyatt GH, et al. Obtaining absolute effect estimates to facilitate
shared decision making in the context of multiple-treatment comparisons.
Abstracts of the Global Evidence Summit, Cape Town, South Africa.
*Cochrane Database of Systematic Reviews* 2017;**9**(Suppl 1):18911.

Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account
for uncertainty due to missing binary outcome data in pairwise
meta-analysis. *Stat Med* 2015;**34**(12):2062–80. doi: 10.1002/sim.6475

White IR, Higgins JP, Wood AM. Allowing for uncertainty due to missing
data in meta-analysis–part 1: two-stage methods. *Stat Med*
2008;**27**(5):711–27. doi: 10.1002/sim.3008

## See also

[`autojags`](https://rdrr.io/pkg/R2jags/man/autojags.html),
[`baseline_model`](https://loukiaspin.github.io/rnmamod/reference/baseline_model.md),
[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md),
[`heterogeneity_param_prior`](https://loukiaspin.github.io/rnmamod/reference/heterogeneity_param_prior.md),
[`jags`](https://rdrr.io/pkg/R2jags/man/jags.html),
[`missingness_param_prior`](https://loukiaspin.github.io/rnmamod/reference/missingness_param_prior.md),
[`prepare_model`](https://loukiaspin.github.io/rnmamod/reference/prepare_model.md)

## Author

Loukia M. Spineli

## Examples

``` r
data("nma.baker2009")

# Show the first six trials of the dataset
head(nma.baker2009)
#>                   study t1 t2 t3 t4 r1 r2 r3 r4 m1 m2 m3 m4  n1  n2 n3 n4
#> 1 Llewellyn-Jones, 1996  1  4 NA NA  3  0 NA NA  1  0 NA NA   8   8 NA NA
#> 2        Paggiaro, 1998  1  4 NA NA 51 45 NA NA 27 19 NA NA 139 142 NA NA
#> 3          Mahler, 1999  1  7 NA NA 47 28 NA NA 23  9 NA NA 143 135 NA NA
#> 4        Casaburi, 2000  1  8 NA NA 41 45 NA NA 18 12 NA NA 191 279 NA NA
#> 5       van Noord, 2000  1  7 NA NA 18 11 NA NA  8  7 NA NA  50  47 NA NA
#> 6         Rennard, 2001  1  7 NA NA 41 38 NA NA 29 22 NA NA 135 132 NA NA

# \donttest{
# Perform a random-effects network meta-analysis
# Note: Ideally, set 'n_iter' to 10000 and 'n_burnin' to 1000
run_model(data = nma.baker2009,
          measure = "OR",
          model = "RE",
          assumption = "IDE-ARM",
          heter_prior = list("halfnormal", 0, 1),
          mean_misspar = c(0, 0),
          var_misspar = 1,
          D = 0,
          ref = 1,
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
#>    Observed stochastic nodes: 100
#>    Unobserved stochastic nodes: 148
#>    Total graph size: 2596
#> 
#> Initializing model
#> 
#> ... Updating the model until convergence
#> $EM
#>                mean        sd         2.5%         25%          50%
#> EM[2,1] -0.94318940 0.4282660 -1.769702800 -1.22834512 -0.953791702
#> EM[3,1] -0.66085480 0.4328575 -1.576573840 -0.93880826 -0.640237486
#> EM[4,1] -0.24350971 0.3247757 -0.797694421 -0.47632380 -0.274695102
#> EM[5,1] -0.51256592 0.2826923 -1.051798760 -0.68606275 -0.534742892
#> EM[6,1] -0.02304669 0.2478971 -0.551750217 -0.19005646 -0.001313911
#> EM[7,1] -0.46436918 0.1818440 -0.795934462 -0.59480594 -0.476224482
#> EM[8,1] -0.53037100 0.1636677 -0.856308808 -0.63603541 -0.528799981
#> EM[3,2]  0.28233460 0.4849089 -0.663818824 -0.05428280  0.281364826
#> EM[4,2]  0.69967969 0.4619547 -0.238505246  0.39258492  0.674835766
#> EM[5,2]  0.43062348 0.4715143 -0.545419238  0.12472308  0.421168841
#> EM[6,2]  0.92014271 0.4492652 -0.009598997  0.64503412  0.938670730
#> EM[7,2]  0.47882022 0.4356892 -0.368431883  0.19361109  0.473565821
#> EM[8,2]  0.41281840 0.4327773 -0.479984957  0.14080569  0.430954138
#> EM[4,3]  0.41734508 0.4636523 -0.470050703  0.11833755  0.402650544
#> EM[5,3]  0.14828888 0.4728140 -0.802514698 -0.16448440  0.149472997
#> EM[6,3]  0.63780811 0.4404039 -0.243781448  0.33197830  0.662458779
#> EM[7,3]  0.19648562 0.4110816 -0.582797497 -0.07387556  0.212225632
#> EM[8,3]  0.13048380 0.4197854 -0.674443022 -0.14976333  0.141303994
#> EM[5,4] -0.26905621 0.3767213 -0.999425235 -0.52295141 -0.270627908
#> EM[6,4]  0.22046302 0.3501291 -0.463590641 -0.01604658  0.228846331
#> EM[7,4] -0.22085947 0.2992480 -0.845049071 -0.41448047 -0.202897004
#> EM[8,4] -0.28686128 0.3202874 -0.980823076 -0.49363073 -0.258193562
#> EM[6,5]  0.48951923 0.3349909 -0.242348071  0.28730669  0.506081192
#> EM[7,5]  0.04819674 0.2954861 -0.567812732 -0.14735367  0.045475539
#> EM[8,5] -0.01780508 0.2889111 -0.630461915 -0.19007078 -0.007521709
#> EM[7,6] -0.44132249 0.2398004 -0.886279100 -0.61049440 -0.448113701
#> EM[8,6] -0.50732431 0.2300243 -0.960203129 -0.65419406 -0.512675598
#> EM[8,7] -0.06600182 0.1640566 -0.387078024 -0.17429735 -0.059609984
#>                  75%       97.5%     Rhat n.eff
#> EM[2,1] -0.671779167 -0.04982525 1.084250    35
#> EM[3,1] -0.375844443  0.13231876 1.063323    52
#> EM[4,1] -0.021187119  0.43804638 1.129526    21
#> EM[5,1] -0.347135242  0.12999656 1.024410   770
#> EM[6,1]  0.165575353  0.39681644 1.131188    23
#> EM[7,1] -0.341474801 -0.08728508 1.039033   120
#> EM[8,1] -0.413478109 -0.21821686 1.023187   350
#> EM[3,2]  0.630363955  1.22180562 1.007434  1300
#> EM[4,2]  1.018310893  1.63247671 1.020274   210
#> EM[5,2]  0.725790286  1.37586034 1.041768    57
#> EM[6,2]  1.220188336  1.73286764 1.048786    58
#> EM[7,2]  0.764643099  1.31379879 1.050465    59
#> EM[8,2]  0.717396997  1.19928243 1.077844    35
#> EM[4,3]  0.693606746  1.39201252 1.018766   290
#> EM[5,3]  0.465053427  1.11915730 1.027751    89
#> EM[6,3]  0.926396441  1.53921184 1.045761    66
#> EM[7,3]  0.463536586  1.01641478 1.037873    87
#> EM[8,3]  0.408280691  0.97618088 1.055155    46
#> EM[5,4] -0.016146626  0.48390863 1.076418    36
#> EM[6,4]  0.467087647  0.82565672 1.135701    20
#> EM[7,4]  0.001926603  0.31349688 1.091873    27
#> EM[8,4] -0.043369786  0.24981603 1.148424    18
#> EM[6,5]  0.717485525  1.11592045 1.046171    50
#> EM[7,5]  0.238850130  0.59664954 1.009634  2000
#> EM[8,5]  0.178830917  0.50179916 1.006052   420
#> EM[7,6] -0.284243842  0.04705296 1.087498    31
#> EM[8,6] -0.352867888 -0.05126907 1.107591    25
#> EM[8,7]  0.046000012  0.24377872 1.030018    81
#> 
#> $dev_o
#>                  mean        sd         2.5%        25%       50%       75%
#> dev.o[1,1]  2.3781592 2.3362406 0.0074928088 0.59623249 1.7125559 3.4186954
#> dev.o[2,1]  0.9677025 1.3462197 0.0010395001 0.10086898 0.4477280 1.2801405
#> dev.o[3,1]  1.1863702 1.5443602 0.0018089309 0.13856063 0.6186585 1.6904974
#> dev.o[4,1]  0.7162283 1.0033377 0.0005864476 0.07738118 0.3260675 0.8990377
#> dev.o[5,1]  0.7391687 1.0215510 0.0009481909 0.08006467 0.3408617 1.0013088
#> dev.o[6,1]  0.9364762 1.2742631 0.0009845114 0.10955964 0.4685048 1.2741555
#> dev.o[7,1]  0.8067370 1.1361442 0.0006366889 0.07994412 0.3810305 1.0738986
#> dev.o[8,1]  0.7356023 1.0821420 0.0007070854 0.07072792 0.3215362 0.9401305
#> dev.o[9,1]  0.7416247 1.0841233 0.0005288078 0.06913609 0.3351249 0.9530078
#> dev.o[10,1] 0.6170400 0.8379483 0.0007142445 0.06873267 0.2945055 0.8323558
#> dev.o[11,1] 0.9792968 1.3215701 0.0013909247 0.10555992 0.4396493 1.3354563
#> dev.o[12,1] 1.2912889 1.5880889 0.0018016100 0.18507071 0.7142733 1.8307768
#> dev.o[13,1] 1.2898297 1.6643804 0.0018884285 0.14097851 0.6510724 1.7938726
#> dev.o[14,1] 0.9132117 1.2461620 0.0010651947 0.10202609 0.4228366 1.2259668
#> dev.o[15,1] 0.8461846 1.1264491 0.0005517014 0.08526399 0.4082004 1.1608607
#> dev.o[16,1] 1.2062596 1.6323540 0.0007069609 0.12275116 0.5572995 1.6330766
#> dev.o[17,1] 1.6418756 1.9446937 0.0030855801 0.22264217 0.9109586 2.4176909
#> dev.o[18,1] 0.9971523 1.4219907 0.0008607413 0.09865149 0.4512387 1.3587211
#> dev.o[19,1] 1.6197726 1.8031357 0.0036483993 0.27063590 1.0220878 2.3785075
#> dev.o[20,1] 0.7505731 1.0928379 0.0006211646 0.07353664 0.3472989 0.9853409
#> dev.o[21,1] 1.0763259 1.4280354 0.0011449648 0.11099901 0.5315871 1.4644610
#> dev.o[1,2]  3.2528238 1.9464467 0.6482235388 1.83865566 2.8756566 4.2717228
#> dev.o[2,2]  0.9695222 1.3081284 0.0006360070 0.09630903 0.4848161 1.3591022
#> dev.o[3,2]  1.0463302 1.3474659 0.0008398151 0.11592869 0.5246374 1.4263828
#> dev.o[4,2]  0.8286589 1.2102489 0.0007865829 0.07227079 0.3619318 1.0909261
#> dev.o[5,2]  0.6089100 0.8868392 0.0005635100 0.06223971 0.2703683 0.8000208
#> dev.o[6,2]  1.0055967 1.3362468 0.0012653573 0.10846500 0.4934576 1.3917448
#> dev.o[7,2]  0.8651850 1.1861516 0.0006449749 0.08699506 0.4236512 1.1507982
#> dev.o[8,2]  0.7133597 1.0124883 0.0004214583 0.07462075 0.3091599 0.9719862
#> dev.o[9,2]  0.6992260 1.0401658 0.0006728956 0.07500464 0.3209577 0.8692467
#> dev.o[10,2] 1.5749895 1.8691974 0.0027171419 0.22142505 0.9032006 2.2526562
#> dev.o[11,2] 0.9784296 1.4006716 0.0012880248 0.10098823 0.4518714 1.2811199
#> dev.o[12,2] 0.8178376 1.1126811 0.0007984376 0.09152524 0.3739715 1.0748010
#> dev.o[13,2] 0.9790166 1.3837251 0.0009325170 0.10405034 0.4375725 1.2925432
#> dev.o[14,2] 0.7602221 1.0879182 0.0008098462 0.07728314 0.3335235 0.9923854
#> dev.o[15,2] 0.9487169 1.2390592 0.0011571847 0.10707651 0.4807316 1.3152069
#> dev.o[16,2] 1.2072451 1.7096792 0.0011490034 0.12533895 0.5605912 1.6005865
#> dev.o[17,2] 1.9280589 1.7895219 0.0096300314 0.56522386 1.4477336 2.8198032
#> dev.o[18,2] 0.9107474 1.2662133 0.0011051204 0.10066204 0.4372561 1.2216145
#> dev.o[19,2] 0.3964745 0.5864798 0.0004450458 0.03959436 0.1746349 0.5064449
#> dev.o[20,2] 0.7371697 1.0015309 0.0009654283 0.07375008 0.3498350 0.9995161
#> dev.o[21,2] 0.9416026 1.1985451 0.0012017964 0.12145811 0.4855115 1.2913916
#> dev.o[9,3]  0.7546427 1.0339791 0.0005993761 0.07772240 0.3533254 1.0514297
#> dev.o[10,3] 0.7055757 1.0015327 0.0007278120 0.07020217 0.3328251 0.9182303
#> dev.o[12,3] 1.1922746 1.5964033 0.0017952114 0.13476229 0.5787751 1.6393526
#> dev.o[13,3] 1.0430266 1.4990127 0.0010257036 0.10869586 0.4698240 1.3769388
#> dev.o[19,3] 1.2945553 1.1420883 0.0091498935 0.43088186 1.0269545 1.8518900
#> dev.o[10,4] 1.0914692 1.3420188 0.0015454724 0.14138617 0.5846418 1.5611129
#> dev.o[12,4] 0.8287740 1.1438309 0.0006723164 0.08920458 0.3742728 1.0836134
#> dev.o[13,4] 1.1223126 1.4919084 0.0008702342 0.12305765 0.5128668 1.5621390
#>                97.5%     Rhat n.eff
#> dev.o[1,1]  8.570850 1.002034  1300
#> dev.o[2,1]  4.869970 1.006262   350
#> dev.o[3,1]  5.560072 1.000623  3000
#> dev.o[4,1]  3.619871 1.002762  1300
#> dev.o[5,1]  3.702596 1.002099  3000
#> dev.o[6,1]  4.291187 1.003137   760
#> dev.o[7,1]  3.919780 1.003337   800
#> dev.o[8,1]  4.020055 1.001699  1700
#> dev.o[9,1]  3.906149 1.000569  3000
#> dev.o[10,1] 3.061910 1.001450  2100
#> dev.o[11,1] 4.906873 1.000831  3000
#> dev.o[12,1] 5.616254 1.009133   270
#> dev.o[13,1] 6.121272 1.000617  3000
#> dev.o[14,1] 4.524863 1.001343  2400
#> dev.o[15,1] 3.935522 1.002165  1600
#> dev.o[16,1] 5.996497 1.004031   630
#> dev.o[17,1] 6.989381 1.007276   330
#> dev.o[18,1] 4.811063 1.001833  2100
#> dev.o[19,1] 6.257214 1.001612  2500
#> dev.o[20,1] 3.888415 1.000683  3000
#> dev.o[21,1] 5.178243 1.000879  3000
#> dev.o[1,2]  8.108033 1.003116   770
#> dev.o[2,2]  4.618727 1.003943   580
#> dev.o[3,2]  4.970592 1.001647  1700
#> dev.o[4,2]  4.075323 1.000942  3000
#> dev.o[5,2]  3.088382 1.000674  3000
#> dev.o[6,2]  4.663466 1.000971  3000
#> dev.o[7,2]  4.182459 1.000711  3000
#> dev.o[8,2]  3.715050 1.006258   350
#> dev.o[9,2]  3.649485 1.004360   960
#> dev.o[10,2] 6.687557 1.006888   540
#> dev.o[11,2] 5.018996 1.007412   860
#> dev.o[12,2] 3.902024 1.000845  3000
#> dev.o[13,2] 5.051275 1.001787  3000
#> dev.o[14,2] 3.955349 1.000967  3000
#> dev.o[15,2] 4.520155 1.002919  1500
#> dev.o[16,2] 5.998561 1.002699   910
#> dev.o[17,2] 6.535571 1.009751   260
#> dev.o[18,2] 4.436317 1.000718  3000
#> dev.o[19,2] 2.076721 1.000718  3000
#> dev.o[20,2] 3.482819 1.000679  3000
#> dev.o[21,2] 4.275800 1.002594   960
#> dev.o[9,3]  3.640744 1.003111  1200
#> dev.o[10,3] 3.587085 1.001177  2900
#> dev.o[12,3] 5.691489 1.004141   550
#> dev.o[13,3] 5.318047 1.001506  2000
#> dev.o[19,3] 4.132635 1.003240  1500
#> dev.o[10,4] 4.627563 1.003228  1100
#> dev.o[12,4] 4.185791 1.002376  1400
#> dev.o[13,4] 5.188098 1.003110   830
#> 
#> $hat_par
#>                     mean         sd        2.5%         25%        50%
#> hat.par[1,1]    1.561427  0.7706284   0.3662796   0.9849925   1.462782
#> hat.par[2,1]   49.203802  4.7988181  39.6905431  45.9969588  49.180796
#> hat.par[3,1]   43.779836  4.6742097  34.9290590  40.4442186  43.669848
#> hat.par[4,1]   42.219089  4.6626781  33.6511434  38.8618484  42.033077
#> hat.par[5,1]   17.095147  2.5365116  12.2529891  15.3430732  17.006142
#> hat.par[6,1]   43.554272  4.1941302  35.3620110  40.7662234  43.571073
#> hat.par[7,1]  157.555705  7.5015190 143.1114568 152.4415751 157.476214
#> hat.par[8,1]   68.456949  5.5462697  57.4165276  64.7962097  68.332694
#> hat.par[9,1]   90.134388  5.1540080  79.9729570  86.5557915  90.276751
#> hat.par[10,1]  78.867790  3.7806326  71.3883614  76.3286471  78.914891
#> hat.par[11,1]  72.142410  5.5169288  61.3780249  68.4094305  72.086337
#> hat.par[12,1]  77.470999  4.2396341  69.1828405  74.6863878  77.435378
#> hat.par[13,1]  49.364823  4.7499948  40.1054279  46.0915804  49.277263
#> hat.par[14,1]  33.892331  4.6516388  25.3103659  30.6737642  33.689718
#> hat.par[15,1]  36.082603  4.5375572  27.4390838  33.0537627  36.017623
#> hat.par[16,1] 303.330535 13.2421766 277.9510559 294.3403075 303.056981
#> hat.par[17,1]  11.238427  2.5899255   6.6798941   9.3327374  11.086406
#> hat.par[18,1]  22.134921  3.2933754  15.9423755  19.8050800  22.124818
#> hat.par[19,1]   4.077727  1.3909076   1.8223568   3.0496352   3.935959
#> hat.par[20,1]  25.275856  3.6616339  18.5829993  22.7199312  25.068166
#> hat.par[21,1]  32.236698  4.3563327  24.1356712  29.2640778  32.074563
#> hat.par[1,2]    1.425590  0.7469654   0.3176340   0.8684717   1.316024
#> hat.par[2,2]   46.857132  4.9823619  37.6070527  43.3596085  46.676111
#> hat.par[3,2]   30.900562  4.0742714  23.1049232  28.0988853  30.882236
#> hat.par[4,2]   43.716867  5.2791948  33.9058211  40.1326388  43.740594
#> hat.par[5,2]   11.858009  2.1029943   7.9066821  10.3831789  11.785630
#> hat.par[6,2]   35.211722  3.9534026  27.8727277  32.5520204  35.114953
#> hat.par[7,2]  197.033154  9.9884803 177.9326681 190.1411063 196.990858
#> hat.par[8,2]   51.374094  5.0527469  41.9388585  48.0062976  51.228983
#> hat.par[9,2]   81.854852  5.7235351  70.8452322  77.9534514  81.768805
#> hat.par[10,2]  72.988705  3.9551303  65.0851469  70.3369025  73.096502
#> hat.par[11,2] 119.842055  8.2097404 104.1621520 114.2698024 119.777866
#> hat.par[12,2]  80.503147  4.8613769  70.9894048  77.0843591  80.602595
#> hat.par[13,2]  26.037456  4.4891101  17.6747573  23.0123477  25.918713
#> hat.par[14,2]  31.721757  4.4013701  23.5255550  28.7105885  31.528781
#> hat.par[15,2]  32.980432  4.5564565  24.8170795  29.7154039  32.697876
#> hat.par[16,2] 247.748870 12.4702544 223.4233340 239.3315416 247.618189
#> hat.par[17,2]   6.906834  1.8403456   3.7175309   5.6106079   6.765455
#> hat.par[18,2]  12.821288  2.4987651   8.2650180  11.0448742  12.742070
#> hat.par[19,2]   2.440677  0.9327663   1.0083055   1.7643408   2.323747
#> hat.par[20,2]  18.776538  3.0975817  13.3002870  16.5601251  18.615084
#> hat.par[21,2]  21.551418  3.4890502  15.0692225  19.1212768  21.498364
#> hat.par[9,3]   78.948658  5.6727163  67.5725853  75.1683037  78.990541
#> hat.par[10,3]  68.289725  4.3305752  59.8793483  65.3552249  68.312066
#> hat.par[12,3]  67.336957  4.8536411  57.7980055  64.0922778  67.402768
#> hat.par[13,3]  35.334099  5.2614737  25.2408058  31.7598985  35.197067
#> hat.par[19,3]   2.443813  0.9042399   0.9882637   1.7941055   2.352753
#> hat.par[10,4]  66.419634  4.3064497  57.4786199  63.6186894  66.522341
#> hat.par[12,4]  62.897842  4.5043942  53.8659665  59.8132395  62.851996
#> hat.par[13,4]  41.131911  4.8086276  31.9125268  37.8365558  40.976520
#>                      75%      97.5%     Rhat n.eff
#> hat.par[1,1]    2.052004   3.309770 1.002595  1000
#> hat.par[2,1]   52.438929  58.585048 1.031741    71
#> hat.par[3,1]   46.892831  53.086099 1.008328   340
#> hat.par[4,1]   45.318859  51.897281 1.000510  3000
#> hat.par[5,1]   18.781952  22.369368 1.000883  3000
#> hat.par[6,1]   46.432291  51.594127 1.006253   350
#> hat.par[7,1]  162.683379 172.089264 1.006981   330
#> hat.par[8,1]   72.090005  79.500692 1.010525   280
#> hat.par[9,1]   93.564816  99.991789 1.008714   250
#> hat.par[10,1]  81.571696  86.032703 1.008390   390
#> hat.par[11,1]  75.741863  83.440718 1.008820   240
#> hat.par[12,1]  80.340190  85.695683 1.024591   100
#> hat.par[13,1]  52.634593  58.648238 1.004430   590
#> hat.par[14,1]  36.957067  43.901549 1.010549   200
#> hat.par[15,1]  39.106583  44.919288 1.012693   190
#> hat.par[16,1] 312.366411 329.839983 1.001908  1500
#> hat.par[17,1]  13.022499  16.637254 1.007958   270
#> hat.par[18,1]  24.229119  28.769954 1.001455  2600
#> hat.par[19,1]   4.962324   7.167159 1.001122  3000
#> hat.par[20,1]  27.733846  33.020230 1.003570   650
#> hat.par[21,1]  35.038419  41.064353 1.001607  3000
#> hat.par[1,2]    1.874509   3.180407 1.003049   790
#> hat.par[2,2]   50.293579  56.638008 1.021942    95
#> hat.par[3,2]   33.549897  39.201511 1.002100  2500
#> hat.par[4,2]   47.104536  54.682661 1.000748  3000
#> hat.par[5,2]   13.224854  16.197280 1.001552  1900
#> hat.par[6,2]   37.811547  43.355891 1.001526  1900
#> hat.par[7,2]  203.918289 216.637310 1.004975   450
#> hat.par[8,2]   54.715578  61.446808 1.002531  3000
#> hat.par[9,2]   85.762117  93.307847 1.007336   340
#> hat.par[10,2]  75.741740  80.504196 1.001808  1500
#> hat.par[11,2] 125.227628 136.454114 1.007865   270
#> hat.par[12,2]  83.925862  89.716494 1.009531   220
#> hat.par[13,2]  29.037634  34.954641 1.005252   420
#> hat.par[14,2]  34.495447  41.137005 1.004508   500
#> hat.par[15,2]  35.890027  42.472098 1.004295   530
#> hat.par[16,2] 256.055134 271.808399 1.005475   400
#> hat.par[17,2]   8.095055  10.877710 1.012113   170
#> hat.par[18,2]  14.389988  18.085557 1.001142  3000
#> hat.par[19,2]   2.972008   4.634375 1.001065  3000
#> hat.par[20,2]  20.779111  25.232063 1.001412  2200
#> hat.par[21,2]  23.838917  28.571842 1.002767   880
#> hat.par[9,3]   82.832588  89.448284 1.008218   290
#> hat.par[10,3]  71.326035  76.579854 1.001721  1600
#> hat.par[12,3]  70.604059  76.802308 1.005674   390
#> hat.par[13,3]  38.839963  46.039759 1.004772   660
#> hat.par[19,3]   2.983798   4.450476 1.002234  1200
#> hat.par[10,4]  69.412314  74.210121 1.006360   350
#> hat.par[12,4]  66.145396  71.417867 1.001808  2900
#> hat.par[13,4]  44.338001  50.573515 1.008332   260
#> 
#> $leverage_o
#>  [1] 0.9178643 0.8510408 0.8171146 0.6693542 0.6588168 0.6805782 0.7736544
#>  [8] 0.7306590 0.6557525 0.6162903 0.7746729 0.6283291 0.8233217 0.7669657
#> [15] 0.6931867 0.9208495 0.8763151 0.7690094 0.7616037 0.6661642 0.7866836
#> [22] 0.2358243 0.8499928 0.6775364 0.7839729 0.5193755 0.6851417 0.8571872
#> [29] 0.7027378 0.6987791 0.7308802 0.8727649 0.7477496 0.9789499 0.6539649
#> [36] 0.8005483 0.9077653 0.3501212 0.5679223 0.3082090 0.6439792 0.5865080
#> [43] 0.6746835 0.7024714 0.7587432 1.0388765 0.1568469 0.6770514 0.6809297
#> [50] 0.7763250
#> 
#> $sign_dev_o
#>  [1]  1  1  1 -1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1  1  1  1 -1  1 -1 -1 -1  1
#> [26] -1  1  1  1  1  1 -1 -1 -1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1
#> 
#> $phi
#>                mean        sd       2.5%        25%          50%         75%
#> phi[1] -0.295801708 0.6039683 -1.5946569 -0.6849308 -0.263458934  0.12010915
#> phi[2]  0.006954976 0.9595987 -1.9020168 -0.6191357 -0.007567471  0.65344284
#> phi[3]  0.069208857 0.9653974 -1.8443533 -0.5688018  0.093231023  0.72894653
#> phi[4] -0.950775923 0.9020072 -2.5273057 -1.5918959 -1.037160043 -0.38462971
#> phi[5] -0.643154747 0.9475764 -2.4771181 -1.2371010 -0.661152659 -0.05066862
#> phi[6]  0.923465106 0.7818165 -0.8132052  0.4560011  0.950357388  1.42622205
#> phi[7] -0.305660490 0.6711014 -1.6982961 -0.7457354 -0.284004741  0.15086061
#> phi[8] -0.271627019 0.8985750 -2.0415904 -0.8880262 -0.269958736  0.32515615
#>            97.5%     Rhat n.eff
#> phi[1] 0.7794479 1.072094    73
#> phi[2] 1.9341151 1.020814   110
#> phi[3] 1.8983370 1.014165   160
#> phi[4] 1.0240063 1.085902    29
#> phi[5] 1.3442365 1.006265  1900
#> phi[6] 2.4332065 1.084614    30
#> phi[7] 0.9541814 1.001147  3000
#> phi[8] 1.5542673 1.012024   180
#> 
#> $model_assessment
#>       DIC       pD      dev n_data
#> 1 88.1037 35.46406 52.63963     50
#> 
#> $data
#>                    study t1 t2 t3 t4  r1  r2 r3 r4  m1 m2 m3 m4  n1  n2  n3  n4
#> 1  Llewellyn-Jones, 1996  1  4 NA NA   3   0 NA NA   1  0 NA NA   8   8  NA  NA
#> 2         Paggiaro, 1998  1  4 NA NA  51  45 NA NA  27 19 NA NA 139 142  NA  NA
#> 3           Mahler, 1999  1  7 NA NA  47  28 NA NA  23  9 NA NA 143 135  NA  NA
#> 4         Casaburi, 2000  1  8 NA NA  41  45 NA NA  18 12 NA NA 191 279  NA  NA
#> 5        van Noord, 2000  1  7 NA NA  18  11 NA NA   8  7 NA NA  50  47  NA  NA
#> 6          Rennard, 2001  1  7 NA NA  41  38 NA NA  29 22 NA NA 135 132  NA  NA
#> 7         Casaburi, 2002  1  8 NA NA 156 198 NA NA  77 66 NA NA 371 550  NA  NA
#> 8          Chapman, 2002  1  7 NA NA  68  52 NA NA  28 20 NA NA 207 201  NA  NA
#> 9          Donohue, 2002  1  7  8 NA  92  82 77 NA  37 20 10 NA 201 213 209  NA
#> 10          Mahler, 2002  1  4  7  5  79  77 63 68  69 68 45 52 181 168 160 165
#> 11           Rossi, 2002  1  6 NA NA  75 117 NA NA  59 92 NA NA 220 425  NA  NA
#> 12         Hanania, 2003  1  4  7  5  73  79 65 71  59 49 57 53 185 183 177 178
#> 13      Szafranski, 2003  1  2  6  3  53  26 38 35  90 62 64 59 205 198 201 208
#> 14          Briggs, 2005  8  7 NA NA  30  36 NA NA  29 41 NA NA 328 325  NA  NA
#> 15        Campbell, 2005  1  6 NA NA  34  35 NA NA  39 30 NA NA 217 215  NA  NA
#> 16      Niewoehner, 2005  1  8 NA NA 296 255 NA NA 111 75 NA NA 915 914  NA  NA
#> 17       van Noord, 2005  8  6 NA NA   4  14 NA NA   1  1 NA NA  70  69  NA  NA
#> 18          Barnes, 2006  1  5 NA NA  24  11 NA NA   4  8 NA NA  73  67  NA  NA
#> 19       O Donnell, 2006  1  7  5 NA   6   1  2 NA   5  1  3 NA  64  59  62  NA
#> 20     Baumgartner, 2007  1  7 NA NA  24  20 NA NA  32 26 NA NA 143 144  NA  NA
#> 21         Freeman, 2007  1  8 NA NA  35  19 NA NA  33 18 NA NA 195 200  NA  NA
#> 
#> $measure
#> [1] "OR"
#> 
#> $model
#> [1] "RE"
#> 
#> $assumption
#> [1] "IDE-ARM"
#> 
#> $mean_misspar
#> [1] 1e-04 1e-04
#> 
#> $var_misspar
#> [1] 1
#> 
#> $D
#> [1] 0
#> 
#> $ref
#> [1] 1
#> 
#> $indic
#>       [,1] [,2] [,3] [,4]
#>  [1,]    1    1   NA   NA
#>  [2,]    1    1   NA   NA
#>  [3,]    1    1   NA   NA
#>  [4,]    1    1   NA   NA
#>  [5,]    1    1   NA   NA
#>  [6,]    1    1   NA   NA
#>  [7,]    1    1   NA   NA
#>  [8,]    1    1   NA   NA
#>  [9,]    1    1    1   NA
#> [10,]    1    1    1    1
#> [11,]    1    1   NA   NA
#> [12,]    1    1    1    1
#> [13,]    1    1    1    1
#> [14,]    1    1   NA   NA
#> [15,]    1    1   NA   NA
#> [16,]    1    1   NA   NA
#> [17,]    1    1   NA   NA
#> [18,]    1    1   NA   NA
#> [19,]    1    1    1   NA
#> [20,]    1    1   NA   NA
#> [21,]    1    1   NA   NA
#> 
#> $jagsfit
#> Inference for Bugs model at "5", fit using jags,
#>  3 chains, each with 1000 iterations (first 0 discarded)
#>  n.sims = 3000 iterations saved. Running time = secs
#>                    mu.vect sd.vect    2.5%     25%     50%     75%   97.5%
#> EM[2,1]             -0.943   0.428  -1.770  -1.228  -0.954  -0.672  -0.050
#> EM[3,1]             -0.661   0.433  -1.577  -0.939  -0.640  -0.376   0.132
#> EM[4,1]             -0.244   0.325  -0.798  -0.476  -0.275  -0.021   0.438
#> EM[5,1]             -0.513   0.283  -1.052  -0.686  -0.535  -0.347   0.130
#> EM[6,1]             -0.023   0.248  -0.552  -0.190  -0.001   0.166   0.397
#> EM[7,1]             -0.464   0.182  -0.796  -0.595  -0.476  -0.341  -0.087
#> EM[8,1]             -0.530   0.164  -0.856  -0.636  -0.529  -0.413  -0.218
#> EM[3,2]              0.282   0.485  -0.664  -0.054   0.281   0.630   1.222
#> EM[4,2]              0.700   0.462  -0.239   0.393   0.675   1.018   1.632
#> EM[5,2]              0.431   0.472  -0.545   0.125   0.421   0.726   1.376
#> EM[6,2]              0.920   0.449  -0.010   0.645   0.939   1.220   1.733
#> EM[7,2]              0.479   0.436  -0.368   0.194   0.474   0.765   1.314
#> EM[8,2]              0.413   0.433  -0.480   0.141   0.431   0.717   1.199
#> EM[4,3]              0.417   0.464  -0.470   0.118   0.403   0.694   1.392
#> EM[5,3]              0.148   0.473  -0.803  -0.164   0.149   0.465   1.119
#> EM[6,3]              0.638   0.440  -0.244   0.332   0.662   0.926   1.539
#> EM[7,3]              0.196   0.411  -0.583  -0.074   0.212   0.464   1.016
#> EM[8,3]              0.130   0.420  -0.674  -0.150   0.141   0.408   0.976
#> EM[5,4]             -0.269   0.377  -0.999  -0.523  -0.271  -0.016   0.484
#> EM[6,4]              0.220   0.350  -0.464  -0.016   0.229   0.467   0.826
#> EM[7,4]             -0.221   0.299  -0.845  -0.414  -0.203   0.002   0.313
#> EM[8,4]             -0.287   0.320  -0.981  -0.494  -0.258  -0.043   0.250
#> EM[6,5]              0.490   0.335  -0.242   0.287   0.506   0.717   1.116
#> EM[7,5]              0.048   0.295  -0.568  -0.147   0.045   0.239   0.597
#> EM[8,5]             -0.018   0.289  -0.630  -0.190  -0.008   0.179   0.502
#> EM[7,6]             -0.441   0.240  -0.886  -0.610  -0.448  -0.284   0.047
#> EM[8,6]             -0.507   0.230  -0.960  -0.654  -0.513  -0.353  -0.051
#> EM[8,7]             -0.066   0.164  -0.387  -0.174  -0.060   0.046   0.244
#> EM.pred[2,1]        -0.934   0.470  -1.861  -1.242  -0.946  -0.630   0.032
#> EM.pred[3,1]        -0.654   0.469  -1.635  -0.939  -0.636  -0.359   0.242
#> EM.pred[4,1]        -0.243   0.378  -0.926  -0.501  -0.256   0.004   0.506
#> EM.pred[5,1]        -0.512   0.345  -1.223  -0.718  -0.518  -0.313   0.206
#> EM.pred[6,1]        -0.025   0.317  -0.702  -0.229   0.006   0.194   0.544
#> EM.pred[7,1]        -0.464   0.261  -0.974  -0.636  -0.464  -0.293   0.028
#> EM.pred[8,1]        -0.534   0.256  -1.083  -0.682  -0.519  -0.373  -0.052
#> EM.pred[3,2]         0.282   0.523  -0.757  -0.077   0.280   0.651   1.286
#> EM.pred[4,2]         0.699   0.503  -0.312   0.359   0.677   1.044   1.681
#> EM.pred[5,2]         0.428   0.511  -0.658   0.090   0.428   0.754   1.439
#> EM.pred[6,2]         0.920   0.490  -0.112   0.611   0.940   1.252   1.827
#> EM.pred[7,2]         0.475   0.477  -0.468   0.185   0.472   0.791   1.381
#> EM.pred[8,2]         0.413   0.475  -0.560   0.111   0.435   0.732   1.299
#> EM.pred[4,3]         0.422   0.500  -0.569   0.105   0.422   0.709   1.512
#> EM.pred[5,3]         0.154   0.505  -0.843  -0.178   0.156   0.485   1.179
#> EM.pred[6,3]         0.638   0.479  -0.312   0.310   0.661   0.946   1.587
#> EM.pred[7,3]         0.195   0.457  -0.704  -0.095   0.206   0.481   1.133
#> EM.pred[8,3]         0.125   0.459  -0.796  -0.164   0.130   0.413   1.052
#> EM.pred[5,4]        -0.261   0.418  -1.082  -0.536  -0.255   0.016   0.569
#> EM.pred[6,4]         0.226   0.397  -0.562  -0.031   0.233   0.504   0.956
#> EM.pred[7,4]        -0.229   0.360  -0.983  -0.453  -0.217   0.020   0.461
#> EM.pred[8,4]        -0.289   0.369  -1.072  -0.521  -0.265  -0.019   0.350
#> EM.pred[6,5]         0.484   0.382  -0.361   0.253   0.507   0.733   1.210
#> EM.pred[7,5]         0.047   0.344  -0.648  -0.171   0.051   0.264   0.717
#> EM.pred[8,5]        -0.013   0.346  -0.729  -0.222   0.002   0.210   0.645
#> EM.pred[7,6]        -0.447   0.308  -1.043  -0.645  -0.458  -0.253   0.186
#> EM.pred[8,6]        -0.510   0.306  -1.135  -0.701  -0.519  -0.307   0.109
#> EM.pred[8,7]        -0.064   0.245  -0.557  -0.204  -0.062   0.085   0.411
#> SUCRA[1]             0.123   0.120   0.000   0.000   0.143   0.143   0.429
#> SUCRA[2]             0.879   0.208   0.286   0.857   1.000   1.000   1.000
#> SUCRA[3]             0.703   0.276   0.000   0.571   0.857   0.857   1.000
#> SUCRA[4]             0.346   0.242   0.000   0.143   0.286   0.429   0.857
#> SUCRA[5]             0.603   0.250   0.143   0.429   0.571   0.857   1.000
#> SUCRA[6]             0.143   0.161   0.000   0.000   0.143   0.286   0.571
#> SUCRA[7]             0.561   0.176   0.286   0.429   0.571   0.714   0.857
#> SUCRA[8]             0.642   0.182   0.286   0.571   0.714   0.714   1.000
#> abs_risk[1]          0.392   0.000   0.392   0.392   0.392   0.392   0.392
#> abs_risk[2]          0.209   0.071   0.099   0.159   0.199   0.247   0.380
#> abs_risk[3]          0.258   0.079   0.117   0.201   0.253   0.307   0.424
#> abs_risk[4]          0.339   0.073   0.225   0.286   0.328   0.387   0.499
#> abs_risk[5]          0.282   0.057   0.184   0.245   0.274   0.313   0.423
#> abs_risk[6]          0.388   0.057   0.271   0.347   0.391   0.432   0.489
#> abs_risk[7]          0.290   0.037   0.225   0.262   0.286   0.314   0.371
#> abs_risk[8]          0.276   0.032   0.215   0.254   0.275   0.299   0.341
#> delta[1,1]           0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[2,1]           0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[3,1]           0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[4,1]           0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[5,1]           0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[6,1]           0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[7,1]           0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[8,1]           0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[9,1]           0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[10,1]          0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[11,1]          0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[12,1]          0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[13,1]          0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[14,1]          0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[15,1]          0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[16,1]          0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[17,1]          0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[18,1]          0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[19,1]          0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[20,1]          0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[21,1]          0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> delta[1,2]          -0.304   0.371  -0.980  -0.560  -0.320  -0.062   0.438
#> delta[2,2]          -0.285   0.302  -0.813  -0.506  -0.306  -0.080   0.337
#> delta[3,2]          -0.547   0.229  -1.007  -0.694  -0.553  -0.383  -0.115
#> delta[4,2]          -0.484   0.182  -0.861  -0.606  -0.478  -0.359  -0.149
#> delta[5,2]          -0.486   0.248  -0.959  -0.653  -0.491  -0.316  -0.021
#> delta[6,2]          -0.381   0.221  -0.787  -0.539  -0.388  -0.241   0.061
#> delta[7,2]          -0.486   0.176  -0.830  -0.607  -0.482  -0.370  -0.145
#> delta[8,2]          -0.437   0.191  -0.819  -0.570  -0.438  -0.309  -0.064
#> delta[9,2]          -0.477   0.201  -0.849  -0.614  -0.484  -0.350  -0.075
#> delta[10,2]         -0.165   0.373  -0.805  -0.429  -0.193   0.069   0.671
#> delta[11,2]         -0.097   0.258  -0.614  -0.276  -0.083   0.102   0.335
#> delta[12,2]         -0.198   0.332  -0.768  -0.437  -0.228   0.029   0.482
#> delta[13,2]         -1.002   0.409  -1.831  -1.262  -1.008  -0.745  -0.166
#> delta[14,2]         -0.114   0.196  -0.526  -0.242  -0.105   0.017   0.251
#> delta[15,2]          0.043   0.244  -0.466  -0.118   0.053   0.210   0.504
#> delta[16,2]         -0.351   0.135  -0.608  -0.442  -0.353  -0.262  -0.080
#> delta[17,2]         -0.620   0.271  -1.203  -0.773  -0.613  -0.439  -0.121
#> delta[18,2]         -0.573   0.302  -1.169  -0.759  -0.573  -0.398   0.062
#> delta[19,2]         -0.554   0.332  -1.178  -0.762  -0.558  -0.356   0.130
#> delta[20,2]         -0.434   0.224  -0.870  -0.592  -0.433  -0.282  -0.009
#> delta[21,2]         -0.587   0.215  -1.038  -0.722  -0.579  -0.439  -0.196
#> delta[9,3]          -0.576   0.191  -0.961  -0.702  -0.574  -0.445  -0.197
#> delta[10,3]         -0.516   0.312  -1.116  -0.708  -0.527  -0.341   0.194
#> delta[12,3]         -0.395   0.302  -0.934  -0.585  -0.421  -0.230   0.277
#> delta[13,3]         -0.718   0.411  -1.628  -0.982  -0.689  -0.445   0.027
#> delta[19,3]         -0.522   0.268  -1.084  -0.679  -0.512  -0.346  -0.044
#> delta[10,4]         -0.516   0.248  -1.033  -0.669  -0.514  -0.343  -0.058
#> delta[12,4]         -0.371   0.235  -0.809  -0.534  -0.386  -0.220   0.140
#> delta[13,4]         -0.135   0.308  -0.784  -0.344  -0.097   0.092   0.364
#> dev.o[1,1]           2.378   2.336   0.007   0.596   1.713   3.419   8.571
#> dev.o[2,1]           0.968   1.346   0.001   0.101   0.448   1.280   4.870
#> dev.o[3,1]           1.186   1.544   0.002   0.139   0.619   1.690   5.560
#> dev.o[4,1]           0.716   1.003   0.001   0.077   0.326   0.899   3.620
#> dev.o[5,1]           0.739   1.022   0.001   0.080   0.341   1.001   3.703
#> dev.o[6,1]           0.936   1.274   0.001   0.110   0.469   1.274   4.291
#> dev.o[7,1]           0.807   1.136   0.001   0.080   0.381   1.074   3.920
#> dev.o[8,1]           0.736   1.082   0.001   0.071   0.322   0.940   4.020
#> dev.o[9,1]           0.742   1.084   0.001   0.069   0.335   0.953   3.906
#> dev.o[10,1]          0.617   0.838   0.001   0.069   0.295   0.832   3.062
#> dev.o[11,1]          0.979   1.322   0.001   0.106   0.440   1.335   4.907
#> dev.o[12,1]          1.291   1.588   0.002   0.185   0.714   1.831   5.616
#> dev.o[13,1]          1.290   1.664   0.002   0.141   0.651   1.794   6.121
#> dev.o[14,1]          0.913   1.246   0.001   0.102   0.423   1.226   4.525
#> dev.o[15,1]          0.846   1.126   0.001   0.085   0.408   1.161   3.936
#> dev.o[16,1]          1.206   1.632   0.001   0.123   0.557   1.633   5.996
#> dev.o[17,1]          1.642   1.945   0.003   0.223   0.911   2.418   6.989
#> dev.o[18,1]          0.997   1.422   0.001   0.099   0.451   1.359   4.811
#> dev.o[19,1]          1.620   1.803   0.004   0.271   1.022   2.379   6.257
#> dev.o[20,1]          0.751   1.093   0.001   0.074   0.347   0.985   3.888
#> dev.o[21,1]          1.076   1.428   0.001   0.111   0.532   1.464   5.178
#> dev.o[1,2]           3.253   1.946   0.648   1.839   2.876   4.272   8.108
#> dev.o[2,2]           0.970   1.308   0.001   0.096   0.485   1.359   4.619
#> dev.o[3,2]           1.046   1.347   0.001   0.116   0.525   1.426   4.971
#> dev.o[4,2]           0.829   1.210   0.001   0.072   0.362   1.091   4.075
#> dev.o[5,2]           0.609   0.887   0.001   0.062   0.270   0.800   3.088
#> dev.o[6,2]           1.006   1.336   0.001   0.108   0.493   1.392   4.663
#> dev.o[7,2]           0.865   1.186   0.001   0.087   0.424   1.151   4.182
#> dev.o[8,2]           0.713   1.012   0.000   0.075   0.309   0.972   3.715
#> dev.o[9,2]           0.699   1.040   0.001   0.075   0.321   0.869   3.650
#> dev.o[10,2]          1.575   1.869   0.003   0.221   0.903   2.253   6.688
#> dev.o[11,2]          0.978   1.401   0.001   0.101   0.452   1.281   5.019
#> dev.o[12,2]          0.818   1.113   0.001   0.092   0.374   1.075   3.902
#> dev.o[13,2]          0.979   1.384   0.001   0.104   0.438   1.293   5.051
#> dev.o[14,2]          0.760   1.088   0.001   0.077   0.334   0.992   3.955
#> dev.o[15,2]          0.949   1.239   0.001   0.107   0.481   1.315   4.520
#> dev.o[16,2]          1.207   1.710   0.001   0.125   0.561   1.601   5.999
#> dev.o[17,2]          1.928   1.790   0.010   0.565   1.448   2.820   6.536
#> dev.o[18,2]          0.911   1.266   0.001   0.101   0.437   1.222   4.436
#> dev.o[19,2]          0.396   0.586   0.000   0.040   0.175   0.506   2.077
#> dev.o[20,2]          0.737   1.002   0.001   0.074   0.350   1.000   3.483
#> dev.o[21,2]          0.942   1.199   0.001   0.121   0.486   1.291   4.276
#> dev.o[9,3]           0.755   1.034   0.001   0.078   0.353   1.051   3.641
#> dev.o[10,3]          0.706   1.002   0.001   0.070   0.333   0.918   3.587
#> dev.o[12,3]          1.192   1.596   0.002   0.135   0.579   1.639   5.691
#> dev.o[13,3]          1.043   1.499   0.001   0.109   0.470   1.377   5.318
#> dev.o[19,3]          1.295   1.142   0.009   0.431   1.027   1.852   4.133
#> dev.o[10,4]          1.091   1.342   0.002   0.141   0.585   1.561   4.628
#> dev.o[12,4]          0.829   1.144   0.001   0.089   0.374   1.084   4.186
#> dev.o[13,4]          1.122   1.492   0.001   0.123   0.513   1.562   5.188
#> effectiveness[1,1]   0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> effectiveness[2,1]   0.611   0.488   0.000   0.000   1.000   1.000   1.000
#> effectiveness[3,1]   0.224   0.417   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,1]   0.008   0.091   0.000   0.000   0.000   0.000   0.000
#> effectiveness[5,1]   0.080   0.272   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,1]   0.001   0.036   0.000   0.000   0.000   0.000   0.000
#> effectiveness[7,1]   0.017   0.129   0.000   0.000   0.000   0.000   0.000
#> effectiveness[8,1]   0.058   0.234   0.000   0.000   0.000   0.000   1.000
#> effectiveness[1,2]   0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> effectiveness[2,2]   0.210   0.407   0.000   0.000   0.000   0.000   1.000
#> effectiveness[3,2]   0.305   0.460   0.000   0.000   0.000   1.000   1.000
#> effectiveness[4,2]   0.044   0.206   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,2]   0.197   0.398   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,2]   0.002   0.041   0.000   0.000   0.000   0.000   0.000
#> effectiveness[7,2]   0.086   0.280   0.000   0.000   0.000   0.000   1.000
#> effectiveness[8,2]   0.157   0.364   0.000   0.000   0.000   0.000   1.000
#> effectiveness[1,3]   0.001   0.032   0.000   0.000   0.000   0.000   0.000
#> effectiveness[2,3]   0.055   0.229   0.000   0.000   0.000   0.000   1.000
#> effectiveness[3,3]   0.122   0.327   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,3]   0.089   0.284   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,3]   0.220   0.414   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,3]   0.007   0.085   0.000   0.000   0.000   0.000   0.000
#> effectiveness[7,3]   0.217   0.412   0.000   0.000   0.000   0.000   1.000
#> effectiveness[8,3]   0.289   0.454   0.000   0.000   0.000   1.000   1.000
#> effectiveness[1,4]   0.003   0.052   0.000   0.000   0.000   0.000   0.000
#> effectiveness[2,4]   0.037   0.190   0.000   0.000   0.000   0.000   1.000
#> effectiveness[3,4]   0.106   0.308   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,4]   0.104   0.305   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,4]   0.145   0.352   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,4]   0.024   0.153   0.000   0.000   0.000   0.000   0.000
#> effectiveness[7,4]   0.305   0.460   0.000   0.000   0.000   1.000   1.000
#> effectiveness[8,4]   0.275   0.447   0.000   0.000   0.000   1.000   1.000
#> effectiveness[1,5]   0.034   0.182   0.000   0.000   0.000   0.000   1.000
#> effectiveness[2,5]   0.047   0.212   0.000   0.000   0.000   0.000   1.000
#> effectiveness[3,5]   0.104   0.305   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,5]   0.175   0.380   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,5]   0.165   0.372   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,5]   0.068   0.251   0.000   0.000   0.000   0.000   1.000
#> effectiveness[7,5]   0.248   0.432   0.000   0.000   0.000   0.000   1.000
#> effectiveness[8,5]   0.159   0.365   0.000   0.000   0.000   0.000   1.000
#> effectiveness[1,6]   0.166   0.372   0.000   0.000   0.000   0.000   1.000
#> effectiveness[2,6]   0.019   0.138   0.000   0.000   0.000   0.000   0.000
#> effectiveness[3,6]   0.069   0.254   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,6]   0.287   0.452   0.000   0.000   0.000   1.000   1.000
#> effectiveness[5,6]   0.122   0.327   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,6]   0.159   0.366   0.000   0.000   0.000   0.000   1.000
#> effectiveness[7,6]   0.119   0.324   0.000   0.000   0.000   0.000   1.000
#> effectiveness[8,6]   0.059   0.235   0.000   0.000   0.000   0.000   1.000
#> effectiveness[1,7]   0.411   0.492   0.000   0.000   0.000   1.000   1.000
#> effectiveness[2,7]   0.013   0.115   0.000   0.000   0.000   0.000   0.000
#> effectiveness[3,7]   0.041   0.198   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,7]   0.141   0.348   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,7]   0.057   0.232   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,7]   0.326   0.469   0.000   0.000   0.000   1.000   1.000
#> effectiveness[7,7]   0.008   0.089   0.000   0.000   0.000   0.000   0.000
#> effectiveness[8,7]   0.003   0.052   0.000   0.000   0.000   0.000   0.000
#> effectiveness[1,8]   0.385   0.487   0.000   0.000   0.000   1.000   1.000
#> effectiveness[2,8]   0.007   0.083   0.000   0.000   0.000   0.000   0.000
#> effectiveness[3,8]   0.029   0.169   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,8]   0.152   0.359   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,8]   0.014   0.116   0.000   0.000   0.000   0.000   0.000
#> effectiveness[6,8]   0.413   0.492   0.000   0.000   0.000   1.000   1.000
#> effectiveness[7,8]   0.000   0.018   0.000   0.000   0.000   0.000   0.000
#> effectiveness[8,8]   0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> hat.par[1,1]         1.561   0.771   0.366   0.985   1.463   2.052   3.310
#> hat.par[2,1]        49.204   4.799  39.691  45.997  49.181  52.439  58.585
#> hat.par[3,1]        43.780   4.674  34.929  40.444  43.670  46.893  53.086
#> hat.par[4,1]        42.219   4.663  33.651  38.862  42.033  45.319  51.897
#> hat.par[5,1]        17.095   2.537  12.253  15.343  17.006  18.782  22.369
#> hat.par[6,1]        43.554   4.194  35.362  40.766  43.571  46.432  51.594
#> hat.par[7,1]       157.556   7.502 143.111 152.442 157.476 162.683 172.089
#> hat.par[8,1]        68.457   5.546  57.417  64.796  68.333  72.090  79.501
#> hat.par[9,1]        90.134   5.154  79.973  86.556  90.277  93.565  99.992
#> hat.par[10,1]       78.868   3.781  71.388  76.329  78.915  81.572  86.033
#> hat.par[11,1]       72.142   5.517  61.378  68.409  72.086  75.742  83.441
#> hat.par[12,1]       77.471   4.240  69.183  74.686  77.435  80.340  85.696
#> hat.par[13,1]       49.365   4.750  40.105  46.092  49.277  52.635  58.648
#> hat.par[14,1]       33.892   4.652  25.310  30.674  33.690  36.957  43.902
#> hat.par[15,1]       36.083   4.538  27.439  33.054  36.018  39.107  44.919
#> hat.par[16,1]      303.331  13.242 277.951 294.340 303.057 312.366 329.840
#> hat.par[17,1]       11.238   2.590   6.680   9.333  11.086  13.022  16.637
#> hat.par[18,1]       22.135   3.293  15.942  19.805  22.125  24.229  28.770
#> hat.par[19,1]        4.078   1.391   1.822   3.050   3.936   4.962   7.167
#> hat.par[20,1]       25.276   3.662  18.583  22.720  25.068  27.734  33.020
#> hat.par[21,1]       32.237   4.356  24.136  29.264  32.075  35.038  41.064
#> hat.par[1,2]         1.426   0.747   0.318   0.868   1.316   1.875   3.180
#> hat.par[2,2]        46.857   4.982  37.607  43.360  46.676  50.294  56.638
#> hat.par[3,2]        30.901   4.074  23.105  28.099  30.882  33.550  39.202
#> hat.par[4,2]        43.717   5.279  33.906  40.133  43.741  47.105  54.683
#> hat.par[5,2]        11.858   2.103   7.907  10.383  11.786  13.225  16.197
#> hat.par[6,2]        35.212   3.953  27.873  32.552  35.115  37.812  43.356
#> hat.par[7,2]       197.033   9.988 177.933 190.141 196.991 203.918 216.637
#> hat.par[8,2]        51.374   5.053  41.939  48.006  51.229  54.716  61.447
#> hat.par[9,2]        81.855   5.724  70.845  77.953  81.769  85.762  93.308
#> hat.par[10,2]       72.989   3.955  65.085  70.337  73.097  75.742  80.504
#> hat.par[11,2]      119.842   8.210 104.162 114.270 119.778 125.228 136.454
#> hat.par[12,2]       80.503   4.861  70.989  77.084  80.603  83.926  89.716
#> hat.par[13,2]       26.037   4.489  17.675  23.012  25.919  29.038  34.955
#> hat.par[14,2]       31.722   4.401  23.526  28.711  31.529  34.495  41.137
#> hat.par[15,2]       32.980   4.556  24.817  29.715  32.698  35.890  42.472
#> hat.par[16,2]      247.749  12.470 223.423 239.332 247.618 256.055 271.808
#> hat.par[17,2]        6.907   1.840   3.718   5.611   6.765   8.095  10.878
#> hat.par[18,2]       12.821   2.499   8.265  11.045  12.742  14.390  18.086
#> hat.par[19,2]        2.441   0.933   1.008   1.764   2.324   2.972   4.634
#> hat.par[20,2]       18.777   3.098  13.300  16.560  18.615  20.779  25.232
#> hat.par[21,2]       21.551   3.489  15.069  19.121  21.498  23.839  28.572
#> hat.par[9,3]        78.949   5.673  67.573  75.168  78.991  82.833  89.448
#> hat.par[10,3]       68.290   4.331  59.879  65.355  68.312  71.326  76.580
#> hat.par[12,3]       67.337   4.854  57.798  64.092  67.403  70.604  76.802
#> hat.par[13,3]       35.334   5.261  25.241  31.760  35.197  38.840  46.040
#> hat.par[19,3]        2.444   0.904   0.988   1.794   2.353   2.984   4.450
#> hat.par[10,4]       66.420   4.306  57.479  63.619  66.522  69.412  74.210
#> hat.par[12,4]       62.898   4.504  53.866  59.813  62.852  66.145  71.418
#> hat.par[13,4]       41.132   4.809  31.913  37.837  40.977  44.338  50.574
#> phi[1]              -0.296   0.604  -1.595  -0.685  -0.263   0.120   0.779
#> phi[2]               0.007   0.960  -1.902  -0.619  -0.008   0.653   1.934
#> phi[3]               0.069   0.965  -1.844  -0.569   0.093   0.729   1.898
#> phi[4]              -0.951   0.902  -2.527  -1.592  -1.037  -0.385   1.024
#> phi[5]              -0.643   0.948  -2.477  -1.237  -0.661  -0.051   1.344
#> phi[6]               0.923   0.782  -0.813   0.456   0.950   1.426   2.433
#> phi[7]              -0.306   0.671  -1.698  -0.746  -0.284   0.151   0.954
#> phi[8]              -0.272   0.899  -2.042  -0.888  -0.270   0.325   1.554
#> tau                  0.172   0.085   0.040   0.108   0.164   0.220   0.367
#> totresdev.o         52.640   8.844  36.561  46.532  52.133  58.226  70.970
#> deviance           580.351  13.448 556.084 570.663 579.747 588.987 608.105
#>                     Rhat n.eff
#> EM[2,1]            1.084    35
#> EM[3,1]            1.063    52
#> EM[4,1]            1.130    21
#> EM[5,1]            1.024   770
#> EM[6,1]            1.131    23
#> EM[7,1]            1.039   120
#> EM[8,1]            1.023   350
#> EM[3,2]            1.007  1300
#> EM[4,2]            1.020   210
#> EM[5,2]            1.042    57
#> EM[6,2]            1.049    58
#> EM[7,2]            1.050    59
#> EM[8,2]            1.078    35
#> EM[4,3]            1.019   290
#> EM[5,3]            1.028    89
#> EM[6,3]            1.046    66
#> EM[7,3]            1.038    87
#> EM[8,3]            1.055    46
#> EM[5,4]            1.076    36
#> EM[6,4]            1.136    20
#> EM[7,4]            1.092    27
#> EM[8,4]            1.148    18
#> EM[6,5]            1.046    50
#> EM[7,5]            1.010  2000
#> EM[8,5]            1.006   420
#> EM[7,6]            1.087    31
#> EM[8,6]            1.108    25
#> EM[8,7]            1.030    81
#> EM.pred[2,1]       1.072    41
#> EM.pred[3,1]       1.055    63
#> EM.pred[4,1]       1.102    26
#> EM.pred[5,1]       1.020  1800
#> EM.pred[6,1]       1.066    46
#> EM.pred[7,1]       1.028   210
#> EM.pred[8,1]       1.015   520
#> EM.pred[3,2]       1.008   980
#> EM.pred[4,2]       1.018   250
#> EM.pred[5,2]       1.036    68
#> EM.pred[6,2]       1.048    55
#> EM.pred[7,2]       1.045    64
#> EM.pred[8,2]       1.060    44
#> EM.pred[4,3]       1.014   330
#> EM.pred[5,3]       1.030    88
#> EM.pred[6,3]       1.043    73
#> EM.pred[7,3]       1.030   110
#> EM.pred[8,3]       1.050    53
#> EM.pred[5,4]       1.056    49
#> EM.pred[6,4]       1.088    28
#> EM.pred[7,4]       1.067    36
#> EM.pred[8,4]       1.103    25
#> EM.pred[6,5]       1.042    58
#> EM.pred[7,5]       1.011  2500
#> EM.pred[8,5]       1.006   550
#> EM.pred[7,6]       1.056    53
#> EM.pred[8,6]       1.066    39
#> EM.pred[8,7]       1.018   200
#> SUCRA[1]           1.071    51
#> SUCRA[2]           1.053    76
#> SUCRA[3]           1.021   130
#> SUCRA[4]           1.072    32
#> SUCRA[5]           1.007   330
#> SUCRA[6]           1.130    24
#> SUCRA[7]           1.020   110
#> SUCRA[8]           1.072    33
#> abs_risk[1]        1.000     1
#> abs_risk[2]        1.086    34
#> abs_risk[3]        1.068    52
#> abs_risk[4]        1.129    21
#> abs_risk[5]        1.024   730
#> abs_risk[6]        1.141    22
#> abs_risk[7]        1.041   120
#> abs_risk[8]        1.025   310
#> delta[1,1]         1.000     1
#> delta[2,1]         1.000     1
#> delta[3,1]         1.000     1
#> delta[4,1]         1.000     1
#> delta[5,1]         1.000     1
#> delta[6,1]         1.000     1
#> delta[7,1]         1.000     1
#> delta[8,1]         1.000     1
#> delta[9,1]         1.000     1
#> delta[10,1]        1.000     1
#> delta[11,1]        1.000     1
#> delta[12,1]        1.000     1
#> delta[13,1]        1.000     1
#> delta[14,1]        1.000     1
#> delta[15,1]        1.000     1
#> delta[16,1]        1.000     1
#> delta[17,1]        1.000     1
#> delta[18,1]        1.000     1
#> delta[19,1]        1.000     1
#> delta[20,1]        1.000     1
#> delta[21,1]        1.000     1
#> delta[1,2]         1.085    32
#> delta[2,2]         1.112    23
#> delta[3,2]         1.018   320
#> delta[4,2]         1.004  2500
#> delta[5,2]         1.029   360
#> delta[6,2]         1.020   210
#> delta[7,2]         1.014  1600
#> delta[8,2]         1.028   140
#> delta[9,2]         1.025   320
#> delta[10,2]        1.151    19
#> delta[11,2]        1.116    25
#> delta[12,2]        1.154    18
#> delta[13,2]        1.096    33
#> delta[14,2]        1.027    81
#> delta[15,2]        1.094    30
#> delta[16,2]        1.006   930
#> delta[17,2]        1.072    40
#> delta[18,2]        1.019   960
#> delta[19,2]        1.020  2600
#> delta[20,2]        1.044    70
#> delta[21,2]        1.012   290
#> delta[9,3]         1.021   280
#> delta[10,3]        1.037   390
#> delta[12,3]        1.037   140
#> delta[13,3]        1.074    51
#> delta[19,3]        1.029   430
#> delta[10,4]        1.055   140
#> delta[12,4]        1.049    59
#> delta[13,4]        1.133    22
#> dev.o[1,1]         1.002  1300
#> dev.o[2,1]         1.006   350
#> dev.o[3,1]         1.001  3000
#> dev.o[4,1]         1.003  1300
#> dev.o[5,1]         1.002  3000
#> dev.o[6,1]         1.003   760
#> dev.o[7,1]         1.003   800
#> dev.o[8,1]         1.002  1700
#> dev.o[9,1]         1.001  3000
#> dev.o[10,1]        1.001  2100
#> dev.o[11,1]        1.001  3000
#> dev.o[12,1]        1.009   270
#> dev.o[13,1]        1.001  3000
#> dev.o[14,1]        1.001  2400
#> dev.o[15,1]        1.002  1600
#> dev.o[16,1]        1.004   630
#> dev.o[17,1]        1.007   330
#> dev.o[18,1]        1.002  2100
#> dev.o[19,1]        1.002  2500
#> dev.o[20,1]        1.001  3000
#> dev.o[21,1]        1.001  3000
#> dev.o[1,2]         1.003   770
#> dev.o[2,2]         1.004   580
#> dev.o[3,2]         1.002  1700
#> dev.o[4,2]         1.001  3000
#> dev.o[5,2]         1.001  3000
#> dev.o[6,2]         1.001  3000
#> dev.o[7,2]         1.001  3000
#> dev.o[8,2]         1.006   350
#> dev.o[9,2]         1.004   960
#> dev.o[10,2]        1.007   540
#> dev.o[11,2]        1.007   860
#> dev.o[12,2]        1.001  3000
#> dev.o[13,2]        1.002  3000
#> dev.o[14,2]        1.001  3000
#> dev.o[15,2]        1.003  1500
#> dev.o[16,2]        1.003   910
#> dev.o[17,2]        1.010   260
#> dev.o[18,2]        1.001  3000
#> dev.o[19,2]        1.001  3000
#> dev.o[20,2]        1.001  3000
#> dev.o[21,2]        1.003   960
#> dev.o[9,3]         1.003  1200
#> dev.o[10,3]        1.001  2900
#> dev.o[12,3]        1.004   550
#> dev.o[13,3]        1.002  2000
#> dev.o[19,3]        1.003  1500
#> dev.o[10,4]        1.003  1100
#> dev.o[12,4]        1.002  1400
#> dev.o[13,4]        1.003   830
#> effectiveness[1,1] 1.000     1
#> effectiveness[2,1] 1.012   170
#> effectiveness[3,1] 1.002  1700
#> effectiveness[4,1] 1.108   470
#> effectiveness[5,1] 1.046   140
#> effectiveness[6,1] 1.293   750
#> effectiveness[7,1] 1.107   240
#> effectiveness[8,1] 1.065   130
#> effectiveness[1,2] 1.000     1
#> effectiveness[2,2] 1.003   840
#> effectiveness[3,2] 1.013   180
#> effectiveness[4,2] 1.029   390
#> effectiveness[5,2] 1.001  3000
#> effectiveness[6,2] 1.189  1200
#> effectiveness[7,2] 1.027   230
#> effectiveness[8,2] 1.063    61
#> effectiveness[1,3] 1.292  1000
#> effectiveness[2,3] 1.047   190
#> effectiveness[3,3] 1.007   620
#> effectiveness[4,3] 1.026   230
#> effectiveness[5,3] 1.002  1100
#> effectiveness[6,3] 1.121   470
#> effectiveness[7,3] 1.002  1500
#> effectiveness[8,3] 1.007   350
#> effectiveness[1,4] 1.295   370
#> effectiveness[2,4] 1.040   330
#> effectiveness[3,4] 1.032   160
#> effectiveness[4,4] 1.038   140
#> effectiveness[5,4] 1.009   440
#> effectiveness[6,4] 1.071   270
#> effectiveness[7,4] 1.002  1100
#> effectiveness[8,4] 1.006   400
#> effectiveness[1,5] 1.129   100
#> effectiveness[2,5] 1.072   150
#> effectiveness[3,5] 1.043   120
#> effectiveness[4,5] 1.006   510
#> effectiveness[5,5] 1.004   900
#> effectiveness[6,5] 1.165    43
#> effectiveness[7,5] 1.015   170
#> effectiveness[8,5] 1.040    96
#> effectiveness[1,6] 1.097    39
#> effectiveness[2,6] 1.040   620
#> effectiveness[3,6] 1.004  1600
#> effectiveness[4,6] 1.002  1700
#> effectiveness[5,6] 1.078    60
#> effectiveness[6,6] 1.013   280
#> effectiveness[7,6] 1.017   260
#> effectiveness[8,6] 1.059   140
#> effectiveness[1,7] 1.068    35
#> effectiveness[2,7] 1.071   480
#> effectiveness[3,7] 1.056   210
#> effectiveness[4,7] 1.042    96
#> effectiveness[5,7] 1.024   370
#> effectiveness[6,7] 1.009   250
#> effectiveness[7,7] 1.115   460
#> effectiveness[8,7] 1.088  1800
#> effectiveness[1,8] 1.054    44
#> effectiveness[2,8] 1.198   260
#> effectiveness[3,8] 1.040   420
#> effectiveness[4,8] 1.063    62
#> effectiveness[5,8] 1.019  1800
#> effectiveness[6,8] 1.077    31
#> effectiveness[7,8] 1.291  3000
#> effectiveness[8,8] 1.000     1
#> hat.par[1,1]       1.003  1000
#> hat.par[2,1]       1.032    71
#> hat.par[3,1]       1.008   340
#> hat.par[4,1]       1.001  3000
#> hat.par[5,1]       1.001  3000
#> hat.par[6,1]       1.006   350
#> hat.par[7,1]       1.007   330
#> hat.par[8,1]       1.011   280
#> hat.par[9,1]       1.009   250
#> hat.par[10,1]      1.008   390
#> hat.par[11,1]      1.009   240
#> hat.par[12,1]      1.025   100
#> hat.par[13,1]      1.004   590
#> hat.par[14,1]      1.011   200
#> hat.par[15,1]      1.013   190
#> hat.par[16,1]      1.002  1500
#> hat.par[17,1]      1.008   270
#> hat.par[18,1]      1.001  2600
#> hat.par[19,1]      1.001  3000
#> hat.par[20,1]      1.004   650
#> hat.par[21,1]      1.002  3000
#> hat.par[1,2]       1.003   790
#> hat.par[2,2]       1.022    95
#> hat.par[3,2]       1.002  2500
#> hat.par[4,2]       1.001  3000
#> hat.par[5,2]       1.002  1900
#> hat.par[6,2]       1.002  1900
#> hat.par[7,2]       1.005   450
#> hat.par[8,2]       1.003  3000
#> hat.par[9,2]       1.007   340
#> hat.par[10,2]      1.002  1500
#> hat.par[11,2]      1.008   270
#> hat.par[12,2]      1.010   220
#> hat.par[13,2]      1.005   420
#> hat.par[14,2]      1.005   500
#> hat.par[15,2]      1.004   530
#> hat.par[16,2]      1.005   400
#> hat.par[17,2]      1.012   170
#> hat.par[18,2]      1.001  3000
#> hat.par[19,2]      1.001  3000
#> hat.par[20,2]      1.001  2200
#> hat.par[21,2]      1.003   880
#> hat.par[9,3]       1.008   290
#> hat.par[10,3]      1.002  1600
#> hat.par[12,3]      1.006   390
#> hat.par[13,3]      1.005   660
#> hat.par[19,3]      1.002  1200
#> hat.par[10,4]      1.006   350
#> hat.par[12,4]      1.002  2900
#> hat.par[13,4]      1.008   260
#> phi[1]             1.072    73
#> phi[2]             1.021   110
#> phi[3]             1.014   160
#> phi[4]             1.086    29
#> phi[5]             1.006  1900
#> phi[6]             1.085    30
#> phi[7]             1.001  3000
#> phi[8]             1.012   180
#> tau                1.103    28
#> totresdev.o        1.008   270
#> deviance           1.004   540
#> 
#> For each parameter, n.eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
#> 
#> DIC info (using the rule: pV = var(deviance)/2)
#> pV = 90.1 and DIC = 670.5
#> DIC is an estimate of expected predictive error (lower deviance is better).
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
#> $EM_pred
#>                     mean        sd       2.5%         25%          50%
#> EM.pred[2,1] -0.93439096 0.4704194 -1.8614694 -1.24168961 -0.945917720
#> EM.pred[3,1] -0.65447817 0.4689424 -1.6348849 -0.93892930 -0.636007958
#> EM.pred[4,1] -0.24266729 0.3775789 -0.9263943 -0.50086641 -0.255940883
#> EM.pred[5,1] -0.51249373 0.3446975 -1.2232883 -0.71768736 -0.517654534
#> EM.pred[6,1] -0.02499362 0.3171654 -0.7019794 -0.22923936  0.006346857
#> EM.pred[7,1] -0.46408359 0.2605853 -0.9741160 -0.63626871 -0.463969222
#> EM.pred[8,1] -0.53440540 0.2559723 -1.0834746 -0.68196860 -0.519379512
#> EM.pred[3,2]  0.28239358 0.5225964 -0.7572123 -0.07682769  0.280168323
#> EM.pred[4,2]  0.69948748 0.5027392 -0.3119523  0.35894421  0.676560974
#> EM.pred[5,2]  0.42752658 0.5114875 -0.6577451  0.09042370  0.428445960
#> EM.pred[6,2]  0.92019240 0.4901500 -0.1117104  0.61132289  0.939923249
#> EM.pred[7,2]  0.47512310 0.4768725 -0.4675294  0.18538482  0.471539352
#> EM.pred[8,2]  0.41340263 0.4747481 -0.5604168  0.11069493  0.434811085
#> EM.pred[4,3]  0.42229359 0.4997595 -0.5691731  0.10486654  0.421767299
#> EM.pred[5,3]  0.15374842 0.5053933 -0.8432533 -0.17764469  0.156086470
#> EM.pred[6,3]  0.63760030 0.4791304 -0.3121988  0.30982717  0.661371142
#> EM.pred[7,3]  0.19457240 0.4567659 -0.7041167 -0.09526075  0.206255204
#> EM.pred[8,3]  0.12483266 0.4594094 -0.7960228 -0.16358152  0.130299035
#> EM.pred[5,4] -0.26064612 0.4181152 -1.0821877 -0.53579344 -0.255345586
#> EM.pred[6,4]  0.22612849 0.3972014 -0.5622487 -0.03102765  0.233178567
#> EM.pred[7,4] -0.22891762 0.3601705 -0.9825723 -0.45269435 -0.216638600
#> EM.pred[8,4] -0.28885883 0.3685978 -1.0722271 -0.52116624 -0.265429240
#> EM.pred[6,5]  0.48355438 0.3820296 -0.3609312  0.25321169  0.506640986
#> EM.pred[7,5]  0.04695220 0.3440098 -0.6480747 -0.17090491  0.051018359
#> EM.pred[8,5] -0.01261039 0.3456008 -0.7292571 -0.22240927  0.002055410
#> EM.pred[7,6] -0.44734297 0.3081427 -1.0432561 -0.64515355 -0.458080275
#> EM.pred[8,6] -0.50989819 0.3057031 -1.1345109 -0.70100857 -0.519196894
#> EM.pred[8,7] -0.06412286 0.2454193 -0.5571835 -0.20392534 -0.062210405
#>                       75%       97.5%     Rhat n.eff
#> EM.pred[2,1] -0.630458325  0.03164389 1.072033    41
#> EM.pred[3,1] -0.359170256  0.24194114 1.054582    63
#> EM.pred[4,1]  0.004391119  0.50559384 1.101752    26
#> EM.pred[5,1] -0.313086964  0.20585867 1.020245  1800
#> EM.pred[6,1]  0.194074138  0.54351576 1.066256    46
#> EM.pred[7,1] -0.293386855  0.02808350 1.028330   210
#> EM.pred[8,1] -0.372937561 -0.05243824 1.014955   520
#> EM.pred[3,2]  0.650558005  1.28555630 1.008459   980
#> EM.pred[4,2]  1.043868560  1.68069533 1.018297   250
#> EM.pred[5,2]  0.753830869  1.43865029 1.035836    68
#> EM.pred[6,2]  1.251872459  1.82671052 1.048476    55
#> EM.pred[7,2]  0.790803087  1.38061819 1.044528    64
#> EM.pred[8,2]  0.732011018  1.29857088 1.060415    44
#> EM.pred[4,3]  0.708619680  1.51178973 1.014048   330
#> EM.pred[5,3]  0.484972329  1.17886513 1.030255    88
#> EM.pred[6,3]  0.945822482  1.58678844 1.043321    73
#> EM.pred[7,3]  0.480999282  1.13309407 1.029974   110
#> EM.pred[8,3]  0.412834899  1.05224282 1.050467    53
#> EM.pred[5,4]  0.015697714  0.56894551 1.056435    49
#> EM.pred[6,4]  0.504158712  0.95553147 1.088144    28
#> EM.pred[7,4]  0.019801957  0.46145679 1.067302    36
#> EM.pred[8,4] -0.018816292  0.34973276 1.103475    25
#> EM.pred[6,5]  0.733370853  1.21049586 1.041712    58
#> EM.pred[7,5]  0.264470553  0.71716542 1.010742  2500
#> EM.pred[8,5]  0.209783535  0.64480703 1.005913   550
#> EM.pred[7,6] -0.253412601  0.18610010 1.055625    53
#> EM.pred[8,6] -0.307304124  0.10911158 1.065884    39
#> EM.pred[8,7]  0.085292484  0.41086400 1.017754   200
#> 
#> $tau
#>        mean          sd        2.5%         25%         50%         75% 
#>  0.17204709  0.08503728  0.04034496  0.10812468  0.16447152  0.21951310 
#>       97.5%        Rhat       n.eff 
#>  0.36724217  1.10256287 28.00000000 
#> 
#> $delta
#>                    mean        sd       2.5%        25%         50%         75%
#> delta[1,2]  -0.30360684 0.3705316 -0.9804842 -0.5597374 -0.31978681 -0.06163109
#> delta[2,2]  -0.28513924 0.3018242 -0.8130273 -0.5055970 -0.30631918 -0.07979417
#> delta[3,2]  -0.54745428 0.2287103 -1.0068736 -0.6940086 -0.55290765 -0.38317315
#> delta[4,2]  -0.48415564 0.1821835 -0.8608569 -0.6056770 -0.47771304 -0.35942222
#> delta[5,2]  -0.48596682 0.2479148 -0.9590630 -0.6530709 -0.49119085 -0.31613492
#> delta[6,2]  -0.38105303 0.2206252 -0.7865092 -0.5385429 -0.38831676 -0.24094187
#> delta[7,2]  -0.48589965 0.1763856 -0.8303081 -0.6072086 -0.48169695 -0.36970840
#> delta[8,2]  -0.43730798 0.1910446 -0.8185150 -0.5698080 -0.43788158 -0.30856391
#> delta[9,2]  -0.47715431 0.2011313 -0.8493187 -0.6139492 -0.48351989 -0.34977536
#> delta[10,2] -0.16530047 0.3733135 -0.8053398 -0.4293337 -0.19267827  0.06923352
#> delta[11,2] -0.09715015 0.2577933 -0.6138716 -0.2764856 -0.08340293  0.10217794
#> delta[12,2] -0.19799537 0.3318056 -0.7679960 -0.4374009 -0.22777140  0.02932860
#> delta[13,2] -1.00174973 0.4087018 -1.8306606 -1.2616452 -1.00849746 -0.74547967
#> delta[14,2] -0.11429727 0.1956355 -0.5263633 -0.2424923 -0.10480327  0.01747769
#> delta[15,2]  0.04260936 0.2440289 -0.4658495 -0.1183496  0.05338254  0.20964446
#> delta[16,2] -0.35114501 0.1350339 -0.6082381 -0.4422329 -0.35289596 -0.26237493
#> delta[17,2] -0.62022687 0.2714952 -1.2030210 -0.7734245 -0.61272807 -0.43873494
#> delta[18,2] -0.57270486 0.3021451 -1.1691907 -0.7586894 -0.57280096 -0.39846036
#> delta[19,2] -0.55376438 0.3319558 -1.1781687 -0.7620615 -0.55757814 -0.35587001
#> delta[20,2] -0.43363439 0.2240522 -0.8701833 -0.5921158 -0.43263091 -0.28222668
#> delta[21,2] -0.58738616 0.2145471 -1.0384818 -0.7216971 -0.57938474 -0.43913582
#> delta[9,3]  -0.57551145 0.1913594 -0.9606819 -0.7016396 -0.57360520 -0.44480359
#> delta[10,3] -0.51557133 0.3121005 -1.1159240 -0.7078375 -0.52737107 -0.34070904
#> delta[12,3] -0.39525223 0.3018513 -0.9338872 -0.5849580 -0.42114774 -0.23019241
#> delta[13,3] -0.71797138 0.4112892 -1.6281570 -0.9815761 -0.68935699 -0.44485897
#> delta[19,3] -0.52232202 0.2681658 -1.0837366 -0.6789361 -0.51155399 -0.34622689
#> delta[10,4] -0.51571395 0.2482143 -1.0333548 -0.6694216 -0.51443932 -0.34303701
#> delta[12,4] -0.37081619 0.2354454 -0.8090300 -0.5339618 -0.38591251 -0.22017771
#> delta[13,4] -0.13511926 0.3082854 -0.7838225 -0.3441567 -0.09722973  0.09196899
#>                   97.5%     Rhat n.eff
#> delta[1,2]   0.43760712 1.085278    32
#> delta[2,2]   0.33651124 1.111735    23
#> delta[3,2]  -0.11458314 1.018275   320
#> delta[4,2]  -0.14854008 1.003847  2500
#> delta[5,2]  -0.02148006 1.029304   360
#> delta[6,2]   0.06067270 1.019850   210
#> delta[7,2]  -0.14514264 1.014450  1600
#> delta[8,2]  -0.06389999 1.027682   140
#> delta[9,2]  -0.07491692 1.024779   320
#> delta[10,2]  0.67052254 1.150931    19
#> delta[11,2]  0.33466107 1.115638    25
#> delta[12,2]  0.48189918 1.153628    18
#> delta[13,2] -0.16566961 1.096418    33
#> delta[14,2]  0.25136207 1.027155    81
#> delta[15,2]  0.50421575 1.093609    30
#> delta[16,2] -0.08046684 1.006233   930
#> delta[17,2] -0.12148618 1.072150    40
#> delta[18,2]  0.06176846 1.018599   960
#> delta[19,2]  0.12997628 1.019847  2600
#> delta[20,2] -0.00868698 1.044186    70
#> delta[21,2] -0.19607312 1.011896   290
#> delta[9,3]  -0.19669213 1.021363   280
#> delta[10,3]  0.19438038 1.037480   390
#> delta[12,3]  0.27744021 1.036949   140
#> delta[13,3]  0.02707925 1.074299    51
#> delta[19,3] -0.04421840 1.029426   430
#> delta[10,4] -0.05803979 1.054922   140
#> delta[12,4]  0.14045257 1.049369    59
#> delta[13,4]  0.36382253 1.132654    22
#> 
#> $heter_prior
#> [1] 0 1 1
#> 
#> $SUCRA
#>               mean        sd      2.5%       25%       50%       75%     97.5%
#> SUCRA[1] 0.1230000 0.1203099 0.0000000 0.0000000 0.1428571 0.1428571 0.4285714
#> SUCRA[2] 0.8791429 0.2075278 0.2857143 0.8571429 1.0000000 1.0000000 1.0000000
#> SUCRA[3] 0.7030000 0.2757550 0.0000000 0.5714286 0.8571429 0.8571429 1.0000000
#> SUCRA[4] 0.3462857 0.2422618 0.0000000 0.1428571 0.2857143 0.4285714 0.8571429
#> SUCRA[5] 0.6027143 0.2501997 0.1428571 0.4285714 0.5714286 0.8571429 1.0000000
#> SUCRA[6] 0.1427619 0.1607225 0.0000000 0.0000000 0.1428571 0.2857143 0.5714286
#> SUCRA[7] 0.5610952 0.1755616 0.2857143 0.4285714 0.5714286 0.7142857 0.8571429
#> SUCRA[8] 0.6420000 0.1815187 0.2857143 0.5714286 0.7142857 0.7142857 1.0000000
#>              Rhat n.eff
#> SUCRA[1] 1.071489    51
#> SUCRA[2] 1.053231    76
#> SUCRA[3] 1.020851   130
#> SUCRA[4] 1.072432    32
#> SUCRA[5] 1.006562   330
#> SUCRA[6] 1.129914    24
#> SUCRA[7] 1.020369   110
#> SUCRA[8] 1.071823    33
#> 
#> $effectiveness
#>                            mean         sd 2.5% 25% 50% 75% 97.5%     Rhat
#> effectiveness[1,1] 0.0000000000 0.00000000    0   0   0   0     0 1.000000
#> effectiveness[2,1] 0.6110000000 0.48760461    0   0   1   1     1 1.012093
#> effectiveness[3,1] 0.2240000000 0.41699156    0   0   0   0     1 1.001654
#> effectiveness[4,1] 0.0083333333 0.09092109    0   0   0   0     0 1.107615
#> effectiveness[5,1] 0.0803333333 0.27185386    0   0   0   0     1 1.046365
#> effectiveness[6,1] 0.0013333333 0.03649657    0   0   0   0     0 1.292577
#> effectiveness[7,1] 0.0170000000 0.12929258    0   0   0   0     0 1.107453
#> effectiveness[8,1] 0.0580000000 0.23378242    0   0   0   0     1 1.064789
#> effectiveness[1,2] 0.0000000000 0.00000000    0   0   0   0     0 1.000000
#> effectiveness[2,2] 0.2096666667 0.40713856    0   0   0   0     1 1.003097
#> effectiveness[3,2] 0.3046666667 0.46034284    0   0   0   1     1 1.012982
#> effectiveness[4,2] 0.0443333333 0.20586893    0   0   0   0     1 1.028807
#> effectiveness[5,2] 0.1966666667 0.39754442    0   0   0   0     1 1.000567
#> effectiveness[6,2] 0.0016666667 0.04079759    0   0   0   0     0 1.189187
#> effectiveness[7,2] 0.0856666667 0.27991786    0   0   0   0     1 1.026792
#> effectiveness[8,2] 0.1573333333 0.36417546    0   0   0   0     1 1.063245
#> effectiveness[1,3] 0.0010000000 0.03161223    0   0   0   0     0 1.292018
#> effectiveness[2,3] 0.0553333333 0.22866785    0   0   0   0     1 1.046989
#> effectiveness[3,3] 0.1216666667 0.32695492    0   0   0   0     1 1.007121
#> effectiveness[4,3] 0.0886666667 0.28430940    0   0   0   0     1 1.026068
#> effectiveness[5,3] 0.2196666667 0.41408982    0   0   0   0     1 1.002310
#> effectiveness[6,3] 0.0073333333 0.08533454    0   0   0   0     0 1.120680
#> effectiveness[7,3] 0.2170000000 0.41227134    0   0   0   0     1 1.001806
#> effectiveness[8,3] 0.2893333333 0.45352852    0   0   0   1     1 1.006567
#> effectiveness[1,4] 0.0026666667 0.05157948    0   0   0   0     0 1.294827
#> effectiveness[2,4] 0.0373333333 0.18960891    0   0   0   0     1 1.039633
#> effectiveness[3,4] 0.1063333333 0.30831517    0   0   0   0     1 1.032325
#> effectiveness[4,4] 0.1040000000 0.30531143    0   0   0   0     1 1.037858
#> effectiveness[5,4] 0.1453333333 0.35249535    0   0   0   0     1 1.008628
#> effectiveness[6,4] 0.0240000000 0.15307453    0   0   0   0     0 1.071314
#> effectiveness[7,4] 0.3050000000 0.46048418    0   0   0   1     1 1.002362
#> effectiveness[8,4] 0.2753333333 0.44675655    0   0   0   1     1 1.005782
#> effectiveness[1,5] 0.0343333333 0.18211428    0   0   0   0     1 1.128849
#> effectiveness[2,5] 0.0470000000 0.21167413    0   0   0   0     1 1.071537
#> effectiveness[3,5] 0.1040000000 0.30531143    0   0   0   0     1 1.042713
#> effectiveness[4,5] 0.1753333333 0.38031535    0   0   0   0     1 1.006315
#> effectiveness[5,5] 0.1653333333 0.37154305    0   0   0   0     1 1.003547
#> effectiveness[6,5] 0.0676666667 0.25121490    0   0   0   0     1 1.165340
#> effectiveness[7,5] 0.2476666667 0.43172910    0   0   0   0     1 1.015329
#> effectiveness[8,5] 0.1586666667 0.36542587    0   0   0   0     1 1.039602
#> effectiveness[1,6] 0.1656666667 0.37184313    0   0   0   0     1 1.096627
#> effectiveness[2,6] 0.0193333333 0.13771666    0   0   0   0     0 1.039942
#> effectiveness[3,6] 0.0693333333 0.25406247    0   0   0   0     1 1.004379
#> effectiveness[4,6] 0.2866666667 0.45227986    0   0   0   1     1 1.001648
#> effectiveness[5,6] 0.1220000000 0.32734037    0   0   0   0     1 1.077560
#> effectiveness[6,6] 0.1590000000 0.36573705    0   0   0   0     1 1.012850
#> effectiveness[7,6] 0.1193333333 0.32423438    0   0   0   0     1 1.017284
#> effectiveness[8,6] 0.0586666667 0.23503894    0   0   0   0     1 1.059141
#> effectiveness[1,7] 0.4110000000 0.49209727    0   0   0   1     1 1.067836
#> effectiveness[2,7] 0.0133333333 0.11471679    0   0   0   0     0 1.070838
#> effectiveness[3,7] 0.0406666667 0.19754973    0   0   0   0     1 1.056379
#> effectiveness[4,7] 0.1410000000 0.34807957    0   0   0   0     1 1.042102
#> effectiveness[5,7] 0.0570000000 0.23188127    0   0   0   0     1 1.023654
#> effectiveness[6,7] 0.3263333333 0.46894903    0   0   0   1     1 1.008905
#> effectiveness[7,7] 0.0080000000 0.08909908    0   0   0   0     0 1.115065
#> effectiveness[8,7] 0.0026666667 0.05157948    0   0   0   0     0 1.088408
#> effectiveness[1,8] 0.3853333333 0.48675511    0   0   0   1     1 1.053983
#> effectiveness[2,8] 0.0070000000 0.08338656    0   0   0   0     0 1.197953
#> effectiveness[3,8] 0.0293333333 0.16876725    0   0   0   0     1 1.039539
#> effectiveness[4,8] 0.1516666667 0.35875729    0   0   0   0     1 1.062882
#> effectiveness[5,8] 0.0136666667 0.11612228    0   0   0   0     0 1.019439
#> effectiveness[6,8] 0.4126666667 0.49239588    0   0   0   1     1 1.076706
#> effectiveness[7,8] 0.0003333333 0.01825742    0   0   0   0     0 1.290904
#> effectiveness[8,8] 0.0000000000 0.00000000    0   0   0   0     0 1.000000
#>                    n.eff
#> effectiveness[1,1]     1
#> effectiveness[2,1]   170
#> effectiveness[3,1]  1700
#> effectiveness[4,1]   470
#> effectiveness[5,1]   140
#> effectiveness[6,1]   750
#> effectiveness[7,1]   240
#> effectiveness[8,1]   130
#> effectiveness[1,2]     1
#> effectiveness[2,2]   840
#> effectiveness[3,2]   180
#> effectiveness[4,2]   390
#> effectiveness[5,2]  3000
#> effectiveness[6,2]  1200
#> effectiveness[7,2]   230
#> effectiveness[8,2]    61
#> effectiveness[1,3]  1000
#> effectiveness[2,3]   190
#> effectiveness[3,3]   620
#> effectiveness[4,3]   230
#> effectiveness[5,3]  1100
#> effectiveness[6,3]   470
#> effectiveness[7,3]  1500
#> effectiveness[8,3]   350
#> effectiveness[1,4]   370
#> effectiveness[2,4]   330
#> effectiveness[3,4]   160
#> effectiveness[4,4]   140
#> effectiveness[5,4]   440
#> effectiveness[6,4]   270
#> effectiveness[7,4]  1100
#> effectiveness[8,4]   400
#> effectiveness[1,5]   100
#> effectiveness[2,5]   150
#> effectiveness[3,5]   120
#> effectiveness[4,5]   510
#> effectiveness[5,5]   900
#> effectiveness[6,5]    43
#> effectiveness[7,5]   170
#> effectiveness[8,5]    96
#> effectiveness[1,6]    39
#> effectiveness[2,6]   620
#> effectiveness[3,6]  1600
#> effectiveness[4,6]  1700
#> effectiveness[5,6]    60
#> effectiveness[6,6]   280
#> effectiveness[7,6]   260
#> effectiveness[8,6]   140
#> effectiveness[1,7]    35
#> effectiveness[2,7]   480
#> effectiveness[3,7]   210
#> effectiveness[4,7]    96
#> effectiveness[5,7]   370
#> effectiveness[6,7]   250
#> effectiveness[7,7]   460
#> effectiveness[8,7]  1800
#> effectiveness[1,8]    44
#> effectiveness[2,8]   260
#> effectiveness[3,8]   420
#> effectiveness[4,8]    62
#> effectiveness[5,8]  1800
#> effectiveness[6,8]    31
#> effectiveness[7,8]  3000
#> effectiveness[8,8]     1
#> 
#> $abs_risk
#>                  mean         sd       2.5%       25%       50%       75%
#> abs_risk[1] 0.3916667 0.00000000 0.39166667 0.3916667 0.3916667 0.3916667
#> abs_risk[2] 0.2088031 0.07097138 0.09885475 0.1586033 0.1987535 0.2474813
#> abs_risk[3] 0.2578731 0.07894550 0.11744140 0.2011504 0.2534021 0.3065802
#> abs_risk[4] 0.3389695 0.07255406 0.22478322 0.2856439 0.3284934 0.3866303
#> abs_risk[5] 0.2816933 0.05742292 0.18360515 0.2448308 0.2738741 0.3127171
#> abs_risk[6] 0.3878941 0.05749216 0.27050497 0.3474270 0.3913537 0.4317447
#> abs_risk[7] 0.2895043 0.03738017 0.22509005 0.2620931 0.2856642 0.3139350
#> abs_risk[8] 0.2759401 0.03236339 0.21473431 0.2541979 0.2750576 0.2986384
#>                 97.5%     Rhat n.eff
#> abs_risk[1] 0.3916667 1.000000     1
#> abs_risk[2] 0.3798613 1.085941    34
#> abs_risk[3] 0.4236047 1.068186    52
#> abs_risk[4] 0.4994336 1.128924    21
#> abs_risk[5] 0.4230378 1.024433   730
#> abs_risk[6] 0.4891279 1.140591    22
#> abs_risk[7] 0.3710775 1.040769   120
#> abs_risk[8] 0.3410702 1.024887   310
#> 
#> $base_risk
#> [1] 0.3916667
#> 
#> attr(,"class")
#> [1] "run_model"
# }
```
