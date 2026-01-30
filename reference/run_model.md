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
#>                mean        sd        2.5%         25%         50%          75%
#> EM[2,1] -1.08356577 0.4134275 -1.85517734 -1.39396693 -1.08301970 -0.771128029
#> EM[3,1] -0.76780323 0.3950076 -1.50458195 -1.04253244 -0.78918749 -0.489777892
#> EM[4,1] -0.35768059 0.2676751 -0.82674772 -0.53422407 -0.37301852 -0.193287387
#> EM[5,1] -0.46810429 0.3036987 -1.04915425 -0.70451723 -0.46158698 -0.226013390
#> EM[6,1] -0.16384500 0.2432555 -0.67836556 -0.31371809 -0.15772696  0.008266464
#> EM[7,1] -0.44795818 0.1631749 -0.76194361 -0.56310252 -0.44954775 -0.328996523
#> EM[8,1] -0.52183475 0.1424579 -0.79622997 -0.62157268 -0.52587720 -0.408946047
#> EM[3,2]  0.31576254 0.4429617 -0.56562623  0.01181802  0.30736298  0.631626576
#> EM[4,2]  0.72588517 0.4525660 -0.05966864  0.39689776  0.68615276  1.030211089
#> EM[5,2]  0.61546148 0.5168517 -0.41430963  0.27142292  0.59479238  0.986640219
#> EM[6,2]  0.91972077 0.4719601  0.01570235  0.57069761  0.92715167  1.230362926
#> EM[7,2]  0.63560759 0.4044169 -0.20867169  0.37524977  0.63196620  0.928151762
#> EM[8,2]  0.56173102 0.4073165 -0.28257036  0.30524133  0.55203085  0.837381271
#> EM[4,3]  0.41012264 0.4194597 -0.36550531  0.08723769  0.43992636  0.702858118
#> EM[5,3]  0.29969895 0.4862564 -0.57838624 -0.06357009  0.26663856  0.625186669
#> EM[6,3]  0.60395824 0.4631669 -0.35878088  0.30214562  0.62284596  0.930455645
#> EM[7,3]  0.31984505 0.3747346 -0.41814273  0.05668522  0.33323152  0.578874077
#> EM[8,3]  0.24596849 0.3814765 -0.51775839 -0.02190225  0.25534068  0.509264951
#> EM[5,4] -0.11042369 0.4102821 -0.83969933 -0.42329027 -0.10951213  0.177771285
#> EM[6,4]  0.19383560 0.3412380 -0.47250122 -0.03800298  0.17914769  0.429708596
#> EM[7,4] -0.09027759 0.2801199 -0.70744887 -0.27142394 -0.05698246  0.113370877
#> EM[8,4] -0.16415415 0.2850610 -0.80669580 -0.35994492 -0.12567745  0.036177831
#> EM[6,5]  0.30425929 0.3481205 -0.33696876  0.05116941  0.31130131  0.553569623
#> EM[7,5]  0.02014610 0.3163163 -0.56203541 -0.22343510  0.03592213  0.241313059
#> EM[8,5] -0.05373046 0.3090578 -0.61209270 -0.29400498 -0.03620529  0.189363754
#> EM[7,6] -0.28411319 0.2876700 -0.85174973 -0.49111331 -0.28483508 -0.075651802
#> EM[8,6] -0.35798975 0.2659676 -0.90587825 -0.53541822 -0.35966606 -0.178580954
#> EM[8,7] -0.07387656 0.1526641 -0.37883065 -0.17498151 -0.07639117  0.019376770
#>               97.5%     Rhat n.eff
#> EM[2,1] -0.29185159 1.010676   200
#> EM[3,1] -0.02457772 1.083741    31
#> EM[4,1]  0.21371604 1.034842   710
#> EM[5,1]  0.05960768 1.111555    23
#> EM[6,1]  0.28060349 1.059770    48
#> EM[7,1] -0.14741517 1.215029    14
#> EM[8,1] -0.25226692 1.084618    30
#> EM[3,2]  1.19450501 1.031636    68
#> EM[4,2]  1.67003405 1.026768   130
#> EM[5,2]  1.63348596 1.038147   110
#> EM[6,2]  1.84883306 1.046505    62
#> EM[7,2]  1.37118786 1.014463   320
#> EM[8,2]  1.34704273 1.006925  1700
#> EM[4,3]  1.17785784 1.088546    28
#> EM[5,3]  1.31495730 1.127278    25
#> EM[6,3]  1.43763597 1.081337    32
#> EM[7,3]  1.02427671 1.015672   140
#> EM[8,3]  0.98981057 1.035783    62
#> EM[5,4]  0.64825338 1.080247    36
#> EM[6,4]  0.83674991 1.042978   120
#> EM[7,4]  0.37837288 1.108166    26
#> EM[8,4]  0.31350617 1.061332    62
#> EM[6,5]  0.95904631 1.218532    13
#> EM[7,5]  0.61201202 1.123837    22
#> EM[8,5]  0.48534353 1.093934    27
#> EM[7,6]  0.22583522 1.142103    19
#> EM[8,6]  0.16079630 1.098612    26
#> EM[8,7]  0.22652460 1.044024    62
#> 
#> $dev_o
#>                  mean        sd         2.5%        25%       50%       75%
#> dev.o[1,1]  2.3595453 2.4879537 0.0061475717 0.53521827 1.6599339 3.4396319
#> dev.o[2,1]  0.7869031 1.1054362 0.0009774747 0.07315485 0.3574974 1.0555491
#> dev.o[3,1]  1.2049459 1.4972431 0.0010541789 0.14027180 0.6403530 1.7443521
#> dev.o[4,1]  0.6943490 1.0689509 0.0008403616 0.06544637 0.2924322 0.8787675
#> dev.o[5,1]  0.7288581 1.0233459 0.0008195145 0.07309802 0.3392500 0.9780301
#> dev.o[6,1]  0.9126697 1.2107064 0.0011193771 0.09776983 0.4469389 1.2407489
#> dev.o[7,1]  0.7232231 1.0471315 0.0005334526 0.06587194 0.3426080 0.9242639
#> dev.o[8,1]  0.6954469 0.9771132 0.0006512589 0.07319712 0.3103599 0.8873605
#> dev.o[9,1]  0.7455268 0.9865272 0.0010541278 0.08778283 0.3594195 1.0140124
#> dev.o[10,1] 0.5864596 0.8831695 0.0005206356 0.05212638 0.2360605 0.7571391
#> dev.o[11,1] 0.9247128 1.2958548 0.0011394954 0.09585769 0.4192176 1.2026921
#> dev.o[12,1] 1.3524077 1.4761989 0.0024856974 0.21767768 0.8411543 2.0167711
#> dev.o[13,1] 1.2210623 1.5235417 0.0014941421 0.15059845 0.6186794 1.6966180
#> dev.o[14,1] 0.8967458 1.1904468 0.0008260939 0.10411622 0.4674457 1.2334310
#> dev.o[15,1] 1.1071601 1.5046946 0.0013659341 0.13145848 0.5337654 1.5079248
#> dev.o[16,1] 1.5236325 1.9505102 0.0021325952 0.19119479 0.7803749 2.0981890
#> dev.o[17,1] 2.2847252 2.4818372 0.0052327117 0.39806748 1.4250930 3.3449272
#> dev.o[18,1] 1.1424774 1.5709462 0.0010359639 0.12183543 0.5307297 1.5394803
#> dev.o[19,1] 1.6851581 1.7416886 0.0045924229 0.35334638 1.1357949 2.4587850
#> dev.o[20,1] 0.7250381 1.0280016 0.0004627470 0.06956439 0.3279362 0.9425122
#> dev.o[21,1] 1.0898010 1.5086706 0.0009043592 0.11622786 0.5130651 1.4539702
#> dev.o[1,2]  3.1887063 1.9334222 0.6307041515 1.73076501 2.8439517 4.1948526
#> dev.o[2,2]  0.7981675 1.1179932 0.0010506216 0.08235215 0.3590675 1.0657616
#> dev.o[3,2]  1.1892850 1.5054337 0.0011917280 0.14815460 0.6372120 1.6692379
#> dev.o[4,2]  0.7587825 1.0487330 0.0007346132 0.07618996 0.3434236 1.0484265
#> dev.o[5,2]  0.6501080 0.9316078 0.0005491074 0.06620975 0.2924598 0.8866004
#> dev.o[6,2]  0.9662009 1.3138392 0.0010614747 0.10823103 0.4626961 1.2914926
#> dev.o[7,2]  0.8390667 1.1983985 0.0007629653 0.08634170 0.3770590 1.1030322
#> dev.o[8,2]  0.6304144 0.8823872 0.0006741771 0.06642348 0.2994503 0.8246915
#> dev.o[9,2]  0.6452567 0.9138424 0.0005190007 0.05898926 0.2845138 0.8344876
#> dev.o[10,2] 1.6223357 1.9203324 0.0023517832 0.20532481 0.9117656 2.3688636
#> dev.o[11,2] 0.9633826 1.3877900 0.0005942917 0.10024965 0.4273175 1.2704482
#> dev.o[12,2] 0.7959089 1.1456262 0.0005503570 0.07998954 0.3491108 0.9958017
#> dev.o[13,2] 1.0070397 1.4250089 0.0012321073 0.10004669 0.4524245 1.3249341
#> dev.o[14,2] 0.8000532 1.1102697 0.0010312894 0.08088683 0.3597813 1.0638137
#> dev.o[15,2] 1.1983059 1.5904842 0.0020142857 0.15295949 0.6160505 1.6242247
#> dev.o[16,2] 1.5323485 1.9679823 0.0019403368 0.18865162 0.7905705 2.1085405
#> dev.o[17,2] 2.6001708 2.2455649 0.0391263675 0.88955415 2.0301707 3.7581939
#> dev.o[18,2] 1.0312410 1.3450017 0.0013665919 0.11218881 0.5465865 1.4493747
#> dev.o[19,2] 0.4278334 0.6151845 0.0003167565 0.03980365 0.1910982 0.5457929
#> dev.o[20,2] 0.6934293 0.9423857 0.0013155656 0.07342778 0.3354219 0.9642190
#> dev.o[21,2] 1.0680006 1.3924807 0.0014959534 0.12969690 0.5440334 1.4957087
#> dev.o[9,3]  0.8031434 1.1110664 0.0009339477 0.08043373 0.3650001 1.0444308
#> dev.o[10,3] 0.6859092 0.9695806 0.0005318669 0.06743659 0.3023108 0.9225433
#> dev.o[12,3] 1.1109131 1.4336065 0.0015326213 0.13284053 0.5602833 1.5227698
#> dev.o[13,3] 0.9877133 1.3973112 0.0006511930 0.08567996 0.4345929 1.3025425
#> dev.o[19,3] 1.3359981 1.2045092 0.0214585313 0.46385165 1.0154903 1.9067441
#> dev.o[10,4] 1.2344793 1.5032964 0.0017410188 0.17507505 0.7081806 1.7551562
#> dev.o[12,4] 0.8037942 1.0842771 0.0011030756 0.09452639 0.3892033 1.0644325
#> dev.o[13,4] 1.2890297 1.6840654 0.0013061326 0.14999724 0.6594104 1.7575939
#>                97.5%     Rhat n.eff
#> dev.o[1,1]  8.528505 1.000931  3000
#> dev.o[2,1]  3.979959 1.001939  1400
#> dev.o[3,1]  5.163079 1.015317   230
#> dev.o[4,1]  3.629947 1.001432  2100
#> dev.o[5,1]  3.725363 1.000577  3000
#> dev.o[6,1]  4.327977 1.004313   530
#> dev.o[7,1]  3.779407 1.000869  3000
#> dev.o[8,1]  3.507480 1.000866  3000
#> dev.o[9,1]  3.568125 1.007101   360
#> dev.o[10,1] 3.195597 1.003261   910
#> dev.o[11,1] 4.906933 1.002744   890
#> dev.o[12,1] 5.154938 1.030220    93
#> dev.o[13,1] 5.520087 1.004064   870
#> dev.o[14,1] 4.227254 1.000742  3000
#> dev.o[15,1] 5.277548 1.002004  1700
#> dev.o[16,1] 7.315382 1.003684  1100
#> dev.o[17,1] 8.983146 1.022521   110
#> dev.o[18,1] 5.557218 1.004387   520
#> dev.o[19,1] 6.371501 1.001711  1700
#> dev.o[20,1] 3.626120 1.001339  2400
#> dev.o[21,1] 5.669149 1.001690  1700
#> dev.o[1,2]  7.972489 1.006945  1400
#> dev.o[2,2]  4.144985 1.004306   660
#> dev.o[3,2]  5.359176 1.009050   240
#> dev.o[4,2]  3.816099 1.003037  1100
#> dev.o[5,2]  3.172401 1.000558  3000
#> dev.o[6,2]  4.879069 1.003060   780
#> dev.o[7,2]  4.038718 1.001945  1400
#> dev.o[8,2]  3.238918 1.002819   860
#> dev.o[9,2]  3.296262 1.002725   900
#> dev.o[10,2] 6.822178 1.001521  2100
#> dev.o[11,2] 5.043810 1.000774  3000
#> dev.o[12,2] 4.187402 1.003197  1200
#> dev.o[13,2] 4.989039 1.001154  3000
#> dev.o[14,2] 3.959716 1.000714  3000
#> dev.o[15,2] 5.440927 1.002556   970
#> dev.o[16,2] 7.143053 1.003424   740
#> dev.o[17,2] 8.393960 1.010732   230
#> dev.o[18,2] 4.738652 1.016503   130
#> dev.o[19,2] 2.258729 1.004934   450
#> dev.o[20,2] 3.228714 1.001538  1900
#> dev.o[21,2] 5.019969 1.004937   450
#> dev.o[9,3]  3.895831 1.002127  1200
#> dev.o[10,3] 3.290821 1.000652  3000
#> dev.o[12,3] 4.992169 1.006601   330
#> dev.o[13,3] 5.104206 1.001215  2800
#> dev.o[19,3] 4.471436 1.012127   170
#> dev.o[10,4] 5.313332 1.009190   370
#> dev.o[12,4] 3.940318 1.000999  3000
#> dev.o[13,4] 6.027898 1.001524  2000
#> 
#> $hat_par
#>                     mean         sd        2.5%         25%        50%
#> hat.par[1,1]    1.588308  0.7842345   0.3690817   0.9806083   1.481083
#> hat.par[2,1]   49.915048  4.5062949  41.2083772  46.9050813  49.984895
#> hat.par[3,1]   43.287557  4.3461301  35.3528515  40.2168180  43.154866
#> hat.par[4,1]   42.264524  4.5809089  33.7539984  39.0444391  42.021593
#> hat.par[5,1]   17.053496  2.4970547  12.1679003  15.3054571  17.014113
#> hat.par[6,1]   43.654943  4.0554506  35.5296481  40.9596246  43.668190
#> hat.par[7,1]  157.126072  7.1697773 142.6176902 152.3786412 156.989648
#> hat.par[8,1]   68.316982  5.4015526  57.7072273  64.5652326  68.319077
#> hat.par[9,1]   89.507471  4.9058160  80.1943389  86.2552140  89.335900
#> hat.par[10,1]  78.518227  3.6827288  71.1176017  76.2411922  78.615760
#> hat.par[11,1]  72.753380  5.6071198  61.4006606  69.0351933  72.706356
#> hat.par[12,1]  77.815889  4.0786463  69.6769072  75.0793957  77.955134
#> hat.par[13,1]  49.260775  4.4660298  40.7316684  46.2708354  49.306784
#> hat.par[14,1]  33.684127  4.5083847  25.7169650  30.6060589  33.369106
#> hat.par[15,1]  37.385899  4.7245596  28.4460978  34.1429893  37.236089
#> hat.par[16,1] 306.296415 13.5842040 280.2877692 297.0199736 306.291243
#> hat.par[17,1]  10.457420  2.5720079   5.9392943   8.6040624  10.328557
#> hat.par[18,1]  21.818390  3.3904760  15.4013457  19.4604804  21.778221
#> hat.par[19,1]   3.988104  1.3583922   1.7984224   3.0116964   3.832037
#> hat.par[20,1]  25.270464  3.5834962  18.5408700  22.8586773  25.141514
#> hat.par[21,1]  32.091813  4.2820987  23.7528732  29.1361995  31.911754
#> hat.par[1,2]    1.399704  0.7471643   0.3092175   0.8202201   1.302766
#> hat.par[2,2]   46.055238  4.6746072  37.0157226  42.8572043  46.009008
#> hat.par[3,2]   31.638046  3.9535004  24.0439777  28.9178211  31.582194
#> hat.par[4,2]   43.908168  5.1110603  34.4835335  40.4151947  43.673062
#> hat.par[5,2]   12.004034  2.1396041   8.0430032  10.4926945  11.892415
#> hat.par[6,2]   35.176926  3.8057335  27.5951282  32.6659991  35.155499
#> hat.par[7,2]  196.819068  9.8092409 177.8769326 190.1657559 196.788064
#> hat.par[8,2]   51.683835  4.7910382  42.5606036  48.3576603  51.562704
#> hat.par[9,2]   81.952913  5.5028786  71.4075376  78.1786553  82.002560
#> hat.par[10,2]  72.865242  3.9595704  64.9573500  70.1982077  73.010466
#> hat.par[11,2] 119.271026  8.3057013 103.6341366 113.6310303 119.041686
#> hat.par[12,2]  80.047824  4.9188298  70.1591275  76.7788921  80.114645
#> hat.par[13,2]  25.673719  4.5146013  17.3398221  22.6040777  25.459968
#> hat.par[14,2]  32.184060  4.3702293  23.7000522  29.1710146  32.090022
#> hat.par[15,2]  31.770697  4.4386184  23.8108025  28.7169746  31.537110
#> hat.par[16,2] 245.281356 12.9702488 220.3935982 236.2485394 245.034312
#> hat.par[17,2]   7.567388  2.0277433   4.1179218   6.0902718   7.368195
#> hat.par[18,2]  13.161623  2.5153262   8.6340149  11.3438325  13.051873
#> hat.par[19,2]   2.524509  0.9560512   1.0639702   1.8259296   2.371632
#> hat.par[20,2]  18.827248  3.0262538  13.5619473  16.7460748  18.670465
#> hat.par[21,2]  22.100714  3.4800623  15.7328843  19.7123403  21.951405
#> hat.par[9,3]   79.940463  5.4715570  69.3906807  76.3333563  79.801474
#> hat.par[10,3]  68.783361  4.1950974  60.3021730  66.0620605  68.870677
#> hat.par[12,3]  67.507023  4.7168826  58.5523131  64.3355205  67.408580
#> hat.par[13,3]  35.106040  5.1293539  25.2206609  31.6617570  34.871076
#> hat.par[19,3]   2.484442  0.9162668   1.1304128   1.8302979   2.343172
#> hat.par[10,4]  67.182071  4.0726492  58.9675825  64.5130301  67.237722
#> hat.par[12,4]  62.641838  4.2915545  54.3071716  59.7770894  62.582106
#> hat.par[13,4]  41.866085  4.8204451  32.4936470  38.5522003  41.765999
#>                      75%      97.5%     Rhat n.eff
#> hat.par[1,1]    2.096001   3.297560 1.010109   640
#> hat.par[2,1]   52.946650  58.544593 1.004563   590
#> hat.par[3,1]   46.270138  51.948100 1.020294   100
#> hat.par[4,1]   45.243676  51.796554 1.007339   330
#> hat.par[5,1]   18.709540  22.057966 1.005875   370
#> hat.par[6,1]   46.333030  51.499535 1.003616  1900
#> hat.par[7,1]  162.142260 171.236697 1.014188   150
#> hat.par[8,1]   71.867322  79.053228 1.006103   360
#> hat.par[9,1]   92.739287  99.542818 1.039857    56
#> hat.par[10,1]  80.962099  85.604137 1.008900   490
#> hat.par[11,1]  76.585175  83.798685 1.002515   990
#> hat.par[12,1]  80.709443  85.185774 1.039413    60
#> hat.par[13,1]  52.230180  57.949802 1.005927   450
#> hat.par[14,1]  36.474856  43.718159 1.001340  2400
#> hat.par[15,1]  40.441786  47.054721 1.010971   200
#> hat.par[16,1] 315.392398 333.548646 1.007115   420
#> hat.par[17,1]  12.136335  15.866731 1.033813    65
#> hat.par[18,1]  24.085935  28.825294 1.017231   120
#> hat.par[19,1]   4.781125   7.208238 1.004289   530
#> hat.par[20,1]  27.581500  32.621074 1.010527   200
#> hat.par[21,1]  34.898877  40.561816 1.004931   450
#> hat.par[1,2]    1.845009   3.139404 1.007857  1200
#> hat.par[2,2]   49.214628  55.637922 1.004544  1100
#> hat.par[3,2]   34.242354  39.691126 1.025505    84
#> hat.par[4,2]   47.162248  54.539999 1.006696   320
#> hat.par[5,2]   13.385703  16.259534 1.003433   680
#> hat.par[6,2]   37.646249  42.853542 1.004682   550
#> hat.par[7,2]  203.343605 217.109228 1.001094  3000
#> hat.par[8,2]   55.074874  61.403326 1.027980    78
#> hat.par[9,2]   85.577101  93.063370 1.011326   200
#> hat.par[10,2]  75.648303  80.212937 1.003658   640
#> hat.par[11,2] 124.691534 136.107836 1.005653   390
#> hat.par[12,2]  83.248159  89.790745 1.004042   570
#> hat.par[13,2]  28.558196  35.053639 1.007580   280
#> hat.par[14,2]  35.021873  41.320498 1.001550  1900
#> hat.par[15,2]  34.542548  40.952776 1.006551   330
#> hat.par[16,2] 254.063149 271.506491 1.003805   620
#> hat.par[17,2]   8.874005  12.067509 1.015047   150
#> hat.par[18,2]  14.811180  18.385537 1.020099   110
#> hat.par[19,2]   3.056149   4.745964 1.015869   130
#> hat.par[20,2]  20.758200  25.289949 1.001974  1400
#> hat.par[21,2]  24.355034  29.452245 1.013336   160
#> hat.par[9,3]   83.475999  90.462181 1.005742   430
#> hat.par[10,3]  71.655088  76.750130 1.001289  2500
#> hat.par[12,3]  70.722452  76.891887 1.023234    90
#> hat.par[13,3]  38.471849  45.758885 1.017833   120
#> hat.par[19,3]   3.022787   4.650562 1.016481   130
#> hat.par[10,4]  69.930866  74.983116 1.009170   290
#> hat.par[12,4]  65.460803  70.974870 1.003301   720
#> hat.par[13,4]  44.991517  51.623390 1.005295   560
#> 
#> $leverage_o
#>  [1] 0.9620279 0.7444184 0.7126572 0.6439405 0.6408963 0.6363404 0.7058946
#>  [8] 0.6930669 0.5924060 0.5765438 0.7983535 0.5819577 0.7273136 0.7195595
#> [15] 0.7099277 0.9620080 0.9782936 0.8288005 0.7331343 0.6413326 0.7682817
#> [22] 0.2341361 0.7593966 0.6154288 0.7265027 0.5280993 0.6375879 0.8271322
#> [29] 0.6277107 0.6452097 0.7272867 0.8958168 0.7619129 1.0019446 0.6305388
#> [36] 0.8121142 0.9922989 0.3625633 0.5546433 0.3055567 0.6079265 0.5513252
#> [43] 0.6216015 0.6631565 0.7165684 0.9872939 0.1478672 0.6130328 0.6178067
#> [50] 0.7654677
#> 
#> $sign_dev_o
#>  [1]  1  1  1 -1  1 -1 -1 -1  1  1  1 -1  1  1 -1 -1  1  1  1 -1  1 -1 -1 -1  1
#> [26] -1  1  1  1  1  1 -1 -1  1 -1  1  1 -1 -1 -1  1 -1 -1 -1  1 -1 -1 -1  1 -1
#> 
#> $phi
#>               mean        sd      2.5%        25%          50%        75%
#> phi[1] -0.13401388 0.5014454 -1.134329 -0.4722309 -0.128715724  0.2189310
#> phi[2] -0.14985614 0.9352908 -2.010423 -0.7616085 -0.120979913  0.4824477
#> phi[3] -0.01772129 0.9600138 -1.897949 -0.6604549 -0.006073649  0.6199807
#> phi[4] -1.20428940 0.8572177 -2.694767 -1.8104522 -1.263304757 -0.6464176
#> phi[5] -0.42105020 0.9623343 -2.157990 -1.0886378 -0.456509078  0.2034791
#> phi[6]  0.68509097 0.8703122 -1.249711  0.1490288  0.738350112  1.2617287
#> phi[7] -0.12473766 0.6623573 -1.465792 -0.5483722 -0.134840617  0.3169678
#> phi[8] -0.19807987 0.9135630 -1.940269 -0.8117067 -0.214115612  0.3992876
#>            97.5%     Rhat n.eff
#> phi[1] 0.8163179 1.065718    48
#> phi[2] 1.6466632 1.002681   920
#> phi[3] 1.8224277 1.034415    64
#> phi[4] 0.6885693 1.038244    89
#> phi[5] 1.6447577 1.088786    30
#> phi[6] 2.2838079 1.072620    35
#> phi[7] 1.1642828 1.050984    45
#> phi[8] 1.6641696 1.006505  1700
#> 
#> $model_assessment
#>        DIC       pD      dev n_data
#> 1 89.34295 34.29508 55.04787     50
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
#> EM[2,1]             -1.084   0.413  -1.855  -1.394  -1.083  -0.771  -0.292
#> EM[3,1]             -0.768   0.395  -1.505  -1.043  -0.789  -0.490  -0.025
#> EM[4,1]             -0.358   0.268  -0.827  -0.534  -0.373  -0.193   0.214
#> EM[5,1]             -0.468   0.304  -1.049  -0.705  -0.462  -0.226   0.060
#> EM[6,1]             -0.164   0.243  -0.678  -0.314  -0.158   0.008   0.281
#> EM[7,1]             -0.448   0.163  -0.762  -0.563  -0.450  -0.329  -0.147
#> EM[8,1]             -0.522   0.142  -0.796  -0.622  -0.526  -0.409  -0.252
#> EM[3,2]              0.316   0.443  -0.566   0.012   0.307   0.632   1.195
#> EM[4,2]              0.726   0.453  -0.060   0.397   0.686   1.030   1.670
#> EM[5,2]              0.615   0.517  -0.414   0.271   0.595   0.987   1.633
#> EM[6,2]              0.920   0.472   0.016   0.571   0.927   1.230   1.849
#> EM[7,2]              0.636   0.404  -0.209   0.375   0.632   0.928   1.371
#> EM[8,2]              0.562   0.407  -0.283   0.305   0.552   0.837   1.347
#> EM[4,3]              0.410   0.419  -0.366   0.087   0.440   0.703   1.178
#> EM[5,3]              0.300   0.486  -0.578  -0.064   0.267   0.625   1.315
#> EM[6,3]              0.604   0.463  -0.359   0.302   0.623   0.930   1.438
#> EM[7,3]              0.320   0.375  -0.418   0.057   0.333   0.579   1.024
#> EM[8,3]              0.246   0.381  -0.518  -0.022   0.255   0.509   0.990
#> EM[5,4]             -0.110   0.410  -0.840  -0.423  -0.110   0.178   0.648
#> EM[6,4]              0.194   0.341  -0.473  -0.038   0.179   0.430   0.837
#> EM[7,4]             -0.090   0.280  -0.707  -0.271  -0.057   0.113   0.378
#> EM[8,4]             -0.164   0.285  -0.807  -0.360  -0.126   0.036   0.314
#> EM[6,5]              0.304   0.348  -0.337   0.051   0.311   0.554   0.959
#> EM[7,5]              0.020   0.316  -0.562  -0.223   0.036   0.241   0.612
#> EM[8,5]             -0.054   0.309  -0.612  -0.294  -0.036   0.189   0.485
#> EM[7,6]             -0.284   0.288  -0.852  -0.491  -0.285  -0.076   0.226
#> EM[8,6]             -0.358   0.266  -0.906  -0.535  -0.360  -0.179   0.161
#> EM[8,7]             -0.074   0.153  -0.379  -0.175  -0.076   0.019   0.227
#> EM.pred[2,1]        -1.085   0.437  -1.926  -1.405  -1.082  -0.763  -0.252
#> EM.pred[3,1]        -0.763   0.425  -1.578  -1.056  -0.776  -0.472   0.039
#> EM.pred[4,1]        -0.358   0.302  -0.899  -0.560  -0.379  -0.181   0.303
#> EM.pred[5,1]        -0.466   0.343  -1.143  -0.712  -0.459  -0.204   0.105
#> EM.pred[6,1]        -0.165   0.290  -0.791  -0.337  -0.158   0.026   0.377
#> EM.pred[7,1]        -0.454   0.222  -0.915  -0.595  -0.445  -0.307  -0.016
#> EM.pred[8,1]        -0.520   0.206  -0.961  -0.643  -0.505  -0.380  -0.154
#> EM.pred[3,2]         0.318   0.467  -0.628   0.006   0.316   0.642   1.266
#> EM.pred[4,2]         0.721   0.479  -0.122   0.385   0.688   1.042   1.717
#> EM.pred[5,2]         0.616   0.541  -0.484   0.260   0.605   0.996   1.663
#> EM.pred[6,2]         0.917   0.493  -0.045   0.564   0.916   1.240   1.874
#> EM.pred[7,2]         0.635   0.433  -0.265   0.364   0.626   0.944   1.441
#> EM.pred[8,2]         0.563   0.428  -0.302   0.297   0.559   0.847   1.385
#> EM.pred[4,3]         0.410   0.445  -0.420   0.067   0.437   0.726   1.257
#> EM.pred[5,3]         0.302   0.508  -0.644  -0.065   0.269   0.644   1.336
#> EM.pred[6,3]         0.598   0.487  -0.392   0.289   0.622   0.937   1.493
#> EM.pred[7,3]         0.322   0.400  -0.449   0.038   0.340   0.590   1.100
#> EM.pred[8,3]         0.249   0.410  -0.588  -0.030   0.273   0.518   1.063
#> EM.pred[5,4]        -0.114   0.439  -0.949  -0.419  -0.109   0.197   0.681
#> EM.pred[6,4]         0.193   0.369  -0.543  -0.055   0.182   0.451   0.892
#> EM.pred[7,4]        -0.089   0.318  -0.793  -0.278  -0.049   0.130   0.439
#> EM.pred[8,4]        -0.166   0.324  -0.894  -0.376  -0.121   0.053   0.385
#> EM.pred[6,5]         0.305   0.383  -0.380   0.019   0.307   0.569   1.025
#> EM.pred[7,5]         0.020   0.352  -0.632  -0.238   0.032   0.265   0.708
#> EM.pred[8,5]        -0.051   0.344  -0.672  -0.314  -0.049   0.203   0.593
#> EM.pred[7,6]        -0.282   0.328  -0.951  -0.513  -0.281  -0.050   0.324
#> EM.pred[8,6]        -0.356   0.307  -0.975  -0.545  -0.360  -0.153   0.229
#> EM.pred[8,7]        -0.071   0.214  -0.543  -0.189  -0.069   0.054   0.346
#> SUCRA[1]             0.063   0.089   0.000   0.000   0.000   0.143   0.286
#> SUCRA[2]             0.924   0.164   0.429   0.857   1.000   1.000   1.000
#> SUCRA[3]             0.733   0.266   0.143   0.571   0.857   0.857   1.000
#> SUCRA[4]             0.422   0.243   0.000   0.286   0.429   0.571   0.857
#> SUCRA[5]             0.536   0.275   0.000   0.286   0.571   0.714   1.000
#> SUCRA[6]             0.235   0.216   0.000   0.143   0.143   0.286   0.857
#> SUCRA[7]             0.488   0.172   0.143   0.429   0.429   0.571   0.857
#> SUCRA[8]             0.599   0.169   0.286   0.429   0.571   0.714   0.857
#> abs_risk[1]          0.392   0.000   0.392   0.392   0.392   0.392   0.392
#> abs_risk[2]          0.187   0.062   0.091   0.138   0.179   0.229   0.325
#> abs_risk[3]          0.237   0.071   0.125   0.185   0.226   0.283   0.386
#> abs_risk[4]          0.313   0.058   0.220   0.274   0.307   0.347   0.444
#> abs_risk[5]          0.291   0.062   0.184   0.241   0.289   0.339   0.406
#> abs_risk[6]          0.355   0.054   0.246   0.320   0.355   0.394   0.460
#> abs_risk[7]          0.293   0.034   0.231   0.268   0.291   0.317   0.357
#> abs_risk[8]          0.277   0.028   0.225   0.257   0.276   0.300   0.333
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
#> delta[1,2]          -0.388   0.290  -0.918  -0.576  -0.395  -0.206   0.228
#> delta[2,2]          -0.389   0.256  -0.848  -0.566  -0.396  -0.221   0.154
#> delta[3,2]          -0.508   0.195  -0.936  -0.631  -0.495  -0.365  -0.190
#> delta[4,2]          -0.490   0.165  -0.818  -0.604  -0.484  -0.381  -0.177
#> delta[5,2]          -0.461   0.204  -0.882  -0.589  -0.448  -0.321  -0.083
#> delta[6,2]          -0.388   0.198  -0.754  -0.520  -0.390  -0.261   0.029
#> delta[7,2]          -0.504   0.152  -0.808  -0.610  -0.501  -0.393  -0.224
#> delta[8,2]          -0.427   0.174  -0.781  -0.538  -0.426  -0.306  -0.089
#> delta[9,2]          -0.471   0.180  -0.848  -0.589  -0.463  -0.341  -0.168
#> delta[10,2]         -0.318   0.308  -0.864  -0.534  -0.334  -0.133   0.348
#> delta[11,2]         -0.206   0.248  -0.748  -0.364  -0.195  -0.030   0.242
#> delta[12,2]         -0.328   0.285  -0.826  -0.527  -0.351  -0.157   0.290
#> delta[13,2]         -1.124   0.392  -1.842  -1.416  -1.137  -0.819  -0.377
#> delta[14,2]         -0.104   0.182  -0.491  -0.212  -0.103   0.012   0.244
#> delta[15,2]         -0.107   0.246  -0.601  -0.263  -0.109   0.063   0.383
#> delta[16,2]         -0.392   0.130  -0.647  -0.480  -0.391  -0.302  -0.136
#> delta[17,2]         -0.419   0.287  -0.999  -0.607  -0.416  -0.226   0.141
#> delta[18,2]         -0.501   0.317  -1.085  -0.737  -0.511  -0.247   0.065
#> delta[19,2]         -0.495   0.333  -1.173  -0.734  -0.477  -0.231   0.064
#> delta[20,2]         -0.427   0.204  -0.835  -0.562  -0.421  -0.287  -0.018
#> delta[21,2]         -0.565   0.192  -0.975  -0.691  -0.553  -0.422  -0.235
#> delta[9,3]          -0.564   0.171  -0.910  -0.679  -0.561  -0.435  -0.255
#> delta[10,3]         -0.486   0.318  -1.113  -0.721  -0.479  -0.235   0.062
#> delta[12,3]         -0.393   0.305  -0.964  -0.626  -0.377  -0.162   0.147
#> delta[13,3]         -0.809   0.373  -1.510  -1.060  -0.834  -0.549  -0.099
#> delta[19,3]         -0.484   0.222  -0.963  -0.612  -0.474  -0.329  -0.092
#> delta[10,4]         -0.491   0.218  -0.971  -0.628  -0.478  -0.339  -0.112
#> delta[12,4]         -0.388   0.203  -0.766  -0.528  -0.395  -0.254   0.046
#> delta[13,4]         -0.245   0.278  -0.890  -0.405  -0.227  -0.060   0.253
#> dev.o[1,1]           2.360   2.488   0.006   0.535   1.660   3.440   8.529
#> dev.o[2,1]           0.787   1.105   0.001   0.073   0.357   1.056   3.980
#> dev.o[3,1]           1.205   1.497   0.001   0.140   0.640   1.744   5.163
#> dev.o[4,1]           0.694   1.069   0.001   0.065   0.292   0.879   3.630
#> dev.o[5,1]           0.729   1.023   0.001   0.073   0.339   0.978   3.725
#> dev.o[6,1]           0.913   1.211   0.001   0.098   0.447   1.241   4.328
#> dev.o[7,1]           0.723   1.047   0.001   0.066   0.343   0.924   3.779
#> dev.o[8,1]           0.695   0.977   0.001   0.073   0.310   0.887   3.507
#> dev.o[9,1]           0.746   0.987   0.001   0.088   0.359   1.014   3.568
#> dev.o[10,1]          0.586   0.883   0.001   0.052   0.236   0.757   3.196
#> dev.o[11,1]          0.925   1.296   0.001   0.096   0.419   1.203   4.907
#> dev.o[12,1]          1.352   1.476   0.002   0.218   0.841   2.017   5.155
#> dev.o[13,1]          1.221   1.524   0.001   0.151   0.619   1.697   5.520
#> dev.o[14,1]          0.897   1.190   0.001   0.104   0.467   1.233   4.227
#> dev.o[15,1]          1.107   1.505   0.001   0.131   0.534   1.508   5.278
#> dev.o[16,1]          1.524   1.951   0.002   0.191   0.780   2.098   7.315
#> dev.o[17,1]          2.285   2.482   0.005   0.398   1.425   3.345   8.983
#> dev.o[18,1]          1.142   1.571   0.001   0.122   0.531   1.539   5.557
#> dev.o[19,1]          1.685   1.742   0.005   0.353   1.136   2.459   6.372
#> dev.o[20,1]          0.725   1.028   0.000   0.070   0.328   0.943   3.626
#> dev.o[21,1]          1.090   1.509   0.001   0.116   0.513   1.454   5.669
#> dev.o[1,2]           3.189   1.933   0.631   1.731   2.844   4.195   7.972
#> dev.o[2,2]           0.798   1.118   0.001   0.082   0.359   1.066   4.145
#> dev.o[3,2]           1.189   1.505   0.001   0.148   0.637   1.669   5.359
#> dev.o[4,2]           0.759   1.049   0.001   0.076   0.343   1.048   3.816
#> dev.o[5,2]           0.650   0.932   0.001   0.066   0.292   0.887   3.172
#> dev.o[6,2]           0.966   1.314   0.001   0.108   0.463   1.291   4.879
#> dev.o[7,2]           0.839   1.198   0.001   0.086   0.377   1.103   4.039
#> dev.o[8,2]           0.630   0.882   0.001   0.066   0.299   0.825   3.239
#> dev.o[9,2]           0.645   0.914   0.001   0.059   0.285   0.834   3.296
#> dev.o[10,2]          1.622   1.920   0.002   0.205   0.912   2.369   6.822
#> dev.o[11,2]          0.963   1.388   0.001   0.100   0.427   1.270   5.044
#> dev.o[12,2]          0.796   1.146   0.001   0.080   0.349   0.996   4.187
#> dev.o[13,2]          1.007   1.425   0.001   0.100   0.452   1.325   4.989
#> dev.o[14,2]          0.800   1.110   0.001   0.081   0.360   1.064   3.960
#> dev.o[15,2]          1.198   1.590   0.002   0.153   0.616   1.624   5.441
#> dev.o[16,2]          1.532   1.968   0.002   0.189   0.791   2.109   7.143
#> dev.o[17,2]          2.600   2.246   0.039   0.890   2.030   3.758   8.394
#> dev.o[18,2]          1.031   1.345   0.001   0.112   0.547   1.449   4.739
#> dev.o[19,2]          0.428   0.615   0.000   0.040   0.191   0.546   2.259
#> dev.o[20,2]          0.693   0.942   0.001   0.073   0.335   0.964   3.229
#> dev.o[21,2]          1.068   1.392   0.001   0.130   0.544   1.496   5.020
#> dev.o[9,3]           0.803   1.111   0.001   0.080   0.365   1.044   3.896
#> dev.o[10,3]          0.686   0.970   0.001   0.067   0.302   0.923   3.291
#> dev.o[12,3]          1.111   1.434   0.002   0.133   0.560   1.523   4.992
#> dev.o[13,3]          0.988   1.397   0.001   0.086   0.435   1.303   5.104
#> dev.o[19,3]          1.336   1.205   0.021   0.464   1.015   1.907   4.471
#> dev.o[10,4]          1.234   1.503   0.002   0.175   0.708   1.755   5.313
#> dev.o[12,4]          0.804   1.084   0.001   0.095   0.389   1.064   3.940
#> dev.o[13,4]          1.289   1.684   0.001   0.150   0.659   1.758   6.028
#> effectiveness[1,1]   0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> effectiveness[2,1]   0.716   0.451   0.000   0.000   1.000   1.000   1.000
#> effectiveness[3,1]   0.206   0.404   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,1]   0.004   0.060   0.000   0.000   0.000   0.000   0.000
#> effectiveness[5,1]   0.046   0.210   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,1]   0.002   0.045   0.000   0.000   0.000   0.000   0.000
#> effectiveness[7,1]   0.005   0.073   0.000   0.000   0.000   0.000   0.000
#> effectiveness[8,1]   0.020   0.141   0.000   0.000   0.000   0.000   0.000
#> effectiveness[1,2]   0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> effectiveness[2,2]   0.181   0.385   0.000   0.000   0.000   0.000   1.000
#> effectiveness[3,2]   0.422   0.494   0.000   0.000   0.000   1.000   1.000
#> effectiveness[4,2]   0.059   0.236   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,2]   0.166   0.372   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,2]   0.034   0.181   0.000   0.000   0.000   0.000   1.000
#> effectiveness[7,2]   0.043   0.204   0.000   0.000   0.000   0.000   1.000
#> effectiveness[8,2]   0.095   0.293   0.000   0.000   0.000   0.000   1.000
#> effectiveness[1,3]   0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> effectiveness[2,3]   0.041   0.198   0.000   0.000   0.000   0.000   1.000
#> effectiveness[3,3]   0.096   0.294   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,3]   0.176   0.381   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,3]   0.228   0.420   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,3]   0.030   0.170   0.000   0.000   0.000   0.000   1.000
#> effectiveness[7,3]   0.125   0.331   0.000   0.000   0.000   0.000   1.000
#> effectiveness[8,3]   0.305   0.460   0.000   0.000   0.000   1.000   1.000
#> effectiveness[1,4]   0.000   0.018   0.000   0.000   0.000   0.000   0.000
#> effectiveness[2,4]   0.017   0.128   0.000   0.000   0.000   0.000   0.000
#> effectiveness[3,4]   0.069   0.253   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,4]   0.146   0.353   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,4]   0.124   0.330   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,4]   0.054   0.226   0.000   0.000   0.000   0.000   1.000
#> effectiveness[7,4]   0.281   0.450   0.000   0.000   0.000   1.000   1.000
#> effectiveness[8,4]   0.309   0.462   0.000   0.000   0.000   1.000   1.000
#> effectiveness[1,5]   0.003   0.055   0.000   0.000   0.000   0.000   0.000
#> effectiveness[2,5]   0.022   0.147   0.000   0.000   0.000   0.000   0.000
#> effectiveness[3,5]   0.079   0.269   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,5]   0.176   0.381   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,5]   0.120   0.325   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,5]   0.087   0.282   0.000   0.000   0.000   0.000   1.000
#> effectiveness[7,5]   0.330   0.470   0.000   0.000   0.000   1.000   1.000
#> effectiveness[8,5]   0.184   0.387   0.000   0.000   0.000   0.000   1.000
#> effectiveness[1,6]   0.062   0.241   0.000   0.000   0.000   0.000   1.000
#> effectiveness[2,6]   0.013   0.115   0.000   0.000   0.000   0.000   0.000
#> effectiveness[3,6]   0.050   0.218   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,6]   0.215   0.411   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,6]   0.175   0.380   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,6]   0.240   0.427   0.000   0.000   0.000   0.000   1.000
#> effectiveness[7,6]   0.166   0.372   0.000   0.000   0.000   0.000   1.000
#> effectiveness[8,6]   0.079   0.270   0.000   0.000   0.000   0.000   1.000
#> effectiveness[1,7]   0.305   0.460   0.000   0.000   0.000   1.000   1.000
#> effectiveness[2,7]   0.006   0.077   0.000   0.000   0.000   0.000   0.000
#> effectiveness[3,7]   0.062   0.241   0.000   0.000   0.000   0.000   1.000
#> effectiveness[4,7]   0.158   0.365   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,7]   0.089   0.285   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,7]   0.323   0.468   0.000   0.000   0.000   1.000   1.000
#> effectiveness[7,7]   0.049   0.216   0.000   0.000   0.000   0.000   1.000
#> effectiveness[8,7]   0.008   0.087   0.000   0.000   0.000   0.000   0.000
#> effectiveness[1,8]   0.630   0.483   0.000   0.000   1.000   1.000   1.000
#> effectiveness[2,8]   0.004   0.063   0.000   0.000   0.000   0.000   0.000
#> effectiveness[3,8]   0.016   0.127   0.000   0.000   0.000   0.000   0.000
#> effectiveness[4,8]   0.067   0.251   0.000   0.000   0.000   0.000   1.000
#> effectiveness[5,8]   0.052   0.223   0.000   0.000   0.000   0.000   1.000
#> effectiveness[6,8]   0.230   0.421   0.000   0.000   0.000   0.000   1.000
#> effectiveness[7,8]   0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> effectiveness[8,8]   0.000   0.000   0.000   0.000   0.000   0.000   0.000
#> hat.par[1,1]         1.588   0.784   0.369   0.981   1.481   2.096   3.298
#> hat.par[2,1]        49.915   4.506  41.208  46.905  49.985  52.947  58.545
#> hat.par[3,1]        43.288   4.346  35.353  40.217  43.155  46.270  51.948
#> hat.par[4,1]        42.265   4.581  33.754  39.044  42.022  45.244  51.797
#> hat.par[5,1]        17.053   2.497  12.168  15.305  17.014  18.710  22.058
#> hat.par[6,1]        43.655   4.055  35.530  40.960  43.668  46.333  51.500
#> hat.par[7,1]       157.126   7.170 142.618 152.379 156.990 162.142 171.237
#> hat.par[8,1]        68.317   5.402  57.707  64.565  68.319  71.867  79.053
#> hat.par[9,1]        89.507   4.906  80.194  86.255  89.336  92.739  99.543
#> hat.par[10,1]       78.518   3.683  71.118  76.241  78.616  80.962  85.604
#> hat.par[11,1]       72.753   5.607  61.401  69.035  72.706  76.585  83.799
#> hat.par[12,1]       77.816   4.079  69.677  75.079  77.955  80.709  85.186
#> hat.par[13,1]       49.261   4.466  40.732  46.271  49.307  52.230  57.950
#> hat.par[14,1]       33.684   4.508  25.717  30.606  33.369  36.475  43.718
#> hat.par[15,1]       37.386   4.725  28.446  34.143  37.236  40.442  47.055
#> hat.par[16,1]      306.296  13.584 280.288 297.020 306.291 315.392 333.549
#> hat.par[17,1]       10.457   2.572   5.939   8.604  10.329  12.136  15.867
#> hat.par[18,1]       21.818   3.390  15.401  19.460  21.778  24.086  28.825
#> hat.par[19,1]        3.988   1.358   1.798   3.012   3.832   4.781   7.208
#> hat.par[20,1]       25.270   3.583  18.541  22.859  25.142  27.581  32.621
#> hat.par[21,1]       32.092   4.282  23.753  29.136  31.912  34.899  40.562
#> hat.par[1,2]         1.400   0.747   0.309   0.820   1.303   1.845   3.139
#> hat.par[2,2]        46.055   4.675  37.016  42.857  46.009  49.215  55.638
#> hat.par[3,2]        31.638   3.954  24.044  28.918  31.582  34.242  39.691
#> hat.par[4,2]        43.908   5.111  34.484  40.415  43.673  47.162  54.540
#> hat.par[5,2]        12.004   2.140   8.043  10.493  11.892  13.386  16.260
#> hat.par[6,2]        35.177   3.806  27.595  32.666  35.155  37.646  42.854
#> hat.par[7,2]       196.819   9.809 177.877 190.166 196.788 203.344 217.109
#> hat.par[8,2]        51.684   4.791  42.561  48.358  51.563  55.075  61.403
#> hat.par[9,2]        81.953   5.503  71.408  78.179  82.003  85.577  93.063
#> hat.par[10,2]       72.865   3.960  64.957  70.198  73.010  75.648  80.213
#> hat.par[11,2]      119.271   8.306 103.634 113.631 119.042 124.692 136.108
#> hat.par[12,2]       80.048   4.919  70.159  76.779  80.115  83.248  89.791
#> hat.par[13,2]       25.674   4.515  17.340  22.604  25.460  28.558  35.054
#> hat.par[14,2]       32.184   4.370  23.700  29.171  32.090  35.022  41.320
#> hat.par[15,2]       31.771   4.439  23.811  28.717  31.537  34.543  40.953
#> hat.par[16,2]      245.281  12.970 220.394 236.249 245.034 254.063 271.506
#> hat.par[17,2]        7.567   2.028   4.118   6.090   7.368   8.874  12.068
#> hat.par[18,2]       13.162   2.515   8.634  11.344  13.052  14.811  18.386
#> hat.par[19,2]        2.525   0.956   1.064   1.826   2.372   3.056   4.746
#> hat.par[20,2]       18.827   3.026  13.562  16.746  18.670  20.758  25.290
#> hat.par[21,2]       22.101   3.480  15.733  19.712  21.951  24.355  29.452
#> hat.par[9,3]        79.940   5.472  69.391  76.333  79.801  83.476  90.462
#> hat.par[10,3]       68.783   4.195  60.302  66.062  68.871  71.655  76.750
#> hat.par[12,3]       67.507   4.717  58.552  64.336  67.409  70.722  76.892
#> hat.par[13,3]       35.106   5.129  25.221  31.662  34.871  38.472  45.759
#> hat.par[19,3]        2.484   0.916   1.130   1.830   2.343   3.023   4.651
#> hat.par[10,4]       67.182   4.073  58.968  64.513  67.238  69.931  74.983
#> hat.par[12,4]       62.642   4.292  54.307  59.777  62.582  65.461  70.975
#> hat.par[13,4]       41.866   4.820  32.494  38.552  41.766  44.992  51.623
#> phi[1]              -0.134   0.501  -1.134  -0.472  -0.129   0.219   0.816
#> phi[2]              -0.150   0.935  -2.010  -0.762  -0.121   0.482   1.647
#> phi[3]              -0.018   0.960  -1.898  -0.660  -0.006   0.620   1.822
#> phi[4]              -1.204   0.857  -2.695  -1.810  -1.263  -0.646   0.689
#> phi[5]              -0.421   0.962  -2.158  -1.089  -0.457   0.203   1.645
#> phi[6]               0.685   0.870  -1.250   0.149   0.738   1.262   2.284
#> phi[7]              -0.125   0.662  -1.466  -0.548  -0.135   0.317   1.164
#> phi[8]              -0.198   0.914  -1.940  -0.812  -0.214   0.399   1.664
#> tau                  0.127   0.080   0.019   0.064   0.110   0.181   0.310
#> totresdev.o         55.048   9.088  38.838  48.819  54.705  60.754  74.948
#> deviance           582.331  13.389 558.342 572.891 581.712 590.860 610.992
#>                     Rhat n.eff
#> EM[2,1]            1.011   200
#> EM[3,1]            1.084    31
#> EM[4,1]            1.035   710
#> EM[5,1]            1.112    23
#> EM[6,1]            1.060    48
#> EM[7,1]            1.215    14
#> EM[8,1]            1.085    30
#> EM[3,2]            1.032    68
#> EM[4,2]            1.027   130
#> EM[5,2]            1.038   110
#> EM[6,2]            1.047    62
#> EM[7,2]            1.014   320
#> EM[8,2]            1.007  1700
#> EM[4,3]            1.089    28
#> EM[5,3]            1.127    25
#> EM[6,3]            1.081    32
#> EM[7,3]            1.016   140
#> EM[8,3]            1.036    62
#> EM[5,4]            1.080    36
#> EM[6,4]            1.043   120
#> EM[7,4]            1.108    26
#> EM[8,4]            1.061    62
#> EM[6,5]            1.219    13
#> EM[7,5]            1.124    22
#> EM[8,5]            1.094    27
#> EM[7,6]            1.142    19
#> EM[8,6]            1.099    26
#> EM[8,7]            1.044    62
#> EM.pred[2,1]       1.010   210
#> EM.pred[3,1]       1.068    37
#> EM.pred[4,1]       1.035   600
#> EM.pred[5,1]       1.093    28
#> EM.pred[6,1]       1.059    58
#> EM.pred[7,1]       1.108    24
#> EM.pred[8,1]       1.047    57
#> EM.pred[3,2]       1.031    71
#> EM.pred[4,2]       1.027   150
#> EM.pred[5,2]       1.038   110
#> EM.pred[6,2]       1.044    69
#> EM.pred[7,2]       1.015   270
#> EM.pred[8,2]       1.007  2400
#> EM.pred[4,3]       1.076    31
#> EM.pred[5,3]       1.121    27
#> EM.pred[6,3]       1.073    35
#> EM.pred[7,3]       1.019   120
#> EM.pred[8,3]       1.030    75
#> EM.pred[5,4]       1.072    43
#> EM.pred[6,4]       1.036   170
#> EM.pred[7,4]       1.091    32
#> EM.pred[8,4]       1.048    90
#> EM.pred[6,5]       1.171    16
#> EM.pred[7,5]       1.108    26
#> EM.pred[8,5]       1.077    34
#> EM.pred[7,6]       1.115    23
#> EM.pred[8,6]       1.078    35
#> EM.pred[8,7]       1.031   140
#> SUCRA[1]           1.019   150
#> SUCRA[2]           1.043   110
#> SUCRA[3]           1.066    44
#> SUCRA[4]           1.026    82
#> SUCRA[5]           1.104    25
#> SUCRA[6]           1.107    28
#> SUCRA[7]           1.035    82
#> SUCRA[8]           1.020   120
#> abs_risk[1]        1.000     1
#> abs_risk[2]        1.011   190
#> abs_risk[3]        1.082    32
#> abs_risk[4]        1.032  1400
#> abs_risk[5]        1.104    25
#> abs_risk[6]        1.066    48
#> abs_risk[7]        1.209    14
#> abs_risk[8]        1.083    30
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
#> delta[1,2]         1.030  3000
#> delta[2,2]         1.021   720
#> delta[3,2]         1.181    16
#> delta[4,2]         1.034    68
#> delta[5,2]         1.146    19
#> delta[6,2]         1.068    36
#> delta[7,2]         1.061    39
#> delta[8,2]         1.125    21
#> delta[9,2]         1.168    16
#> delta[10,2]        1.041   420
#> delta[11,2]        1.049    70
#> delta[12,2]        1.037   360
#> delta[13,2]        1.023   100
#> delta[14,2]        1.024   110
#> delta[15,2]        1.056    46
#> delta[16,2]        1.017   200
#> delta[17,2]        1.130    20
#> delta[18,2]        1.109    24
#> delta[19,2]        1.076    33
#> delta[20,2]        1.109    23
#> delta[21,2]        1.056    41
#> delta[9,3]         1.104    25
#> delta[10,3]        1.100    26
#> delta[12,3]        1.147    18
#> delta[13,3]        1.134    21
#> delta[19,3]        1.161    17
#> delta[10,4]        1.178    16
#> delta[12,4]        1.067    36
#> delta[13,4]        1.077    51
#> dev.o[1,1]         1.001  3000
#> dev.o[2,1]         1.002  1400
#> dev.o[3,1]         1.015   230
#> dev.o[4,1]         1.001  2100
#> dev.o[5,1]         1.001  3000
#> dev.o[6,1]         1.004   530
#> dev.o[7,1]         1.001  3000
#> dev.o[8,1]         1.001  3000
#> dev.o[9,1]         1.007   360
#> dev.o[10,1]        1.003   910
#> dev.o[11,1]        1.003   890
#> dev.o[12,1]        1.030    93
#> dev.o[13,1]        1.004   870
#> dev.o[14,1]        1.001  3000
#> dev.o[15,1]        1.002  1700
#> dev.o[16,1]        1.004  1100
#> dev.o[17,1]        1.023   110
#> dev.o[18,1]        1.004   520
#> dev.o[19,1]        1.002  1700
#> dev.o[20,1]        1.001  2400
#> dev.o[21,1]        1.002  1700
#> dev.o[1,2]         1.007  1400
#> dev.o[2,2]         1.004   660
#> dev.o[3,2]         1.009   240
#> dev.o[4,2]         1.003  1100
#> dev.o[5,2]         1.001  3000
#> dev.o[6,2]         1.003   780
#> dev.o[7,2]         1.002  1400
#> dev.o[8,2]         1.003   860
#> dev.o[9,2]         1.003   900
#> dev.o[10,2]        1.002  2100
#> dev.o[11,2]        1.001  3000
#> dev.o[12,2]        1.003  1200
#> dev.o[13,2]        1.001  3000
#> dev.o[14,2]        1.001  3000
#> dev.o[15,2]        1.003   970
#> dev.o[16,2]        1.003   740
#> dev.o[17,2]        1.011   230
#> dev.o[18,2]        1.017   130
#> dev.o[19,2]        1.005   450
#> dev.o[20,2]        1.002  1900
#> dev.o[21,2]        1.005   450
#> dev.o[9,3]         1.002  1200
#> dev.o[10,3]        1.001  3000
#> dev.o[12,3]        1.007   330
#> dev.o[13,3]        1.001  2800
#> dev.o[19,3]        1.012   170
#> dev.o[10,4]        1.009   370
#> dev.o[12,4]        1.001  3000
#> dev.o[13,4]        1.002  2000
#> effectiveness[1,1] 1.000     1
#> effectiveness[2,1] 1.048    56
#> effectiveness[3,1] 1.046    70
#> effectiveness[4,1] 1.070  1700
#> effectiveness[5,1] 1.021   510
#> effectiveness[6,1] 1.294   500
#> effectiveness[7,1] 1.143   530
#> effectiveness[8,1] 1.036   650
#> effectiveness[1,2] 1.000     1
#> effectiveness[2,2] 1.064    55
#> effectiveness[3,2] 1.012   180
#> effectiveness[4,2] 1.069   120
#> effectiveness[5,2] 1.029   120
#> effectiveness[6,2] 1.141    91
#> effectiveness[7,2] 1.005  2400
#> effectiveness[8,2] 1.037   160
#> effectiveness[1,3] 1.000     1
#> effectiveness[2,3] 1.003  3000
#> effectiveness[3,3] 1.008   710
#> effectiveness[4,3] 1.009   360
#> effectiveness[5,3] 1.052    57
#> effectiveness[6,3] 1.105   140
#> effectiveness[7,3] 1.046   100
#> effectiveness[8,3] 1.001  2700
#> effectiveness[1,4] 1.291  3000
#> effectiveness[2,4] 1.018  1600
#> effectiveness[3,4] 1.009   840
#> effectiveness[4,4] 1.012   320
#> effectiveness[5,4] 1.038   120
#> effectiveness[6,4] 1.064   140
#> effectiveness[7,4] 1.011   220
#> effectiveness[8,4] 1.004   570
#> effectiveness[1,5] 1.109  1300
#> effectiveness[2,5] 1.009  2300
#> effectiveness[3,5] 1.001  3000
#> effectiveness[4,5] 1.007   480
#> effectiveness[5,5] 1.013   360
#> effectiveness[6,5] 1.017   360
#> effectiveness[7,5] 1.012   180
#> effectiveness[8,5] 1.016   200
#> effectiveness[1,6] 1.028   300
#> effectiveness[2,6] 1.034  1100
#> effectiveness[3,6] 1.018   550
#> effectiveness[4,6] 1.023   130
#> effectiveness[5,6] 1.154    25
#> effectiveness[6,6] 1.006   400
#> effectiveness[7,6] 1.077    48
#> effectiveness[8,6] 1.091    72
#> effectiveness[1,7] 1.004   590
#> effectiveness[2,7] 1.136   500
#> effectiveness[3,7] 1.140    56
#> effectiveness[4,7] 1.003  1100
#> effectiveness[5,7] 1.067    89
#> effectiveness[6,7] 1.018   130
#> effectiveness[7,7] 1.051   200
#> effectiveness[8,7] 1.063   940
#> effectiveness[1,8] 1.010   210
#> effectiveness[2,8] 1.135   750
#> effectiveness[3,8] 1.042   700
#> effectiveness[4,8] 1.118    62
#> effectiveness[5,8] 1.105    86
#> effectiveness[6,8] 1.030    96
#> effectiveness[7,8] 1.000     1
#> effectiveness[8,8] 1.000     1
#> hat.par[1,1]       1.010   640
#> hat.par[2,1]       1.005   590
#> hat.par[3,1]       1.020   100
#> hat.par[4,1]       1.007   330
#> hat.par[5,1]       1.006   370
#> hat.par[6,1]       1.004  1900
#> hat.par[7,1]       1.014   150
#> hat.par[8,1]       1.006   360
#> hat.par[9,1]       1.040    56
#> hat.par[10,1]      1.009   490
#> hat.par[11,1]      1.003   990
#> hat.par[12,1]      1.039    60
#> hat.par[13,1]      1.006   450
#> hat.par[14,1]      1.001  2400
#> hat.par[15,1]      1.011   200
#> hat.par[16,1]      1.007   420
#> hat.par[17,1]      1.034    65
#> hat.par[18,1]      1.017   120
#> hat.par[19,1]      1.004   530
#> hat.par[20,1]      1.011   200
#> hat.par[21,1]      1.005   450
#> hat.par[1,2]       1.008  1200
#> hat.par[2,2]       1.005  1100
#> hat.par[3,2]       1.026    84
#> hat.par[4,2]       1.007   320
#> hat.par[5,2]       1.003   680
#> hat.par[6,2]       1.005   550
#> hat.par[7,2]       1.001  3000
#> hat.par[8,2]       1.028    78
#> hat.par[9,2]       1.011   200
#> hat.par[10,2]      1.004   640
#> hat.par[11,2]      1.006   390
#> hat.par[12,2]      1.004   570
#> hat.par[13,2]      1.008   280
#> hat.par[14,2]      1.002  1900
#> hat.par[15,2]      1.007   330
#> hat.par[16,2]      1.004   620
#> hat.par[17,2]      1.015   150
#> hat.par[18,2]      1.020   110
#> hat.par[19,2]      1.016   130
#> hat.par[20,2]      1.002  1400
#> hat.par[21,2]      1.013   160
#> hat.par[9,3]       1.006   430
#> hat.par[10,3]      1.001  2500
#> hat.par[12,3]      1.023    90
#> hat.par[13,3]      1.018   120
#> hat.par[19,3]      1.016   130
#> hat.par[10,4]      1.009   290
#> hat.par[12,4]      1.003   720
#> hat.par[13,4]      1.005   560
#> phi[1]             1.066    48
#> phi[2]             1.003   920
#> phi[3]             1.034    64
#> phi[4]             1.038    89
#> phi[5]             1.089    30
#> phi[6]             1.073    35
#> phi[7]             1.051    45
#> phi[8]             1.007  1700
#> tau                1.180    16
#> totresdev.o        1.010   220
#> deviance           1.006   380
#> 
#> For each parameter, n.eff is a crude measure of effective sample size,
#> and Rhat is the potential scale reduction factor (at convergence, Rhat=1).
#> 
#> DIC info (using the rule: pV = var(deviance)/2)
#> pV = 89.2 and DIC = 671.5
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
#>                     mean        sd        2.5%          25%         50%
#> EM.pred[2,1] -1.08467173 0.4366460 -1.92567716 -1.405290358 -1.08221752
#> EM.pred[3,1] -0.76278577 0.4249281 -1.57842817 -1.056405141 -0.77602661
#> EM.pred[4,1] -0.35826332 0.3019428 -0.89884304 -0.560428866 -0.37909568
#> EM.pred[5,1] -0.46641888 0.3426449 -1.14321954 -0.712028142 -0.45863031
#> EM.pred[6,1] -0.16537910 0.2902573 -0.79103642 -0.336880557 -0.15794130
#> EM.pred[7,1] -0.45376710 0.2218684 -0.91463228 -0.594623706 -0.44531576
#> EM.pred[8,1] -0.52021399 0.2058627 -0.96130865 -0.643472122 -0.50472958
#> EM.pred[3,2]  0.31810794 0.4671492 -0.62832166  0.006128761  0.31612475
#> EM.pred[4,2]  0.72112138 0.4790516 -0.12195748  0.385315184  0.68813385
#> EM.pred[5,2]  0.61555866 0.5412715 -0.48389661  0.259852549  0.60510467
#> EM.pred[6,2]  0.91737637 0.4932709 -0.04509473  0.564048326  0.91560703
#> EM.pred[7,2]  0.63519538 0.4329326 -0.26466195  0.363978589  0.62617047
#> EM.pred[8,2]  0.56330776 0.4280803 -0.30179360  0.296987244  0.55948814
#> EM.pred[4,3]  0.41029245 0.4452615 -0.41956300  0.066621969  0.43698346
#> EM.pred[5,3]  0.30168024 0.5083531 -0.64354181 -0.065460513  0.26919065
#> EM.pred[6,3]  0.59835994 0.4865012 -0.39238862  0.288616334  0.62161357
#> EM.pred[7,3]  0.32196913 0.4003096 -0.44864106  0.037704480  0.34031578
#> EM.pred[8,3]  0.24868898 0.4095823 -0.58848844 -0.030403113  0.27269293
#> EM.pred[5,4] -0.11363478 0.4393768 -0.94902286 -0.418703833 -0.10864098
#> EM.pred[6,4]  0.19349742 0.3692715 -0.54268308 -0.054882707  0.18228677
#> EM.pred[7,4] -0.08935595 0.3182274 -0.79336467 -0.278362863 -0.04909552
#> EM.pred[8,4] -0.16639389 0.3235570 -0.89445146 -0.376124766 -0.12083654
#> EM.pred[6,5]  0.30508463 0.3827299 -0.38028753  0.019129320  0.30676358
#> EM.pred[7,5]  0.02036938 0.3520192 -0.63172581 -0.237951815  0.03209049
#> EM.pred[8,5] -0.05149945 0.3435155 -0.67230549 -0.313587677 -0.04903220
#> EM.pred[7,6] -0.28229929 0.3283443 -0.95137103 -0.513045364 -0.28108661
#> EM.pred[8,6] -0.35550837 0.3066574 -0.97507459 -0.545037994 -0.35989144
#> EM.pred[8,7] -0.07121551 0.2144172 -0.54317925 -0.189018213 -0.06945708
#>                      75%       97.5%     Rhat n.eff
#> EM.pred[2,1] -0.76264963 -0.25200316 1.009974   210
#> EM.pred[3,1] -0.47222992  0.03859115 1.068102    37
#> EM.pred[4,1] -0.18099786  0.30265944 1.034531   600
#> EM.pred[5,1] -0.20385173  0.10493648 1.092979    28
#> EM.pred[6,1]  0.02600996  0.37698612 1.058613    58
#> EM.pred[7,1] -0.30699800 -0.01640321 1.107870    24
#> EM.pred[8,1] -0.37966309 -0.15445076 1.046884    57
#> EM.pred[3,2]  0.64167750  1.26642116 1.031439    71
#> EM.pred[4,2]  1.04201849  1.71745765 1.026601   150
#> EM.pred[5,2]  0.99596997  1.66283983 1.037571   110
#> EM.pred[6,2]  1.24029732  1.87386915 1.044227    69
#> EM.pred[7,2]  0.94386127  1.44075421 1.014778   270
#> EM.pred[8,2]  0.84650087  1.38491449 1.007386  2400
#> EM.pred[4,3]  0.72619881  1.25663849 1.075774    31
#> EM.pred[5,3]  0.64446784  1.33632323 1.121491    27
#> EM.pred[6,3]  0.93650911  1.49309381 1.073461    35
#> EM.pred[7,3]  0.58975913  1.10034479 1.019383   120
#> EM.pred[8,3]  0.51814764  1.06348730 1.029922    75
#> EM.pred[5,4]  0.19687621  0.68112413 1.072329    43
#> EM.pred[6,4]  0.45124821  0.89192630 1.035902   170
#> EM.pred[7,4]  0.13005653  0.43897014 1.091393    32
#> EM.pred[8,4]  0.05297870  0.38466718 1.047550    90
#> EM.pred[6,5]  0.56940195  1.02459451 1.171388    16
#> EM.pred[7,5]  0.26526662  0.70770698 1.108290    26
#> EM.pred[8,5]  0.20292219  0.59273309 1.077345    34
#> EM.pred[7,6] -0.05042549  0.32394290 1.114790    23
#> EM.pred[8,6] -0.15261386  0.22850555 1.078127    35
#> EM.pred[8,7]  0.05412859  0.34588736 1.030687   140
#> 
#> $tau
#>        mean          sd        2.5%         25%         50%         75% 
#>  0.12695209  0.07984776  0.01913408  0.06389442  0.11046393  0.18133506 
#>       97.5%        Rhat       n.eff 
#>  0.31041361  1.17967052 16.00000000 
#> 
#> $delta
#>                   mean        sd       2.5%        25%        50%         75%
#> delta[1,2]  -0.3877233 0.2903228 -0.9181648 -0.5763986 -0.3947342 -0.20605910
#> delta[2,2]  -0.3892433 0.2559630 -0.8477256 -0.5657029 -0.3964507 -0.22058204
#> delta[3,2]  -0.5081054 0.1946517 -0.9360039 -0.6310886 -0.4953308 -0.36545817
#> delta[4,2]  -0.4898420 0.1646070 -0.8179129 -0.6036128 -0.4840714 -0.38078336
#> delta[5,2]  -0.4606613 0.2041363 -0.8819322 -0.5891083 -0.4480068 -0.32080119
#> delta[6,2]  -0.3877150 0.1979125 -0.7537991 -0.5199508 -0.3897765 -0.26083552
#> delta[7,2]  -0.5038366 0.1522428 -0.8075554 -0.6097587 -0.5006116 -0.39262137
#> delta[8,2]  -0.4265632 0.1739327 -0.7810351 -0.5380534 -0.4260225 -0.30573578
#> delta[9,2]  -0.4713464 0.1797038 -0.8484650 -0.5887129 -0.4634319 -0.34070994
#> delta[10,2] -0.3183508 0.3077716 -0.8639165 -0.5344079 -0.3337421 -0.13291966
#> delta[11,2] -0.2064307 0.2476471 -0.7475996 -0.3641260 -0.1945881 -0.03001666
#> delta[12,2] -0.3284730 0.2847451 -0.8260801 -0.5267727 -0.3513912 -0.15743797
#> delta[13,2] -1.1237872 0.3921329 -1.8423701 -1.4155614 -1.1371032 -0.81949432
#> delta[14,2] -0.1041357 0.1816044 -0.4911982 -0.2115361 -0.1025628  0.01237067
#> delta[15,2] -0.1072727 0.2458582 -0.6010313 -0.2626005 -0.1093913  0.06349389
#> delta[16,2] -0.3917362 0.1303968 -0.6465950 -0.4799467 -0.3910737 -0.30158047
#> delta[17,2] -0.4186279 0.2869045 -0.9992253 -0.6074110 -0.4164897 -0.22551540
#> delta[18,2] -0.5007266 0.3168685 -1.0845513 -0.7368813 -0.5111483 -0.24694264
#> delta[19,2] -0.4951453 0.3328552 -1.1734328 -0.7342970 -0.4768681 -0.23089139
#> delta[20,2] -0.4273595 0.2038348 -0.8348447 -0.5621759 -0.4211233 -0.28666202
#> delta[21,2] -0.5653500 0.1922115 -0.9754544 -0.6914354 -0.5531994 -0.42221487
#> delta[9,3]  -0.5635279 0.1708487 -0.9103961 -0.6794193 -0.5607780 -0.43481908
#> delta[10,3] -0.4864293 0.3183214 -1.1127655 -0.7205162 -0.4786876 -0.23522605
#> delta[12,3] -0.3933413 0.3046332 -0.9635469 -0.6264132 -0.3772520 -0.16237328
#> delta[13,3] -0.8092069 0.3729214 -1.5099691 -1.0601985 -0.8339221 -0.54894197
#> delta[19,3] -0.4840963 0.2216821 -0.9630369 -0.6115137 -0.4736482 -0.32896763
#> delta[10,4] -0.4908771 0.2182532 -0.9706643 -0.6276205 -0.4776945 -0.33902357
#> delta[12,4] -0.3878595 0.2034742 -0.7660479 -0.5280040 -0.3949774 -0.25398051
#> delta[13,4] -0.2454372 0.2776055 -0.8898658 -0.4051891 -0.2268704 -0.05974512
#>                   97.5%     Rhat n.eff
#> delta[1,2]   0.22763262 1.030054  3000
#> delta[2,2]   0.15388589 1.020500   720
#> delta[3,2]  -0.19037443 1.180641    16
#> delta[4,2]  -0.17737594 1.033684    68
#> delta[5,2]  -0.08315796 1.145707    19
#> delta[6,2]   0.02907250 1.067833    36
#> delta[7,2]  -0.22391250 1.061376    39
#> delta[8,2]  -0.08906307 1.125327    21
#> delta[9,2]  -0.16798835 1.167531    16
#> delta[10,2]  0.34777418 1.041403   420
#> delta[11,2]  0.24234019 1.048723    70
#> delta[12,2]  0.28950900 1.037389   360
#> delta[13,2] -0.37713035 1.022734   100
#> delta[14,2]  0.24447726 1.023948   110
#> delta[15,2]  0.38319219 1.056110    46
#> delta[16,2] -0.13555833 1.016550   200
#> delta[17,2]  0.14094266 1.130036    20
#> delta[18,2]  0.06509985 1.108979    24
#> delta[19,2]  0.06357894 1.076074    33
#> delta[20,2] -0.01780548 1.108869    23
#> delta[21,2] -0.23503792 1.055865    41
#> delta[9,3]  -0.25475241 1.103667    25
#> delta[10,3]  0.06233988 1.099587    26
#> delta[12,3]  0.14739978 1.147086    18
#> delta[13,3] -0.09850386 1.134353    21
#> delta[19,3] -0.09167523 1.161473    17
#> delta[10,4] -0.11196418 1.178308    16
#> delta[12,4]  0.04643525 1.066615    36
#> delta[13,4]  0.25281161 1.076931    51
#> 
#> $heter_prior
#> [1] 0 1 1
#> 
#> $SUCRA
#>                mean         sd      2.5%       25%       50%       75%
#> SUCRA[1] 0.06261905 0.08940405 0.0000000 0.0000000 0.0000000 0.1428571
#> SUCRA[2] 0.92409524 0.16351506 0.4285714 0.8571429 1.0000000 1.0000000
#> SUCRA[3] 0.73261905 0.26592230 0.1428571 0.5714286 0.8571429 0.8571429
#> SUCRA[4] 0.42219048 0.24289783 0.0000000 0.2857143 0.4285714 0.5714286
#> SUCRA[5] 0.53604762 0.27499770 0.0000000 0.2857143 0.5714286 0.7142857
#> SUCRA[6] 0.23533333 0.21582244 0.0000000 0.1428571 0.1428571 0.2857143
#> SUCRA[7] 0.48847619 0.17166497 0.1428571 0.4285714 0.4285714 0.5714286
#> SUCRA[8] 0.59861905 0.16914595 0.2857143 0.4285714 0.5714286 0.7142857
#>              97.5%     Rhat n.eff
#> SUCRA[1] 0.2857143 1.019417   150
#> SUCRA[2] 1.0000000 1.043263   110
#> SUCRA[3] 1.0000000 1.066187    44
#> SUCRA[4] 0.8571429 1.025810    82
#> SUCRA[5] 1.0000000 1.103648    25
#> SUCRA[6] 0.8571429 1.106835    28
#> SUCRA[7] 0.8571429 1.035246    82
#> SUCRA[8] 0.8571429 1.020052   120
#> 
#> $effectiveness
#>                            mean         sd 2.5% 25% 50% 75% 97.5%     Rhat
#> effectiveness[1,1] 0.0000000000 0.00000000    0   0   0   0     0 1.000000
#> effectiveness[2,1] 0.7163333333 0.45085213    0   0   1   1     1 1.047623
#> effectiveness[3,1] 0.2060000000 0.40449789    0   0   0   0     1 1.045550
#> effectiveness[4,1] 0.0036666667 0.06045197    0   0   0   0     0 1.070399
#> effectiveness[5,1] 0.0463333333 0.21024103    0   0   0   0     1 1.021224
#> effectiveness[6,1] 0.0020000000 0.04468406    0   0   0   0     0 1.293699
#> effectiveness[7,1] 0.0053333333 0.07284681    0   0   0   0     0 1.142756
#> effectiveness[8,1] 0.0203333333 0.14116137    0   0   0   0     0 1.036332
#> effectiveness[1,2] 0.0000000000 0.00000000    0   0   0   0     0 1.000000
#> effectiveness[2,2] 0.1806666667 0.38480590    0   0   0   0     1 1.064111
#> effectiveness[3,2] 0.4223333333 0.49401340    0   0   0   1     1 1.011777
#> effectiveness[4,2] 0.0590000000 0.23566398    0   0   0   0     1 1.069140
#> effectiveness[5,2] 0.1656666667 0.37184313    0   0   0   0     1 1.028889
#> effectiveness[6,2] 0.0340000000 0.18125935    0   0   0   0     1 1.140966
#> effectiveness[7,2] 0.0433333333 0.20364032    0   0   0   0     1 1.004532
#> effectiveness[8,2] 0.0950000000 0.29326382    0   0   0   0     1 1.036533
#> effectiveness[1,3] 0.0000000000 0.00000000    0   0   0   0     0 1.000000
#> effectiveness[2,3] 0.0410000000 0.19832325    0   0   0   0     1 1.003477
#> effectiveness[3,3] 0.0956666667 0.29418260    0   0   0   0     1 1.007530
#> effectiveness[4,3] 0.1756666667 0.38059976    0   0   0   0     1 1.008984
#> effectiveness[5,3] 0.2280000000 0.41961255    0   0   0   0     1 1.052437
#> effectiveness[6,3] 0.0296666667 0.16969430    0   0   0   0     1 1.105048
#> effectiveness[7,3] 0.1253333333 0.33115169    0   0   0   0     1 1.045651
#> effectiveness[8,3] 0.3046666667 0.46034284    0   0   0   1     1 1.001233
#> effectiveness[1,4] 0.0003333333 0.01825742    0   0   0   0     0 1.290904
#> effectiveness[2,4] 0.0166666667 0.12804044    0   0   0   0     0 1.017840
#> effectiveness[3,4] 0.0690000000 0.25349639    0   0   0   0     1 1.008636
#> effectiveness[4,4] 0.1456666667 0.35283053    0   0   0   0     1 1.012346
#> effectiveness[5,4] 0.1240000000 0.32963650    0   0   0   0     1 1.038237
#> effectiveness[6,4] 0.0540000000 0.22605538    0   0   0   0     1 1.063766
#> effectiveness[7,4] 0.2810000000 0.44956242    0   0   0   1     1 1.010854
#> effectiveness[8,4] 0.3093333333 0.46229586    0   0   0   1     1 1.004052
#> effectiveness[1,5] 0.0030000000 0.05469915    0   0   0   0     0 1.109096
#> effectiveness[2,5] 0.0220000000 0.14670779    0   0   0   0     0 1.009428
#> effectiveness[3,5] 0.0786666667 0.26926268    0   0   0   0     1 1.000814
#> effectiveness[4,5] 0.1756666667 0.38059976    0   0   0   0     1 1.006823
#> effectiveness[5,5] 0.1196666667 0.32462545    0   0   0   0     1 1.012857
#> effectiveness[6,5] 0.0870000000 0.28188204    0   0   0   0     1 1.016640
#> effectiveness[7,5] 0.3303333333 0.47041151    0   0   0   1     1 1.012379
#> effectiveness[8,5] 0.1836666667 0.38727667    0   0   0   0     1 1.016475
#> effectiveness[1,6] 0.0616666667 0.24058924    0   0   0   0     1 1.027579
#> effectiveness[2,6] 0.0133333333 0.11471679    0   0   0   0     0 1.033780
#> effectiveness[3,6] 0.0500000000 0.21798128    0   0   0   0     1 1.018216
#> effectiveness[4,6] 0.2146666667 0.41065935    0   0   0   0     1 1.022580
#> effectiveness[5,6] 0.1750000000 0.38003045    0   0   0   0     1 1.153686
#> effectiveness[6,6] 0.2403333333 0.42735711    0   0   0   0     1 1.006371
#> effectiveness[7,6] 0.1656666667 0.37184313    0   0   0   0     1 1.077274
#> effectiveness[8,6] 0.0793333333 0.27030337    0   0   0   0     1 1.090666
#> effectiveness[1,7] 0.3046666667 0.46034284    0   0   0   1     1 1.003880
#> effectiveness[2,7] 0.0060000000 0.07723981    0   0   0   0     0 1.135746
#> effectiveness[3,7] 0.0620000000 0.24119575    0   0   0   0     1 1.139778
#> effectiveness[4,7] 0.1583333333 0.36511413    0   0   0   0     1 1.002903
#> effectiveness[5,7] 0.0890000000 0.28479121    0   0   0   0     1 1.066507
#> effectiveness[6,7] 0.3233333333 0.46782672    0   0   0   1     1 1.017587
#> effectiveness[7,7] 0.0490000000 0.21590400    0   0   0   0     1 1.050630
#> effectiveness[8,7] 0.0076666667 0.08723775    0   0   0   0     0 1.063109
#> effectiveness[1,8] 0.6303333333 0.48279490    0   0   1   1     1 1.010085
#> effectiveness[2,8] 0.0040000000 0.06312946    0   0   0   0     0 1.134929
#> effectiveness[3,8] 0.0163333333 0.12677505    0   0   0   0     0 1.041609
#> effectiveness[4,8] 0.0673333333 0.25064017    0   0   0   0     1 1.118431
#> effectiveness[5,8] 0.0523333333 0.22273548    0   0   0   0     1 1.105253
#> effectiveness[6,8] 0.2296666667 0.42068858    0   0   0   0     1 1.030221
#> effectiveness[7,8] 0.0000000000 0.00000000    0   0   0   0     0 1.000000
#> effectiveness[8,8] 0.0000000000 0.00000000    0   0   0   0     0 1.000000
#>                    n.eff
#> effectiveness[1,1]     1
#> effectiveness[2,1]    56
#> effectiveness[3,1]    70
#> effectiveness[4,1]  1700
#> effectiveness[5,1]   510
#> effectiveness[6,1]   500
#> effectiveness[7,1]   530
#> effectiveness[8,1]   650
#> effectiveness[1,2]     1
#> effectiveness[2,2]    55
#> effectiveness[3,2]   180
#> effectiveness[4,2]   120
#> effectiveness[5,2]   120
#> effectiveness[6,2]    91
#> effectiveness[7,2]  2400
#> effectiveness[8,2]   160
#> effectiveness[1,3]     1
#> effectiveness[2,3]  3000
#> effectiveness[3,3]   710
#> effectiveness[4,3]   360
#> effectiveness[5,3]    57
#> effectiveness[6,3]   140
#> effectiveness[7,3]   100
#> effectiveness[8,3]  2700
#> effectiveness[1,4]  3000
#> effectiveness[2,4]  1600
#> effectiveness[3,4]   840
#> effectiveness[4,4]   320
#> effectiveness[5,4]   120
#> effectiveness[6,4]   140
#> effectiveness[7,4]   220
#> effectiveness[8,4]   570
#> effectiveness[1,5]  1300
#> effectiveness[2,5]  2300
#> effectiveness[3,5]  3000
#> effectiveness[4,5]   480
#> effectiveness[5,5]   360
#> effectiveness[6,5]   360
#> effectiveness[7,5]   180
#> effectiveness[8,5]   200
#> effectiveness[1,6]   300
#> effectiveness[2,6]  1100
#> effectiveness[3,6]   550
#> effectiveness[4,6]   130
#> effectiveness[5,6]    25
#> effectiveness[6,6]   400
#> effectiveness[7,6]    48
#> effectiveness[8,6]    72
#> effectiveness[1,7]   590
#> effectiveness[2,7]   500
#> effectiveness[3,7]    56
#> effectiveness[4,7]  1100
#> effectiveness[5,7]    89
#> effectiveness[6,7]   130
#> effectiveness[7,7]   200
#> effectiveness[8,7]   940
#> effectiveness[1,8]   210
#> effectiveness[2,8]   750
#> effectiveness[3,8]   700
#> effectiveness[4,8]    62
#> effectiveness[5,8]    86
#> effectiveness[6,8]    96
#> effectiveness[7,8]     1
#> effectiveness[8,8]     1
#> 
#> $abs_risk
#>                  mean         sd       2.5%       25%       50%       75%
#> abs_risk[1] 0.3916667 0.00000000 0.39166667 0.3916667 0.3916667 0.3916667
#> abs_risk[2] 0.1867084 0.06238818 0.09149723 0.1377293 0.1789715 0.2294464
#> abs_risk[3] 0.2372068 0.07061976 0.12511121 0.1849982 0.2262691 0.2829065
#> abs_risk[4] 0.3132276 0.05757604 0.21976102 0.2739773 0.3071813 0.3466949
#> abs_risk[5] 0.2912603 0.06150182 0.18400188 0.2414348 0.2886604 0.3393202
#> abs_risk[6] 0.3553794 0.05445212 0.24625669 0.3199438 0.3547925 0.3936380
#> abs_risk[7] 0.2926059 0.03356584 0.23107423 0.2682706 0.2911388 0.3166287
#> abs_risk[8] 0.2773508 0.02843270 0.22503851 0.2569495 0.2756408 0.2995886
#>                 97.5%     Rhat n.eff
#> abs_risk[1] 0.3916667 1.000000     1
#> abs_risk[2] 0.3247202 1.011044   190
#> abs_risk[3] 0.3858265 1.081666    32
#> abs_risk[4] 0.4435922 1.032370  1400
#> abs_risk[5] 0.4059571 1.103975    25
#> abs_risk[6] 0.4601576 1.066158    48
#> abs_risk[7] 0.3571566 1.209327    14
#> abs_risk[8] 0.3334597 1.082997    30
#> 
#> $base_risk
#> [1] 0.3916667
#> 
#> attr(,"class")
#> [1] "run_model"
# }
```
