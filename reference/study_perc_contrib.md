# Calculate study percentage contributions to summary treatment effects or regression coefficients

A data-frame on the percentage contributions of each study to every
possible pairwise comparison in the investigated network. Study
percentage contributions are based on the proposed methodology of
Donegan and colleagues (2018).

## Usage

``` r
study_perc_contrib(
  study_name,
  base_t,
  exp_t,
  ref_t,
  obs_se,
  obs_cov = NULL,
  covar,
  covar_assum,
  model,
  tau = NULL,
  tau_beta = NULL
)
```

## Arguments

- study_name:

  A vector of labels with the names of the studies included in the
  investigated network. For multi-arm studies, the study name should
  appear as many times as the number of possible comparisons among the
  compared treatments.

- base_t:

  A vector of numbers referring to the treatment identifier for the
  baseline arm (comparator) of each study.

- exp_t:

  A vector of numbers referring to the treatment identifier for the
  experimental arm of each study.

- ref_t:

  A scalar for the selected reference treatment in the network.

- obs_se:

  A vector of numbers referring to the estimated standard error of the
  treatment effect of each study. For multi-arm studies, the standard
  error of the treatment effect of each possible comparison among the
  compared treatments should be included.

- obs_cov:

  A vector of numbers referring to the covariance in the block
  variance-covariance matrix of the estimated treatments effects for the
  multi-arm studies only. This argument should be left unspecified if
  there are no multi-arm studies in the network.

- covar:

  A vector of numbers referring to a continuous covariate that indicates
  a study characteristic or summary patient characteristic.

- covar_assum:

  Character string indicating the structure of the
  treatment-by-covariate interaction, as described in Cooper et
  al. (2009) if interest also lies on the study percentage contributions
  to the estimated regression coefficients. Set `covar_assumption` equal
  to `"no"`, `"exchangeable"`, `"independent"`, or `"common"`. When
  `covar_assum = "no"`, only the study percentage contributions to the
  summary treatment effects will be calculated. There is no default
  argument.

- model:

  Character string indicating the analysis model with values `"RE"`, or
  `"FE"` for the random-effects and fixed-effect model, respectively.
  There is no default argument.

- tau:

  A scalar referring to the estimated between-study standard deviation
  obtained from network meta-analysis, if `covar_assum = "no"`, or
  network meta-regression for a specific treatment-by-covariate
  interaction assumption. This argument should be left unspecified when
  `model = "FE"`.

- tau_beta:

  A scalar referring to the estimated standard deviation of the
  exchangeable regression coefficients obtained from network
  meta-regression with exchangeable treatment-by-covariate interaction.
  This argument should be left unspecified when `covar_assum` is not
  "exchangeable". There is no default argument.

## Value

A list of the following two elements:

- perc_contribute:

  A data-frame with four columns referring to the study name, baseline
  and experimental treatment arm, and the covariate and as many columns
  as the number of possible comparisons with the study percentage
  contributions to summary treatment effects referring to the basic
  parameters, and followed by the functional parameters. If interest
  lies also on the regression coefficients, extra columns appear
  referring to the study percentage contributions to the regression
  coefficients specified from the argument `covar_assum`.

- covar_assumption:

  The estimated summary odd ratio in the logarithmic scale when
  `measure = "RR"` or `measure = "RD"`.

## Details

Note that the columns referring to the study percentage contributions to
summary treatment effects are indicated by the letter 'd' with two
numbers in decreasing order for the comparison: the first number refers
to the comparator and the second number refers to the experimental
treatment of the comparison. If interest lies also on the regression
coefficients, the correspoding columns are indicated by 'beta'.

The function centers the covariate to the mean but presents the original
version of the covariate.

## References

Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing
between-study heterogeneity and inconsistency in mixed treatment
comparisons: Application to stroke prevention treatments in individuals
with non-rheumatic atrial fibrillation. *Stat Med*
2009;**28**(14):1861–81. doi: 10.1002/sim.3594

Donegan S, Dias S, Tudur-Smith C, Marinho V, Welton NJ. Graphs of study
contributions and covariate distributions for network meta-regression.
*Res Synth Methods* 2018;**9**(2):243–60. doi: 10.1002/jrsm.1292

## Author

Loukia M. Spineli

## Examples

``` r
data("nma.malaria.donegan2018")

# Get study contributions to fixed-effect network meta-regression
# results under the assumption of independent treatment-by-covariate
# interaction
study_perc_contrib(study_name = nma.malaria.donegan2018$s,
                   base_t = nma.malaria.donegan2018$t1,
                   exp_t = nma.malaria.donegan2018$t2,
                   ref_t = 1,
                   obs_se = nma.malaria.donegan2018$se,
                   covar = nma.malaria.donegan2018$x,
                   covar_assum = "independent",
                   model = "FE")
#> $perc_contribute
#>          study_name comparator_arm experimental_arm    d12    d13    d23 beta12
#> study 1           1              1                2  0.389  0.001  0.245  0.320
#> study 2           2              1                2  2.970  0.018  1.872  2.532
#> study 3           3              1                2 20.874  3.708 10.953 24.132
#> study 4           4              1                2  2.259  0.033  1.397  1.415
#> study 5           5              1                2  0.848  0.156  0.442  1.041
#> study 6           6              1                2  4.286  0.667  2.303  3.923
#> study 7           7              1                2  1.336  0.026  0.822  0.771
#> study 8           8              1                2  5.075  0.077  3.224  4.829
#> study 9           9              1                2  0.665  0.002  0.418  0.543
#> study 10         10              1                2  5.007  0.044  3.163  4.416
#> study 11         11              1                2  0.265  0.005  0.163  0.148
#> study 12         12              1                2  4.957  0.046  3.133  4.403
#> study 13         13              1                2 25.697  0.072 16.144 20.996
#> study 14         14              1                2  1.820  0.018  1.151  1.628
#> study 15         15              1                3  0.296  0.851  0.303  0.363
#> study 16         16              1                3  0.978  3.061  1.143  1.201
#> study 17         17              1                3  0.020  0.865  0.484  0.022
#> study 18         18              1                3  9.455 32.222 12.559 11.605
#> study 19         19              1                3  1.140 49.865 29.320  1.545
#> study 20         20              1                3  0.002  0.280  0.159  0.001
#> study 21         21              1                3  0.446  1.493  0.577  0.547
#> study 22         22              1                3  0.580  2.174  0.884  0.711
#> study 23         23              2                3  7.652  3.107  6.577 10.132
#> study 24         24              2                3  2.984  1.211  2.565  2.774
#>          beta13 beta23 covariate
#> study 1   0.002  0.213  3.845000
#> study 2   0.025  1.682  3.500000
#> study 3   4.480 13.418 30.000000
#> study 4   0.038  0.972  6.200000
#> study 5   0.189  0.581 31.050000
#> study 6   0.805  2.134 26.500000
#> study 7   0.030  0.535  6.800000
#> study 8   0.099  3.177  2.333333
#> study 9   0.003  0.362  3.930000
#> study 10  0.059  2.924  3.150000
#> study 11  0.006  0.103  7.000000
#> study 12  0.062  2.914  3.075000
#> study 13  0.117 14.004  3.916667
#> study 14  0.024  1.077  3.000000
#> study 15  1.171  0.478 34.470000
#> study 16  3.681  1.462 30.950000
#> study 17  0.557  0.358  5.500000
#> study 18 33.519 12.869 27.900000
#> study 19 46.264 27.464  2.850000
#> study 20  0.209  0.130  4.500000
#> study 21  1.602  0.620 28.500000
#> study 22  1.901  0.694 25.000000
#> study 23  4.047  9.285 32.250000
#> study 24  1.110  2.543 27.333333
#> 
#> $covar_assumption
#> [1] "independent"
#> 
#> attr(,"class")
#> [1] "study_perc_contrib"
```
