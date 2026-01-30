# Visualising study percentage contributions against a covariate

A scatter plot of the study percentage contributions against the values
of a continuous study-level covariate for the treatment effects of
comparisons referring to the basic parameters, functional parameters or
both. Contributions on the estimated regression coefficients are also
presented. Study percentage contributions are based on the proposed
methodology of Donegan and colleagues (2018).

## Usage

``` r
covar_contribution_plot(
  contr_res,
  comparisons = "basic",
  drug_names,
  upper_limit = 100,
  name_x_axis = NULL,
  axis_title_size = 14,
  axis_text_size = 14,
  strip_text_size = 14,
  subtitle_size = 14,
  label_size = 4,
  seq_by = 0.1,
  percentage = FALSE
)
```

## Arguments

- contr_res:

  An object of S3 class
  [`study_perc_contrib`](https://loukiaspin.github.io/rnmamod/reference/study_perc_contrib.md).
  This object contains the study percentage contributions to the
  treatment effects (or regression coefficients, if relevant) of all
  possible comparisons in the network. See 'Value' in
  [`study_perc_contrib`](https://loukiaspin.github.io/rnmamod/reference/study_perc_contrib.md).

- comparisons:

  Character string indicating the type of comparisons to plot, with
  possible values: `"basic"`, `"functional"`, or `"all"` to consider
  only the basic parameters, only the functional parameters, or both,
  respectively. The default argument is `"basic"`.

- drug_names:

  A vector of labels with the name of the interventions in the order
  they appear in the argument `contr_res`. If `drug_names` is not
  defined, the order of the interventions as they appear in `contr_res`
  is used, instead.

- upper_limit:

  A positive number to define the upper bound of range of percentage
  values for the y-axis. The default argument is 100.

- name_x_axis:

  Text for the x axis title through the `labs` function found in the
  R-package [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- axis_title_size:

  A positive integer for the font size of x axis title.
  `axis_title_size` determines the axis.title (and legend.title)
  arguments found in the theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- axis_text_size:

  A positive integer for the font size of axis text (both axes).
  `axis_text_size` determines the axis.text (and legend.text) arguments
  found in the theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- strip_text_size:

  A positive integer for the font size of strip text in facets.
  `strip_text_size` determines the strip.text argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- subtitle_size:

  A positive integer for the font size of subtitle. `subtitle_size`
  determines the plot.subtitle argument found in the theme's properties
  in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- label_size:

  A positive integer for the font size of labels appearing on each data
  point. `label_size` determines the size argument found in the geom's
  aesthetic properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- seq_by:

  A positive integer for the sequence of values in the x-axis. `seq_by`
  appears in the arguments breaks and labels found in the
  scale_x_continuous aesthetic properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- percentage:

  Logical with values `TRUE` if the covariate is measured in per cent
  and `FALSE` otherwise. The default argument is `FALSE`.

## Value

If interest lies only on the study percentage contributions to the
summary treatment effects of all possible pairwise comparisons, the
function returns one plot named 'plot_treat'. If interest lies also on
the study percentage contributions to the regression coefficient(s), the
function returns also the plot named 'plot_reg'.

## Details

A panel of scatter plots is returned on the study percentage
contributions to the treatment effects (and also regression
coefficients, if relevant) against a continuous covariate for each
comparison defined by the argument `comparisons`; namely, only those
referring to the basic or functional parameters or all possible pairwise
comparisons. Blue and red points indicate the studies investigating the
corresponding comparisons directly and indirectly, respectively. Each
point displays the number of the corresponding study in the dataset.

If interest also lies on the study percentage contributions to the
regression coefficients, the regression coefficients can be determined
to be common across the comparisons, independent or exchangeable and
this assumption is specified in the
[`study_perc_contrib`](https://loukiaspin.github.io/rnmamod/reference/study_perc_contrib.md)
function.

## References

Donegan S, Dias S, Tudur-Smith C, Marinho V, Welton NJ. Graphs of study
contributions and covariate distributions for network meta-regression.
*Res Synth Methods* 2018;**9**(2):243–60. doi: 10.1002/jrsm.1292

## See also

[`study_perc_contrib`](https://loukiaspin.github.io/rnmamod/reference/study_perc_contrib.md)

## Author

Loukia M. Spineli

## Examples

``` r
if (FALSE) { # \dontrun{
data("nma.fluoride.donegan2018")

# Get study contributions to random-effects network meta-regression
# results under the assumption of independent treatment-by-covariate
# interaction
res <- study_perc_contrib(study_name = nma.fluoride.donegan2018$study,
                          base_t = nma.fluoride.donegan2018$t1,
                          exp_t = nma.fluoride.donegan2018$t2,
                          ref_t = 1,
                          obs_se = nma.fluoride.donegan2018$SE,
                          obs_cov = nma.fluoride.donegan2018$Cov,
                          covar = nma.fluoride.donegan2018$year,
                          covar_assum = "independent",
                          model = "RE",
                          tau = sqrt(0.03))

# Covariate-contribution plot on the basic parameters only
covar_contribution_plot(contr_res = res,
                        comparisons = "basic",
                        drug_names = c("NT", "PL", "DE", "RI", "GE", "VA"),
                        upper_limit = 15,
                        name_x_axis = "Randomisation year",
                        seq_by = 10)
} # }
```
