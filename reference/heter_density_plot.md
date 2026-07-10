# Visualising the density of two prior distributions for the heterogeneity parameter

Creating the density plot of two prior distributions for the
between-study variance (log-normal and location-scale t distributions)
or between-study standard deviation (half-normal distribution).

## Usage

``` r
heter_density_plot(
  distr,
  heter_prior1,
  heter_prior2,
  heter1 = "tau",
  heter2 = "tau",
  caption = FALSE,
  x_axis_name = TRUE,
  y_axis_name = TRUE,
  title_name = NULL,
  axis_title_size = 13,
  axis_text_size = 13,
  legend_title_size = 13,
  legend_text_size = 13
)
```

## Arguments

- distr:

  Character string indicating the prior distribution. Set `distr` equal
  to one of the following: `"lognormal"`, `"logt"`, or `"halfnormal"`,
  which refers to a log-normal, location-scale, or half-normal
  distribution, respectively.

- heter_prior1:

  A numeric vector with two values for the first prior distribution: 1)
  the mean value and 2) the standard deviation. When
  `distr = "halfnormal"`, the first value should zero and the second a
  non-negative value referring to the scale parameter of the
  distribution.

- heter_prior2:

  A numeric vector with two values for the second prior distribution: 1)
  the mean value and 2) the standard deviation. When
  `distr = "halfnormal"`, the first value should zero and the second a
  non-negative value referring to the scale parameter of the
  distribution.

- heter1:

  Character string indicating the heterogeneity parameter for
  `heter_prior1`. Set `heter1` equal to one of the following: `"tau"`,
  or `"tau_omega"`, which refers to a between-study heterogeneity or
  between-design heterogeneity (inconsistency), respectively. This
  argument is relevant only when `distr = "lognormal"` or
  `distr = "logt"`. The default is `"tau"`.

- heter2:

  Character string indicating the heterogeneity parameter for
  `heter_prior2`. Set `heter2` equal to one of the following: `"tau"`,
  or `"tau_omega"`, which refers to a between-study heterogeneity or
  between-design heterogeneity (inconsistency), respectively. This
  argument is relevant only when `distr = "lognormal"` or
  `distr = "logt"`. The default is `"tau"`.

- caption:

  Logical to indicate whether to report a caption at the bottom right of
  the plot. It is relevant only when `distr = "lognormal"` and
  `distr = "logt"`. The default is `FALSE` (do not report).

- x_axis_name:

  Logical to indicate whether to present the title of x-axis
  ('Between-study standard deviation'). The default is `TRUE` (report).

- y_axis_name:

  Logical to indicate whether to present the title of y-axis
  ('Density'). The default is `TRUE` (report).

- title_name:

  Text for the title of the plot. `title_name` determines the labs
  argument of the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- axis_title_size:

  A positive integer for the font size of axis title. `axis_title_size`
  determines the axis.title argument found in the theme's properties in
  the R-package [ggplot2](https://CRAN.R-project.org/package=ggplot2).
  The default option is 13.

- axis_text_size:

  A positive integer for the font size of axis text. `axis_text_size`
  determines the axis.text argument found in the theme's properties in
  the R-package [ggplot2](https://CRAN.R-project.org/package=ggplot2).
  The default option is 13.

- legend_title_size:

  A positive integer for the font size of legend title.
  `legend_text_size` determines the legend.text argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2). The default
  option is 13.

- legend_text_size:

  A positive integer for the font size of legend text.
  `legend_text_size` determines the legend.text argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2). The default
  option is 13.

## Value

A plot with the density of two selected prior distributions for the
heterogeneity parameter. Two different colours are used to discern the
distributions. A legend is also created with the name and
hyper-parameters of the selected prior distributions. The filled area
under each curved indicates the values up to the median of the
corresponding distribution. The x-axis present the 0.1

`heter_density_plot` also returns a table with the percentiles of each
distribution.

## Details

Use this function to inspect the shape of the distribution and the range
of between-study variance or standard deviation values before you define
the argument `heter_prior` in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md))
to run random-effects network meta-analysis.

Turner et al. (2012), Turner et al. (2015), and Rhodes et al. (2016)
provide predictive prior distributions for the between-study variance
for a binary outcome, measured in the log-odds ratio scale, and a
continuous outcome, measured in the standardised mean difference scale,
respectively.

## References

Rhodes KM, Turner RM, Higgins JP. Predictive distributions were
developed for the extent of heterogeneity in meta-analyses of continuous
outcome data. *J Clin Epidemiol* 2015;**68**(1):52–60. doi:
10.1016/j.jclinepi.2014.08.012

Turner RM, Jackson D, Wei Y, Thompson SG, Higgins JP. Predictive
distributions for between-study heterogeneity and simple methods for
their application in Bayesian meta-analysis. *Stat Med*
2015;**34**(6):984–98. doi: 10.1002/sim.6381

Turner RM, Davey J, Clarke MJ, Thompson SG, Higgins JP. Predicting the
extent of heterogeneity in meta-analysis, using empirical data from the
Cochrane Database of Systematic Reviews. *Int J Epidemiol*
2012;**41**(3):818–27. doi: 10.1093/ije/dys041

## See also

[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)

## Author

Loukia M. Spineli

## Examples

``` r

if (FALSE) { # \dontrun{
## Two empirical priors for between-study variance of log odds ratio.
heter_density_plot(distr = "lognormal",
                   heter_prior1 = c(-2.56, 1.74),  # General healthcare setting
                   heter_prior2 = c(-1.83, 1.52))  # Pain and pharma vs. placebo/ctrl

## Two empirical priors for between-study variance of standardised mean
## difference.
heter_density_plot(distr = "logt",
                   heter_prior1 = c(-3.44, 2.59),  # General healthcare setting
                   heter_prior2 = c(-0.60, 2.61))  # Pain and pharma vs. placebo/ctrl for cancer

## Two half-normal prior distributions for between-study standard deviation
heter_density_plot(distr = "halfnormal",
                   heter_prior1 = c(0, 1),
                   heter_prior2 = c(0, 0.5))
} # }
```
