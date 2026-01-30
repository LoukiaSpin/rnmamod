# Leverage plot

Plots the leverage against the square root of the posterior mean of
residual deviance of the trial-arms under the model of interest.

## Usage

``` r
leverage_plot(net, drug_names, title)
```

## Arguments

- net:

  An object of S3 class
  [`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md),
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
  or
  [`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md).
  See 'Value' in
  [`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md),
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
  or
  [`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md).

- drug_names:

  A vector of labels with the name of the interventions in the order
  they appear in the argument `data` of
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  If `drug_names` is not defined, the order of the interventions as they
  appear in `data` is used, instead.

- title:

  A title to indicate the model (consistency model, network
  meta-regression or unrelated mean effects model).

## Value

A scatterplot of the leverage against the square root of the posterior
mean of residual deviance of the trial-arms under the model of interest.
The green, yellow, and red curves correspond to the parabola \\x^2 + y =
k\\ with \\k\\ = 1, 2, and 3, respectively. The data points correspond
to trial-arms. Data points found outside the yellow parabola are linked
with a pair of numbers. The first number refers to the position of the
trial in the dataset, and the second number refers to the corresponding
trial-arm (see 'Arguments' and 'Value' in
[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md)).
These trial-arms contribute more than 1 to the deviance information
criterion and, hence, the model's poor fit.

## Details

`leverage_plot` is integrated in the
[`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md)
function to create the leverage plot for the consistency model and the
unrelated mean effects model. These plots appear side-by-side in the
output of
[`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md).
Dias et al. (2010) used leverage plots to investigate the fit of the
consistency and inconsistency models–the latter through the
node-splitting approach.

## References

Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
treatment comparison meta-analysis. *Stat Med* 2010;**29**(7-8):932–44.
doi: 10.1002/sim.3767

## See also

[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md),
[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md),
[`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md)

## Author

Loukia M. Spineli
