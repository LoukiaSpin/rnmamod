# Scatterplot of SUCRA values

Creates a scatterplot of the SUCRA values from the network meta-analysis
and the network meta-regression for a specified level or value of the
investigated covariate.

## Usage

``` r
scatterplot_sucra(full, reg, cov_name = "covariate value", drug_names)
```

## Arguments

- full:

  An object of S3 class
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  See 'Value' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

- reg:

  An object of S3 class
  [`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md).
  See 'Value' in
  [`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md).

- cov_name:

  A character or text to indicate the name of the covariate.

- drug_names:

  A vector of labels with the name of the interventions in the order
  they appear in the argument `data` of
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  If `drug_names` is not defined, the order of the interventions as they
  appear in `data` is used, instead.

## Value

A scatterplot of the SUCRA values under the network meta-analysis
(y-axis) against the SUCRA values under the network meta-regression
(x-axis) for a specified level or value of the investigated covariate.

## Details

The names of the interventions appear above each point in the plot.
Three coloured rectangles are drawn in the scatterplot: a red rectangle
for SUCRA values up to 50%, a yellow rectangular for SUCRA values
between 50% and 80%, and a green rectangle for SUCRA values over 80%.
Interventions falling at the green area are considered as the highest
ranked interventions, whilst interventions falling at the red area are
considered as the lowest ranked interventions.

`scatterplot_sucra` is integrated in
[`metareg_plot`](https://loukiaspin.github.io/rnmamod/reference/metareg_plot.md).

`scatterplot_sucra` can be used only for a network of interventions.
Otherwise, the execution of the function will be stopped and an error
message will be printed on the R console.

## References

Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical
summaries for presenting results from multiple-treatment meta-analysis:
an overview and tutorial. *J Clin Epidemiol* 2011;**64**(2):163–71. doi:
10.1016/j.jclinepi.2010.03.016

## See also

[`metareg_plot`](https://loukiaspin.github.io/rnmamod/reference/metareg_plot.md),
[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)

## Author

Loukia M. Spineli
