# Comparator-specific forest plot for network meta-regression

Provides a forest plot with the posterior median and 95% credible and
prediction intervals for comparisons with the selected intervention
(comparator) in the network under the network meta-analysis *and*
network meta-regression for a specified level or value of the
investigated covariate, and a forest plot with the corresponding SUCRA
values.

## Usage

``` r
forestplot_metareg(full, reg, compar, cov_name = "covariate value", drug_names)
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

- compar:

  A character to indicate the comparator intervention. It must be any
  name found in `drug_names`.

- cov_name:

  A character or text to indicate the name of the covariate.

- drug_names:

  A vector of labels with the name of the interventions in the order
  they appear in the argument `data` of
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  If `drug_names` is not defined, the order of the interventions as they
  appear in `data` is used, instead.

## Value

A panel of two forest plots: (1) a forest plot on the effect estimates
and predictions of comparisons with the selected intervention in the
network under the network meta-analysis and network meta-regression for
a specified level or value of the investigated covariate, and (2) a
forest plot on the posterior mean and 95% credible interval of SUCRA
values of the interventions (Salanti et al., 2011).

## Details

The y-axis of the forest plot on **effect sizes** displays the labels of
the interventions in the network; the selected intervention that
comprises the `compar` argument is annotated in the plot with the label
'Comparator intervention'. For each comparison with the selected
intervention, the 95% credible and prediction intervals are displayed as
overlapping lines. Black lines refer to estimation under both analyses.
Green and red lines refer to prediction under network meta-analysis and
network meta-regression, respectively. The corresponding numerical
results are displayed above each line: 95% credible intervals are found
in parentheses, and 95% predictive intervals are found in brackets. Odds
ratios, relative risks, and ratio of means are reported in the original
scale after exponentiation of the logarithmic scale.

The y-axis for the forest plot on **SUCRA** values displays the labels
of the interventions in the network. The corresponding numerical results
are displayed above each line. Three coloured rectangles appear in the
forest plot: a red rectangle for SUCRA values up to 50%, a yellow
rectangular for SUCRA values between 50% and 80%, and a green rectangle
for SUCRA values over 80%. Interventions falling at the green area are
considered as the highest ranked interventions, whilst interventions
falling at the red area are considered as the lowest ranked
interventions.

In both plots, the interventions are sorted in the descending order of
their SUCRA values based on the network meta-analysis.

`forestplot_metareg` is integrated in
[`metareg_plot`](https://loukiaspin.github.io/rnmamod/reference/metareg_plot.md).

`forestplot_metareg` can be used only for a network of interventions. In
the case of two interventions, the execution of the function will be
stopped and an error message will be printed on the R console.

## References

Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical
summaries for presenting results from multiple-treatment meta-analysis:
an overview and tutorial. *J Clin Epidemiol* 2011;**64**(2):163–71.
[doi:10.1016/j.jclinepi.2010.03.016](https://doi.org/10.1016/j.jclinepi.2010.03.016)

## See also

[`metareg_plot`](https://loukiaspin.github.io/rnmamod/reference/metareg_plot.md),
[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)

## Author

Loukia M. Spineli
