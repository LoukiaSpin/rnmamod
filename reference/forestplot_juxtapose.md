# Forest plot of juxtaposing several network meta-analysis models

Provides a forest plot with the posterior median and 95% credible and
prediction intervals for comparisons with the selected intervention
(comparator) in the network under several network meta-analyses models,
as well as a forest plot with the corresponding SUCRA values.

## Usage

``` r
forestplot_juxtapose(
  results,
  compar,
  name,
  drug_names,
  axis_title_size = 12,
  axis_text_size = 12,
  caption_text_size = 9,
  label_size = 3.5,
  position_width = 0.8
)
```

## Arguments

- results:

  A list of at least two objects of S3 class
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
  or
  [`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md).
  See 'Value' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
  and
  [`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md).

- compar:

  A character to indicate the comparator intervention. It must be any
  name found in `drug_names`.

- name:

  A vector of characters referring to the juxtaposed models. If the
  argument is left unspecified, the names of models appear as 'Model X'
  with 'X' being the order/position of each model in the argument
  `results`.

- drug_names:

  A vector of labels with the name of the interventions in the order
  they appear in the argument `data` of
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  If `drug_names` is not defined, the order of the interventions as they
  appear in `data` is used, instead.

- axis_title_size:

  A positive integer for the font size of x axis title.
  `axis_title_size` determines the axis.title argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- axis_text_size:

  A positive integer for the font size of axis text (both axes).
  `axis_text_size` determines the axis.text argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- caption_text_size:

  A positive integer for the font size of caption text.
  `caption_text_size` determines the plot.caption argument found in the
  theme's properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- label_size:

  A positive integer for the font size of labels appearing on each
  interval. `label_size` determines the size argument found in the
  geom's aesthetic properties in the R-package
  [ggplot2](https://CRAN.R-project.org/package=ggplot2).

- position_width:

  A positive integer specifying the vertical position of the intervals.
  `position_width` is found in the geom's aesthetic properties in the
  R-package [ggplot2](https://CRAN.R-project.org/package=ggplot2).

## Value

A list of the following two figures:

- forest_plots:

  A panel of two forest plots: (1) a forest plot on the posterior median
  and 95% credible and prediction intervals for comparisons with the
  selected comparator treatment (specified with `compar`), and (2) a
  forest plot on the posterior mean and 95% credible interval of SUCRA
  values of the treatments (Salanti et al., 2011).

- tau_plot:

  A forest plot on the posterior median and 95% credible interval of the
  between-study standard deviation.

## Details

The y-axis of the forest plot on **forest_plots** displays the labels of
the treatments in the network; the selected treatment that comprises the
`compar` argument is annotated in the plot with the label 'Comparator
intervention'. For each comparison with the selected treatment, the 95%
credible and prediction intervals are displayed as overlapping lines.
Black lines refer to estimation under both analyses. Coloured lines
refer to prediction under each model, respectively. The corresponding
numerical results are displayed above each line: 95% credible intervals
are found in parentheses, and 95% predictive intervals are found in
brackets. Odds ratios, relative risks, and ratio of means are reported
in the original scale after exponentiation of the logarithmic scale.

If one of the models refer to network meta-regression
([`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md))
the results on treatment effects (estimation and prediction) and SUCRA
values refer to the covariate value selected when employing
[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md).

The y-axis for the forest plot on **SUCRA** values displays the labels
of the treatments in the network. The corresponding numerical results
are displayed above each line.

In **forest_plots** and **tau_plot**, the treatments are sorted in the
descending order of their SUCRA values based on the first model
specified in `results`.

**Important note:** `forestplot_juxtapose` should be used to compare the
results from several network meta-analysis models that contain the same
treatments, have the same meta-analysis model (fixed-effect or
random-effects) and the same effect measure; otherwise the execution of
the function will be stopped and an error message will be printed on the
R console.

`forestplot_juxtapose` is used only for a network of treatments. In the
case of two treatments, the execution of the function will be stopped
and an error message will be printed on the R console.

## References

Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical
summaries for presenting results from multiple-treatment meta-analysis:
an overview and tutorial. *J Clin Epidemiol* 2011;**64**(2):163–71.
[doi:10.1016/j.jclinepi.2010.03.016](https://doi.org/10.1016/j.jclinepi.2010.03.016)

## See also

[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)

## Author

Loukia M. Spineli
