# A panel of interval plots for the unrelated mean effects model

Creates a panel of interval plots on the summary effect sizes under the
consistency model and the unrelated mean effects model. The number of
interval plots equals the number of pairwise comparisons observed in the
network.

## Usage

``` r
intervalplot_panel_ume(full, ume, drug_names)
```

## Arguments

- full:

  An object of S3 class
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  See 'Value' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

- ume:

  An object of S3 class
  [`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md).
  See 'Value' in
  [`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md).

- drug_names:

  A vector of labels with the name of the interventions in the order
  they appear in the argument `data` of
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  If the argument `drug_names` is not defined, the order of the
  interventions as they appear in `data` is used, instead.

## Value

A panel of interval plots on the posterior mean and 95% credible
interval of the summary effect size under the consistency model and the
improved unrelated mean effects model (Spineli, 2021) of all pairwise
comparisons observed in the network.

## Details

`intervalplot_panel_ume` is integrated in the
[`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md)
function. The consistency model and the unrelated mean effects model are
abbreviated in the y-axis as 'NMA model' and 'UME model', respectively.
The intervals are highlighted with green, when the corresponding summary
effect sizes do not cross the vertical line of no difference, and red
otherwise. Grey panels refer to the frail comparisons as detected by the
[`improved_ume`](https://loukiaspin.github.io/rnmamod/reference/improved_ume.md)
function (see 'Details' in
[`improved_ume`](https://loukiaspin.github.io/rnmamod/reference/improved_ume.md)).

For a binary outcome, when `measure` is "RR" (relative risk) or "RD"
(risk difference) in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
`intervalplot_panel_ume` currently presents the results in the odds
ratio scale.

## References

Spineli LM. A revised framework to evaluate the consistency assumption
globally in a network of interventions. *Med Decis Making* 2021. doi:
10.1177/0272989X211068005

## See also

[`improved_ume`](https://loukiaspin.github.io/rnmamod/reference/improved_ume.md)
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md),
[`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md)

## Author

Loukia M. Spineli
