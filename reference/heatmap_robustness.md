# Heatmap of robustness

Facilitates the detection of comparisons that are associated with a lack
of robustness in the context of a sensitivity analysis.

## Usage

``` r
heatmap_robustness(robust, drug_names)
```

## Arguments

- robust:

  An object of S3 class
  [`robustness_index`](https://loukiaspin.github.io/rnmamod/reference/robustness_index.md)
  and
  [`robustness_index_user`](https://loukiaspin.github.io/rnmamod/reference/robustness_index_user.md).
  See 'Value' in
  [`robustness_index`](https://loukiaspin.github.io/rnmamod/reference/robustness_index.md)
  and
  [`robustness_index_user`](https://loukiaspin.github.io/rnmamod/reference/robustness_index_user.md).

- drug_names:

  A vector of labels with the name of the interventions in the order
  they appear in the argument `data` of
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  If `drug_names` is not defined, the order of the interventions as they
  appear in `data` is used, instead.

## Value

`heatmap_robustness` first prints on the R console a message on the
threshold of robustness determined by the user in
[`robustness_index`](https://loukiaspin.github.io/rnmamod/reference/robustness_index.md)
and
[`robustness_index_user`](https://loukiaspin.github.io/rnmamod/reference/robustness_index_user.md).
Then, it returns a lower triangular heatmap matrix with the robustness
index value of all possible pairwise comparisons.

## Details

The heatmap illustrates the robustness index for each possible pairwise
comparison in the network. The pairwise comparisons are read from left
to right. Comparisons highlighted with green or red colour imply robust
or frail conclusions for the primary analysis, respectively. This
corresponds to robustness index below or at least the selected threshold
of robustness. `heatmap_robustness` inherits the threshold of robustness
selected in the
[`robustness_index`](https://loukiaspin.github.io/rnmamod/reference/robustness_index.md)
or
[`robustness_index_user`](https://loukiaspin.github.io/rnmamod/reference/robustness_index_user.md)
function. The robustness index of each pairwise comparison also appears
in the corresponding cell. When there is at least one comparison with
frail conclusions, the primary analysis results may be questionable for
the whole network (Spineli et al., 2021).

`heatmap_robustness` is *not* restricted to the sensitivity analysis
concerning the impact of missing participant outcome data.

`heatmap_robustness` can be used only for a network of interventions.
Otherwise, the execution of the function will be stopped and an error
message will be printed on the R console.

## References

Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness
of primary analysis results: A case study on missing outcome data in
pairwise and network meta-analysis. *Res Synth Methods*
2021;**12**(4):475–90. doi: 10.1002/jrsm.1478

## See also

[`robustness_index`](https://loukiaspin.github.io/rnmamod/reference/robustness_index.md),
[`robustness_index_user`](https://loukiaspin.github.io/rnmamod/reference/robustness_index_user.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)

## Author

Loukia M. Spineli

## Examples

``` r
data("nma.baker2009")

# Read results from 'run_sensitivity' (using the default arguments)
res_sens <- readRDS(system.file('extdata/res_sens_baker.rds',
                    package = 'rnmamod'))

# Calculate the robustness index
robust <- robustness_index(sens = res_sens,
                           threshold = 0.28)
#> The value 0.28 was assigned as 'threshold' for odds ratio.

# The names of the interventions in the order they appear in the dataset
interv_names <- c("placebo", "budesonide", "budesonide plus formoterol",
                  "fluticasone", "fluticasone plus salmeterol",
                  "formoterol", "salmeterol", "tiotropium")

# Create the heatmap of robustness
heatmap_robustness(robust = robust,
                   drug_names = interv_names)
#> The value 0.28 was assigned as 'threshold' for odds ratio.

```
