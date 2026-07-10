# Robustness index when 'metafor' or 'netmeta' are used

Calculates the robustness index for a sensitivity analysis (Spineli et
al., 2021) performed using the results of the analysis performed via the
R-package [netmeta](https://CRAN.R-project.org/package=netmeta) or
[metafor](https://CRAN.R-project.org/package=metafor). The user defines
the input and the function returns the robustness index.

## Usage

``` r
robustness_index_user(sens, pkg, attribute, threshold)
```

## Arguments

- sens:

  A list of R objects of class
  [`netmeta`](https://rdrr.io/pkg/netmeta/man/netmeta.html),
  [`netmetabin`](https://rdrr.io/pkg/netmeta/man/netmetabin.html) (see
  [netmeta](https://CRAN.R-project.org/package=netmeta)) or
  [`rma`](https://wviechtb.github.io/metafor/reference/rma.uni.html),
  [`rma.glmm`](https://wviechtb.github.io/metafor/reference/rma.glmm.html),
  [`rma.mh`](https://wviechtb.github.io/metafor/reference/rma.mh.html),
  [`rma.mv`](https://wviechtb.github.io/metafor/reference/rma.mv.html),
  [`rma.peto`](https://wviechtb.github.io/metafor/reference/rma.peto.html),
  and
  [`rma.uni`](https://wviechtb.github.io/metafor/reference/rma.uni.html)
  (see [metafor](https://CRAN.R-project.org/package=metafor)). The
  number of elements equals the number of analyses using the same
  dataset and the same R-package. The first element should refer to the
  primary analysis. Hence, the list should include at least two elements
  (see 'Details').

- pkg:

  Character string indicating the R-package with values `"netmeta"`, or
  `"metafor"`.

- attribute:

  This is relevant only for
  [netmeta](https://CRAN.R-project.org/package=netmeta). A vector of at
  least two characters with values `"TE.common"` or `"TE.random"`. See
  'Values' in [`netmeta`](https://rdrr.io/pkg/netmeta/man/netmeta.html)
  or [`netmetabin`](https://rdrr.io/pkg/netmeta/man/netmetabin.html).

- threshold:

  A number indicating the threshold of robustness, that is, the
  minimally allowed deviation between the primary analysis (the first
  element in `sens`) and re-analysis results. See 'Details' below.

## Value

`robustness_index_user` prints on the R console a message in red text on
the threshold of robustness determined by the user. Then, the function
returns the following list of elements:

- robust_index:

  A numeric scalar or vector on the robustness index values. In the case
  of a pairwise meta-analysis, `robust_index` is scalar as only one
  summary effect size is obtained. In the case of network meta-analysis,
  `robust_index` is a vector with length equal to the number of possible
  pairwise comparisons; one robustness index per pairwise comparison.

- robust:

  A character or character vector (of same length with `robust_index`)
  on whether the primary analysis results are *robust* or *frail* to the
  different re-analyses.

- kld:

  A vector or matrix on the Kullback-Leibler divergence (KLD) measure in
  the summary effect size from a subsequent re-analysis to the primary
  analysis. In the case of a pairwise meta-analysis, `kld` is a vector
  with length equal to the number of total analyses (one KLD value is
  obtained per analysis). The number of total analyses equals the length
  of `sens`. In the case of network meta-analysis, `robust_index` is a
  matrix with number of rows equal to the number of total analyses and
  number of columns equal to the number of possible pairwise
  comparisons; one KLD value per analysis and possible comparison.

- attribute:

  The attributes considered.

- threshold:

  The threshold used to be inherited by the
  [`heatmap_robustness`](https://loukiaspin.github.io/rnmamod/reference/heatmap_robustness.md)
  function. See 'Details'.

## Details

Thresholds of robustness have been proposed only for the odds ratio and
standardised mean difference (Spineli et al., 2021). The user may
consider the values 0.28 and 0.17 in the argument `threshold` for the
odds ratio and standardised mean difference effect measures (the default
values), respectively, or consider other plausible values. When the
argument `threshold` has not been defined, `robustness_index` considers
the default values 0.28 and 0.17 as threshold for robustness for binary
and continuous outcome, respectively, regardless of the effect measure
(the default thresholds may not be proper choices for other effect
measures; hence, use these threshold with great caution in this case).
Spineli et al. (2021) offers a discussion on specifying the `threshold`
of robustness.

When other effect measure is used (other than odds ratio or standardised
mean difference) or the elements in `sens` refer to different effect
measures, the execution of the function will be stopped and an error
message will be printed in the R console.

In `robust`, the value `"robust"` appears when the calculated
`robust_index` is less than `threshold`; otherwise, the value `"frail"`
appears.

## References

Kullback S, Leibler RA. On information and sufficiency. *Ann Math Stat*
1951;**22**(1):79–86. doi: 10.1214/aoms/1177729694

Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness
of primary analysis results: A case study on missing outcome data in
pairwise and network meta-analysis. *Res Synth Methods*
2021;**12**(4):475–90. doi: 10.1002/jrsm.1478

## See also

[`kld_measure`](https://loukiaspin.github.io/rnmamod/reference/kld_measure.md),
[`rma`](https://wviechtb.github.io/metafor/reference/rma.uni.html),
[`rma.glmm`](https://wviechtb.github.io/metafor/reference/rma.glmm.html),
[`rma.mh`](https://wviechtb.github.io/metafor/reference/rma.mh.html),
[`rma.mv`](https://wviechtb.github.io/metafor/reference/rma.mv.html),
[`rma.peto`](https://wviechtb.github.io/metafor/reference/rma.peto.html),
[`rma.uni`](https://wviechtb.github.io/metafor/reference/rma.uni.html),
[`netmeta`](https://rdrr.io/pkg/netmeta/man/netmeta.html),
[`netmetabin`](https://rdrr.io/pkg/netmeta/man/netmetabin.html),
[`heatmap_robustness`](https://loukiaspin.github.io/rnmamod/reference/heatmap_robustness.md)

## Author

Loukia M. Spineli

## Examples

``` r

if (FALSE) { # \dontrun{
library(netmeta)

data(Baker2009)

# Transform from arm-based to contrast-based format
p1 <- pairwise(treatment, exac, total, studlab = paste(study, year),
data = Baker2009, sm = "OR")

# Conduct standard network meta-analysis
net1 <- netmeta(p1, ref = "Placebo")

# Calculate the robustness index (random-effects versus fixed-effect)
robustness_index_user(sens = list(net1, net1),
                      pkg = "netmeta",
                      attribute = c("TE.random", "TE.common"),
                      threshold = 0.28)
} # }
```
