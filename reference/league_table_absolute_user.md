# League table for relative and absolute effects (user defined)

In line with
[`league_table_absolute`](https://loukiaspin.github.io/rnmamod/reference/league_table_absolute.md),
provides a league table of the estimated odds ratio, and risk difference
per 1000 participants for all possible comparisons of interventions in
the network. The main diagonal of the table presents the absolute risk
for each intervention in the network. `league_table_absolute_user`
requires users to input the summary effect and 95% credible or
confidence interval of the basic parameters in the reported effect
measure. This function should be used when the user has access to the
results of a published systematic review rather than the raw trial-level
data. In the latter case, the user should consider the function
[`league_table_absolute`](https://loukiaspin.github.io/rnmamod/reference/league_table_absolute.md).
`league_table_absolute_user` is applied for one binary outcome only.

## Usage

``` r
league_table_absolute_user(
  data,
  measure,
  base_risk,
  drug_names,
  show = NULL,
  save_xls
)
```

## Arguments

- data:

  A data-frame with the summary effects of comparisons with the
  reference intervention of the network, known as basic parameters. The
  data-frame has `T` rows (`T` is the number of interventions in the
  network) and four columns that contain the point estimate, the lower
  and upper bound of the 95% (confidence or credible) interval of the
  corresponding basic parameters, and a ranking measure to indicate the
  order of the interventions in the hierarchy from the best to the worst
  with possible choices a non-zero positive integer for the rank, the
  SUCRA value (Salanti et al., 2011) or p-score value (Ruecker and
  Schwarzer, 2015). The first row of the data-frame refers to the
  selected reference intervention and should include (1) the null value
  three times at the investigated effect measure (i.e. 1 for odds ratio
  and relative risk, and 0 for risk difference), and (2) the value of
  the ranking measure.

- measure:

  Character string indicating the effect measure of `data`. For a binary
  outcome, the following can be considered: `"OR"`, `"RR"` or `"RD"` for
  the odds ratio, relative risk, and risk difference, respectively.

- base_risk:

  A number in the interval (0, 1) that indicates the baseline risk for
  the selected reference intervention.

- drug_names:

  A vector of labels with the name of the interventions in the order
  they appear in the argument `data`. The first intervention should be
  the selected reference intervention.

- show:

  A vector of at least three character strings that refer to the names
  of the interventions *exactly* as defined in `drug_names`. Then, the
  league table will be created for these interventions only. If `show`
  is not defined, the league table will present all interventions as
  defined in `drug_names`.

- save_xls:

  Logical to indicate whether to export the tabulated results to an
  'xlsx' file (via the
  [`write_xlsx`](https://docs.ropensci.org/writexl//reference/write_xlsx.html)
  function of the R-package
  [writexl](https://CRAN.R-project.org/package=writexl)) to the working
  directory of the user. The default is `FALSE` (do not export).

## Value

A league table showing the estimate and 95% confidence interval of the
odds ratio (upper off-diagonals), risk difference per 1000 participants
(lower off-diagonals), and absolute risks per 1000 participants (main
diagonal).

## Details

When the published results are reported in the relative risk scale
(i.e., `measure = "RR"`), the function calculates odds ratios and risk
differences (point estimate and 95% confidence interval) for all
possible pairwise comparisons in the network based on the obtained
absolute risks and the selected baseline risk. Likewise, when the
published results are in the odds ratio or risk difference scale (i.e.,
`measure = "OR"` or `measure = "RD"`, respectively), the function
calculates risk differences or odds ratios (point estimate and 95%
confidence interval), respectively, for all possible pairwise
comparisons in the network based on the obtained absolute risks and the
selected baseline risk.

The rows and columns of the league table display the names of the
interventions sorted by decreasing order from the best to the worst
based on the ranking measure in the fourth column of the argument
`data`. The upper off-diagonals contain the estimate and 95% confidence
interval of the odds ratio, the lower off-diagonals contain the estimate
and 95% confidence interval of the risk difference (per 1000
participants), and the main diagonal comprises the absolute risks and
their 95% confidence interval (per 1000 participants) of the
corresponding non-reference interventions. The reference intervention of
the network (which the baseline risk has been selected for) is indicated
in the main diagonal with a black, thick frame.

Comparisons between interventions should be read from left to right.
Results that indicate strong evidence in favour of the row-defining
intervention (i.e. the respective 95% confidence interval does not
include the null value) are indicated in bold.

Furthermore, `league_table_absolute_user` exports
`table_relative_absolute_effect`, a table with the relative and absolute
effects of the basic parameters, as an 'xlsx' file (via the
[`write_xlsx`](https://docs.ropensci.org/writexl//reference/write_xlsx.html)
function) to the working directory of the user.

To obtain unique absolute risks for each intervention, we have
considered the transitive risks framework, namely, an intervention has
the same absolute risk regardless of the comparator intervention(s) in a
trial (Spineli et al., 2017).

`league_table_absolute_user` can be used only for a network of
interventions. In the case of two interventions, the execution of the
function will be stopped and an error message will be printed in the R
console.

## References

Ruecker G, Schwarzer G. Ranking treatments in frequentist network
meta-analysis works without resampling methods. *BMC Med Res Methodol*
2015;**15**:58. doi: 10.1186/s12874-015-0060-8

Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical
summaries for presenting results from multiple-treatment meta-analysis:
an overview and tutorial. *J Clin Epidemiol* 2011;**64**(2):163–71. doi:
10.1016/j.jclinepi.2010.03.016

Spineli LM, Brignardello-Petersen R, Heen AF, Achille F, Brandt L,
Guyatt GH, et al. Obtaining absolute effect estimates to facilitate
shared decision making in the context of multiple-treatment comparisons.
Abstracts of the Global Evidence Summit, Cape Town, South Africa.
*Cochrane Database of Systematic Reviews* 2017;**9**(Suppl 1):1891.

## See also

[`league_table_absolute`](https://loukiaspin.github.io/rnmamod/reference/league_table_absolute.md),
[`write_xlsx`](https://docs.ropensci.org/writexl//reference/write_xlsx.html)

## Author

Loukia M. Spineli
