# League table for relative and absolute effects

Provides a league table of the estimated odds ratio, and risk difference
per 1000 participants for all possible comparisons of interventions in
the network. The main diagonal of the table presents the absolute risk
for each intervention in the network. `league_table_absolute` can be
used for a random-effects or fixed-effect network meta-analysis. This
function should be used when the user has access to the raw trial-level
data (one-trial-per-row format with arm-level data).
`league_table_absolute` is applied for one binary outcome only.

## Usage

``` r
league_table_absolute(full, drug_names, show = NULL)
```

## Arguments

- full:

  An object of S3 class
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
  See 'Value' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

- drug_names:

  A vector of labels with the name of the interventions in the order
  they appear in the argument `data` of
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

- show:

  A vector of at least three character strings that refer to the names
  of the interventions *exactly* as defined in `drug_names`. Then, the
  league table will be created for these interventions only. If `show`
  is not defined, the league table will present all interventions as
  defined in `drug_names`.

## Value

A league table showing the posterior estimate and 95% credible interval
of the odds ratio (upper off-diagonals), risk difference per 1000
participants (lower off-diagonals), and absolute risks per 1000
participants (main diagonal).

## Details

The user must define the argument `measure = "RD"` in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md);
otherwise, the function will be stopped and an error message will be
printed in the R console.

The rows and columns of the league table display the names of the
interventions sorted by decreasing order from the best to the worst
based on their SUCRA value (Salanti et al., 2011) for the odds ratio.
The upper off-diagonals contain the posterior median and 95% credible
interval of the odds ratio, the lower off-diagonals contain the
posterior median and 95% credible interval of the risk difference (per
1000 participants), and the main diagonal comprises the posterior median
and 95% credible interval of the absolute risks (per 1000 participants)
of the corresponding interventions. The reference intervention of the
network (which the baseline risk has been selected for) is indicated in
the main diagonal with a black, thick frame.

Comparisons between interventions should be read from left to right.
Results that indicate strong evidence in favor of the row-defining
intervention (i.e. the respective 95% credible interval does not include
the null value) are indicated in bold.

To obtain unique absolute risks for each intervention, the network
meta-analysis model has been extended to incorporate the transitive
risks framework, namely, an intervention has the same absolute risk
regardless of the comparator intervention(s) in a trial (Spineli et al.,
2017). The absolute risks are a function of the odds ratio (the
**base-case** effect measure for a binary outcome) and the selected
baseline risk for the reference intervention (Appendix in Dias et al.,
2013). See 'Arguments' in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).
We advocate using the odds ratio as an effect measure for its desired
mathematical properties. Then, the risk difference can be obtained as a
function of the absolute risks of the corresponding interventions in the
comparison of interest.

`league_table_absolute` can be used only for a network of interventions.
In the case of two interventions, the execution of the function will be
stopped and an error message will be printed in the R console.

## References

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

[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)

## Author

Loukia M. Spineli
