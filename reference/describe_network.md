# A function to describe the evidence base

Calculates the necessary elements to describe the evidence base for an
outcome across the network, the interventions, and observed comparisons.

## Usage

``` r
describe_network(data, drug_names, measure, save_xls)
```

## Arguments

- data:

  A data-frame of a one-trial-per-row format containing arm-level data
  of each trial. See 'Format' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

- drug_names:

  A vector of labels with the name of the interventions in the order
  they appear in the argument `data`.

- measure:

  Character string indicating the effect measure. For a binary outcome,
  the following can be considered: `"OR"`, `"RR"` or `"RD"` for the odds
  ratio, relative risk, and risk difference, respectively. For a
  continuous outcome, the following can be considered: `"MD"`, `"SMD"`,
  or `"ROM"` for mean difference, standardised mean difference and ratio
  of means, respectively.

- save_xls:

  Logical to indicate whether to export the tabulated results to an
  'xlsx' file (via the
  [`write_xlsx`](https://docs.ropensci.org/writexl//reference/write_xlsx.html)
  function of the R-package
  [writexl](https://CRAN.R-project.org/package=writexl)) at the working
  directory of the user. The default is `FALSE` (do not export).

## Value

`describe_network` returns the following data-frames that describe the
evidence base:

- network_description:

  The number of: interventions, possible comparisons, direct and
  indirect comparisons, number of trials in total, number of two-arm and
  multi-arm trials, number of randomised participants, and proportion of
  participants completing the trial (completers). When the outcome is
  binary, the number of trials with at least one zero event, and the
  number of trials with all zero events are also presented.

- table_interventions:

  For each intervention, the number of trials, number of randomised
  participants, and proportion of completers. When the outcome is
  binary, the data-frame presents also the corresponding proportion of
  total observed events, the minimum, median and maximum proportion of
  observed events across the corresponding trials.

- table_comparisons:

  Identical structure to `table_interventions` but for each observed
  comparison in the network.

## Details

`describe_network` calls
[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md)
to facilitate the calculations.

Furthermore, `describe_network` exports the data-frames to separate
'xlsx' files (via the
[`write_xlsx`](https://docs.ropensci.org/writexl//reference/write_xlsx.html)
function of the R-package
[writexl](https://CRAN.R-project.org/package=writexl)) at the working
directory of the user.

## See also

[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
[`write_xlsx`](https://docs.ropensci.org/writexl//reference/write_xlsx.html)

## Author

Loukia M. Spineli
