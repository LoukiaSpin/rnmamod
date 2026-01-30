# Prepare the dataset in the proper format for R2jags

`data_preparation` prepares the dataset in the proper format for R2jags
and returns a list of elements that
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)
inherits via the argument `data`.

## Usage

``` r
data_preparation(data, measure)
```

## Arguments

- data:

  A data-frame of the one-trial-per-row format with arm-level data. See
  'Format' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md).

- measure:

  Character string indicating the effect measure. For a binary outcome,
  the following can be considered: `"OR"`, `"RR"` or `"RD"` for the odds
  ratio, relative risk, and risk difference, respectively. For a
  continuous outcome, the following can be considered: `"MD"`, `"SMD"`,
  or `"ROM"` for mean difference, standardised mean difference and ratio
  of means, respectively.

## Value

A list of data-frames on the following elements to be passed to
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md):

- pseudo_m:

  A pseudo-data-frame with values -1 and **m** for the corresponding
  trial-arms with unreported and reported missing participant outcome
  data, respectively (see 'Details').

- m:

  The number of missing participant outcome data in each trial-arm (see
  'Details').

- N:

  The number of randomised participants in each trial-arm.

- t:

  The intervention identifier in each trial-arm.

- I:

  A pseudo-data-frame that indicates whether missing participant outcome
  data have been reported or not for each observed trial-arm (see
  'Details').

- measure:

  The effect measure for the analysed outcome.

- y0:

  The observed mean value of the outcome in each trial-arm, when the
  outcome is continuous.

- se0:

  The observed standard deviation of the outcome in each trial-arm, when
  the outcome is continuous.

- r:

  The number of observed events of the outcome in each trial-arm, when
  the outcome is binary.

## Details

`data_preparation` prepares the data for the Bayesian analysis (See
'Format' in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).
`data_preparation` creates the pseudo-data-frames `m_new`, `I`, and
`m_pseudo` that have the same dimensions with the element `N`. `m_new`
takes the zero value for the observed trial-arms with unreported missing
participant outcome data (i.e., **m** equals `NA` for the corresponding
trial-arms), the same value with **m** for the observed trial-arms with
reported missing participant outcome data, and `NA` for the unobserved
trial-arms. `I` is a dummy data-frame and takes the value one for the
observed trial-arms with reported missing participant outcome data, the
zero value for the observed trial-arms with unreported missing
participant outcome data (i.e., `m_new` equals zero for the
corresponding trial-arms), and `NA` for the unobserved trial-arms. Thus,
`I` indicates whether missing participant outcome data have been
collected for the observed trial-arms. If the user has not defined the
element **m** in `data_preparation`, `m_new` and `I` take the zero value
for all observed trial-arms to indicate that no missing participant
outcome data have been collected for the analysed outcome. `I` and
`m_new` are used from the following functions of the package:
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md),
[`prepare_model`](https://loukiaspin.github.io/rnmamod/reference/prepare_model.md),
[`run_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/run_nodesplit.md),
[`prepare_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/prepare_nodesplit.md),
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md),
[`prepare_ume`](https://loukiaspin.github.io/rnmamod/reference/prepare_ume.md),
and
[`run_sensitivity`](https://loukiaspin.github.io/rnmamod/reference/run_sensitivity.md).
Lastly, `m_pseudo` is a variant of `m_new`: it takes the value -1 for
the observed trial-arms with unreported missing participant outcome data
(i.e., **m** equals `NA` for the corresponding trial-arms), the same
value with **m** for the observed trial-arms with reported missing
participant outcome data, and `NA` for the unobserved trial-arms. It is
used in function
[`heatmap_missing_network`](https://loukiaspin.github.io/rnmamod/reference/heatmap_missing_network.md)
to calculate and illustrate the percentage of missing participant
outcome data across the observed comparisons and interventions of the
network and the function
[`heatmap_missing_dataset`](https://loukiaspin.github.io/rnmamod/reference/heatmap_missing_dataset.md)
to illustrate the trial-arms with unreported missing participant outcome
data. All pseudo-data-frames aim to retain the trials without
information on missing participant outcome data.

Furthermore, `data_preparation` sorts the interventions across the arms
of each trial in an ascending order and correspondingly the remaining
elements in `data` (See 'Format' in
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).
`data_preparation` considers the first column in **t** as being the
control arm for every trial. Thus, this sorting ensures that
interventions with a lower identifier are consistently treated as the
control arm in each trial. This case is relevant in non-star-shaped
networks.

## See also

[`heatmap_missing_dataset`](https://loukiaspin.github.io/rnmamod/reference/heatmap_missing_dataset.md),
[`heatmap_missing_network`](https://loukiaspin.github.io/rnmamod/reference/heatmap_missing_network.md),
[R2jags](https://CRAN.R-project.org/package=R2jags),
[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`run_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/run_nodesplit.md),
[`run_sensitivity`](https://loukiaspin.github.io/rnmamod/reference/run_sensitivity.md),
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md),
[`prepare_model`](https://loukiaspin.github.io/rnmamod/reference/prepare_model.md),
[`prepare_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/prepare_nodesplit.md),
[`prepare_ume`](https://loukiaspin.github.io/rnmamod/reference/prepare_ume.md)

## Author

Loukia M. Spineli
