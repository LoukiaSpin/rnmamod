# Detect the frail comparisons in multi-arm trials

Detects the frail comparisons in multi-arm trials, that is, comparisons
between non-baseline interventions not investigated in any two-arm trial
in the network (Spineli, 2021). The 'original' model of Dias et al.
(2013) omits the frail comparisons from the estimation process of the
unrelated mean effects model. Consequently, their posterior distribution
coincides with the prior distribution yielding implausible posterior
standard deviations.

## Usage

``` r
improved_ume(t, N, ns, na)
```

## Arguments

- t:

  A data-frame of the one-trial-per-row format containing the
  intervention identifier in each arm of every trial (see 'Details'
  below, and 'Format' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).

- N:

  A data-frame of the one-trial-per-row format containing the number of
  participants randomised to the assigned intervention in each arm of
  every trial (see 'Details' below, and 'Format' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).

- ns:

  A scale parameter on the number trials.

- na:

  A vector of length equal to `ns` with the number of arms in each
  trial.

## Value

The output of `improved_ume` is a list of elements that are inherited by
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md):

- nbase_multi:

  A scalar parameter on the number of frail comparisons.

- t1_bn:

  A vector with numeric values referring to the first arm of each frail
  comparison.

- t2_bn:

  A vector with numeric values referring to the second arm of each frail
  comparison.

- ref_base:

  A scalar referring to the reference intervention for the subnetwork of
  interventions in frail comparisons.

- base:

  A vector with numeric values referring to the baseline intervention of
  the multi-arm trials that contain the frail comparisons.

- obs_comp:

  A data-frame that indicates how many two-arm and multi-arm trials have
  included each pairwise comparison observed in the network.

## Details

`improved_ume` is integrated in
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md)
and calls the output of
[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md)
after sorting the rows so that multi-arm trials appear at the bottom of
the dataset. When there are no multi-arm trials or no frail comparisons
in the network, `improved_ume` returns only the element `obs_comp` (see,
'Value').

## References

Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence
synthesis for decision making 4: inconsistency in networks of evidence
based on randomized controlled trials. *Med Decis Making*
2013;**33**(5):641–56. doi: 10.1177/0272989X12455847

Spineli LM. A revised framework to evaluate the consistency assumption
globally in a network of interventions. *Med Decis Making* 2021. doi:
10.1177/0272989X211068005

## See also

[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md)

## Author

Loukia M. Spineli
