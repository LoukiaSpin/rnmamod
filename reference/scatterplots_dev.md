# Deviance scatterplots

Illustrates the posterior mean of deviance contribution of the
individual data points under the unrelated mean effects model (via
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md))
against the posterior mean of deviance contribution under the
consistency model (via
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).

## Usage

``` r
scatterplots_dev(full, ume, colour)
```

## Arguments

- full:

  A numeric vector with the posterior mean of deviance obtained using
  the consistency model (see 'Value' in
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md)).

- ume:

  A numeric vector with the posterior mean of deviance obtained using
  the unrelated mean effects model (see 'Value' in
  [`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md)).

- colour:

  A string to define the colour of the points in the plot.

## Value

A scatterplot of the posterior mean deviance contribution of the
individual data points from the unrelated mean effects model against
those from the consistency model. Each data point corresponds to a
trial-arm indicated by a pair of numbers. The first number refers to the
position of the trial in the dataset, and the second number refers to
the corresponding trial-arm (see 'Arguments' and 'Value' in
[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md)).

## Details

`scatterplots_dev` is integrated in the
[`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md)
function to compare the models regarding the posterior mean of deviance.
This scatterplot has also been considered by Dias et al. (2013). When
the majority of data points are scattered across the diagonal line, we
may conclude that the compared models have a good agreement. Data points
systematically scattered above or below the diagonal line may contribute
more to the poor fit of the unrelated mean effects model and the
consistency model, respectively.

## References

Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence
synthesis for decision making 4: inconsistency in networks of evidence
based on randomized controlled trials. *Med Decis Making*
2013a;**33**(5):641–56. doi: 10.1177/0272989X12455847

## See also

[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md),
[`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md)

## Author

Loukia M. Spineli
