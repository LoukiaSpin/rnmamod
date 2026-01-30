# The Bland-Altman plot

This function facilitates creating the Bland-Altman plot on the
posterior mean deviance contribution for two models using only three
arguments.

## Usage

``` r
bland_altman_plot(model1, model2, colour)
```

## Arguments

- model1:

  A vector with the numeric values of the target model (for instance,
  the consistency model).

- model2:

  A vector with the numeric values of the reference model (for instance,
  the unrelated mean effects model).

- colour:

  A string to define the colour of the data points in the plot.

## Value

Bland-Altman plot on the posterior mean deviance contribution of the
individual data points under model 1 and model 2. Each data point
corresponds to a trial-arm indicated by a pair of numbers. The first
number refers to the position of the trial in the dataset, and the
second arm refers to the corresponding trial-arm (see 'Arguments' and
'Value' in
[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md)).
The plot also displays the average bias and the 95% limits of agreement
with horizontal solid black lines.

## Details

`bland_altman_plot` is integrated in
[`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md)
to create the Bland-Altman plot on the posterior mean of deviance under
the consistency model (via
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md))
and the unrelated mean effects model (via
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md)).

A uniform scattering of the data points within the 95% limits of
agreement and average bias close to 0 indicate that the compared models
have a good agreement. Data points positioned above or below the 95%
limits of agreement correspond to trials that contribute to the poor fit
of the consistency model or unrelated mean effects model, respectively.

`bland_altman_plot` can be used to compare the following models
regarding deviance contribution:

- the consistency model (via
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md))
  with the unrelated effect means model (via
  [`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md));

- the network meta-analysis model (via
  [`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md))
  with the network meta-analysis model (via
  [`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md)).

## References

Bland JM, Altman DG. Measuring agreement in method comparison studies.
*Stat Methods Med Res* 1999;**8**:135–60. doi:
10.1177/096228029900800204

## See also

[`data_preparation`](https://loukiaspin.github.io/rnmamod/reference/data_preparation.md),
[`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md),
[`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md),
[`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md),
[`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md)

## Author

Loukia M. Spineli
