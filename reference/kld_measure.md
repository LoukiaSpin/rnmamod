# Function for the Kullback-Leibler Divergence of two normally distributed treatment effects for the same pairwise comparison

The user specify the (posterior) mean and standard error (or posterior
standard deviation) of two estimated treatment effects, X and Y, that
refer to the **same** pairwise comparison and are assumed to follow a
normal distribution. The function returns the Kullback-Leibler
Divergence (KLD) measure of 1) approximating X with Y, 2) approximating
Y with X, and 3) their average.

## Usage

``` r
kld_measure(mean_y, sd_y, mean_x, sd_x)
```

## Arguments

- mean_y:

  A real number that refers to the mean of the estimated treatment
  effect Y on the scale of the selected effect measure (in logarithmic
  scale for relative effect measures).

- sd_y:

  A positive integer that refers to the posterior standard deviation or
  the standard error of the estimated treatment effect Y on the scale of
  the selected effect measure (in logarithmic scale for relative effect
  measures).

- mean_x:

  A real number that refers to the mean of the estimated treatment
  effect X on the scale of the selected effect measure (in logarithmic
  scale for relative effect measures).

- sd_x:

  A positive integer that refers to the posterior standard deviation or
  the standard error of the estimated treatment effect X on the scale of
  the selected effect measure (in logarithmic scale for relative effect
  measures).

## Value

The function return the following numeric results:

|                |                                                             |
|----------------|-------------------------------------------------------------|
| **kld_sym**    | The symmetric KLD value as the average of two KLD values .  |
|                |                                                             |
| **kld_x_true** | The KLD value when approximating X by Y (X is the 'truth'). |
|                |                                                             |
| **kld_y_true** | The KLD value when approximating Y by X (Y is the 'truth'). |

## References

Kullback S, Leibler RA. On information and sufficiency. *Ann Math Stat*
1951;**22**(1):79–86. doi: 10.1214/aoms/1177729694

## See also

[`kld_inconsistency`](https://loukiaspin.github.io/rnmamod/reference/kld_inconsistency.md),
[`kld_inconsistency_user`](https://loukiaspin.github.io/rnmamod/reference/kld_inconsistency_user.md),
[`robustness_index`](https://loukiaspin.github.io/rnmamod/reference/robustness_index.md),
[`robustness_index_user`](https://loukiaspin.github.io/rnmamod/reference/robustness_index_user.md)
