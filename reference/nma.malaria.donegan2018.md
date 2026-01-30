# Artemether, artesunate and quinine for severe malaria

A dataset of 24 two-arm trials comparing artemether, artesunate and
quinine with each other in patients with severe malaria. The analysed
binary outcome is death.

## Usage

``` r
data(nma.malaria.donegan2018)
```

## Format

A data frame with 24 rows of arm-based data and 11 columns referring to
the trial number, the treatment identifier in the compared arms, the
odds ratio in the logarithmic scale and its standard error for each
trial, the average age in years and its centered version, the number of
events and number randomised for each arm and trial.

## Source

Donegan S, Dias S, Tudur-Smith C, Marinho V, Welton NJ. Graphs of study
contributions and covariate distributions for network meta-regression.
*Res Synth Methods* 2018;**9**(2):243–60. doi: 10.1002/jrsm.1292

## Details

The interventions have been coded as follows: 1, quinine; 2, artemether;
and 3, artesunate
