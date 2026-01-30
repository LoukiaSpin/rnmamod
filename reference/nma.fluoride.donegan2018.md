# Topical fluoride interventions for preventing dental caries

A dataset of 130 trials comparing different forms of topical fluoride
interventions for preventing dental caries in children or adolescents
with at least 1 year or school year. The analysed continuous outcome is
the change from baseline in decayed, missing, and filled tooth surfaces.
The dataset contains also multi-arm trials.

## Usage

``` r
data(nma.fluoride.donegan2018)
```

## Format

A data frame with 140 rows of arm-based data and 16 columns referring to
the trial number, the treatment identifier in the compared arms, the
standardised mean difference and its standard error for each trial and
possible comparison (in the case of multi-arm trial), the randomisation
year, the standard deviation and number randomised for each arm and
trial, the pooled standard deviation and the within-study covariance in
multi-arm trials.

## Source

Donegan S, Dias S, Tudur-Smith C, Marinho V, Welton NJ. Graphs of study
contributions and covariate distributions for network meta-regression.
*Res Synth Methods* 2018;**9**(2):243–60. doi: 10.1002/jrsm.1292

## Details

The interventions have been coded as follows: 1, no treatment; 2,
placebo; 3, dentifrice; 4, rinse; 5, gel; and 6, varnish
