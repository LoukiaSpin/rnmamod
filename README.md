*** 

# rnmamod: package to perform Bayesian network meta-analysis methods

**rnmamod** is an R package to perform one-stage Bayesian fixed-effect or random-effects network meta-analysis while adjusting for *missing participant outcome data* using the pattern-mixture model. In the case of two inteventions, rnmamod performs one-stage Bayesian pairwise meta-analysis. The package handles a data-frame of binary or continuous outcome data in the arm-based format. The odds ratio, mean difference, standardised mean difference, and ratio of means are currently considered. The pattern-mixture model allows the incorporation of the informative missingness odds ratio for binary outcomes, whilst the informative missingness difference of means and the informative missingness ratio of means for continuous outcomes. The package comprises a suite of all necessary models for estimation and prediction of the intervention effect, and evaluation of the consistency assumption locally and globally. Missing participant outcome data are addressed in all models of the rnmamod package. The rnmamod package also includes a rich suite of visualisation tools that aid the interpretation and accommodation of the results in the submitted research work for publication. 

The rnmamod package is currently in development version.

## Installation

Run the following code to install rnmamod:

    devtools::install_github("LoukiaSpin/rnmamod")
    library(rnmamod)

## Example

We will use the dataset of [Baker et al. (2009)](https://pubmed.ncbi.nlm.nih.gov/19637942/) that includes 21 trials comparing seven pharmacologic interventions with each other and placebo in chronic obstructive pulmonary disease (COPD) patients. The exacerbation of COPD is the analysed binary outcome.

``` r
head(nma.baker2009)
#>                 study t1 t2 t3 t4  r1  r2 r3 r4 m1 m2 m3 m4  n1  n2 n3 n4
#> Llewellyn-Jones, 1996  3  8 NA NA   8   4 NA NA  0  1 NA NA   8   8 NA NA
#>        Paggiaro, 1998  3  8 NA NA  78  61 NA NA 19 27 NA NA 142 139 NA NA
#>          Mahler, 1999  6  8 NA NA  98  73 NA NA  9 23 NA NA 135 143 NA NA
#>        Casaburi, 2000  7  8 NA NA 222 132 NA NA 12 18 NA NA 279 191 NA NA
#>       van Noord, 2000  6  8 NA NA  29  24 NA NA  7  8 NA NA  47  50 NA NA
#>         Rennard, 2001  6  8 NA NA  72  65 NA NA 22 29 NA NA 132 135 NA NA
```

Create the network plot:

``` r
# The names of the interventions in the order they appear in the dataset
interv.names <- c("budesodine", "budesodine plus formoterol", "fluticasone", "fluticasone plus
                   salmeterol", "formoterol", "salmeterol", "tiotropium", "placebo")

netplot(data = nma.baker2009, drug.names = interv.names, text.cex = 1.5)
```

<div style="text-align: center"> <img src="figures/Network Baker.png" width="650" height="500" align="center"></div>

The following code performs a Bayesian random-effects network meta-analysis under the missing at random assumption and using intervention-specific informative missingness odds ratio (`IDE-ARM`) in the logarithmic scale:

``` r
res <- run.model(data = nma.baker2009,
                 measure = "OR",
                 model = "RE",
                 assumption = "IDE-ARM",
                 heter.prior = list("halfnormal", 0, 1),
                 mean.misspar = 0,
                 var.misspar = 1,
                 D = 1,
                 n.chains = 3,
                 n.iter = 10000,
                 n.burnin = 1000,
                 n.thin = 1)
```
<br/>
Illustrate all possible pairwise comparisons of the interventions using a league heatmap:

``` r
league.heatmap(full = res, drug.names = interv.names)
```

<div style="text-align: center"> <img src="figures/League Baker.png" width="750" height="600" align="center"></div>
<br/>
The following code presents the hierarchy of the interventions in the network using integrated rankograms and SUCRA (surfacw under the cumulative ranking) curves:

``` r
rankosucra.plot(full = res, drug.names = interv.names)
```

<div style="text-align: center"> <img src="figures/Sucra Baker.png" width="750" height="600" align="center"></div>
<br/>
<div style="text-align: right"> <img src="figures/dfg_logo_schriftzug_blau_foerderung_en.png" width="300" height="200" align="right"></div>
