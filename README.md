*** 

<!---# rnmamod: package to perform Bayesian network meta-analysis methods--->
# rnmamod: A Bayesian network meta-analysis package for R 


  <!-- badges: start -->
  [![R-CMD-check](https://github.com/LoukiaSpin/rnmamod/workflows/R-CMD-check/badge.svg)](https://github.com/LoukiaSpin/rnmamod/actions)
  [!(https://img.shields.io/github/license/LoukiaSpin/rnmamod.svg)](https://github.com/LoukiaSpin/rnmamod/actions)
  <!-- badges: end -->
  

  
## Description 

The `rnmamod` package performs one-stage fixed-effect or random-effects Bayesian network meta-analysis (NMA) while adjusting for *missing participant outcome data* using the pattern-mixture model. In the case of two interventions, the package performs one-stage Bayesian pairwise meta-analysis. 

The package handles data inputs in *arm-based format*:
- binary data (effect size: odds ratio)  
- continuous data (effect size: raw and standardised mean differences, ratio of means).

<!--- The pattern-mixture model allows the incorporation of the informative missingness odds ratio for binary outcomes, whilst the informative missingness difference of means and the informative missingness ratio of means for continuous outcomes. --->

The `rnmamod` package comprises a suite of all necessary models for estimation and prediction of the intervention effect, and evaluation of the consistency assumption locally and globally. Also includes a rich suite of visualisation tools to aid in interpretation of the results and preparation of NMA manuscript submission. 
Missing participant outcome data are addressed in all models of the rnmamod package.

<!---that aid the interpretation and accommodation of the results in the submitted research work for publication. --->

The `rnmamod` package is currently in development version.

## Getting started

Run the following code to install the development version of the package:

    install.packages("devtools")
    devtools::install_github("LoukiaSpin/rnmamod")

## Example

<!--- We will use the dataset of --->
[Baker et al. (2009)](https://pubmed.ncbi.nlm.nih.gov/19637942/) comprising 21 trials comparing seven pharmacologic interventions with each other and placebo in chronic obstructive pulmonary disease (COPD) patients. The exacerbation of COPD (harmful outcome) is the analysed binary outcome.

``` r
head(nma.baker2009)
#>                 study t1 t2 t3 t4 r1 r2 r3 r4 m1 m2 m3 m4  n1  n2 n3 n4
#> Llewellyn-Jones, 1996  1  4 NA NA  3  0 NA NA  1  0 NA NA   8   8 NA NA
#>        Paggiaro, 1998  1  4 NA NA 51 45 NA NA 27 19 NA NA 139 142 NA NA
#>          Mahler, 1999  1  7 NA NA 47 28 NA NA 23  9 NA NA 143 135 NA NA
#>        Casaburi, 2000  1  8 NA NA 41 45 NA NA 18 12 NA NA 191 279 NA NA
#>       van Noord, 2000  1  7 NA NA 18 11 NA NA  8  7 NA NA  50  47 NA NA
#>         Rennard, 2001  1  7 NA NA 41 38 NA NA 29 22 NA NA 135 132 NA NA
```

Create the network plot using the *nma.networkplot* function of the [pcnetmeta](https://cran.r-project.org/web/packages/pcnetmeta/index.html) R package:

``` r
# The names of the interventions in the order they appear in the dataset
interv.names <- c("placebo", "budesonide", "budesonide plus formoterol", "fluticasone", "fluticasone plus salmeterol", "formoterol", "salmeterol", "tiotropium")

netplot(data = nma.baker2009, drug.names = interv.names, text.cex = 1.5)
```

<div style="text-align: center"> 
<img src="man/figures/Network_Baker.png" width="800" height="551" align="center">
</div>

The following code performs a Bayesian random-effects network meta-analysis under the *missing at random assumption* while using an intervention-specific informative missingness odds ratio (`assumption = "IDE-ARM"`) in the logarithmic scale:

``` r
res <- run.model(data = nma.baker2009,
                 measure = "OR",
                 model = "RE",
                 assumption = "IDE-ARM",
                 heter.prior = list("halfnormal", 0, 1),
                 mean.misspar = c(0, 0), 
                 var.misspar = 1,
                 D = 0,
                 n.chains = 3,
                 n.iter = 10000,
                 n.burnin = 1000,
                 n.thin = 1)
```

Illustrate all possible pairwise comparisons of the interventions using a league heatmap. Interventions are sorted in decreasing order by their posterior mean SUCRA (surface under the cumulative ranking) value in the main diagonal:

``` r
league.heatmap(full = res, drug.names = interv.names)
```

<div style="text-align: center"> 
<!--- <img src="man/figures/League_Baker.png" width="950" height="700" align="center"> --->
<!--- <img src="man/figures/League_Baker.png" width="70%" height="70%" align="center"> --->  
<img src="man/figures/League_Baker.png" width="70%" height="90%" align="center"> 
</div>
<br/>

The following code presents the hierarchy of the interventions in the network using integrated rankograms and SUCRA curves:

``` r
rankosucra.plot(full = res, drug.names = interv.names)
```

<div style="text-align: center"> 
<!--- <img src="man/figures/Sucra_Baker.png" width="800" height="600" align="center"> --->
<!--- <img src="man/figures/Sucra_Baker.png" width="70%" height="70%" align="center"> --->
<img src="man/figures/Sucra_Baker.png" width="70%" height="90%" align="center">
</div>

##
The development of the **rnmamod** package is funded by the **German Research Foundation** (Deutsche Forschungsgemeinschaft) (grant no. [SP 1664/1-3](https://gepris.dfg.de/gepris/projekt/339420617)) <div style="text-align: right"> <img src="man/figures/dfg_logo_schriftzug_blau_foerderung_en.png" width="340" height="150" align="right"></div>
