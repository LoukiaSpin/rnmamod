*** 

# rnmamod: package to perform Bayesian network meta-analysis methods

rnmamod is an R package to perform one-stage Bayesian fixed-effect or random-effects network meta-analysis while adjusting for missing participant outcome data using the pattern-mixture model. In the case of two inteventions, rnmamod performs one-stage Bayesian pairwise meta-analysis. The package handles a data-frame of binary and continuous outcome data in the arm-based format. The odds ratio, mean difference, standardised mean difference, and ratio of means are currently considered. The pattern-mixture model allows the incorporation of the informative missingness odds ratio for binary outcomes, whilst the informative missingness difference of means and the informative missingness ratio of means for continuous outcomes. The package comprises a suite of all necessary models for estimation and prediction of the intervention effect, and evaluation of the consistency assumption locally and globally. Missing participant outcome data are addressed in all models of the rnmamod package. The rnmamod package also includes a rich suite of visualisation tools that aid the interpretation and accommodation of the results in the submitted research work for publication. 

The rnmamod package is currently in development version.

## Installation

Run the following code to install rnmamod:

    devtools::install_github("LoukiaSpin/rnmamod")
    library(rnmamod)

## Example

The following code performs a Bayesian random-effects network meta-analysis under the missing at random assumption and using intervention-specific informative missingness odds ratio in the logarithmic scale:

