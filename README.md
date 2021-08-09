*** 

# rnmamod: package to perform Bayesian network meta-analysis methods

rnmamod is an R package to perform one-stage Bayesian network meta-analysis while adjusting for missing participant outcome data using the pattern-mixture model. The package handles a data-frame of binary and continuous outcome data in the arm-based format. The folloiwing effect measures are currently considered: the odds ratio, mean difference, standardised mean difference, and ratio of means. The pattern-mixture models allows the accommodation of the informative missingness odds ration for binary outcomes, whilst the informative missingeness difference of means and the informative missingness ratio of means for continuous outcomes. The package comprises a suite of all necessary models for estimation and prediction of the intervention effect, and evaluation of the consistency assumption locally and globally.

