# Changelog

## rnmamod, version 0.5.0 (2025-06-13)

CRAN release: 2025-06-13

- Replaced **mcmcplots** with **coda** R package.
- Function **plot_study_dissimilarities**:
  - Presents the range of Gower’s dissimilarity values for each study in
    the network, as well as their between- and within-comparison
    dissimilarities (based on the function **comp_clustering**).
- Function **study_perc_contrib**:
  - Calculates the percentage contributions of each study to every
    possible pairwise comparison in the investigated network and returns
    a data-frame. Study percentage contributions are based on Donegan et
    al.,

  2018. <doi:10.1002/jrsm.1292>
- Function **covar_contribution_plot**:
  - Returns a scatter plot of the study percentage contributions against
    the values of a continuous study-level covariate for the treatment
    effects of the basic parameters, functional parameters or both
    (based on the function **study_perc_contrib**).
- Function **forestplot_juxtapose**:
  - Provides a forest plot juxtaposing several NMA models (via the
    functions **run_model** and **run_metareg**) based on posterior
    treatment effects (including predictions) of all treatments versus a
    selected comparator and a forest plot with the corresponding SUCRA
    values.
- Function **heter_density_plot** :
  - Creates the density plot of two prior distributions for the
    between-study variance (log-normal and location-scale t
    distributions) or between-study standard deviation (half-normal
    distribution). This plot aids in deciding how to define the argument
    *heter_prior* in **run_model** to run random-effects network
    meta-analysis.
- Function **inconsistency_variance_prior**:
  - Calculates the hyperparameters of the log-normal distribution and
    location-scale t-distribution of the inconsistency variance in the
    log-odds ratio and standardised mean difference scales,
    respectively, based on selected empirical distributions for the
    between-study variance proposed by Turner et al. (2015)
    <doi:10.1002/sim.6381> and Rhodes et al. (2015)
    <doi:10.1016/j.jclinepi.2014.08.012>. Calculations are based on Law
    et al.,

  2016. <doi:10.1186/s12874-016-0184-5>.
- Function **table_tau2_prior**:
  - Returns a table with the hyperparameters of the predictive
    distributions for the between-study variance developed by Turner et
    al. (2015) <doi:10.1002/sim.6381> and Rhodes et al. (2015)
    <doi:10.1016/j.jclinepi.2014.08.012>. This table aids in selecting
    the hyperparameters for the function **heterogeneity_param_prior**
    when considering an informative prior distribution for the
    between-study variance parameter to conduct random-effects network
    meta-analysis.

## rnmamod, version 0.4.0 (2024-03-24)

CRAN release: 2024-03-24

- Function **comp_clustering**:
  - Performs quantitative evaluation of the transitivity assumption
    using inter-trial dissimilarities for various trial-level aggregate
    participant and methodological characteristics that may act as
    effect modifiers.
- Function **dendro_heatmap**:
  - Returns the dendrogram with integrated heatmap of the clustered
    comparisons and trials based on hierarchical agglomerative
    clustering (performed using the function **comp_clustering**). The R
    packages *heatmaply* and *dendextend* have been used.
- Function **distr_characteristics**:
  - It returns violin plots with integrated box plots and dots for
    quantitative characteristics, and stacked barplots for qualitative
    characteristics across the observed treatment comparisons. The
    function can also be used to illustrate the distribution of the
    characteristics across the clusters defined from
    **comp_clustering**.
- Function **miss_characteristics**:
  - It returns various plots to visualise the missing cases in the
    characteristics across trials and treatment comparisons.
- Function **gower_distance**:
  - It returns the N-by-N matrix on Gower’s dissimilarity coefficient
    for all pairs of N trials in a network.
- Function **mcmc_diagnostics**:
  - returns a bar plot on the ratio of MCMC error to the posterior
    standard deviation and a bar plot on the Gelman-Rubin R diagnostic.
    Green bars indicate ratio less than 0.05 and R less than 1.10;
    otherwise, the bars are red.
- Functions **baseline_model**, **run_metareg**, **run_model**,
  **run_nodesplit**, **run_sensitivity**, **run_series_meta**, and
  **run_ume**:
  - The corresponding models are updated until convergence via the
    *autojags* function of the R package *R2jags*.
  - The argument *inits* has been added to allow the user define the
    initial values for the parameters, following the documentation of
    the *jags* function in the R package *R2jags*.
- Function **describe_network**:
  - It reports only the tables with the evidence base information on one
    outcome. The network plot is not reported (see and use **netplot**,
    instead).
- Function **netplot**:
  - Self-created function using the R package *igraph*. This function
    creates the network plot.

## rnmamod, version 0.3.0 (2022-11-01)

CRAN release: 2022-11-01

- Function **baseline_model**:
  - processes the elements in the argument *base_risk* for a fixed,
    random or predicted baseline model and passes the output to
    run_model or run_metareg to obtain the absolute risks for all
    interventions.
  - when a predicted baseline model is conducted, it returns a forest
    plot with the trial-specific and summary probability of an event for
    the selected reference intervention.
- Function **forestplot_metareg**:
  - upgraded plot that presents two forest plots side-by-side: (i) one
    on estimation and prediction from network meta-analysis and network
    meta-regression for a selected comparator intervention (allows
    comparison of these two analyses), and (ii) one on SUCRA values from
    both analyses. Both forest plots present results from network
    meta-regression for a selected value of the investigated covariate.
- Function **league_table_absolute_user**:
  - (only for binary outcome) yields the same graph with
    league_table_absolute, but the input is not *rnmamod* object: the
    user defines the input and it includes the summary effect and
    corresponding (credible or confidence) interval for comparisons with
    a reference intervention. These results stem from a network
    meta-analysis conducted using another R-package or statistical
    software.
- Function **robustness_index_user**:
  - calculates the robustness index for a sensitivity analysis performed
    using the R-package *netmeta* or *metafor*. The user defines the
    input and the function returns the robustness index. This function
    returns the same output with the **robustness_index** function.
- Function **trans_quality**:
  - classifies a systematic review with multiple interventions as having
    low, unclear or high quality regarding the transitivity evaluation
    based on five quality criteria.

## rnmamod, version 0.2.0 (2022-04-06)

CRAN release: 2022-04-06

- Typos and links (for functions and packages) in the documentation are
  corrected.
- Function **run_model**:
  - allows the user to define the reference intervention of the network
    via the argument *ref*;
  - (only for binary outcome) estimates the absolute risks for all
    non-reference interventions using a selected baseline risk for the
    reference intervention (argument *base_risk*);
  - (only for binary outcome) estimates the relative risks and risk
    difference as functions of the estimated absolute risks.
- Function **league_table_absolute**:
  - (only for binary outcome) it presents the absolute risks per 1000
    participants in main diagonal, the odds ratio on the upper
    off-diagonals, and the risk difference per 1000 participants in the
    lower off-diagonals;
  - allow the user to select the interventions to present via the
    argument *show* (ideal for very large networks that make the league
    table cluttered).
- Functions **league_heatmap** and **league_heatmap_pred**:
  - allow the user to select the interventions to present via the
    argument *show* (ideal for very large networks that make the league
    table cluttered);
  - allow the user to illustrate the results of two outcomes for the
    same model (i.e. via run_model or run_metareg) or the results of two
    models on the same outcome (applicable for: (i) run_model versus
    run_metareg, and (ii) run_model versus run_series_meta).
- Functions **series_meta_plot** and **nodesplit_plot**:
  - present the extent of heterogeneity in the forest plot of tau using
    different colours for low, reasonable, fairly high, and fairly
    extreme tau (the classification has been suggested by Spiegelhalter
    et al., 2004; ISBN 0471499757).

## rnmamod, version 0.1.0 (2021-11-21)

CRAN release: 2021-11-29

- First CRAN submission.
