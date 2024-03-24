# rnmamod, version 0.4.0 (2024-03-24)

 * Function __comp_clustering__: 
   - Performs quantitative evaluation of the transitivity assumption using
   inter-trial dissimilarities for various trial-level aggregate participant 
   and methodological characteristics that may act as effect modifiers. 
 * Function __dendro_heatmap__:  
   - Returns the dendrogram with integrated heatmap of the clustered comparisons 
   and trials based on hierarchical agglomerative clustering (performed using 
   the function __comp_clustering__). The R packages _heatmaply_ and 
   _dendextend_ have been used.
 * Function __distr_characteristics__:  
   - It returns violin plots with integrated box plots and dots for quantitative 
   characteristics, and stacked barplots for qualitative characteristics across
   the observed treatment comparisons. The function can also be used to 
   illustrate the distribution of the characteristics across the clusters 
   defined from __comp_clustering__.
 * Function __miss_characteristics__:  
   - It returns various plots to visualise the missing cases in the 
   characteristics across trials and treatment comparisons.
 * Function __gower_distance__: 
   - It returns the N-by-N matrix on Gower's dissimilarity coefficient for all
   pairs of N trials in a network.
 * Function __mcmc_diagnostics__: 
   - returns a bar plot on the ratio of MCMC error to the posterior standard
   deviation and a bar plot on the Gelman-Rubin R diagnostic. Green bars
   indicate ratio less than 0.05 and R less than 1.10; otherwise, the bars are
   red.
 * Functions __baseline_model__, __run_metareg__, __run_model__, 
   __run_nodesplit__, __run_sensitivity__, __run_series_meta__, and __run_ume__: 
   - The corresponding models are updated until convergence via the _autojags_
   function of the R package _R2jags_.
   - The argument _inits_ has been added to allow the user define the initial
   values for the parameters, following the documentation of the _jags_ function
   in the R package _R2jags_. 
 * Function __describe_network__:
   - It reports only the tables with the evidence base information on one 
   outcome. The network plot is not reported (see and use __netplot__, instead).
 * Function __netplot__:
   - Self-created function using the R package _igraph_. This function creates
   the network plot.

# rnmamod, version 0.3.0 (2022-11-01)

 * Function __baseline_model__:
   - processes the elements in the argument _base_risk_ for a fixed, random or 
   predicted baseline model and passes the output to run_model or run_metareg to 
   obtain the absolute risks for all interventions.
   - when a predicted baseline model is conducted, it returns a forest plot with
   the trial-specific and summary probability of an event for the selected
   reference intervention.
 * Function __forestplot_metareg__:
   - upgraded plot that presents two forest plots side-by-side: (i) one on 
   estimation and prediction from network meta-analysis and network 
   meta-regression for a selected comparator intervention (allows comparison of 
   these two analyses), and (ii) one on SUCRA values from both analyses. 
   Both forest plots present results from network meta-regression for a selected 
   value of the investigated covariate.
 * Function __league_table_absolute_user__:
   - (only for binary outcome) yields the same graph with league_table_absolute,
   but the input is not _rnmamod_ object: the user defines the input and it
   includes the summary effect and corresponding (credible or confidence) 
   interval for comparisons with a reference intervention. These results stem 
   from a network meta-analysis conducted using another R-package or statistical 
   software.
 * Function __robustness_index_user__:
   - calculates the robustness index for a sensitivity analysis performed using 
   the R-package _netmeta_ or _metafor_. The user defines the input and the 
   function returns the robustness index. This function returns the same output 
   with the __robustness_index__ function.
 * Function __trans_quality__:
   - classifies a systematic review with multiple interventions as having low, 
   unclear or high quality regarding the transitivity evaluation based on five
   quality criteria.

# rnmamod, version 0.2.0 (2022-04-06)

 * Typos and links (for functions and packages) in the documentation are 
 corrected.
 * Function __run_model__:
   - allows the user to define the reference intervention of the network via the
   argument _ref_;
   - (only for binary outcome) estimates the absolute risks for all 
   non-reference interventions using a selected baseline risk for the reference 
   intervention (argument _base_risk_); 
   - (only for binary outcome) estimates the relative risks and risk difference 
   as functions of the estimated absolute risks.
 * Function __league_table_absolute__:
   - (only for binary outcome) it presents the absolute risks per 1000 
   participants in main diagonal, the odds ratio on the upper off-diagonals, and 
   the risk difference per 1000 participants in the lower off-diagonals;
   - allow the user to select the interventions to present via the argument 
   _show_ (ideal for very large networks that make the league table cluttered).
 * Functions __league_heatmap__ and __league_heatmap_pred__:
   - allow the user to select the interventions to present via the argument 
   _show_ (ideal for very large networks that make the league table cluttered);
   - allow the user to illustrate the results of two outcomes for the same model
   (i.e. via run_model or run_metareg) or the results of two models on the same 
   outcome (applicable for: (i) run_model versus run_metareg, and (ii) run_model 
   versus run_series_meta).
 * Functions __series_meta_plot__ and __nodesplit_plot__:
   - present the extent of heterogeneity in the forest plot of tau using 
   different colours for low, reasonable, fairly high, and fairly extreme tau 
   (the classification has been suggested by Spiegelhalter et al., 2004; 
   ISBN 0471499757).

# rnmamod, version 0.1.0 (2021-11-21)

 - First CRAN submission.
