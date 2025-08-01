#' @name rnmamod-package
#' @aliases rnmamod rnmamod-package
#'
#' @title rnmamod: Bayesian Network Meta-analysis with Missing Participants
#'
#' @description
#'   An R package for performing Bayesian network meta-analysis while handling
#'   missing participant outcome data properly, assessing the robustness of the
#'   primary analysis results, and exploring the transitivity assumption.
#'
#' @details
#'   R-package \bold{rnmamod} is built upon the WinBUGS program code found in
#'   the series of tutorial papers on evidence synthesis methods for decision
#'   making (Dias et al., 2013a; Dias et al., 2013b; Dias et al., 2013c) and
#'   Dias et al. (2010) that introduces the node-splitting approach.
#'   All models comprise Bayesian hierarchical models for one-stage network
#'   meta-analysis and they are implemented in JAGS through the R-package
#'   \bold{R2jags}.
#'
#'   \bold{rnmamod} comprises a suite of core models implemented in a
#'   systematic review with multiple interventions:
#'   \itemize{
#'    \item fixed-effect and random-effects network meta-analysis
#'    (\code{\link{run_model}}) based on Dias et al. (2013c);
#'    \item fixed-effect and random-effects network meta-regression
#'    (\code{\link{run_metareg}}) based on Cooper et al. (2009), and Dias et
#'    al. (2013b);
#'    \item fixed-effect and random-effects separate pairwise meta-analyses for
#'    comparisons with at least two trials (\code{\link{run_series_meta}});
#'    \item local evaluation of the consistency assumption using the
#'    fixed-effect or random-effects node-splitting approach
#'    (\code{\link{run_nodesplit}}) based on Dias et al. (2010), and
#'    van Valkenhoef et al. (2016);
#'    \item global evaluation of the consistency assumption using the
#'    fixed-effect or random-effects unrelated mean effects model
#'    (\code{\link{run_ume}}) based on Dias et al. (2013a) and Spineli (2021);
#'    \item comprehensive sensitivity analysis for the impact of aggregate
#'    binary and continuous missing participant outcome data
#'    (\code{\link{run_sensitivity}}) based on Spineli et al. (2021a);
#'    \item hierarchical baseline model for the selected reference intervention
#'    (\code{\link{baseline_model}}) based in Dias et al. (2013d).
#'   }
#'
#'   \bold{rnmamod} also includes a rich suite of visualisation tools to aid in
#'   the interpretation of the results and preparation of the manuscript for
#'   submission:
#'   \itemize{
#'    \item network plot and description of the evidence base
#'    (\code{\link{netplot}} and \code{\link{describe_network}},
#'    respectively) following the  PRISMA statement for systematic reviews with
#'    network meta-analysis (Hutton et al., 2015);
#'    \item illustration of the R-hat (Gelman and Rubin, 1992) and MCMC error
#'    for all monitored nodes and creation of an HTML file with a panel of
#'    diagnostic plots for each monitored parameter
#'    (\code{\link{mcmc_diagnostics}});
#'    \item heatmap on the proportion of missing participants across the network
#'    (\code{\link{heatmap_missing_network}}) and across the intervention arms
#'    of each trial in the dataset (\code{\link{heatmap_missing_dataset}});
#'    \item league heatmap with the estimated and predicted summary effects of
#'    all possible pairwise comparisons in the network and integrated SUCRA
#'    (Salanti et al., 2011) or P-scores (Ruecker and Schwarzer, 2015)
#'    (\code{\link{league_heatmap}} and
#'    \code{\link{league_heatmap_pred}}, respectively) after performing network
#'    meta-analysis or network meta-regression;
#'    \item league table for relative and absolute effects for all pairwise
#'    comparisons and interventions when conducting network meta-analysis anew
#'    via the package (\code{\link{league_table_absolute}}) or using the results
#'    of a published systematic review with network meta-analysis
#'    (\code{\link{league_table_absolute_user}});
#'    \item forest plot with the trial-specific and summary absolute risks when
#'    employing the hierarchical baseline model for the selected reference
#'    intervention (\code{\link{baseline_model}}) as described in
#'    Dias et al. (2013d);
#'    \item rankograms with integrated SUCRA values for each intervention in
#'    the network (\code{\link{rankosucra_plot}}) after performing network
#'    meta-analysis (Salanti et al., 2011);
#'    \item forest plot with the estimated and predicted summary effects of all
#'    comparisons with a selected intervention (\code{\link{forestplot}}) as
#'    obtained from the network meta-analysis model, and a forest plot with the
#'    corresponding SUCRA values (Salanti et al., 2011);
#'    \item tabulation of the estimated regression coefficient(s), the estimated
#'    and predicted summary effects, measures of model fit and estimated
#'    between-trial standard deviation before and after adjusting for a
#'    trial-specific covariate (\code{\link{metareg_plot}}), and visualisation
#'    of the summary effects and SUCRA values from both models
#'    (\code{\link{forestplot_metareg}}, and \code{\link{scatterplot_sucra}},
#'    respectively--both found in \code{\link{metareg_plot}});
#'    \item tabulation of the estimated direct and indirect effects of the split
#'    nodes and corresponding inconsistency factors, measures of model fit and
#'    estimated between-trial standard deviation after each split node, and
#'    visualisation of these results (\code{\link{nodesplit_plot}});
#'    \item calculation of the mean and standard deviation for the prior
#'    distribution for the inconsistency variance, and the median inconsistency
#'    standard deviation (\code{\link{inconsistency_variance_prior}}) based on a
#'    selected empirical prior distribution for the between-study variance
#'    (proposed by Turner et al., 2015 and Rhodes et al., 2015);
#'    \item a panel of density plots for each target comparison (based on the
#'    back-calculation or the node-splitting approach) illustrating the
#'    posterior distribution of direct and indirect estimates, the inconsistency
#'    parameter estimate and 95% interval and the Kullback-Leibler divergence
#'    value and a percent stacked bar plot on the percentage contribution of
#'    approximating direct with indirect estimate (and vice-versa) to the total
#'    information loss for each split node (\code{\link{kld_inconsistency}} and
#'    \code{\link{kld_inconsistency_user}} when using the results of a published
#'    systematic review with network meta-analysis) (Spineli, 2024);
#'    \item tabulation of the estimated summary effects of all comparisons
#'    observed in the network, measures of model fit and estimated between-trial
#'    standard deviation under the unrelated mean effects model and network
#'    meta-analysis, as well as visualisation of the summary effects from both
#'    models (\code{\link{intervalplot_panel_ume}}) and the goodness of fit of
#'    each model using a series of complementary plots
#'    (\code{\link{scatterplots_dev}} (Dias et al., 2013a),
#'    \code{\link{bland_altman_plot}} (Bland and Altman, 1999), and
#'    \code{\link{leverage_plot}} (Dias et al., 2010)--all found in
#'    \code{\link{ume_plot}});
#'    \item tabulation of the estimated summary effects and corresponding
#'    between-trial standard deviation for comparisons with at least two trials
#'    under pairwise and network meta-analysis, as well as visualisation of
#'    these results (\code{\link{series_meta_plot}});
#'    \item calculation and visualisation of the robustness index for all
#'    possible comparisons in the network (\code{\link{robustness_index}},
#'    \code{\link{robustness_index_user}} and \code{\link{heatmap_robustness}})
#'    (Spineli et al., 2021a);
#'    \item enhanced balloon plot with the summary effects and between-trial
#'    standard deviation for a selected pairwise comparison under several
#'    scenarios about the missingness parameter (\code{\link{balloon_plot}})
#'    (Spineli et al., 2021a);
#'    \item barplot with the Kullback-Leibler divergence measure from each
#'    informative scenario to the missing-at-random assumption about the
#'    missingness parameter for a selected pairwise comparison
#'    (\code{\link{kld_barplot}}) (Spineli et al., 2021a).
#'   }
#'
#'   \bold{rnmamod} also assists the researcher in assessing the transitivity
#'   assumption quantitatively based on trial dissimilarities for various
#'   trial-level aggregate participant and methodological characteristics
#'   calculated using the Gower's dissimilarity coefficient (Gower, 1971)
#'   (\code{\link{gower_distance}} and \code{\link{comp_clustering}}) (Spineli
#'   et al., 2025; Spineli, 2024). Results on the clustered comparisons based on
#'   hierarchical agglomerative clustering are illustrated using a dendrogram
#'   with integrated heatmap (\code{\link{dendro_heatmap}}) (Spineli et al., 2025).
#'   The distribution of the characteristics is presented using violin plots
#'   with integrated box plots and dots, and stacked bar plots across the
#'   observed treatment comparisons (\code{\link{distr_characteristics}}).
#'   Missing data in the characteristics across the trials and observed
#'   comparisons are visualised using bar plots and tile plot
#'   (\code{\link{miss_characteristics}}).
#'
#'   Missing participant outcome data are addressed in all models of the package
#'   after extending the code to incorporate the pattern-mixture model
#'   (Spineli et al., 2021b; Spineli, 2019).
#'
#'   Type \code{citation("rnmamod")} on how to cite \bold{rnmamod}.
#'
#'   To report possible bugs and errors, send an email to Loukia Spineli
#'   (\email{Spineli.Loukia@mh-hannove.de}).
#'
#'   The development version of \bold{rnmamod} is available on
#'   \href{https://github.com/LoukiaSpin/rnmamod}{GitHub} under the
#'   GPL-3.0 License.
#'
#' @author {Loukia M. Spineli}
#'
#' @references
#' Bland JM, Altman DG. Measuring agreement in method comparison studies.
#' \emph{Stat Methods Med Res} 1999;\bold{8}(2):135--60.
#' doi: 10.1177/096228029900800204
#'
#' Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing between-study
#' heterogeneity and inconsistency in mixed treatment comparisons: Application
#' to stroke prevention treatments in individuals with non-rheumatic atrial
#' fibrillation. \emph{Stat Med} 2009;\bold{28}(14):1861--81.
#' doi: 10.1002/sim.3594
#'
#' Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence synthesis
#' for decision making 4: inconsistency in networks of evidence based on
#' randomized controlled trials.
#' \emph{Med Decis Making} 2013a;\bold{33}(5):641--56.
#' doi: 10.1177/0272989X12455847
#'
#' Dias S, Sutton AJ, Welton NJ, Ades AE. Evidence synthesis for decision making
#' 3: heterogeneity--subgroups, meta-regression, bias, and bias-adjustment.
#' \emph{Med Decis Making} 2013b;\bold{33}(5):618--40.
#' doi: 10.1177/0272989X13485157
#'
#' Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
#' making 2: a generalized linear modeling framework for pairwise and network
#' meta-analysis of randomized controlled trials. \emph{Med Decis Making}
#' 2013c;\bold{33}(5):607--17. doi: 10.1177/0272989X12458724
#'
#' Dias S, Welton NJ, Sutton AJ, Ades AE. Evidence synthesis for decision
#' making 5: the baseline natural history model. \emph{Med Decis Making}
#' 2013d;\bold{33}(5):657--70. doi: 10.1177/0272989X13485155
#'
#' Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
#' treatment comparison meta-analysis.
#' \emph{Stat Med} 2010;\bold{29}(7-8):932--44.
#' doi: 10.1002/sim.3767
#'
#' Gelman, A, Rubin, DB. Inference from iterative simulation using multiple
#' sequences. \emph{Stat Sci} 1992;\bold{7}(4):457--72.
#' doi: 10.1214/ss/1177011136
#'
#' Gower JC. A General Coefficient of Similarity and Some of Its Properties.
#' \emph{Biometrics} 1971;\bold{27}(4):857--71.
#' http://dx.doi.org/10.2307/2528823
#'
#' Hutton B, Salanti G, Caldwell DM, Chaimani A, Schmid CH, Cameron C, et al.
#' The PRISMA extension statement for reporting of systematic reviews
#' incorporating network meta-analyses of health care interventions: checklist
#' and explanations. \emph{Ann Intern Med} 2015;\bold{162}(11):777--84.
#' doi: 10.7326/M14-2385
#'
#' Rhodes KM, Turner RM, Higgins JP. Predictive distributions were developed for
#' the extent of heterogeneity in meta-analyses of continuous outcome data.
#' \emph{J Clin Epidemiol} 2015;\bold{68}(1):52--60.
#' doi: 10.1016/j.jclinepi.2014.08.012
#'
#' Ruecker G, Schwarzer G. Ranking treatments in frequentist network
#' meta-analysis works without resampling methods.
#' \emph{BMC Med Res Methodol} 2015;\bold{15}:58.
#' doi: 10.1186/s12874-015-0060-8
#'
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries
#' for presenting results from multiple-treatment meta-analysis: an overview and
#' tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71.
#' doi: 10.1016/j.jclinepi.2010.03.016
#'
#' Spineli LM, Papadimitropoulou K, Kalyvas C. Exploring the Transitivity
#' Assumption in Network Meta-Analysis: A Novel Approach and Its Implications.
#' \emph{Stat Med} 2025;\bold{44}(7):e70068.
#' doi: 10.1002/sim.70068.
#'
#' Spineli LM. An empirical study on 209 networks of treatments revealed
#' intransitivity to be common and multiple statistical tests suboptimal to
#' assess transitivity. \emph{BMC Med Res Methodol} 2024;\bold{24}(1):301.
#' doi: 10.1186/s12874-024-02436-7.
#'
#' Spineli LM. Local inconsistency detection using the Kullback-Leibler
#' divergence measure. \emph{Syst Rev} 2024;\bold{13}(1):261.
#' doi: 10.1186/s13643-024-02680-4.
#'
#' Spineli LM. A revised framework to evaluate the consistency assumption
#' globally in a network of interventions. \emph{Med Decis Making} 2021.
#' doi: 10.1177/0272989X211068005
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of
#' primary analysis results: A case study on missing outcome data in pairwise
#' and network meta-analysis.
#' \emph{Res Synth Methods} 2021a;\bold{12}(4):475--90. doi: 10.1002/jrsm.1478
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing outcome
#' data in network meta-analysis: a one-stage pattern-mixture model approach.
#' \emph{Stat Methods Med Res} 2021b;\bold{30}(4):958--75.
#' doi: 10.1177/0962280220983544
#'
#' Spineli LM. An empirical comparison of Bayesian modelling strategies for
#' missing binary outcome data in network meta-analysis.
#' \emph{BMC Med Res Methodol} 2019;\bold{19}(1):86.
#' doi: 10.1186/s12874-019-0731-y
#'
#' Turner RM, Jackson D, Wei Y, Thompson SG, Higgins JP. Predictive distributions
#' for between-study heterogeneity and simple methods for their application in
#' Bayesian meta-analysis. \emph{Stat Med} 2015;\bold{34}(6):984--98.
#' doi: 10.1002/sim.6381
#'
#' van Valkenhoef G, Dias S, Ades AE, Welton NJ. Automated generation of
#' node-splitting models for assessment of inconsistency in network
#' meta-analysis. \emph{Res Synth Methods} 2016;\bold{7}(1):80--93.
#' doi: 10.1002/jrsm.1167
#'
#' @keywords package
NULL
