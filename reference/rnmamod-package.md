# rnmamod: Bayesian Network Meta-analysis with Missing Participants

An R package for performing Bayesian network meta-analysis while
handling missing participant outcome data properly, assessing the
robustness of the primary analysis results, and exploring the
transitivity assumption.

## Details

R-package **rnmamod** is built upon the WinBUGS program code found in
the series of tutorial papers on evidence synthesis methods for decision
making (Dias et al., 2013a; Dias et al., 2013b; Dias et al., 2013c) and
Dias et al. (2010) that introduces the node-splitting approach. All
models comprise Bayesian hierarchical models for one-stage network
meta-analysis and they are implemented in JAGS through the R-package
**R2jags**.

**rnmamod** comprises a suite of core models implemented in a systematic
review with multiple interventions:

- fixed-effect and random-effects network meta-analysis
  ([`run_model`](https://loukiaspin.github.io/rnmamod/reference/run_model.md))
  based on Dias et al. (2013c);

- fixed-effect and random-effects network meta-regression
  ([`run_metareg`](https://loukiaspin.github.io/rnmamod/reference/run_metareg.md))
  based on Cooper et al. (2009), and Dias et al. (2013b);

- fixed-effect and random-effects separate pairwise meta-analyses for
  comparisons with at least two trials
  ([`run_series_meta`](https://loukiaspin.github.io/rnmamod/reference/run_series_meta.md));

- local evaluation of the consistency assumption using the fixed-effect
  or random-effects node-splitting approach
  ([`run_nodesplit`](https://loukiaspin.github.io/rnmamod/reference/run_nodesplit.md))
  based on Dias et al. (2010), and van Valkenhoef et al. (2016);

- global evaluation of the consistency assumption using the fixed-effect
  or random-effects unrelated mean effects model
  ([`run_ume`](https://loukiaspin.github.io/rnmamod/reference/run_ume.md))
  based on Dias et al. (2013a) and Spineli (2021);

- comprehensive sensitivity analysis for the impact of aggregate binary
  and continuous missing participant outcome data
  ([`run_sensitivity`](https://loukiaspin.github.io/rnmamod/reference/run_sensitivity.md))
  based on Spineli et al. (2021a);

- hierarchical baseline model for the selected reference intervention
  ([`baseline_model`](https://loukiaspin.github.io/rnmamod/reference/baseline_model.md))
  based in Dias et al. (2013d).

**rnmamod** also includes a rich suite of visualisation tools to aid in
the interpretation of the results and preparation of the manuscript for
submission:

- network plot and description of the evidence base
  ([`netplot`](https://loukiaspin.github.io/rnmamod/reference/netplot.md)
  and
  [`describe_network`](https://loukiaspin.github.io/rnmamod/reference/describe_network.md),
  respectively) following the PRISMA statement for systematic reviews
  with network meta-analysis (Hutton et al., 2015);

- illustration of the R-hat (Gelman and Rubin, 1992) and MCMC error for
  all monitored nodes and creation of an HTML file with a panel of
  diagnostic plots for each monitored parameter
  ([`mcmc_diagnostics`](https://loukiaspin.github.io/rnmamod/reference/mcmc_diagnostics.md));

- heatmap on the proportion of missing participants across the network
  ([`heatmap_missing_network`](https://loukiaspin.github.io/rnmamod/reference/heatmap_missing_network.md))
  and across the intervention arms of each trial in the dataset
  ([`heatmap_missing_dataset`](https://loukiaspin.github.io/rnmamod/reference/heatmap_missing_dataset.md));

- league heatmap with the estimated and predicted summary effects of all
  possible pairwise comparisons in the network and integrated SUCRA
  (Salanti et al., 2011) or P-scores (Ruecker and Schwarzer, 2015)
  ([`league_heatmap`](https://loukiaspin.github.io/rnmamod/reference/league_heatmap.md)
  and
  [`league_heatmap_pred`](https://loukiaspin.github.io/rnmamod/reference/league_heatmap_pred.md),
  respectively) after performing network meta-analysis or network
  meta-regression;

- league table for relative and absolute effects for all pairwise
  comparisons and interventions when conducting network meta-analysis
  anew via the package
  ([`league_table_absolute`](https://loukiaspin.github.io/rnmamod/reference/league_table_absolute.md))
  or using the results of a published systematic review with network
  meta-analysis
  ([`league_table_absolute_user`](https://loukiaspin.github.io/rnmamod/reference/league_table_absolute_user.md));

- forest plot with the trial-specific and summary absolute risks when
  employing the hierarchical baseline model for the selected reference
  intervention
  ([`baseline_model`](https://loukiaspin.github.io/rnmamod/reference/baseline_model.md))
  as described in Dias et al. (2013d);

- rankograms with integrated SUCRA values for each intervention in the
  network
  ([`rankosucra_plot`](https://loukiaspin.github.io/rnmamod/reference/rankosucra_plot.md))
  after performing network meta-analysis (Salanti et al., 2011);

- forest plot with the estimated and predicted summary effects of all
  comparisons with a selected intervention
  ([`forestplot`](https://loukiaspin.github.io/rnmamod/reference/forestplot.md))
  as obtained from the network meta-analysis model, and a forest plot
  with the corresponding SUCRA values (Salanti et al., 2011);

- tabulation of the estimated regression coefficient(s), the estimated
  and predicted summary effects, measures of model fit and estimated
  between-trial standard deviation before and after adjusting for a
  trial-specific covariate
  ([`metareg_plot`](https://loukiaspin.github.io/rnmamod/reference/metareg_plot.md)),
  and visualisation of the summary effects and SUCRA values from both
  models
  ([`forestplot_metareg`](https://loukiaspin.github.io/rnmamod/reference/forestplot_metareg.md),
  and
  [`scatterplot_sucra`](https://loukiaspin.github.io/rnmamod/reference/scatterplot_sucra.md),
  respectively–both found in
  [`metareg_plot`](https://loukiaspin.github.io/rnmamod/reference/metareg_plot.md));

- tabulation of the estimated direct and indirect effects of the split
  nodes and corresponding inconsistency factors, measures of model fit
  and estimated between-trial standard deviation after each split node,
  and visualisation of these results
  ([`nodesplit_plot`](https://loukiaspin.github.io/rnmamod/reference/nodesplit_plot.md));

- calculation of the mean and standard deviation for the prior
  distribution for the inconsistency variance, and the median
  inconsistency standard deviation
  ([`inconsistency_variance_prior`](https://loukiaspin.github.io/rnmamod/reference/inconsistency_variance_prior.md))
  based on a selected empirical prior distribution for the between-study
  variance (proposed by Turner et al., 2015 and Rhodes et al., 2015);

- a panel of density plots for each target comparison (based on the
  back-calculation or the node-splitting approach) illustrating the
  posterior distribution of direct and indirect estimates, the
  inconsistency parameter estimate and 95 value and a percent stacked
  bar plot on the percentage contribution of approximating direct with
  indirect estimate (and vice-versa) to the total information loss for
  each split node
  ([`kld_inconsistency`](https://loukiaspin.github.io/rnmamod/reference/kld_inconsistency.md)
  and
  [`kld_inconsistency_user`](https://loukiaspin.github.io/rnmamod/reference/kld_inconsistency_user.md)
  when using the results of a published systematic review with network
  meta-analysis) (Spineli, 2024);

- tabulation of the estimated summary effects of all comparisons
  observed in the network, measures of model fit and estimated
  between-trial standard deviation under the unrelated mean effects
  model and network meta-analysis, as well as visualisation of the
  summary effects from both models
  ([`intervalplot_panel_ume`](https://loukiaspin.github.io/rnmamod/reference/intervalplot_panel_ume.md))
  and the goodness of fit of each model using a series of complementary
  plots
  ([`scatterplots_dev`](https://loukiaspin.github.io/rnmamod/reference/scatterplots_dev.md)
  (Dias et al., 2013a),
  [`bland_altman_plot`](https://loukiaspin.github.io/rnmamod/reference/bland_altman_plot.md)
  (Bland and Altman, 1999), and
  [`leverage_plot`](https://loukiaspin.github.io/rnmamod/reference/leverage_plot.md)
  (Dias et al., 2010)–all found in
  [`ume_plot`](https://loukiaspin.github.io/rnmamod/reference/ume_plot.md));

- tabulation of the estimated summary effects and corresponding
  between-trial standard deviation for comparisons with at least two
  trials under pairwise and network meta-analysis, as well as
  visualisation of these results
  ([`series_meta_plot`](https://loukiaspin.github.io/rnmamod/reference/series_meta_plot.md));

- calculation and visualisation of the robustness index for all possible
  comparisons in the network
  ([`robustness_index`](https://loukiaspin.github.io/rnmamod/reference/robustness_index.md),
  [`robustness_index_user`](https://loukiaspin.github.io/rnmamod/reference/robustness_index_user.md)
  and
  [`heatmap_robustness`](https://loukiaspin.github.io/rnmamod/reference/heatmap_robustness.md))
  (Spineli et al., 2021a);

- enhanced balloon plot with the summary effects and between-trial
  standard deviation for a selected pairwise comparison under several
  scenarios about the missingness parameter
  ([`balloon_plot`](https://loukiaspin.github.io/rnmamod/reference/balloon_plot.md))
  (Spineli et al., 2021a);

- barplot with the Kullback-Leibler divergence measure from each
  informative scenario to the missing-at-random assumption about the
  missingness parameter for a selected pairwise comparison
  ([`kld_barplot`](https://loukiaspin.github.io/rnmamod/reference/kld_barplot.md))
  (Spineli et al., 2021a).

**rnmamod** also assists the researcher in assessing the transitivity
assumption quantitatively based on trial dissimilarities for various
trial-level aggregate participant and methodological characteristics
calculated using the Gower's dissimilarity coefficient (Gower, 1971)
([`gower_distance`](https://loukiaspin.github.io/rnmamod/reference/gower_distance.md)
and
[`comp_clustering`](https://loukiaspin.github.io/rnmamod/reference/comp_clustering.md))
(Spineli et al., 2025; Spineli, 2024). Results on the clustered
comparisons based on hierarchical agglomerative clustering are
illustrated using a dendrogram with integrated heatmap
([`dendro_heatmap`](https://loukiaspin.github.io/rnmamod/reference/dendro_heatmap.md))
(Spineli et al., 2025). The distribution of the characteristics is
presented using violin plots with integrated box plots and dots, and
stacked bar plots across the observed treatment comparisons
([`distr_characteristics`](https://loukiaspin.github.io/rnmamod/reference/distr_characteristics.md)).
Missing data in the characteristics across the trials and observed
comparisons are visualised using bar plots and tile plot
([`miss_characteristics`](https://loukiaspin.github.io/rnmamod/reference/miss_characteristics.md)).

Missing participant outcome data are addressed in all models of the
package after extending the code to incorporate the pattern-mixture
model (Spineli et al., 2021b; Spineli, 2019).

Type `citation("rnmamod")` on how to cite **rnmamod**.

To report possible bugs and errors, send an email to Loukia Spineli
(<Spineli.Loukia@mh-hannove.de>).

The development version of **rnmamod** is available on
[GitHub](https://github.com/LoukiaSpin/rnmamod) under the GPL-3.0
License.

## References

Bland JM, Altman DG. Measuring agreement in method comparison studies.
*Stat Methods Med Res* 1999;**8**(2):135–60. doi:
10.1177/096228029900800204

Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing
between-study heterogeneity and inconsistency in mixed treatment
comparisons: Application to stroke prevention treatments in individuals
with non-rheumatic atrial fibrillation. *Stat Med*
2009;**28**(14):1861–81. doi: 10.1002/sim.3594

Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence
synthesis for decision making 4: inconsistency in networks of evidence
based on randomized controlled trials. *Med Decis Making*
2013a;**33**(5):641–56. doi: 10.1177/0272989X12455847

Dias S, Sutton AJ, Welton NJ, Ades AE. Evidence synthesis for decision
making 3: heterogeneity–subgroups, meta-regression, bias, and
bias-adjustment. *Med Decis Making* 2013b;**33**(5):618–40. doi:
10.1177/0272989X13485157

Dias S, Sutton AJ, Ades AE, Welton NJ. Evidence synthesis for decision
making 2: a generalized linear modeling framework for pairwise and
network meta-analysis of randomized controlled trials. *Med Decis
Making* 2013c;**33**(5):607–17. doi: 10.1177/0272989X12458724

Dias S, Welton NJ, Sutton AJ, Ades AE. Evidence synthesis for decision
making 5: the baseline natural history model. *Med Decis Making*
2013d;**33**(5):657–70. doi: 10.1177/0272989X13485155

Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
treatment comparison meta-analysis. *Stat Med* 2010;**29**(7-8):932–44.
doi: 10.1002/sim.3767

Gelman, A, Rubin, DB. Inference from iterative simulation using multiple
sequences. *Stat Sci* 1992;**7**(4):457–72. doi: 10.1214/ss/1177011136

Gower JC. A General Coefficient of Similarity and Some of Its
Properties. *Biometrics* 1971;**27**(4):857–71.
http://dx.doi.org/10.2307/2528823

Hutton B, Salanti G, Caldwell DM, Chaimani A, Schmid CH, Cameron C, et
al. The PRISMA extension statement for reporting of systematic reviews
incorporating network meta-analyses of health care interventions:
checklist and explanations. *Ann Intern Med* 2015;**162**(11):777–84.
doi: 10.7326/M14-2385

Rhodes KM, Turner RM, Higgins JP. Predictive distributions were
developed for the extent of heterogeneity in meta-analyses of continuous
outcome data. *J Clin Epidemiol* 2015;**68**(1):52–60. doi:
10.1016/j.jclinepi.2014.08.012

Ruecker G, Schwarzer G. Ranking treatments in frequentist network
meta-analysis works without resampling methods. *BMC Med Res Methodol*
2015;**15**:58. doi: 10.1186/s12874-015-0060-8

Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical
summaries for presenting results from multiple-treatment meta-analysis:
an overview and tutorial. *J Clin Epidemiol* 2011;**64**(2):163–71. doi:
10.1016/j.jclinepi.2010.03.016

Spineli LM, Papadimitropoulou K, Kalyvas C. Exploring the Transitivity
Assumption in Network Meta-Analysis: A Novel Approach and Its
Implications. *Stat Med* 2025;**44**(7):e70068. doi: 10.1002/sim.70068.

Spineli LM. An empirical study on 209 networks of treatments revealed
intransitivity to be common and multiple statistical tests suboptimal to
assess transitivity. *BMC Med Res Methodol* 2024;**24**(1):301. doi:
10.1186/s12874-024-02436-7.

Spineli LM. Local inconsistency detection using the Kullback-Leibler
divergence measure. *Syst Rev* 2024;**13**(1):261. doi:
10.1186/s13643-024-02680-4.

Spineli LM. A revised framework to evaluate the consistency assumption
globally in a network of interventions. *Med Decis Making* 2021. doi:
10.1177/0272989X211068005

Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness
of primary analysis results: A case study on missing outcome data in
pairwise and network meta-analysis. *Res Synth Methods*
2021a;**12**(4):475–90. doi: 10.1002/jrsm.1478

Spineli LM, Kalyvas C, Papadimitropoulou K. Continuous(ly) missing
outcome data in network meta-analysis: a one-stage pattern-mixture model
approach. *Stat Methods Med Res* 2021b;**30**(4):958–75. doi:
10.1177/0962280220983544

Spineli LM. An empirical comparison of Bayesian modelling strategies for
missing binary outcome data in network meta-analysis. *BMC Med Res
Methodol* 2019;**19**(1):86. doi: 10.1186/s12874-019-0731-y

Turner RM, Jackson D, Wei Y, Thompson SG, Higgins JP. Predictive
distributions for between-study heterogeneity and simple methods for
their application in Bayesian meta-analysis. *Stat Med*
2015;**34**(6):984–98. doi: 10.1002/sim.6381

van Valkenhoef G, Dias S, Ades AE, Welton NJ. Automated generation of
node-splitting models for assessment of inconsistency in network
meta-analysis. *Res Synth Methods* 2016;**7**(1):80–93. doi:
10.1002/jrsm.1167

## Author

Loukia M. Spineli
