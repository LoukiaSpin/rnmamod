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
