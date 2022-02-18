# rnmamod, version 0.1.0 (2021-11-21)

 - First CRAN submission.

# rnmamod, 0.1.X (2021-03-XX)

 * Typos and links (for functions and packages) in the documentation are 
 corrected.
 * Function mcmc_diagnostics:
   - a bubble plot that indicates which monitored parameters have converged 
   (using the green colour; otherwise, the red colour) based on the Rhat
   (Gelman-Rubin convergence statistic) and MCMC error (less than 1.1 and 5\%, 
   respectively).
 * Function run_model:
   - allows the user to define the reference intervention of the network via the
   argument _ref_;
   - (only for binary outcome) estimates the absolute effects for all 
   non-reference interventions using a selected baseline risk for the reference 
   intervention (argument _base_risk_); 
   - (only for binary outcome) estimates the relative risks and risk difference 
   as functions of the estimated absolute risks.
 * Function league_heatmap_absolute:
   - (only for binary outcome) it has the same output with the league_heatmap 
   function; however, it presents the absolute risks per 1000 participants in 
   main diagonal, the odds ratio on the upper off-diagonals, and the relative 
   risks (or risk difference per 1000 participants) in the lower 
   off-diagonals---the user can define the latter.
 * Functions league_heatmap and league_heatmap_pred:
   - allow the user to select the interventions to present via the argument 
   _show_ (ideal for very large networks that make the league table cluttered);
   - allow the user to illustrate the results of two outcomes for the same model
   (i.e. via run_model or run_metareg) or the results of two models on the same 
   outcome (plausible pairwise comparisons include: (i) run_model versus 
   run_metareg, and (ii) run_model versus run_series_meta).
 * Function series_meta_plot:
   - presents the extent of heterogeneity in the forest plot of tau using 
   different colours for low, reasonable, fairly high, and fairly extreme tau 
   (the classification has been suggested by Spiegelhalter et al., 2004; 
   ISBN 0471499757).
