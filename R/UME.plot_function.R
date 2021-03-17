#' Plot the results from the unrelated mean effects model
#'
#' @export
UME.plot <- function(full, ume, drug.names) {


  ## The results on the following parameters will be used:
  # Posterior results on the effect estimates under consistency model
  EM.full <- full$EM

  # Posterior results on the effect estimates under UME model
  EM.ume <- ume$EM

  # Posterior results on between-trial standard deviation under consistency model
  tau.full <- full$tau

  # Posterior results on between-trial standard deviation under UME model
  tau.ume <- ume$tau

  # Posterior mean on deviance contribution for missing outcomes under consistency model
  dev.m.full <- full$dev.m

  # Posterior mean on deviance contribution for missing outcomes under UME model
  dev.m.ume <- ume$dev.m

  # Posterior mean on deviance contribution for observed outcomes under consistency model
  dev.o.full <- full$dev.o

  # Posterior mean on deviance contribution for observed outcomes under UME model
  dev.o.ume <- ume$dev.o

  # Measures of model assessment: DIC, pD, total residual deviance (as obtained from R2jags) under consistency model
  model.assess.full <- full$model.assessment

  # Measures of model assessment: DIC, pD, total residual deviance (as obtained from R2jags) under UME model
  model.assess.ume <- ume$model.assessment



  ## Obtain all unique pairwise comparisons using the 'combn' functions
  poss.pair.comp <- data.frame(t(combn(1:length(drug.names), 2))[, 2], t(combn(1:length(drug.names), 2))[, 1])
  colnames(poss.pair.comp) <- c("treat1", "treat2")


  poss.pair.comp$comp <- paste(poss.pair.comp[, 1], "vs", poss.pair.comp[, 2])
  poss.pair.comp.clean <- poss.pair.comp[is.element(poss.pair.comp$comp, ume$obs.comp), ]


  ## Replace intervention id with their original name
  # For treat1 (non-baseline arm)
  for(i in sort(unique(unlist(poss.pair.comp.clean[, 1])))) {
    poss.pair.comp.clean[poss.pair.comp.clean$treat1 == i, 1] <- drug.names[i]
  }

  # For treat2 (baseline arm)
  for(i in sort(unique(unlist(poss.pair.comp.clean[, 2])))) {
    poss.pair.comp.clean[poss.pair.comp.clean$treat2 == i, 2] <- drug.names[i]
  }



  ## Create a vector with the comparisons (non-baseline versus beseline) using the 'paste' function
  comparison <- paste(poss.pair.comp.clean[, 1], "vs", poss.pair.comp.clean[, 2])



  ## Keep only the effect estimates according to the 'poss.pair.comp.clean' - Consistency model
  EM.full.clean <- format(round(EM.full[is.element(poss.pair.comp$comp, ume$obs.comp), 1:4], 2), nsmall = 2)



  ## Keep only the effect estimates according to the 'poss.pair.comp.clean' - UME model
  EM.ume.clean <- format(round(EM.ume[, 1:4], 2), nsmall = 2)



  ## Keep only the 95% credible intervals (CrI) according to the 'poss.pair.comp.clean' - Consistency model
  CrI.full.clean <- paste0("c(", EM.full.clean[, 3], ",", " ", EM.full.clean[, 4], ")", ifelse(as.numeric(EM.full.clean[, 3]) > 0 | as.numeric(EM.full.clean[, 4]) < 0, "*", " "))



  ## Keep only the 95% credible intervals (CrI) according to the 'poss.pair.comp.clean' - UME model
  CrI.ume.clean <- paste0("c(", EM.ume.clean[, 3], ",", " ", EM.ume.clean[, 4], ")", ifelse(as.numeric(EM.ume.clean[, 3]) > 0 | as.numeric(EM.ume.clean[, 4]) < 0, "*", " "))



  ## Create a data-frame with effect estimates on both models
  EM.both.models <- data.frame(comparison, EM.full.clean[, 1:2], CrI.full.clean, EM.ume.clean[, 1:2], CrI.ume.clean)
  colnames(EM.both.models) <- c("Comparison", "Posterior mean NMA", "Posterior SD NMA", "95% CrI NMA", "Posterior mean UME",
                                "Posterior SD UME", "95% CrI UME  ")
  rownames(EM.both.models) <- NULL



  ## Center the columns of the 'data.frame'
  EM.both.models <- format(EM.both.models, width = max(sapply(names(EM.both.models), nchar)), justify = "centre")



  ## A data-frame with the measures on model assessment
  model.assessment <- data.frame(t(model.assess.full), t(model.assess.ume))
  colnames(model.assessment) <- c("Full NMA", "UME model")
  message(ifelse(model.assess.full[1] - model.assess.ume[1] > 5, "The UME model may be preferred when accounting for model fit and complexity",
                 ifelse(model.assess.full[1] - model.assess.ume[1] < -5, "The consistency model may be preferred when accounting for model fit and complexity",
                        "There is little to choose between the two models")))



  ## A data-frame with the posterior median and 95% CrII on between-trial standard deviation
  between.trial.SD <- rbind(tau.full[c(1, 3:4)], tau.ume[c(1, 3:4)])
  colnames(between.trial.SD) <- c("Posterior median", "Lower 95% CrI", "Upper 95% CrI")
  rownames(between.trial.SD) <- c("Full NMA", "UME model")



  ## Scatterplot on the deviance contribution of consistency versus UME models
  scatterplots <- scatterplots.dev(full, ume, drug.names)



  ## Leverage plots
  # Consistency model for observed outcomes
  lever.full.o <- leverage.plot(full, drug.names, title.o = "Observed outcomes under consistency model", title.m = "Missing outcome data under consistency model")$leverage.plot.observed

  # Consistency model for missing outcomes
  lever.full.m <- leverage.plot(full, drug.names, title.o = "Observed outcomes under consistency model", title.m = "Missing outcome data under consistency model")$leverage.plot.missing

  # UME model for observed outcomes
  lever.ume.o <- leverage.plot(ume, drug.names, title.o = "Observed outcomes under UME model", title.m = "Missing outcome data under UME model")$leverage.plot.observed

  # UME model for missing outcomes
  lever.ume.m <- leverage.plot(ume, drug.names, title.o = "Observed outcomes under UME model", title.m = "Missing outcome data under UME model")$leverage.plot.missing



  ## Bring together all four leverage plots
  lev.plots <- ggarrange(lever.full.o, lever.ume.o, lever.full.m, lever.ume.m, nrow = 2, ncol = 2, labels = c("A)", "", "B)", ""))



  return(list(EM.both.models = EM.both.models, model.assessment = model.assessment, between.trial.SD = between.trial.SD,
              scatterplots = scatterplots, leverage.plots = lev.plots))

}
