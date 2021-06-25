#' Plot the results from the unrelated mean effects model
#'
#' @export
UME.plot <- function(full, ume, drug.names, threshold) {


  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug.names
  }


  if(length(drug.names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }


  if (missing(threshold) & is.element(full$measure, "OR")) {
    threshold <- 0.28
    #message("The value 0.28 was assigned on 'threshold' by default")
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The value 0.28 was assigned on 'threshold' by default", "\033[0m", "\n")))
  } else if (missing(threshold) & is.element(full$measure, c("MD", "SMD", "ROM"))) {
    threshold <- 0.17
    #message("The value 0.17 was assigned on 'threshold' by default")
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The value 0.17 was assigned on 'threshold' by default", "\033[0m", "\n")))
  } else {
    threshold <- threshold
    #message(paste("The value", threshold, "was assigned on 'threshold' for", effect.measure.name(full$measure)))
    message(cat(paste0("\033[0;", col = 32, "m", txt = paste("The value", threshold, "was assigned on 'threshold' for", effect.measure.name(full$measure)), "\033[0m", "\n")))
  }

  ## The results on the following parameters will be used:
  # Analysis model
  model <- if (full$model != ume$model) {
    stop("The argument 'model' differs in 'run.model' and 'run.UME'. Specify the same 'model' and run the analysis again")
  } else {
    full$model
  }

  # Effect measure
  measure <- if (full$measure != ume$measure) {
    stop("The argument 'measure' differs in 'run.model' and 'run.UME'. Specify the same 'measure' and run the analysis again")
  }

  # Posterior results on the effect estimates under consistency model
  EM.full <- full$EM

  # Posterior results on the effect estimates under UME model
  EM.ume <- ume$EM

  # Posterior results on between-trial standard deviation under consistency model
  tau.full <- if (model == "RE") {
    full$tau
  } else {
    NA
  }

  # Posterior results on between-trial standard deviation under UME model
  tau.ume <- if (model == "RE") {
    ume$tau
  } else {
    NA
  }

  # Posterior mean on deviance contribution for observed outcomes under consistency model
  dev.o.full <- full$dev.o

  # Posterior mean on deviance contribution for observed outcomes under UME model
  dev.o.ume <- ume$dev.o

  # Measures of model assessment: DIC, pD, total residual deviance (as obtained from R2jags) under consistency model
  model.assess.full <- full$model.assessment

  # Measures of model assessment: DIC, pD, total residual deviance (as obtained from R2jags) under UME model
  model.assess.ume <- ume$model.assessment

  # Observed comparisons in the network
  obs.comp <- ume$obs.comp


  ## Possible and observed comparisons (with names)
  possible.comp <- possible.observed.frail.comparisons(drug.names, obs.comp)


  ## Keep only the effect estimates according to the 'poss.pair.comp.clean' - Consistency model
  EM.full.clean <- format(round(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), c(1:3, 7)], 2), nsmall = 2)



  ## Keep only the effect estimates according to the 'poss.pair.comp.clean' - UME model
  EM.ume.clean <- format(round(EM.ume[, 1:4], 2), nsmall = 2)



  ## Keep only the 95% credible intervals (CrI) according to the 'poss.pair.comp.clean' - Consistency model
  CrI.full.clean <- paste0("(", EM.full.clean[, 3], ",", " ", EM.full.clean[, 4], ")", ifelse(as.numeric(EM.full.clean[, 3]) > 0 | as.numeric(EM.full.clean[, 4]) < 0, "*", " "))



  ## Keep only the 95% credible intervals (CrI) according to the 'poss.pair.comp.clean' - UME model
  CrI.ume.clean <- paste0("(", EM.ume.clean[, 3], ",", " ", EM.ume.clean[, 4], ")", ifelse(as.numeric(EM.ume.clean[, 3]) > 0 | as.numeric(EM.ume.clean[, 4]) < 0, "*", " "))



  ## Create a data-frame with effect estimates on both models
  EM.both.models <- data.frame(possible.comp$obs.comp[, 4], EM.full.clean[, 1:2], CrI.full.clean, EM.ume.clean[, 1:2], CrI.ume.clean)
  colnames(EM.both.models) <- c("Comparison", "Posterior mean NMA", "Posterior SD NMA", "95% CrI NMA", "Posterior mean UME",
                                "Posterior SD UME", "95% CrI UME")
  rownames(EM.both.models) <- NULL



  ## Center the columns of the 'data.frame'
  EM.both.models <- format(EM.both.models, width = max(sapply(names(EM.both.models), nchar)), justify = "centre")



  ## A data-frame with the measures on model assessment
  model.assessment <- data.frame(t(model.assess.full), t(model.assess.ume))
  colnames(model.assessment) <- c("Full NMA", "UME model")
  message(ifelse(model.assess.full[1] - model.assess.ume[1] > 5, "The UME model may be preferred when accounting for model fit and complexity",
                 ifelse(model.assess.full[1] - model.assess.ume[1] < -5, "The consistency model may be preferred when accounting for model fit and complexity",
                        "There is little to choose between the two models")))



  ## A data-frame with the posterior median and 95% CrI on between-trial standard deviation
  if (model == "RE") {
    between.trial.SD <- rbind(tau.full[c(5, 3, 7)], tau.ume[c(5, 3, 7)])
    colnames(between.trial.SD) <- c("Posterior median", "Lower 95% CrI", "Upper 95% CrI")
    rownames(between.trial.SD) <- c("Full NMA", "UME model")
  } else {
    between.trial.SD <- NA
  }



  ## Scatterplot on the deviance contribution of consistency versus UME models
  scatterplot.o <- scatterplots.dev(dev.o.full[, 1], dev.o.ume[, 1], colour = "#D55E00")


  ## Bland-Altman plot on the deviance contribution of consistency versus UME models
  BA.observed <- BlandAltman.plot(dev.o.full[, 1], dev.o.ume[, 1], colour = "#D55E00")


  ## Bring together all four plots
  scatterplots <- ggarrange(scatterplot.o, BA.observed, nrow = 1, ncol = 2)


  ## Leverage plots
  # Consistency model for observed outcomes
  lever.full.o <- leverage.plot(full, drug.names, title.o = "Observed outcomes under consistency model")

    # UME model for observed outcomes
  lever.ume.o <- leverage.plot(ume, drug.names, title.o = "Observed outcomes under UME model")



  ## Bring together the leverage plots for observed outcome
  lev.plots <- ggarrange(lever.full.o, lever.ume.o, nrow = 1, ncol = 2, labels = c("A)",  "B)"))


  forestplots <- forestplot.panel.UME(full, ume, drug.names)


  heatmap <- heatmap.similarity.UME(full, ume, drug.names, threshold)


  ## Write the table with the EMs from both models as .xlsx
  write_xlsx(EM.both.models, paste0(getwd(),"Table NMA vs UME.xlsx"))


  ## Return results
  results <- if (model == "RE") {
    list(EM.both.models = EM.both.models,
         model.assessment = model.assessment,
         between.trial.SD = between.trial.SD,
         scatterplots = scatterplots,
         levarage.plots = lev.plots,
         forestplots.panel = forestplots,
         heatmap.similarity = heatmap,
         threshold = threshold)
  } else {
    list(EM.both.models = EM.both.models,
         model.assessment = model.assessment,
         scatterplots = scatterplots,
         levarage.plots = lev.plots,
         forestplots.panel = forestplots,
         heatmap.similarity = heatmap,
         threshold = threshold)
  }

  return(results)
}
