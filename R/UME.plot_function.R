#' End-user-ready results: consistency model versus unrelated mean effects model
#'
#' @description This function hosts a toolkit of functions that facilitates the comparison of the consistency model (via \code{run.model}) with the unrelated mean effects model (via \code{run.UME}) regarding the posterior summaries of the summary effect size for
#'   the pairwise comparisons observed in the network, the between-trial standard deviation (\eqn{\tau}) and model assessment parameters.
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param ume An object of S3 class \code{\link{run.UME}}. See 'Value' in \code{\link{run.UME}}.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#' @param save.xls Logical to indicate whether to export the tabulated results to an Excel 'xlsx' format (via the \code{\link[writexl]{write_xlsx}} function) to the working directory of the user.
#'   The default is \code{FALSE} (do not export to an Excel format).
#'
#' @return \code{UME.plot} prints on the R console a message on the most parsimonious model (if any) based on the deviance information criterion (DIC; in red text).
#' Then, the function returns the following list of elements:
#' \tabular{ll}{
#'  \code{Table.effect.size} \tab The posterior mean, posterior standard deviation, and 95\% credible interval of the summary effect size for each pairwise comparison observed in the network under the consistency model and the unrelated mean effects model.\cr
#'  \tab \cr
#'  \code{Table.model.assessment} \tab The DIC, number of effective parameters, and total residual deviance under the consistency model and the unrelated mean effects model (Spiegelhalter et al. (2002)).\cr
#'  \tab \cr
#'  \code{Table.tau} \tab The posterior median and 95\% credible interval of \eqn{\tau} under the consistency model and the unrelated mean effects model.
#'   When a fixed-effect model has been performed, \code{UME.plot} does not return this element.\cr
#'  \tab \cr
#'  \code{Scatterplots} \tab The scatterplot and the Bland-Altman plot on the posterior mean deviance contribution of the individual data points under the consistency model and the unrelated mean effects model.
#'   See 'Details' and 'Value' in \code{\link{scatterplots.dev}} and \code{\link{BlandAltman.plot}}.\cr
#'  \tab \cr
#'  \code{Levarage.plots} \tab The leverage plot on the posterior mean of deviance of the individual data points under the consistency model and the unrelated mean effects model, separately.
#'   See 'Details' and 'Value' in \code{\link{leverage.plot}}.\cr
#'  \tab \cr
#'  \code{Intervalplots} \tab A panel of interval plots on the summary effect size under the consistency model and the unrelated mean effects model for each pairwise comparison observed in the network.
#'   See 'Details' and 'Value' in \code{\link{intervalplot.panel.UME}}.\cr
#' }
#'
#' @details The DIC of the consistency model is compared with the DIC of the unrelated mean effects model (Dias et al. (2013)). If the difference in DIC exceeds 5, the unrelated mean effects model is preferred; if the difference in DIC is less than -5,
#'   the consistency is preferred; otherwise, there is little to choose between the compared models.
#'
#'   Furthermore, \code{UME.plot} exports \code{EM.both.models} to an Excel 'xlsx' format (via the \code{\link[writexl]{write_xlsx}} function) to the working directory of the user.
#'
#'   \code{UME.plot} can be used only for a network of interventions. In the case of two interventions, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{run.UME}}, \code{\link{BlandAltman.plot}}, \code{link{BlandAltman.plot}}, \code{\link{leverage.plot}}, \code{\link{intervalplot.panel.UME}}
#'
#' @references
#' Spineli LM. A novel framework to evaluate the consistency assumption globally in a network of interventions. \emph{submitted} 2021.
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of primary analysis results: A case study on missing outcome data in pairwise and network meta-analysis. \emph{Res Synth Methods} 2021;\bold{12}(4):475--490. [\doi{10.1002/jrsm.1478}]
#'
#' Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence synthesis for decision making 4: inconsistency in networks of evidence based on randomized controlled trials. \emph{Med Decis Making} 2013a;\bold{33}(5):641--56. [\doi{10.1177/0272989X12455847}]
#'
#' Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A. Bayesian measures of model complexity and fit. \emph{J R Stat Soc B} 2002;\bold{64}:583--616. [\doi{10.1111/1467-9868.00353}]
#'
#' @examples
#' data("nma.liu2013")
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis (consistency model)
#' res <- run.model(data = nma.liu2013,
#'                  measure = "OR",
#'                  model = "RE",
#'                  assumption = "IDE-ARM",
#'                  heter.prior = list("halfnormal", 0, 1),
#'                  mean.misspar = c(0, 0),
#'                  var.misspar = 1,
#'                  D = 1,
#'                  n.chains = 3,
#'                  n.iter = 10000,
#'                  n.burnin = 1000,
#'                  n.thin = 1)
#'
#' # Run random-effects network meta-analysis with node-splitting approach
#' ume <- run.UME(full = res, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("placebo", "pramipexole", "serotonin–norepinephrine reuptake inhibitor",
#'                   "serotonin reuptake inhibitor", "tricyclic antidepressant", "pergolide")
#'
#' # Plot the results from the consistency model and the node-splitting approach
#' UME.plot(full = res, ume = ume, drug.names = interv.names)
#' }
#'
#' @export
UME.plot <- function(full, ume, drug.names, save.xls) {


  save.xls <- if (missing(save.xls)) {
    FALSE
  } else {
    save.xls
  }


  model <- full$model
  measure <- full$measure


  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug.names
  }


  if(length(drug.names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }

  # Posterior results on the effect estimates under consistency model
  EM.full <- full$EM[, c(1:3, 7)]
  EM.full[, c(1, 3, 4)] <- if (is.element(measure, c("OR", "ROM"))) {
    exp(EM.full[, c(1, 3, 4)])
  } else {
    EM.full[, c(1, 3, 4)]
  }

  # Posterior results on the effect estimates under UME model
  EM.ume <- ume$EM[, c(1:3, 7)]
  EM.ume[, c(1, 3, 4)] <- if (is.element(measure, c("OR", "ROM"))) {
    exp(EM.ume[, c(1, 3, 4)])
  } else {
    EM.ume[, c(1, 3, 4)]
  }

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
  possible.comp <- possible.observed.comparisons(drug.names, obs.comp)


  ## Keep only the effect estimates according to the 'poss.pair.comp.clean' - Consistency model
  EM.full.clean <- format(round(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), ], 2), nsmall = 2)


  ## Keep only the effect estimates according to the 'poss.pair.comp.clean' - UME model
  EM.ume.clean <- format(round(EM.ume[, ], 2), nsmall = 2)


  ## Keep only the 95% credible intervals (CrI) according to the 'poss.pair.comp.clean' - Consistency model
  CrI.full.clean <- if (is.element(measure, c("OR", "ROM"))) {
    paste0("(", EM.full.clean[, 3], ",", " ", EM.full.clean[, 4], ")", ifelse(as.numeric(EM.full.clean[, 3]) > 1 | as.numeric(EM.full.clean[, 4]) < 1, "*", " "))
  } else {
    paste0("(", EM.full.clean[, 3], ",", " ", EM.full.clean[, 4], ")", ifelse(as.numeric(EM.full.clean[, 3]) > 0 | as.numeric(EM.full.clean[, 4]) < 0, "*", " "))
  }


  ## Keep only the 95% credible intervals (CrI) according to the 'poss.pair.comp.clean' - UME model
  CrI.ume.clean <- if (is.element(measure, c("OR", "ROM"))) {
    paste0("(", EM.ume.clean[, 3], ",", " ", EM.ume.clean[, 4], ")", ifelse(as.numeric(EM.ume.clean[, 3]) > 1 | as.numeric(EM.ume.clean[, 4]) < 1, "*", " "))
  } else {
    paste0("(", EM.ume.clean[, 3], ",", " ", EM.ume.clean[, 4], ")", ifelse(as.numeric(EM.ume.clean[, 3]) > 0 | as.numeric(EM.ume.clean[, 4]) < 0, "*", " "))
  }


  ## Create a data-frame with effect estimates on both models
  EM.both.models <- data.frame(possible.comp$obs.comp[, 4], EM.full.clean[, 1:2], CrI.full.clean, EM.ume.clean[, 1:2], CrI.ume.clean)
  colnames(EM.both.models) <- c("Comparison", "Mean NMA", "SD NMA", "95% CrI NMA", "Mean UME", "SD UME", "95% CrI UME")
  rownames(EM.both.models) <- NULL


  ## Center the columns of the 'data.frame'
  #EM.both.models <- format(EM.both.models, width = max(sapply(names(EM.both.models), nchar)), justify = "centre")


  ## A data-frame with the measures on model assessment
  model.assessment <- data.frame(t(model.assess.full), t(model.assess.ume))
  colnames(model.assessment) <- c("Full NMA", "UME model")
  message(ifelse(model.assess.full[1] - model.assess.ume[1] > 5, "The UME model may be preferred when accounting for model fit and complexity",
                 ifelse(model.assess.full[1] - model.assess.ume[1] < -5, "The consistency model may be preferred when accounting for model fit and complexity",
                        "There is little to choose between the two models")))


  ## A data-frame with the posterior median and 95% CrI on between-trial standard deviation
  if (model == "RE") {
    between.trial.SD <- rbind(tau.full[c(5, 3, 7)], tau.ume[c(5, 3, 7)])
    colnames(between.trial.SD) <- c("Median", "Lower 95% CrI", "Upper 95% CrI")
    rownames(between.trial.SD) <- c("Full NMA", "UME model")
  } else {
    between.trial.SD <- NA
  }


  ## Scatterplot on the deviance contribution of consistency versus UME models
  scatterplot.o <- scatterplots.dev(dev.o.full[, 1], dev.o.ume[, 1], colour = "#D55E00")


  ## Bland-Altman plot on the deviance contribution of consistency versus UME models
  BA.observed <- BlandAltman.plot(dev.o.full[, 1], dev.o.ume[, 1], colour = "#D55E00")


  ## Bring together all four plots
  scatterplots <- ggpubr::ggarrange(scatterplot.o, BA.observed, nrow = 1, ncol = 2)


  ## Leverage plots
  # Consistency model for observed outcomes
  lever.full.o <- leverage.plot(full, drug.names, title = "Consistency model")

    # UME model for observed outcomes
  lever.ume.o <- leverage.plot(ume, drug.names, title = "Unrelated mean effects model")



  ## Bring together the leverage plots for observed outcome
  lev.plots <- ggpubr::ggarrange(lever.full.o, lever.ume.o, nrow = 1, ncol = 2, labels = c("A)",  "B)"))


  intervalplots <- intervalplot.panel.UME(full, ume, drug.names)


  ## Write the table with the EMs from both models as .xlsx
  if (save.xls == TRUE) {
    write_xlsx(EM.both.models, paste0("Table NMA vs UME", ".xlsx"))
    write_xlsx(model.assessment, paste0("Table Model Assessment_NMA vs UME", ".xlsx"))
  }


  ## Return results
  results <- if (model == "RE") {
    list(Table.effect.size = knitr::kable(EM.both.models),
         Table.model.assessment = knitr::kable(model.assessment),
         Table.tau = knitr::kable(between.trial.SD),
         Scatterplots = scatterplots,
         Levarage.plots = lev.plots,
         Intervalplots = intervalplots)
  } else {
    list(Table.effect.size = knitr::kable(EM.both.models),
         Table.model.assessment = knitr::kable(model.assessment),
         Scatterplots = scatterplots,
         Levarage.plots = lev.plots,
         Intervalplots = intervalplots)
  }

  return(results)
}
