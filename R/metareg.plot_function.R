#' End-user-ready results: network meta-analysis versus network meta-regression analysis
#'
#' @description XX
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param reg An object of S3 class \code{\link{run.metareg}}. See 'Value' in \code{\link{run.metareg}}.
#' @param compar A character to indicate the comparator intervention. It must be any name found in \code{drug.names}.
#' @param cov.value A vector of two elements in the following order: a number that corresponds to a value of the covariate considered in \code{\link{run.metareg}},
#'   and a character object to indicate the name of the covariate.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#' @param save.xls Logical to indicate whether to export the tabulated results to an Excel 'xlsx' format (via the \code{\link[writexl]{write_xlsx}} function) to the working directory of the user.
#'   The default is \code{FALSE} (do not export to an Excel format).
#'
#' @return \code{metareg.plot} prints on the R console a message on the most parsimonious model (if any) based on the deviance information criterion (DIC; in red text).
#' Furthermore, the function returns the following list of elements:
#' \tabular{ll}{
#'  \code{Table.estimates} \tab The posterior mean, and 95\% credible interval of the summary effect size for each comparison with the selected intervention under network meta-analysis and meta-regression.\cr
#'  \tab \cr
#'  \code{Table.predictions} \tab The posterior mean, and 95\% predictive interval of the summary effect size for each comparison with the selected intervention under network meta-analysis and meta-regression.\cr
#'  \tab \cr
#'  \code{Table.model.assessment} \tab The DIC, total residual deviance, number of effective parameters, and the posterior mean and 95% credible interval of between-trial standard deviation (\eqn{\tau}) under each model (Spiegelhalter et al. (2002)).
#'   When a fixed-effect model has been performed, \code{metareg.plot} does not return results on \eqn{\tau}.\cr
#'  \tab \cr
#'  \code{Table.regression.coeffients} \tab The posterior mean and 95\% credible interval of the regression coefficient(s).\cr
#'  \tab \cr
#'  \code{Interval.plots} \tab The panel of forest-plots on estimated and predicted effect sizes of comparisons with the selected intervention under network meta-analysis and meta-regression.
#'   See 'Details' and 'Value' in \code{\link{forestplot.metareg}}.\cr
#' }
#'
#' @details The DIC of the network meta-analysis model is compared with the DIC of the network meta-regression model. If the difference in DIC exceeds 5, the network meta-regression model is preferred;
#'   if the difference in DIC is less than -5, the network meta-analysis model is preferred; otherwise, there is little to choose between the compared models.
#'
#'   Furthermore, \code{metareg.plot} exports all tabulated results to separate Excel 'xlsx' formats (via the \code{\link[writexl]{write_xlsx}} function) to the working directory of the user.
#'
#'   \code{metareg.plot} can be used only for a network of interventions. In the case of two interventions, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{run.metareg}}, \code{\link{forestplot.metareg}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: an overview and tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71. [\doi{10.1016/j.jclinepi.2010.03.016}]
#'
#' Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A. Bayesian measures of model complexity and fit. \emph{J R Stat Soc B} 2002;\bold{64}:583--616. [\doi{10.1111/1467-9868.00353}]
#'
#' @examples
#' data("nma.baker2009")
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' res <- run.model(data = nma.baker2009,
#'                  measure = "OR",
#'                  model = "RE",
#'                  assumption = "IDE-ARM",
#'                  heter.prior = list("halfnormal", 0, 1),
#'                  mean.misspar = c(0, 0),
#'                  var.misspar = 1,
#'                  D = 0,
#'                  n.chains = 3,
#'                  n.iter = 10000,
#'                  n.burnin = 1000,
#'                  n.thin = 1)
#'
#' # Publicatiom year
#' pub.year <- c(1996, 1998, 1999, 2000, 2000, 2001, rep(2002, 5), 2003, 2003,
#'               rep(2005, 4), 2006, 2006, 2007, 2007)
#'
#' # Perform a random-effects network meta-regression (exchangeable structure)
#' reg <- run.metareg(full = res,
#'                    covariate = pub.year,
#'                    covar.assumption = "exchangeable",
#'                    n.chains = 3,
#'                    n.iter = 10000,
#'                    n.burnin = 1000,
#'                    n.thin = 1)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("placebo", "budesonide", "budesonide plus formoterol",
#'                   "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium")
#'
#' # Plot the results from the network meta-analysis and meta-regression publication year 2000
#' # For comparisons with salmeterol
#' metareg.plot(full = res,
#'              reg = reg,
#'              compar = "salmeterol",
#'              cov.value = c(2000, "publication year"),
#'              drug.names = interv.names)
#' }
#'
#' @export
metareg.plot <- function(full, reg, compar, cov.value, drug.names, save.xls) {


  options(warn = -1)

  if(length(drug.names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }


 if (length(unique(reg$covariate)) < 3 & !is.element(cov.value[1], reg$covariate)) {
   stop("The first element of the argument 'cov.value' is out of the value range of the analysed covariate", call. = F)
 } else if (length(unique(reg$covariate)) > 2 & (cov.value[1] < min(reg$covariate) | cov.value[1] > max(reg$covariate))) {
   stop("The first element of the argument 'cov.value' is out of the value range of the analysed covariate", call. = F)
 }

  save.xls <- if (missing(save.xls)) {
    FALSE
  } else {
    save.xls
  }

  compar <- if(missing(compar)) {
    stop("The argument 'compar' has not been defined", call. = F)
  } else if(!is.element(compar, drug.names)) {
    stop("The value of the argument 'compar' is not found in the 'drug.names'", call. = F)
  } else if(is.element(compar, drug.names)) {
    compar
  }

  cov.value <- if (missing(cov.value)) {
    stop("The argument 'cov.value' has not been defined", call. = F)
  } else if (length(cov.value) < 2) {
    stop("The argument 'cov.value' must be a vector with elements a number and a character", call. = F)
  } else if (length(cov.value) == 2) {
    cov.value
  }

  covariate <- if (length(unique(reg$covariate)) < 3) {
    unique(reg$covariate)
  } else {
    reg$covariate
  }

  cov.val <- ifelse(length(unique(covariate)) < 3, as.numeric(cov.value[1]), as.numeric(cov.value[1]) - mean(covariate))

  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- as.factor(1:length(full$SUCRA[, 1]))
  } else {
    drug.names
  }

  model <- full$model
  measure <- effect.measure.name(full$measure)

  # Posterior results on the SUCRA value under NMA
  sucra.full <- round(full$SUCRA, 2)
  sucra.full.order <- round(full$SUCRA, 2)[order(sucra.full[, 1], decreasing = T), ]

  # Sort the drugs by their NMA-SUCRA in decreasing order
  drug.names.sorted <- drug.names[order(sucra.full[, 1], decreasing = T)]

  ## A matrix with all possible comparisons in the network
  poss.pair.comp1 <- data.frame(exp = t(combn(drug.names, 2))[, 2], comp = t(combn(drug.names, 2))[, 1])
  poss.pair.comp2 <- data.frame(exp = t(combn(drug.names, 2))[, 1], comp = t(combn(drug.names, 2))[, 2])
  poss.pair.comp <- rbind(poss.pair.comp1, poss.pair.comp2)

  ## Posterior results on effect size for comparisons with the selected intervention (NMA)
  EM.ref00.nma <- cbind(rbind(data.frame(mean = full$EM[, 1], lower = full$EM[, 3], upper = full$EM[, 7]),
                              data.frame(mean = full$EM[, 1]*(-1), lower = full$EM[, 7]*(-1), upper = full$EM[, 3]*(-1))),
                        poss.pair.comp)

  EM.subset.nma <- subset(EM.ref00.nma, EM.ref00.nma[5] == compar)
  EM.ref0.nma <- rbind(EM.subset.nma[, 1:3], c(rep(NA, 3)))
  sucra.full.new <- data.frame(sucra.full[, 1], drug.names)[order(match(data.frame(sucra.full[, 1], drug.names)[, 2], EM.subset.nma[, 4])), 1]
  EM.ref.nma <- EM.ref0.nma[order(sucra.full.new, decreasing = T), ]


  # Posterior mean of regression coefficients for all unique pairwise comparisons
  if (is.element(reg$covar.assumption, c("exchangeable", "independent"))) {
    beta00 <- cbind(rbind(data.frame(mean = reg$beta.all[, 1], lower = reg$beta.all[, 3], upper = reg$beta.all[, 7]),
                          data.frame(mean = reg$beta.all[, 1]*(-1), lower = reg$beta.all[, 7]*(-1), upper = reg$beta.all[, 3]*(-1))),
                    poss.pair.comp)

    beta.all.subset <- subset(beta00, beta00[5] == compar)

    beta0 <- rbind(beta.all.subset[, 1:3], c(rep(NA, 3)))

    beta <- beta0[order(sucra.full.new, decreasing = T), ]
    rownames(beta) <- NULL
  } else {
    beta <- reg$beta[1, c(1, 3, 7)]
  }


  ## Effect size of all possible pairwise comparisons (NMR)
  par.mean <- as.vector(c(reg$EM[, 1] + reg$beta.all[, 1]*cov.val,
                          (reg$EM[, 1]*(-1)) + (reg$beta.all[, 1]*(-1)*cov.val)))
  par.sd <- as.vector(c(sqrt(((reg$EM[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2)),
                        sqrt(((reg$EM[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2))))
  #if (is.element(reg$covar.assumption, c("exchangeable", "independent"))) {
  #  par.mean <- as.vector(c(reg$EM[, 1] + reg$beta.all[, 1]*cov.val,
  #                          (reg$EM[, 1]*(-1)) + (reg$beta.all[, 1]*(-1)*cov.val)))
  #  par.sd <- as.vector(c(sqrt(((reg$EM[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2)),
  #                        sqrt(((reg$EM[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2))))
  #} else {
  #  par.mean <- as.vector(c(reg$EM[, 1], reg$EM[, 1]*(-1)))
  #  par.sd <- as.vector(c(reg$EM[, 2], reg$EM[, 2]))
  #
  #  # Correcting for comparisons with the reference intervention of the network
  #  par.mean[1:(length(drug.names) - 1)] <- as.vector(c(reg$EM[1:(length(drug.names) - 1), 1] + reg$beta[1]*cov.val,
  #                                                      (reg$EM[1:(length(drug.names) - 1), 1]*(-1)) + (reg$beta[1]*(-1)*cov.val)))
  #  par.sd[1:(length(drug.names) - 1)] <- as.vector(c(sqrt(((reg$EM[1:(length(drug.names) - 1), 2])^2) + ((reg$beta[2]*cov.val)^2)),
  #                                                    sqrt(((reg$EM[1:(length(drug.names) - 1), 2])^2) + ((reg$beta[2]*cov.val)^2))))
  #}

  EM.ref00.nmr <- cbind(mean = par.mean, lower = par.mean - 1.96*par.sd, upper = par.mean + 1.96*par.sd,
                        poss.pair.comp)
  EM.subset.nmr <- subset(EM.ref00.nmr, EM.ref00.nmr[5] == compar)
  EM.ref0.nmr <- rbind(EM.subset.nmr[, 1:3], c(rep(NA, 3)))
  EM.ref.nmr <- EM.ref0.nmr[order(sucra.full.new, decreasing = T), ]
  rownames(EM.ref.nma) <- rownames(EM.ref.nmr) <- NULL


  # Posterior results on the predicted estimates of comparisons with the selected comparator as reference
  if (model == "RE") {
    pred.ref00.nma <- cbind(rbind(data.frame(mean = full$EM.pred[, 1], lower = full$EM.pred[, 3], upper = full$EM.pred[, 7]),
                                  data.frame(mean = full$EM.pred[, 1]*(-1), lower = full$EM.pred[, 7]*(-1), upper = full$EM.pred[, 3]*(-1))),
                            poss.pair.comp)
    pred.subset.nma <- subset(pred.ref00.nma, pred.ref00.nma[5] == compar)
    pred.ref0.nma <- rbind(pred.subset.nma[, 1:3], c(rep(NA, 3)))

    par.mean <- as.vector(c(reg$EM.pred[, 1] + reg$beta.all[, 1]*cov.val,
                            (reg$EM.pred[, 1]*(-1)) + (reg$beta.all[, 1]*(-1)*cov.val)))
    par.sd <- as.vector(c(sqrt(((reg$EM.pred[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2)),
                          sqrt(((reg$EM.pred[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2))))

    #if (is.element(reg$covar.assumption, c("exchangeable", "independent"))) {
    #  par.mean <- as.vector(c(reg$EM.pred[, 1] + reg$beta.all[, 1]*cov.val,
    #                          (reg$EM.pred[, 1]*(-1)) + (reg$beta.all[, 1]*(-1)*cov.val)))
    #  par.sd <- as.vector(c(sqrt(((reg$EM.pred[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2)),
    #                        sqrt(((reg$EM.pred[, 2])^2) + ((reg$beta.all[, 2]*cov.val)^2))))
    #} else {
    #  par.mean <- as.vector(c(reg$EM.pred[, 1], reg$EM.pred[, 1]*(-1)))
    #  par.sd <- as.vector(c(reg$EM.pred[, 2], reg$EM.pred[, 2]))
    #
    #  # Correcting for comparisons with the reference intervention of the network
    #  par.mean[1:(length(drug.names) - 1)] <- as.vector(c(reg$EM[1:(length(drug.names) - 1), 1] + reg$beta[1]*cov.val,
    #                                                      (reg$EM[1:(length(drug.names) - 1), 1]*(-1)) + (reg$beta[1]*(-1)*cov.val)))
    #  par.sd[1:(length(drug.names) - 1)] <- as.vector(c(sqrt(((reg$EM[1:(length(drug.names) - 1), 2])^2) + ((reg$beta[2]*cov.val)^2)),
    #                                                    sqrt(((reg$EM[1:(length(drug.names) - 1), 2])^2) + ((reg$beta[2]*cov.val)^2))))
    #}

    pred.ref00.nmr <-  cbind(data.frame(mean = par.mean, lower = par.mean - 1.96*par.sd, upper = par.mean + 1.96*par.sd),
                             poss.pair.comp)
    pred.subset.nmr <- subset(pred.ref00.nmr, pred.ref00.nmr[5] == compar)
    pred.ref0.nmr <- rbind(pred.subset.nmr[, 1:3], c(rep(NA, 3)))
  } else if (model != "RE") {
    pred.ref00.nma <- NA
    pred.ref00.nmr <- NA
  }


  # Sort by SUCRA in decreasing order and remove the reference intervention (number 1)
  if (model == "RE") {
    pred.ref.nma <- pred.ref0.nma[order(sucra.full.new, decreasing = T), ]
    pred.ref.nmr <- pred.ref0.nmr[order(sucra.full.new, decreasing = T), ]
  } else {
    NA
  }
  rownames(pred.ref.nma) <- rownames(pred.ref.nmr) <- NULL


  if (!is.element(measure, c("Odds ratio", "Ratio of means"))) {
    EM.ref.nma <- EM.ref.nma
    EM.ref.nmr <- EM.ref.nmr
    pred.ref.nma <- pred.ref.nma
    pred.ref.nmr <- pred.ref.nmr
    beta <- beta
  } else {
    EM.ref.nma <- exp(EM.ref.nma)
    EM.ref.nmr <- exp(EM.ref.nmr)
    pred.ref.nma <- exp(pred.ref.nma)
    pred.ref.nmr <- exp(pred.ref.nmr)
    beta <- exp(beta)
  }


  # Posterior results on between-trial standard deviation under NMA
  tau.full <- if (model == "RE") {
    round(full$tau, 2)
  } else {
    NA
  }

  # Posterior results on between-trial standard deviation under meta-regression
  tau.meta <- if (model == "RE") {
    round(reg$tau, 2)
  } else {
    NA
  }

  # Posterior mean of model assessment measures under NMA
  model.assess.NMA <- round(full$model.assessment, 2)

  # Posterior mean of model assessment measures under meta-regression
  model.assess.meta <- round(reg$model.assessment, 2)


  ## The 95% CrIs of the between-trial standard deviation under NMA and meta-regression
  if (model == "RE") {
    CrI.tau.full <- paste0("(", tau.full[, 3], ",", " ", tau.full[, 7], ")")
    CrI.tau.meta <- paste0("(", tau.meta[, 3], ",", " ", tau.meta[, 7], ")")
  } else {
    CrI.tau.full <- NA
    CrI.tau.meta <- NA
  }


  ## A data-frame on the model assessment results and between-trial standard deviation under NMA & meta-regression
  if (model == "RE") {
    table.model.assess <- data.frame(c("Network meta-analysis", "Meta-regression"),
                                     rbind(model.assess.NMA[c(1, 3, 2)],
                                           model.assess.meta[c(1, 3, 2)]),
                                     rbind(cbind(tau.full[, 5], tau.full[, 2], CrI.tau.full),
                                           cbind(tau.meta[, 5], tau.meta[, 2], CrI.tau.meta)))
    colnames(table.model.assess) <- c("Analysis", "DIC", "Mean deviance", "pD", "Median tau", "SD tau", "95% CrI tau")
  } else {
    table.model.assess <- data.frame(c("Network meta-analysis", "Meta-regression"),
                                     rbind(model.assess.NMA[c(1, 3, 2)],
                                           model.assess.meta[c(1, 3, 2)]))
    colnames(table.model.assess) <- c("Analysis", "DIC", "Mean deviance", "pD")
  }
  message(ifelse(model.assess.NMA[1] - model.assess.meta[1] > 5, "The network meta-regression model may be preferred when accounting for model fit and complexity",
                 ifelse(model.assess.NMA[1] - model.assess.meta[1] < -5, "The network meta-analysis model may be preferred when accounting for model fit and complexity",
                        "There is little to choose between the two models")))


  ## A data-frame with the effect estimates and regression coefficients of reference-comparisons from both analyses (Sort by NMA-SUCRA in decreasing order)
  if (!is.element(measure, c("Odds ratio", "Ratio of means"))) {
    CrI.est.nma <- paste0("(", round(EM.ref.nma[, 2], 2), ",", " ", round(EM.ref.nma[, 3], 2), ")", ifelse(EM.ref.nma[, 2] > 0 | EM.ref.nma[, 3] < 0, "*", " "))
    CrI.est.nmr <- paste0("(", round(EM.ref.nmr[, 2], 2), ",", " ", round(EM.ref.nmr[, 3], 2), ")", ifelse(EM.ref.nmr[, 2] > 0 | EM.ref.nmr[, 3] < 0, "*", " "))
    CrI.pred.nma <- paste0("(", round(pred.ref.nma[, 2], 2), ",", " ", round(pred.ref.nma[, 3], 2), ")", ifelse(pred.ref.nma[, 2] > 0 | pred.ref.nma[, 3] < 0, "*", " "))
    CrI.pred.nmr <- paste0("(", round(pred.ref.nmr[, 2], 2), ",", " ", round(pred.ref.nmr[, 3], 2), ")", ifelse(pred.ref.nmr[, 2] > 0 | pred.ref.nmr[, 3] < 0, "*", " "))
    CrI.beta <- if (is.element(reg$covar.assumption, c("exchangeable", "independent"))) {
      paste0("(", round(beta[, 2], 2), ",", " ", round(beta[, 3], 2), ")", ifelse(beta[, 2] > 0 | beta[, 3] < 0, "*", " "))
    } else {
      paste0("(", round(reg$beta[2], 2), ",", " ", round(reg$beta[3], 2), ")", ifelse(reg$beta[2] > 0 | reg$beta[3] < 0, "*", " "))
    }
  } else {
    CrI.est.nma <- paste0("(", round(EM.ref.nma[, 2], 2), ",", " ", round(EM.ref.nma[, 3], 2), ")", ifelse(EM.ref.nma[, 2] > 1 | EM.ref.nma[, 3] < 1, "*", " "))
    CrI.est.nmr <- paste0("(", round(EM.ref.nmr[, 2], 2), ",", " ", round(EM.ref.nmr[, 3], 2), ")", ifelse(EM.ref.nmr[, 2] > 1 | EM.ref.nmr[, 3] < 1, "*", " "))
    CrI.pred.nma <- paste0("(", round(pred.ref.nma[, 2], 2), ",", " ", round(pred.ref.nma[, 3], 2), ")", ifelse(pred.ref.nma[, 2] > 1 | pred.ref.nma[, 3] < 1, "*", " "))
    CrI.pred.nmr <- paste0("(", round(pred.ref.nmr[, 2], 2), ",", " ", round(pred.ref.nmr[, 3], 2), ")", ifelse(pred.ref.nmr[, 2] > 1 | pred.ref.nmr[, 3] < 1, "*", " "))
    CrI.beta <- if (is.element(reg$covar.assumption, c("exchangeable", "independent"))) {
      paste0("(", round(beta[, 2], 2), ",", " ", round(beta[, 3], 2), ")", ifelse(beta[, 2] > 1 | beta[, 3] < 1, "*", " "))
    } else {
      paste0("(", round(beta[2], 2), ",", " ", round(beta[3], 2), ")", ifelse(beta[2] > 1 | beta[3] < 1, "*", " "))
    }
  }


  # Tabulate results on comparisons with the reference (both models)
  Est.both.models <- na.omit(data.frame(drug.names.sorted,
                                        round(EM.ref.nma[, 1], 2),
                                        CrI.est.nma,
                                        round(EM.ref.nmr[, 1], 2),
                                        CrI.est.nmr))
  Pred.both.models <- na.omit(data.frame(drug.names.sorted,
                                         round(pred.ref.nma[, 1], 2),
                                         CrI.pred.nma,
                                         round(pred.ref.nmr[, 1], 2),
                                         CrI.pred.nmr))
  colnames(Est.both.models) <- colnames(Pred.both.models) <- c(paste("versus", compar), "Mean NMA", "95% CrI NMA", "Mean NMR", "95% CrI NMR")
  rownames(Est.both.models) <- rownames(Pred.both.models) <- NULL


  # Results on the regression coefficient
  if (is.element(reg$covar.assumption, c("exchangeable", "independent")) ) {
    reg.coeff <- na.omit(data.frame(drug.names.sorted, round(beta[, 1], 2), CrI.beta))
    colnames(reg.coeff) <- c(paste("versus", compar), "Mean beta", "95% CrI beta")
  } else {
    reg.coeff <- data.frame(round(beta[1], 2), CrI.beta)
    colnames(reg.coeff) <- c("Mean beta", "95% CrI beta")
  }
  rownames(reg.coeff) <- NULL


  ## Forest plots of reference-comparisons on effect estimate
  forest.plots <- forestplot.metareg(full, reg, compar, cov.value, drug.names)
  SUCRA.scatterplot <- scatteplot.sucra(full, reg, cov.value, drug.names)

  ## Write all tables as .xlsx
  if (save.xls == TRUE) {
    write_xlsx(Est.both.models, paste0("Table NMA vs NMR_Estimation", ".xlsx"))
    write_xlsx(Pred.both.models, paste0("Table NMA vs NMR_Prediction", ".xlsx"))
    write_xlsx(table.model.assess, paste0("Table Model Assessment_NMA vs NMR", ".xlsx"))

    if (is.element(reg$covar.assumption, c("exchangeable", "independent")) ) {
      write_xlsx(reg.coeff, paste0("Table NMA vs NMR_Coefficient", ".xlsx"))
    }
  }


  ## Return results
  return(list(Table.estimates = knitr::kable(Est.both.models),
              Table.predictions = knitr::kable(Pred.both.models),
              Table.model.assessment = knitr::kable(table.model.assess),
              Table.regression.coeffients = knitr::kable(reg.coeff),
              Interval.plots = forest.plots,
              SUCRA.scatterplot = SUCRA.scatterplot))
}
