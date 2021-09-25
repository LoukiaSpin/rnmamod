#' Forest-plot of comparisons with the selected intervention
#'
#' @description This function illustrates a forest plot of the posterior mean and 95\% credible and predictive interval of comparisons with the selected intervention of the network.
#'
#' @param full An object of S3 class \code{\link{run.model}} or \code{\link{run.metareg}}. See 'Value' in \code{\link{run.model}} and \code{\link{run.metareg}}.
#' @param compar A character to indicate the comparator intervention. it must be any name found in \code{drug.names}.
#' @param cov.value A vector of two elements in the following order: a number that corresponds to a value of the covariate considered in \code{\link{run.metareg}},
#'   and a character object to indicate the name of the covariate.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @return A panel of two forest plots: (1) a forest plot on the effect estimates and predictions of comparisons with the selected intervention of the network, and
#' (2) a forest plot on the posterior mean and 95\% credible interval of SUCRA values of the interventions (Salanti et al., 2011).
#'
#' @details The x-axis in the forest plot of effect sizes displays all interventions in the network; the selected intervention that comprises the \code{compar} is indicated in the plot with a homonymous label.
#'   For each comparison with the selected intervention, the 95\% credible and predictive intervals are displayed as overlapping lines with different colours. When the between-trial variance is very low, these two intervals become indiscernible.
#'   Furthermore, the corresponding numerical results are displayed above each line: 95\% credible intervals are found in parentheses, and 95\% predictive intervals are found in brackets.
#'   Odds ratio and ratio of means are reported in the original scale after exponentiation of the logarithmic scale.
#'
#'   The interventions are sorted in the descending order of their SUCRA values.
#'
#'   \code{forestplot} can be used only for a network of interventions. In the case of two interventions, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{run.metareg}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: an overview and tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71. [\doi{10.1016/j.jclinepi.2010.03.016}]
#'
#' @examples
#' data("nma.liu2013")
#'
#' # Show the first six trials of the dataset (one-trial-per-row format)
#' head(nma.liu2013)
#' #            study t1 t2 t3 r1 r2 r3 m1 m2 m3  n1  n2 n3
#' #    Richard, 2012  1  3  4 15 16 23  6  8  4  39  42 34
#' #     Barone, 2010  1  2 NA 27 38 NA 19 20 NA 152 144 NA
#' # Weinbtraub, 2010  1  3 NA  2  5 NA  6  6 NA  27  28 NA
#' #      Menza, 2009  1  4  5  4  2  9  6  7  5  17  18 17
#' #      Devos, 2008  1  4  5  4  8 11  0  2  1  16  15 17
#' #   Antonini, 2006  4  5 NA 10  8 NA  4  4 NA  16  15 NA
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
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
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("placebo", "pramipexole", "serotonin–norepinephrine reuptake inhibitor",
#'                   "serotonin reuptake inhibitor", "tricyclic antidepressant", "pergolide")
#'
#' # Create the league heatmap
#' forestplot(full = res, compar = "placebo", drug.names = interv.names)
#' }
#'
#' @export
forestplot.metareg <- function(full, reg, compar, cov.value = NULL, drug.names) {


  options(warn = -1)

  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug.names
  }

  compar <- if(missing(compar)) {
    stop("The argument 'compar' has not been defined", call. = F)
  } else if(!is.element(compar, drug.names)) {
    stop("The value of 'compar' is not found in the 'drug.names'", call. = F)
  } else if(is.element(compar, drug.names)) {
    compar
  }

  cov.value <- if (!is.null(reg$beta.all) & missing(cov.value)) {
    stop("The argument 'cov.value' has not been defined", call. = F)
  } else if (!is.null(reg$beta.all) & length(cov.value) < 2) {
    stop("The argument 'cov.value' must be a vector with elements a number and a character", call. = F)
  } else if (!is.null(reg$beta.all) & length(cov.value) == 2) {
    cov.value
  }

  measure <- effect.measure.name(full$measure)
  model <- full$model
  sucra <- full$SUCRA
  cov.val <- ifelse(length(unique(reg$covariate)) < 3, as.numeric(cov.value[1]), as.numeric(cov.value[1]) - mean(reg$covariate))


  ## A matrix with all possible comparisons in the network
  poss.pair.comp1 <- data.frame(exp = t(combn(drug.names, 2))[, 2], comp = t(combn(drug.names, 2))[, 1])
  poss.pair.comp2 <- data.frame(exp = t(combn(drug.names, 2))[, 1], comp = t(combn(drug.names, 2))[, 2])
  poss.pair.comp <- rbind(poss.pair.comp1, poss.pair.comp2)


  ## Effect size of all possible pairwise comparisons (NMA)
  EM.ref00.nma <- cbind(rbind(data.frame(mean = full$EM[, 1], lower = full$EM[, 3], upper = full$EM[, 7]),
                              data.frame(mean = full$EM[, 1]*(-1), lower = full$EM[, 7]*(-1), upper = full$EM[, 3]*(-1))),
                        poss.pair.comp)

  EM.subset.nma <- subset(EM.ref00.nma, EM.ref00.nma[5] == compar)

  EM.ref0.nma <- rbind(EM.subset.nma[, 1:3], c(rep(NA, 3)))

  sucra.new <- data.frame(sucra[, 1], drug.names)[order(match(data.frame(sucra[, 1], drug.names)[, 2], EM.subset.nma[, 4])), 1]

  EM.ref.nma <- EM.ref0.nma[order(sucra.new, decreasing = T), ]

  ## Effect size of all possible pairwise comparisons (NMR)
  EM.ref00.nmr <- cbind(rbind(data.frame(mean = reg$EM[, 1], lower = reg$EM[, 3], upper = reg$EM[, 7]) +
                             (data.frame(mean = reg$beta.all[, 1], lower = reg$beta.all[, 3], upper = reg$beta.all[, 7])*cov.val),
                              data.frame(mean = reg$EM[, 1]*(-1), lower = reg$EM[, 7]*(-1), upper = reg$EM[, 3]*(-1)) +
                             (data.frame(mean = reg$beta.all[, 1]*(-1), lower = reg$beta.all[, 7]*(-1), upper = reg$beta.all[, 3]*(-1))*cov.val)),
                        poss.pair.comp)

  EM.subset.nmr <- subset(EM.ref00.nmr, EM.ref00.nmr[5] == compar)

  EM.ref0.nmr <- rbind(EM.subset.nmr[, 1:3], c(rep(NA, 3)))

  EM.ref.nmr <- EM.ref0.nmr[order(sucra.new, decreasing = T), ]
  rownames(EM.ref.nma) <- rownames(EM.ref.nmr) <- NULL


  # Posterior results on the predicted estimates of comparisons with the selected comparator as reference
  if (model == "RE") {
    pred.ref00.nma <- cbind(rbind(data.frame(mean = full$EM.pred[, 1], lower = full$EM.pred[, 3], upper = full$EM.pred[, 7]),
                                  data.frame(mean = full$EM.pred[, 1]*(-1), lower = full$EM.pred[, 7]*(-1), upper = full$EM.pred[, 3]*(-1))),
                            poss.pair.comp)

    pred.subset.nma <- subset(pred.ref00.nma, pred.ref00.nma[5] == compar)

    pred.ref0.nma <- rbind(pred.subset.nma[, 1:3], c(rep(NA, 3)))

    pred.ref00.nmr <- cbind(rbind(data.frame(mean = reg$EM.pred[, 1], lower = reg$EM.pred[, 3], upper = reg$EM.pred[, 7]) +
                                  (data.frame(mean = reg$beta.all[, 1], lower = reg$beta.all[, 3], upper = reg$beta.all[, 7])*cov.val),
                                  data.frame(mean = reg$EM.pred[, 1]*(-1), lower = reg$EM.pred[, 7]*(-1), upper = reg$EM.pred[, 3]*(-1)) +
                                  (data.frame(mean = reg$beta.all[, 1]*(-1), lower = reg$beta.all[, 7]*(-1), upper = reg$beta.all[, 3]*(-1))*cov.val)),
                            poss.pair.comp)

    pred.subset.nmr <- subset(pred.ref00.nmr, pred.ref00.nmr[5] == compar)

    pred.ref0.nmr <- rbind(pred.subset.nmr[, 1:3], c(rep(NA, 3)))
  } else if (model != "RE") {
    pred.ref00.nma <- NA
    pred.ref00.nmr <- NA
  }


  # Sort by SUCRA in decreasing order and remove the reference intervention (number 1)
  if (model == "RE") {
    pred.ref.nma <- pred.ref0.nma[order(sucra.new, decreasing = T), ]
    pred.ref.nmr <- pred.ref0.nmr[order(sucra.new, decreasing = T), ]
  } else {
    NA
  }
  rownames(pred.ref.nma) <- rownames(pred.ref.nmr) <- NULL


  # Sort the drugs by their SUCRA in decreasing order and remove the reference intervention (number 1)
  drug.names.sorted <- drug.names[order(sucra[, 1], decreasing = T)]


  ## Create a data-frame with credible and predictive intervals of comparisons with the reference intervention
  if (!is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    prepare.EM.nma <- data.frame(as.factor(rep(length(drug.names):1, 2)),
                                 rep(drug.names.sorted, 2),
                                 round(rbind(EM.ref.nma, pred.ref.nma), 2),
                                 rep(c("Credible interval", "Predictive interval"), each = length(drug.names)))
    prepare.EM.nmr <- data.frame(as.factor(rep(length(drug.names):1, 2)),
                                 rep(drug.names.sorted, 2),
                                 round(rbind(EM.ref.nmr, pred.ref.nmr), 2),
                                 rep(c("Credible interval", "Predictive interval"), each = length(drug.names)))
    colnames(prepare.EM.nma) <- colnames(prepare.EM.nmr) <- c("order", "comparison", "mean", "lower", "upper", "interval")
  } else if (is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    prepare.EM.nma <- data.frame(as.factor(rep(length(drug.names):1, 2)),
                                 rep(drug.names.sorted, 2),
                                 round(rbind(exp(EM.ref.nma), exp(pred.ref.nma)), 2),
                                 rep(c("Credible interval", "Predictive interval"), each = length(drug.names)))
    prepare.EM.nmr <- data.frame(as.factor(rep(length(drug.names):1, 2)),
                                 rep(drug.names.sorted, 2),
                                 round(rbind(exp(EM.ref.nmr), exp(pred.ref.nmr)), 2),
                                 rep(c("Credible interval", "Predictive interval"), each = length(drug.names)))
    colnames(prepare.EM.nma) <- colnames(prepare.EM.nmr) <- c("order", "comparison", "mean", "lower", "upper", "interval")
  } else if (!is.element(measure, c("Odds ratio", "Ratio of means")) & model == "FE") {
    prepare.EM.nma <- data.frame(as.factor(length(drug.names):1),
                                 drug.names.sorted,
                                 round(EM.ref.nma, 2))
    prepare.EM.nmr <- data.frame(as.factor(length(drug.names):1),
                                 drug.names.sorted,
                                 round(EM.ref.nmr, 2))
    colnames(prepare.EM.nma) <- colnames(prepare.EM.nmr) <- c("order", "comparison", "mean", "lower", "upper")
  } else if (is.element(measure, c("Odds ratio", "Ratio of means")) & model == "RE") {
    prepare.EM.nma <- data.frame(as.factor(length(drug.names):1),
                                 drug.names.sorted,
                                 round(exp(EM.ref.nma), 2))
    prepare.EM.nmr <- data.frame(as.factor(length(drug.names):1),
                                 drug.names.sorted,
                                 round(exp(EM.ref.nmr), 2))
    colnames(prepare.EM.nma) <- colnames(prepare.EM.nmr) <- c("order", "comparison", "mean", "lower", "upper")
  }

  prepare.EM <- cbind(rbind(prepare.EM.nma, prepare.EM.nmr), analysis = rep(c("NMA", "NMR"), each = length(prepare.EM.nma[, 1])))


  ## Forest plots on credible/predictive intervals of comparisons with the selected comparator


  return(Forest.plots = forest.plots)
}
