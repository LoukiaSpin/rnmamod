#' End-user-ready results: consistency model versus series of pairwise meta-analyses
#'
#' @description \code{series.meta.plot} facilitates the comparison of the estimated effects and between-trial standard deviation (where applicable) from .
#'   Hence, the user can detect any inconsistencies in the estimated effects from the compared models and explore the gains in their precision by applying the consistency model.
#'   Furthermore, the user can investigate the plausibility of the common between-trial heterogeneity assumption which is typically considered in the consistency model.
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param meta An object of S3 class \code{\link{run.series.meta}}. See 'Value' in \code{\link{run.series.meta}}.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @return \code{series.meta.plot} returns a panel of two forest plots: (1) a forest plot on the posterior mean and 95\% credible interval of the effect of the observed comparisons under the consistency model and the corresponding pairwise meta-analyses,
#'   and (2) a forest plot on the posterior median and 95\% credible interval of the between-trial standard deviation (\eqn{\tau}) for the observed comparisons. The estimated \eqn{\tau} from the consistency model appears as a rectangle in the forest plot.
#'   When a fixed-effect model has been performed, only the forest plot on the estimated effects is shown.
#'
#'   The R console prints the data-frame with the estimated effects of the observed comparisons under both models and the data-frame on the estimated (\eqn{\tau})s from the pairwise meta-analyses.
#'   Furthermore, \code{series.meta.plot} exports both data-frames to an Excel 'xlsx' format (via the \code{\link[writexl]{write_xlsx}} function) to the working directory of the user.
#'
#' @details \code{series.meta.plot} can be used only for a network of interventions. In the case of two interventions, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{run.series.meta}}
#'
#' @examples
#' data("nma.baker2009.RData")
#'
#' # Perform a random-effects network meta-analysis
#' res1 <- run.model(data = nma.baker2009, measure = "OR", model = "RE", assumption = "IDE-ARM", heter.prior = list("halfnormal", 0, 1), mean.misspar = 0, var.misspar = 1, D = 1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' # Run separate random-effects pairwise meta-analyses
#' meta1 <- run.series.meta(full = res1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("budesodine", "budesodine plus formoterol", "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium", "placebo")
#'
#' # Plot the results from both models
#' series.meta.plot(full = res1, meta = meta1, drug.names = interv.names)
#'
#' @export
series.meta.plot <- function(full, meta, drug.names) {


  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug.names
  }


  if (length(drug.names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }

  ## The results on the following parameters will be used:
  # Posterior results on the effect estimates under NMA
  EM.full <- full$EM

  # Posterior results on the effect estimates under separate random-effect pairwise meta-analysis (RE-MAs)
  EM.meta <- meta$EM

  # Posterior results on between-trial standard deviation under NMA
  tau.full <- full$tau

  # Posterior results on between-trial standard deviation under RE-MAs
  tau.meta <- meta$tau

  # Possible and observed comparisons
  possible.comp <- possible.observed.comparisons(drug.names, obs.comp = paste0(meta$EM[, "t2"], "vs", meta$EM[, "t1"]))

  # Observed comparisons
  obs.comp <- possible.comp$obs.comp[, 3]

  model <- full$model

  measure <- effect.measure.name(full$measure)

  ## Keep only the effect estimates according to the 'poss.pair.comp.clean' - Consistency model
  EM.full.clean <- format(round(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), c(1:3, 7)], 2), nsmall = 2)



  ## Effect estimate of separate RE-MAs
  EM.meta.clean <- format(round(EM.meta[, c(3:5, 9)], 2), nsmall = 2)



  ## Between-trial standard deviation of separate RE-MAs
  if (model == "RE") {
    tau.meta.clean <- format(round(tau.meta[, c(3:5, 9)], 2), nsmall = 2)
  } else {
    tau.meta.clean <- NA
  }



  ## Keep only the 95% credible intervals (CrI) according to the 'poss.pair.comp.clean' - Consistency model
  CrI.full.clean <- paste0("(", EM.full.clean[, 3], ",", " ", EM.full.clean[, 4], ")", ifelse(as.numeric(EM.full.clean[, 3]) > 0 | as.numeric(EM.full.clean[, 4]) < 0, "*", " "))



  ## The 95% CrIs of the effect estimate of separate RE-MAs
  CrI.meta.clean <- paste0("(", EM.meta.clean[, 3], ",", " ", EM.meta.clean[, 4], ")", ifelse(as.numeric(EM.meta.clean[, 3]) > 0 | as.numeric(EM.meta.clean[, 4]) < 0, "*", " "))



  ## The 95% CrIs of the between-trial standard deviation of separate RE-MAs
  CrI.tau.meta <- if (model == "RE") {
    paste0("(", tau.meta.clean[, 3], ",", " ", tau.meta.clean[, 4], ")")
  } else {
    NA
  }



  ## Create a data-frame with effect estimates on both models
  if (model == "RE") {
    EM.both.models <- data.frame(possible.comp$obs.comp[, 4], EM.full.clean[, 1:2], CrI.full.clean, EM.meta.clean[, 1:2], CrI.meta.clean, tau.meta.clean[, 1:2], CrI.tau.meta)
    colnames(EM.both.models) <- c("Comparison", "Posterior mean NMA", "Posterior SD NMA", "95% CrI NMA", "Posterior mean MA",
                                  "Posterior SD MA", "95% CrI MA", "Posterior median tau", "Posterior SD tau", "95% CrI tau")
  } else {
    EM.both.models <- data.frame(possible.comp$obs.comp[, 4], EM.full.clean[, 1:2], CrI.full.clean, EM.meta.clean[, 1:2], CrI.meta.clean)
    colnames(EM.both.models) <- c("Comparison", "Posterior mean NMA", "Posterior SD NMA", "95% CrI NMA", "Posterior mean MA",
                                  "Posterior SD MA", "95% CrI MA")
  }
  rownames(EM.both.models) <- NULL



  ## Prepare dataset for ggplot2
  # Effect estimate
  prepare <- data.frame(rep(1:length(obs.comp), 2), rep(possible.comp$obs.comp[, 4], 2), rbind(apply(EM.full.clean[, -2], 2, as.numeric), apply(EM.meta.clean[, -2], 2, as.numeric)), rep(c("Network meta-analysis", "Paiwise meta-analysis"), each = length(obs.comp)))
  colnames(prepare) <- c("order", "comparison", "mean", "lower", "upper", "analysis")
  rownames(prepare) <- NULL

  # Between-trial standard deviation
  if (model == "RE") {
    prepare.tau <- data.frame(possible.comp$obs.comp[, 4], tau.meta.clean[, -2])
    colnames(prepare.tau) <- c("comparison", "median", "lower", "upper")
  } else {
    prepare.tau <- NA
  }




  ## Forest plots of comparisons on effect estimate
  p1 <- ggplot(data = prepare, aes(x = as.factor(order), y = mean, ymin = lower, ymax = upper, colour = analysis, group = analysis)) +
          geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
          #geom_hline(yintercept = ifelse(measure != "Odds ratio" & measure != "Ratio of means", 0, 1), lty = 2, size = 1.3, col = "grey53") +
          geom_hline(yintercept = 0, lty = 1, size = 1, col = "grey53") +
          geom_point(size = 1.5,  colour = "black", stroke = 0.3, position = position_dodge(width = 0.5)) +
          geom_text(aes(x = as.factor(order), y = round(as.numeric(mean), 2), label = round(as.numeric(mean), 2)),
                    color = "black", hjust = -0.1, vjust = -0.5, size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.8), inherit.aes = T) +
          labs(x = "", y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), paste(measure, "(in logarithmic scale)"), measure), colour = "Analysis") +
          scale_x_discrete(breaks = as.factor(1:length(obs.comp)), labels = prepare$comparison[1:length(obs.comp)]) +
          #scale_y_continuous(trans = ifelse(measure != "Odds ratio" & measure != "Ratio of means", "identity", "log10")) +
          scale_y_continuous(trans = "identity") +
          scale_color_manual(breaks = c("Network meta-analysis", "Paiwise meta-analysis"), values = c("#009E73", "#D55E00")) +
          coord_flip() +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "none", legend.justification = c(0.13, 0),
                legend.text =  element_text(color = "black", size = 12), legend.title =  element_text(color = "black", face = "bold", size = 12))



  ## Forest plots of comparisons-specific between-trial standard deviation
  p2 <- if (model == "RE") {
    ggplot(data = prepare.tau, aes(x = as.factor(1:length(obs.comp)), y = as.numeric(median), ymin = as.numeric(lower), ymax = as.numeric(upper))) +
      geom_rect(aes(xmin = 0, xmax = Inf, ymin = tau.full[3], ymax = tau.full[7]), fill = "#1f93ff", alpha = 0.1) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = tau.full[5], lty = 1, size = 1, col = "#006CD1") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = as.factor(1:length(obs.comp)), y = round(as.numeric(median), 2), label = round(as.numeric(median), 2)),
                color = "black", hjust = -0.1, vjust = -0.5, size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.8), inherit.aes = T) +
      scale_x_discrete(breaks = as.factor(1:length(obs.comp)), labels = prepare.tau$comparison[1:length(obs.comp)]) +
      labs(x = "", y = "Between-trial standard deviation") +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold", size = 12))
  } else {
    NA
  }



  ## Bring together both forest-plots
  forest.plots <- if (model == "RE") {
    ggarrange(p1, p2, nrow = 1, ncol = 2, labels = c("A)", "B)"), common.legend = T, legend = "bottom")
  } else {
    p1
  }


  ## Write the table as .xlsx
  write_xlsx(EM.both.models, paste0(getwd(),"Table NMA vs MA.xlsx"))


  return(list(EM.both.models = EM.both.models, forest.plots = forest.plots))
}
