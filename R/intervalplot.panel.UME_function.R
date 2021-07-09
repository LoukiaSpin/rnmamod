#' A panel of interval plots: Consistency model versus unrelated mean effects model
#'
#' @description This function creates a panel of interval plots on the summary effect sizes under the consistency model and the unrelated mean effects model.
#'   The number of interval plots equals the number of pairwise comparisons observed in the network.
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param ume An object of S3 class \code{\link{run.UME}}. See 'Value' in \code{\link{run.UME}}.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}.
#'   If the argument \code{drug.names} is not defined, the order of the interventions as they appear in \code{data} is used, instead.
#'
#' @return A panel of interval plots on the posterior mean and 95\% credible interval of the summary effect size under the consistency model and the improved
#'   unrelated mean effects model (Spineli, 2021) of all pairwise comparisons observed in the network.
#'
#' @details \code{intervalplot.panel.UME} is integrated in the \code{UME.plot} function. The consistency model and the unrelated mean effects model are
#'   abbreviated in the y-axis as 'NMA model' and 'UME model'. The intervals are highlighted with green, when the corresponding summary effect sizes do not cross the vertical line of no difference, and red otherwise.
#'   Grey panels refer to the frail comparisons as detected by the \code{run.UME} function (see 'Details' in \code{run.UME}).
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{run.UME}}, \code{\link{UME.plot}}
#'
#' @references
#' Spineli LM. A novel framework to evaluate the consistency assumption globally in a network of interventions. \emph{submitted} 2021.
#'
#' @export
intervalplot.panel.UME <- function(full, ume, drug.names) {


  EM.full <- full$EM
  EM.ume <- ume$EM
  obs.comp <- ume$obs.comp
  frail.comp <- ume$frail.comp
  measure <- effect.measure.name(full$measure)

  ## Possible and observed comparisons
  possible.comp <- possible.observed.comparisons(drug.names, obs.comp)


  # Indicate frail comparisons
  frail.comp.ind <- ifelse(is.element(possible.comp$obs.comp[, 3], frail.comp), "yes", "no")


  # Calculate bounds of 95% CrI
  ume.lower <- EM.ume[, 1] - 1.96*EM.ume[, 2]
  ume.upper <- EM.ume[, 1] + 1.96*EM.ume[, 2]


  ## Keep only the effect estimates according to the 'poss.pair.comp.clean' - Consistency model
  ume.mean <- round(EM.ume[, 1], 2)
  nma.mean <- round(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), 1], 2)
  nma.lower <- round(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), 3], 2)
  nma.upper <- round(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), 7], 2)


  # Indicate statistical significance
  nma.stat.signif <- ifelse(nma.lower > 0 | nma.upper < 0, "strong", "weak")
  ume.stat.signif <- ifelse(ume.lower > 0 | ume.upper < 0, "strong", "weak")


  # Create the dataframe
  data.set <- data.frame(c(nma.mean, ume.mean), c(nma.lower, ume.lower), c(nma.upper, ume.upper), c(nma.stat.signif, ume.stat.signif),
                         rep(possible.comp$obs.comp[, 4], 2), rep(c("NMA model", "UME model"), each = length(obs.comp)), frail.comp.ind)
  colnames(data.set) <- c("mean", "lower", "upper", "stat.sign", "comp", "analysis", "frail")


  ## Obtain forestplot
  ggplot(data = data.set, aes(x = as.factor(analysis), y = mean, ymin = lower, ymax = upper, colour = stat.sign)) +
    geom_rect(aes(fill = frail),xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf, alpha = 0.2) +
    geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
    geom_hline(yintercept = 0, lty = 1, size = 1, col = "black") +
    geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
    geom_text(aes(x = as.factor(analysis), y = mean, label = paste0(sprintf("%.2f", mean), " ", "(",
              sprintf("%.2f", lower), ",", " ", sprintf("%.2f", upper), ")")),
              color = "black", hjust = 0, vjust = -0.5, size = 3.3, check_overlap = F, parse = F,
              position = position_dodge(width = 0.8), inherit.aes = T) +
    facet_wrap(vars(factor(comp, levels = unique(data.set$comp))), scales = "fixed") +
    scale_fill_manual(breaks = c("yes", "no"), values = c("grey53", "white")) +
    scale_color_manual(breaks = c("strong", "weak"), values = c("#009E73", "#D55E00")) +
    #scale_y_continuous(trans = ifelse(!is.element(measure, c("Odds ratio", "Ratio of means")), "identity", "log10")) +
    scale_y_continuous(trans = "identity") +
    labs(x = "", y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), paste(measure, "(in logarithmic scale)"), measure), colour = "Evidence", fill = "") +
    coord_flip() +
    theme_classic() +
    guides(fill = FALSE) +
    theme(axis.text.x = element_text(color = "black", size = 11), axis.text.y = element_text(color = "black", size = 11),
          axis.title.x = element_text(color = "black", face = "bold", size = 11), strip.text = element_text(size = 11),
          legend.position = "bottom", legend.text = element_text(color = "black", size = 11),
          legend.title = element_text(color = "black", face = "bold", size = 11))

}
