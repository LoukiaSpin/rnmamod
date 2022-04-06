#' Leverage plot
#'
#' @description Plots the leverage against the square root of the
#'   posterior mean of residual deviance of the trial-arms under the model of
#'   interest.
#'
#' @param net An object of S3 class \code{\link{run_metareg}},
#'   \code{\link{run_model}}, or \code{\link{run_ume}}.
#'   See 'Value' in \code{\link{run_metareg}}, \code{\link{run_model}}, or
#'   \code{\link{run_ume}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}.
#'   If \code{drug_names} is not defined, the order of the
#'   interventions as they appear in \code{data} is used, instead.
#' @param title A title to indicate the model (consistency model, network
#'   meta-regression or unrelated mean effects model).
#'
#' @return A scatterplot of the leverage against the square root of the
#'   posterior mean of residual deviance of the trial-arms under the model of
#'   interest. The green, yellow, and red curves correspond to the parabola
#'   \eqn{x^2 + y = k} with \eqn{k} = 1, 2, and 3, respectively. The data points
#'   correspond to trial-arms. Data points found outside the yellow parabola are
#'   linked with a pair of numbers. The first number refers to the position
#'   of the trial in the dataset, and the second number refers to the
#'   corresponding trial-arm (see 'Arguments' and 'Value' in
#'   \code{\link{data_preparation}}). These trial-arms contribute more than
#'   1 to the deviance information criterion and, hence, the model's poor fit.
#'
#' @details \code{leverage_plot} is integrated in the \code{\link{ume_plot}}
#'   function to create the leverage plot for the consistency model and the
#'   unrelated mean effects model. These plots appear side-by-side in the output
#'   of \code{\link{ume_plot}}. Dias et al. (2010) used leverage plots to
#'   investigate the fit of the consistency and inconsistency models--the
#'   latter through the node-splitting approach.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{data_preparation}}, \code{\link{run_metareg}},
#'   \code{\link{run_model}}, \code{\link{run_ume}}, \code{\link{ume_plot}}
#'
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
#' treatment comparison meta-analysis.
#' \emph{Stat Med} 2010;\bold{29}(7-8):932--44. doi: 10.1002/sim.3767
#'
#' @export
leverage_plot <- function(net, drug_names, title) {

  # Posterior mean of deviance contribution (observed outcomes)
  dev_o <- net$dev_o

  # Leverage (observed outcomes)
  lev_o <- net$leverage_o

  # The sign of the difference between observed and fitted outcome
  sign_o <- net$sign_dev_o

  # A function to extract numbers from a character.
  # Source: http://stla.github.io/stlapblog/posts/numextract.html
  numextract <- function(string) {
    unlist(regmatches(string, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", string)))
  }

  # A matriX with the trial-arms of the analysed network:
  # first column is trial id, second column is arm id (i.e. 2, 3, and so on)
  trial_arm0 <- matrix(as.numeric(numextract(names(dev_o[, 1]))),
                       nrow = length(dev_o[, 1]),
                       ncol = 2,
                       byrow = TRUE)

  # Create a vector of trial-arms using the 'paste0' function
  trial_arm <- paste0(trial_arm0[, 1], ",", trial_arm0[, 2])

  # Prepare the data-frame for the leverage plot (using ggplot2)
  prepare_lev <- round(data.frame(sign_o * sqrt(as.vector(dev_o[, 1])), lev_o),
                       2)
  colnames(prepare_lev) <- c("signed_dev_o", "lev_o")

  # Keep only trial-arms that exceed the parabola y + (x)^2 = c at c = 2.
  poor_o <- ifelse(prepare_lev$lev_o > 2 - (prepare_lev$signed_dev_o^2) |
                     prepare_lev$lev_o < - (2 - (prepare_lev$signed_dev_o^2)),
                   trial_arm, NA)
  poor_fit_o <- data.frame(prepare_lev[!is.na(poor_o), 1:2],
                           poor_o[!is.na(poor_o)])
  colnames(poor_fit_o) <- c("signed_dev", "leverage", "poor")

  # Leverage plot for observed outcomes
  observed <- ggplot(data = prepare_lev, aes(x = signed_dev_o, y = lev_o)) +
               geom_point(size = 2, colour = "black") +
               geom_smooth(aes(x = signed_dev_o,
                               y = 1 - (signed_dev_o^2)),
                           method = "loess",
                           formula = "y ~ x",
                           colour = "#009E73",
                           linetype = 2) +
               geom_smooth(aes(x = signed_dev_o,
                               y = 2 - (signed_dev_o^2)),
                           method = "loess",
                           formula = "y ~ x",
                           colour = "orange",
                           linetype = 2) +
               geom_smooth(aes(x = signed_dev_o,
                               y = 3 - (signed_dev_o^2)),
                           method = "loess",
                           formula = "y ~ x",
                           colour = "#D55E00",
                           linetype = 2) +
               geom_text_repel(data = poor_fit_o,
                               aes(x = signed_dev,
                                   y = leverage,
                                   label = poor),
                         color = "blue",
                         fontface = "bold",
                         hjust = "right",
                         size = 3.8,
                         max.overlaps = Inf,
                         nudge_x = -0.1,
                         direction = "y") +
               labs(x = expression(
                 "" %+-% sqrt("Posterior mean of the residual deviance")),
                    y = "Leverage (each data point)") +
               coord_cartesian(xlim = c(min(prepare_lev$signed_dev_o),
                                        max(prepare_lev$signed_dev_o)),
                               ylim = c(0,
                                        max(3 - (prepare_lev$signed_dev_o^2))),
                               expand = TRUE) +
               ggtitle(title) +
               theme_classic() +
               theme(axis.title.x = element_text(color = "black", size = 12),
                     axis.title.y = element_text(color = "black", size = 12),
                     axis.text.x = element_text(color = "black", size = 12),
                     axis.text.y = element_text(color = "black", size = 12),
                     plot.title = element_text(color = "black", size = 11,
                                               face = "bold"))

  return(observed)
}
