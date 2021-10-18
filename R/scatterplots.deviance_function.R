#' Create the deviance scatterplot
#'
#' @description Illustrates the posterior mean of deviance of the individual
#'   data points under the unrelated mean effects model (via \code{run_ume})
#'   against the posterior mean of deviance under the consistency model
#'   (via \code{run_model}).
#'
#' @param full A numeric vector with the posterior mean of deviance obtained
#'   using the consistency model.
#' @param ume A numeric vector with the posterior mean of deviance obtained
#'   using the unrelated mean effects model.
#' @param colour A string to define the colour of the points in the plot.
#'
#' @return A scatterplot on the posterior mean deviance contribution of the
#'   individual data points under the consistency model and the unrelated mean
#'   effects model. Each data point corresponds to a trial-arm indicated by a
#'   pair of numbers. The first number refers to the trial position in the
#'   dataset, and the second arm refers to the corresponding trial-arm
#'   (see 'Arguments' and 'Value' in \code{data_preparation}).
#'
#' @details \code{scatterplots_dev} is integrated in the \code{ume_plot}
#'   function to create the scatterplot on the posterior mean of deviance under
#'   the compared models, also considered by Dias et al. (2013). When the
#'   majority of data points are scattered across the diagonal dotted line, we
#'   may conclude that the compared models have a good agreement. Data points
#'   scattered above and below the diagonal dotted line may contribute more to
#'   the poor fit of the unrelated mean effects model or the consistency model,
#'   respectively. However, using the \code{scatterplots_dev}, it is difficult
#'   to declare these data points as possible 'outliers'. Therefore,
#'   \code{scatterplots_dev} is accompanied by \code{bland_altman_plot} in the
#'   \code{ume_plot} function.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{ume_plot}}, \code{\link{run_model}},
#'   \code{\link{run_ume}}
#'
#' @references
#' Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence synthesis
#' for decision making 4: inconsistency in networks of evidence based on
#' randomized controlled trials.
#' \emph{Med Decis Making} 2013a;\bold{33}(5):641--56.
#' [\doi{10.1177/0272989X12455847}]
#'
#' @export
scatterplots_dev <- function(full, ume, colour) {

  # A function to extract numbers from a character.
  # Source: http://stla.github.io/stlapblog/posts/numextract.html
  numextract <- function(string) {
    unlist(regmatches(string, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", string)))
  }

  # A matriX with the trial-arms of the analysed network:
  # first column trial id, second column arm id (i.e. 2, 3, and so on)
  trial_arm0 <- matrix(as.numeric(numextract(names(full))),
                       nrow = length(full),
                       ncol = 2,
                       byrow = T)

  # Create a vector of trial-arms using the 'paste0' function
  trial_arm <- paste0(trial_arm0[, 1], ",", trial_arm0[, 2])

  # Prepare the data-frame for the scatterplot (using ggplot2)
  prepare_dev <- data.frame(trial_arm, full, ume)

  # Round to the second decimal
  scale_fun <- function(x) sprintf("%.2f", x)

  # Scatterplot on the observed outcomes
  ggplot(data = prepare_dev, aes(x = full, y = ume)) +
    geom_point(size = 2, colour = colour) +
    geom_abline(slope = 1, lty = 2, size = 1) +
    geom_text(aes(x = full,
                  y = ume,
                  label = trial_arm),
              color = "black",
              fontface = "bold",
              hjust = -0.2,
              vjust = -0.3,
              size = 4.0,
              check_overlap = T,
              inherit.aes = T) +
    labs(x = "Network meta-analysis",
         y = "Unrelated mean effects model") +
    scale_x_continuous(limits = c(min(prepare_dev$full, prepare_dev$ume),
                                  max(prepare_dev$full, prepare_dev$ume)),
                       labels = scale_fun) +
    scale_y_continuous(limits = c(min(prepare_dev$full, prepare_dev$ume),
                                  max(prepare_dev$full, prepare_dev$ume)),
                       labels = scale_fun) +
    theme_classic() +
    theme(axis.title.x = element_text(color = "black", size = 12),
          axis.title.y = element_text(color = "black", size = 12),
          axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12))
}
