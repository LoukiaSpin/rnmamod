#' Deviance scatterplot
#'
#' @description Illustrates the posterior mean of deviance contribution of the
#'   individual data points under the unrelated mean effects model (via
#'   \code{\link{run_ume}}) against the posterior mean of deviance contribution
#'   under the consistency model (via \code{\link{run_model}}).
#'
#' @param full A numeric vector with the posterior mean of deviance obtained
#'   using the consistency model (see 'Value' in \code{\link{run_model}}).
#' @param ume A numeric vector with the posterior mean of deviance obtained
#'   using the unrelated mean effects model (see 'Value' in
#'   \code{\link{run_ume}}).
#' @param colour A string to define the colour of the points in the plot.
#'
#' @return A scatterplot of the posterior mean deviance contribution of the
#'   individual data points from the unrelated mean effects model against those
#'   from the consistency model. Each data point corresponds to a trial-arm
#'   indicated by a pair of numbers. The first number refers to the position
#'   of the trial in the dataset, and the second number refers to the
#'   corresponding trial-arm (see 'Arguments' and 'Value' in
#'   \code{\link{data_preparation}}).
#'
#' @details \code{scatterplots_dev} is integrated in the \code{\link{ume_plot}}
#'   function to compare the models regarding the posterior mean of deviance.
#'   This scatterplot has also been considered by Dias et al., (2013).
#'   When the majority of data points are scattered across the diagonal line,
#'   we may conclude that the compared models have a good agreement. Data points
#'   systematically scattered above and below the diagonal line may
#'   contribute more to the poor fit of the unrelated mean effects model and the
#'   consistency model, respectively.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{data_preparation}}, \code{\link{run_model}},
#'   \code{\link{run_ume}}, \code{\link{ume_plot}}
#'
#' @references
#' Dias S, Welton NJ, Sutton AJ, Caldwell DM, Lu G, Ades AE. Evidence synthesis
#' for decision making 4: inconsistency in networks of evidence based on
#' randomized controlled trials.
#' \emph{Med Decis Making} 2013a;\bold{33}(5):641--56.
#' \doi{10.1177/0272989X12455847}
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
                       byrow = TRUE)

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
              check_overlap = TRUE,
              inherit.aes = TRUE) +
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
