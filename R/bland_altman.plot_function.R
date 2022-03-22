#' The Bland-Altman plot
#'
#' @description This function facilitates creating the Bland-Altman plot on the
#'  posterior mean deviance contribution for two models using only three
#'  arguments.
#'
#' @param model1 A vector with the numeric values of the target model (for
#'   instance, the consistency model).
#' @param model2 A vector with the numeric values of the reference model (for
#'   instance, the unrelated mean effects model).
#' @param colour A string to define the colour of the data points in the plot.
#'
#' @return Bland-Altman plot on the posterior mean deviance contribution of the
#'  individual data points under model 1 and model 2.
#'  Each data point corresponds to a trial-arm indicated by a pair of numbers.
#'  The first number refers to the position of the trial in the dataset,
#'  and the second arm refers to the corresponding trial-arm (see 'Arguments'
#'  and 'Value' in \code{\link{data_preparation}}).
#'  The plot also displays the average bias and the 95\% limits of agreement
#'  with horizontal solid black lines.
#'
#' @details \code{bland_altman_plot} is integrated in \code{\link{ume_plot}}
#'   to create the Bland-Altman plot on the posterior mean of deviance
#'   under the consistency model (via \code{\link{run_model}}) and the
#'   unrelated mean effects model (via \code{\link{run_ume}}).
#'
#'   A uniform scattering of the data points within the 95\% limits of agreement
#'   and average bias close to 0 indicate that the compared models have a good
#'   agreement. Data points positioned above or below the 95\% limits of
#'   agreement correspond to trials that contribute to the poor fit of the
#'   consistency model or unrelated mean effects model, respectively.
#'
#'   \code{bland_altman_plot} can be used to compare the following models
#'   regarding deviance contribution:
#'   \itemize{
#'    \item the consistency model (via \code{\link{run_model}}) with the
#'    unrelated effect means model (via \code{\link{run_ume}});
#'    \item the network meta-analysis model (via \code{\link{run_model}}) with
#'    the network meta-analysis model (via \code{\link{run_metareg}}).
#'    }
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{data_preparation}}, \code{\link{run_metareg}},
#'   \code{\link{run_model}}, \code{\link{run_ume}}, \code{\link{ume_plot}}
#'
#' @references
#' Bland JM, Altman DG. Measuring agreement in method comparison studies.
#' \emph{Stat Methods Med Res} 1999;\bold{8}:135--60.
#'
#' @export
bland_altman_plot <- function(model1, model2, colour) {

  # A function to extract numbers from a character.
  #Source: http://stla.github.io/stlapblog/posts/numextract.html
  numextract <- function(string) {
    unlist(regmatches(string, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", string)))
  }

  # A matriX with the trial-arms of the analysed network:
  # first column is trial id, second column is arm id (i.e. 2, 3, and so on)
  trial_arm0 <- matrix(as.numeric(numextract(names(model1))),
                        nrow = length(model1),
                        ncol = 2,
                        byrow = TRUE)

  # Use the matrix above to create a vector of trial-arms using the
  # 'paste0' function: 'trial id', 'arm id'
  trial_arm <- paste0(trial_arm0[, 1], ",", trial_arm0[, 2])

  # Round to the second decimal
  scale_fun <- function(x) sprintf("%.2f", x)

  difference <- model1 - model2
  average <- apply(cbind(model1, model2), 1, mean)

  # Dataset to be plotted in ggplot2
  data <- data.frame(round(average, 4), round(difference, 4), trial_arm)

  # Bias and Limits of Agreement (LoA)
  bias <- mean(difference)
  lower_lim <- bias - 1.96 * sd(difference)
  upper_lim  <- bias + 1.96 * sd(difference)

  ggplot(data, aes(x = data[, 1], y = data[, 2])) +
    geom_point(size = 2, colour = colour) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_hline(yintercept = bias, linetype = 1) +
    geom_hline(yintercept = lower_lim, linetype = 1) +
    geom_hline(yintercept = upper_lim, linetype = 1) +
    geom_text(aes(x = data[, 1],
                  y = data[, 2],
                  label = data[, 3]),
              color = "black",
              fontface = "bold",
              hjust = -0.2,
              vjust = -0.3,
              size = 4.0,
              check_overlap = TRUE,
              inherit.aes = TRUE) +
    scale_y_continuous(labels = scale_fun) +
    scale_x_continuous(labels = scale_fun) +
    theme_classic() +
    labs(x = "Average posterior mean deviance contribution",
    y = "Difference in posterior mean deviance contribution") +
    theme(axis.text = element_text(size = 12),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 12))
}
