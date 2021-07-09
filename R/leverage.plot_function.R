#' Create the leverage plot
#'
#' @description \code{leverage.plot} illustrates the posterior mean of deviance of the individual data points under the unrelated mean effects model (via \code{run.UME})
#'   against the posterior mean of deviance under the consistency model (via \code{run.model}).
#'
#' @param net An object of S3 class \code{\link{run.model}} or \code{\link{run.UME}}. See 'Value' in \code{\link{run.model}} and \code{\link{run.UME}}.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}.
#'   If the argument \code{drug.names} is not defined, the order of the interventions as they appear in \code{data} is used, instead.
#' @param title A title to indicate the model (consistency model or unrelated mean effects model).
#'
#' @return A scatterplot of the leverage against the square root of the posterior mean of residual deviance of the trial-arms under the consistency models or
#'   the unrelated mean effects model. The green, yellow, and red curves correspond to the parabola \eqn{x^2 + y = k} with \eqn{k} = 1, 2, and 3, respectively.
#'   Data points correspond to trial-arms. Data points found outside the yellow parabola are linked with a pair of numbers.
#'   The first number refers to the trial position in the dataset, and the second arm refers to the corresponding trial-arm (see 'Arguments' and 'Value' in \code{data.preparation}).
#'   These trial-arms contribute more than 1 to the posterior mean deviance and, hence, the model's poor fit.
#'
#' @details \code{leverage.plot} is integrated in the \code{UME.plot} function to create the leverage on the posterior mean of deviance under the consistency model and
#'   the unrelated mean effects model. These plots appear side-by-side in the output of \code{UME.plot}.
#'   Dias et al. (2010) used leverage plots to investigate the fit of the consistency and inconsistency model - the latter through
#'   the node-splitting approach.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{UME.plot}}, \code{\link{run.model}}, \code{\link{run.UME}}
#'
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed treatment comparison meta-analysis. \emph{Stat Med} 2010;\bold{29}(7-8):932--44. [\doi{10.1002/sim.3767}]
#'
#' @export
leverage.plot <- function(net, drug.names, title) {



  ## The results on the following parameters will be used:

  # Posterior mean of deviance contribution (observed outcomes)
  dev.o <- net$dev.o

  # Leverage (observed outcomes)
  lev.o <- net$leverage.o

  # The sign of the difference between observed and fitted number of observed outcome
  sign.o <- net$sign.dev.o



  ## A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
  Numextract <- function(string){
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
  }



  ## A matriX with the trial-arms of the analysed network: first column is the trial id, second column is the arm id (i.e. 2, 3, and so on)
  (trial.arm0 <- matrix(as.numeric(Numextract(names(dev.o[, 1]))), nrow = length(dev.o[, 1]), ncol = 2, byrow = T))



  ## use the matrix above to create a vector of trial-arms using the 'paste0' function: 'trial id', 'arm id'
  (trial.arm <- paste0(trial.arm0[, 1], ",", trial.arm0[, 2]))



  ## Prepare the data-frame for the leverage plot (using ggplot2)
  prepare.lev <- round(data.frame(sign.o*sqrt(as.vector(dev.o[, 1])), lev.o), 2)
  colnames(prepare.lev) <- c("signed.dev.o", "lev.o")



  ## Keep only trial-arms that exceed the parabola y + (x)^2 = c at c = 2.
  # Observed outcomes
  poor.o <- ifelse(prepare.lev$lev.o > 2 - (prepare.lev$signed.dev.o^2) | prepare.lev$lev.o < -(2 - (prepare.lev$signed.dev.o^2)), trial.arm, NA)
  poor.fit.o <- data.frame(prepare.lev[!is.na(poor.o), 1:2], poor.o[!is.na(poor.o)])
  colnames(poor.fit.o) <- c("signed.dev", "leverage", "poor")


  ## Leverage plot for observed outcomes
  observed <- ggplot(data = prepare.lev, aes(x = signed.dev.o, y = lev.o)) +
               geom_point(size = 2, colour = "black") +
               geom_smooth(aes(x = signed.dev.o, y = 1 - (signed.dev.o^2)), method = 'loess', formula = 'y ~ x',  colour = "#009E73", linetype = 2) +
               geom_smooth(aes(x = signed.dev.o, y = 2 - (signed.dev.o^2)), method = 'loess', formula = 'y ~ x',  colour = "orange", linetype = 2) +
               geom_smooth(aes(x = signed.dev.o, y = 3 - (signed.dev.o^2)), method = 'loess', formula = 'y ~ x',  colour = "#D55E00", linetype = 2) +
               geom_text_repel(data = poor.fit.o, aes(x = signed.dev, y = leverage, label = poor),
                         color = "blue", fontface = "bold", hjust = "right", size = 3.8, max.overlaps = Inf, nudge_x = -0.1, direction = "y") +
               labs(x = expression(""%+-% sqrt("Posterior mean of the residual deviance")), y = "Leverage (each data point)") +
               coord_cartesian(xlim = c(min(prepare.lev$signed.dev.o), max(prepare.lev$signed.dev.o)),
                               ylim = c(0, max(3 - (prepare.lev$signed.dev.o^2))), expand = T) +
               ggtitle(title) +
               theme_classic() +
               theme(axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12),
                     axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                     plot.title = element_text(color = "black", size = 11, face = "bold"))



  return(observed)

}
