#' Visualising the density of two prior distributions for the heterogeneity parameter
#'
#' @description
#' Creating the density plot of two prior distributions for the between-study
#' variance (log-normal and location-scale t distributions) or between-study
#' standard deviation (half-normal distribution).
#'
#' @param distr Character string indicating the prior distribution.
#'   Set \code{distr} equal to one of the following: \code{"lognormal"},
#'   \code{"logt"}, or \code{"halfnormal"}, which refers to a log-normal,
#'   location-scale, or half-normal distribution.
#' @param heter_prior1 A numeric vector with two values for the first prior
#'   distribution: 1) the mean value and 2) the standard deviation. When
#'   \code{distr = "halfnormal"}, the first value should zero and the second a
#'   non-negative value referring to the scale parameter of the distribution.
#' @param heter_prior2 A numeric vector with two values for the second prior
#'   distribution: 1) the mean value and 2) the standard deviation. When
#'   \code{distr = "halfnormal"}, the first value should zero and the second a
#'   non-negative value referring to the scale parameter of the distribution.
#' @param caption Logical to indicate whether to report a caption at the bottom
#'   right of the plot. It is relevant only when \code{distr = "lognormal"} and
#'   \code{distr = "logt"}. The default is \code{FALSE} (do not report).
#'
#' @return A plot with the density of two selected prior distributions for the
#' heterogeneity parameter. Two different colours are used to discern the
#' distributions. A legend is also created with the name and hyper-parameters of
#' the selected prior distributions. The filled area under each curved indicates
#' the values up to the median of the corresponding distribution. The x-axis
#' present the 0.1%, 50%, 99% percentile of each distribution.
#'
#' \code{heter_density_plot} also returns a table with the percentiles of each
#' distribution.
#'
#' @details
#' Use this function to inspect the shape of the distribution and the range of
#' between-study variance or standard deviation values before you define the
#' argument \code{heter_prior} in \code{\link{run_model}}) to run random-effects
#' network meta-analysis.
#'
#' Turner et al. (2012), Turner et al. (2015), and Rhodes et al. (2016) provide
#' predictive prior distributions for the between-study variance for a binary
#' outcome, measured in the log-odds ratio scale, and a continuous outcome,
#' measured in the standardised mean difference scale, respectively.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_model}}
#'
#' @references
#' Rhodes KM, Turner RM, Higgins JP. Predictive distributions were developed
#' for the extent of heterogeneity in meta-analyses of continuous outcome data.
#' \emph{J Clin Epidemiol} 2015;\bold{68}(1):52--60.
#' doi: 10.1016/j.jclinepi.2014.08.012
#'
#' Turner RM, Jackson D, Wei Y, Thompson SG, Higgins JP. Predictive
#' distributions for between-study heterogeneity and simple methods for their
#' application in Bayesian meta-analysis.
#' \emph{Stat Med} 2015;\bold{34}(6):984--98. doi: 10.1002/sim.6381
#'
#' Turner RM, Davey J, Clarke MJ, Thompson SG, Higgins JP. Predicting the extent
#' of heterogeneity in meta-analysis, using empirical data from the Cochrane
#' Database of Systematic Reviews.
#' \emph{Int J Epidemiol} 2012;\bold{41}(3):818--27. doi: 10.1093/ije/dys041
#'
#' @examples
#'
#' ## Two empirical priors for between-study variance of log odds ratio.
#' heter_density_plot(distr = "lognormal",
#'                    heter_prior1 = c(-3.50, 1.26),
#'                    heter_prior2 = c(-4.79, 1.67))
#'
#' ## Two empirical priors for between-study variance of standardised mean
#' ## difference.
#' heter_density_plot(distr = "logt",
#'                    heter_prior1 = c(-2.76, 2.58),
#'                    heter_prior2 = c(-2.43, 2.50))
#'
#' ## Two half-normal prior distributions for between-study standard deviation
#' heter_density_plot(distr = "halfnormal",
#'                    heter_prior1 = c(0, 1),
#'                    heter_prior2 = c(0, 0.5))
#'
#' @export
heter_density_plot <- function (distr,
                                heter_prior1,
                                heter_prior2,
                                caption = TRUE) {


  ## Default arguments
  distr <- if (missing(heter_prior1)) {
    stop("The argument 'distr' must be defined", call. = FALSE)
  } else if (!is.element(distr, c("halfnormal", "lognormal", "logt"))) {
    stop("Insert 'halfnormal', 'lognormal', or 'logt'",
         call. = FALSE)
  } else {
    distr
  }
  heter_prior1 <- if (missing(heter_prior1)) {
    stop("The argument 'heter_prior1' must be defined", call. = FALSE)
  } else if (distr == "halfnormal" & heter_prior1[1] < 0) {
    stop("The first element must be a non-negative number.", call. = FALSE)
  } else if (distr == "halfnormal" & heter_prior1[1] > heter_prior1[[2]]) {
    stop("The second element must be larger than the first element.",
         call. = FALSE)
  } else if (heter_prior1[2] < 0) {
    stop("The second element' must be a positive number", call. = FALSE)
  } else {
    heter_prior1
  }
  heter_prior2 <- if (missing(heter_prior2)) {
    stop("The argument 'heter_prior1' must be defined", call. = FALSE)
  } else if (distr == "halfnormal" & heter_prior2[1] < 0) {
    stop("The first element must be a non-negative number.", call. = FALSE)
  } else if (distr == "halfnormal" & heter_prior2[1] > heter_prior2[[2]]) {
    stop("The second element must be larger than the first element.",
         call. = FALSE)
  } else if (heter_prior2[2] < 0) {
    stop("The second element' must be a positive number", call. = FALSE)
  } else {
    heter_prior2
  }
  caption <- if (caption == TRUE) {
    "*Note: Distributions are plotted on the logarithmic scale."
  } else {
    ""
  }

  ## Function for log-transformed values (relevant for 'lognormal' and 'logt')
  label_value <- function (x) {
    unlist(lapply(x, function(x) if (x < 0.0001)
      format(x, scientific = TRUE, digits = 2) else round(x, 2)))
  }

  value_x <- prob_dens_y <- NULL
  if (distr == "lognormal") {

    ## Give a distinct name to each function
    # For heter_prior1
    name1 <- paste0("LN(", heter_prior1[1], ", ", heter_prior1[2], "\u00b2)")

    # For heter_prior2
    name2 <- paste0("LN(", heter_prior2[1], ", ", heter_prior2[2], "\u00b2)")


    ## Bring into a data-frame
    dataset_dist <-
      data.frame(distr = c(name1, name2),
                 mean = c(heter_prior1[1], heter_prior2[1]),
                 sd = c(heter_prior1[2], heter_prior2[2]))


    ## Define tau2 range
    values <- seq(0.00001, 3, 0.01)


    ## Log-transform tau2 range
    values_tau2_log <- log(values)


    ## Obtain the pdf
    prob_dens <-
      lapply(1:dim(dataset_dist)[1],
             function(x)
               dnorm(values_tau2_log, dataset_dist[x, 2], dataset_dist[x, 3]))


    ## Bring together in a data-frame
    dataset_pdf <-
      data.frame(value_x = rep(values_tau2_log, 2),
                 prob_dens_y = unlist(prob_dens),
                 distr = rep(dataset_dist$distr, length(values)))
    dataset_pdf$distr <- factor(dataset_pdf$distr, levels = c(name1, name2))


    ## Define breaks (for ggplot2)
    breaks_tau2 <-
      unique(sort(unlist(lapply(1:dim(dataset_dist)[1],
                                function(x)
                                  qlnorm(c(0.001, 0.50,  0.99),
                                         dataset_dist[x, 2],
                                         dataset_dist[x, 3])))))


    ## Get table with quartiles
    tab0 <- lapply(1:dim(dataset_dist)[1],
                   function(x) qlnorm(c(0.025, 0.25, 0.50, 0.75, 0.975),
                                      dataset_dist[x, 2],
                                      dataset_dist[x, 3]))
    tab <- matrix(label_value(unlist(tab0)),
                  nrow = 5, ncol = 2, byrow = FALSE)
    colnames(tab) <- c(name1, name2)
    rownames(tab) <- c("2.5%", "25%", "50%", "75%", "97.5%")


    ## Get density plots
    plot <-
      ggplot(dataset_pdf,
             aes(x = value_x,
                 y = prob_dens_y)) +
      stat_function(fun = function(z) {
        dnorm(z, dataset_dist[1, 2], dataset_dist[1, 3])},
        xlim = c(min(log(breaks_tau2)),
                 log(qlnorm(0.5, dataset_dist[1, 2], dataset_dist[1, 3]))),
        geom = "area",
        fill = "#0072B2",
        alpha = 0.2) +
      stat_function(fun = function(z) {
        dnorm(z, dataset_dist[1, 2], dataset_dist[1, 3])},
        col = "#0072B2",
        linewidth = 1.3) +
      stat_function(fun = function(z) {
        dnorm(z, dataset_dist[2, 2], dataset_dist[2, 3])},
        xlim = c(min(log(breaks_tau2)),
                 log(qlnorm(0.5, dataset_dist[2, 2], dataset_dist[2, 3]))),
        geom = "area",
        fill = "grey20",
        alpha = 0.2) +
      stat_function(fun = function(z) {
        dnorm(z, dataset_dist[2, 2], dataset_dist[2, 3])},
        col = "grey40",
        linewidth = 1.3) +
      geom_point(data = dataset_pdf,
                 aes(x = value_x,
                     y = prob_dens_y,
                     fill = distr),
                 alpha = 0) +
      scale_x_continuous(breaks = log(breaks_tau2),
                         labels = label_value(breaks_tau2),
                         #limits = c(min(log(breaks_tau2)),
                         #           max(log(breaks_tau2))),
                         guide = guide_axis(check.overlap = TRUE)) +
      labs(x = "Between-study variance",
           y = "Density",
           fill = "Distribution",
           caption = caption) +
      guides(colour = "none",
             fill = guide_legend(override.aes = list(size = 3,
                                                     alpha = 1,
                                                     colour = c("#0072B2",
                                                                "grey40")))) +
      theme_classic() +
      theme(axis.text = element_text(size = 13),
            axis.title = element_text(size = 13, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 13, face = "bold"))

  } else if (distr == "logt") {

    ## Give a distinct name to each function
    # For heter_prior1
    name1 <- paste0("t(", heter_prior1[1], ", ", heter_prior1[2], "\u00b2, 5)")

    # For heter_prior2
    name2 <- paste0("t(", heter_prior2[1], ", ", heter_prior2[2], "\u00b2, 5)")


    ## Bring into a data-frame
    dataset_dist <-
      data.frame(distr = c(name1, name2),
                 mean = c(heter_prior1[1], heter_prior2[1]),
                 sd = c(heter_prior1[2], heter_prior2[2]))


    ## Define tau2 range
    values <- seq(0.00001, 3, 0.01)


    ## Log-transform tau2 range
    values_tau2_log <- log(values)


    ## Obtain the pdf
    prob_dens <-
      lapply(1:dim(dataset_dist)[1],
             function(x)
               (1 / dataset_dist[x, 3]) *
               dt((values_tau2_log - dataset_dist[x, 2]) / dataset_dist[x, 3],
                  5))


    ## Bring together in a data-frame
    dataset_pdf <-
      data.frame(value_x = rep(values_tau2_log, 2),
                 prob_dens_y = unlist(prob_dens),
                 distr = rep(dataset_dist$distr, length(values)))
    dataset_pdf$distr <- factor(dataset_pdf$distr, levels = c(name1, name2))


    ## Define breaks (for ggplot2)
    breaks_tau2 <-
      unique(sort(unlist(
        lapply(1:dim(dataset_dist)[1],
               function(x)
                 exp((qt(c(0.001, 0.50,  0.99), 5) * dataset_dist[x, 3]) +
                       dataset_dist[x, 2])))))


    ## Get table with quartiles
    tab0 <- lapply(1:dim(dataset_dist)[1],
                   function(x)
                     exp((qt(c(0.025, 0.25, 0.50, 0.75, 0.975), 5) *
                            dataset_dist[x, 3]) + dataset_dist[x, 2]))
    tab <- matrix(label_value(unlist(tab0)),
                  nrow = 5, ncol = 2, byrow = FALSE)
    colnames(tab) <- c(name1, name2)
    rownames(tab) <- c("2.5%", "25%", "50%", "75%", "97.5%")


    ## Get density plots
    plot <-
      ggplot(dataset_pdf,
             aes(x = value_x,
                 y = prob_dens_y)) +
      stat_function(fun = function(z) {(1 / dataset_dist[1, 3]) *
          dt((z - dataset_dist[1, 2]) / dataset_dist[1, 3], 5)},
        xlim = c(min(log(breaks_tau2)),
                 (qt(0.5, 5) * dataset_dist[1, 3]) + dataset_dist[1, 2]),
        geom = "area",
        fill = "#0072B2",
        alpha = 0.2) +
      stat_function(fun = function(z) {(1 / dataset_dist[1, 3]) *
          dt((z - dataset_dist[1, 2]) / dataset_dist[1, 3], 5)},
        col = "#0072B2",
        linewidth = 1.3) +
      stat_function(fun = function(z) {(1 / dataset_dist[2, 3]) *
          dt((z - dataset_dist[2, 2]) / dataset_dist[2, 3], 5)},
        xlim = c(min(log(breaks_tau2)),
                 (qt(0.5, 5) * dataset_dist[2, 3]) + dataset_dist[2, 2]),
        geom = "area",
        fill = "grey20",
        alpha = 0.2) +
      stat_function(fun = function(z) {(1 / dataset_dist[2, 3]) *
          dt((z - dataset_dist[2, 2]) / dataset_dist[2, 3], 5)},
        col = "grey40",
        linewidth = 1.3) +
      geom_point(data = dataset_pdf,
                 aes(x = value_x,
                     y = prob_dens_y,
                     fill = distr),
                 alpha = 0) +
      scale_x_continuous(breaks = log(breaks_tau2),
                         labels = label_value(breaks_tau2),
                         #limits = c(min(log(breaks_tau2)),
                         #           max(log(breaks_tau2))),
                         guide = guide_axis(check.overlap = TRUE)) +
      labs(x = "Between-study variance",
           y = "Density",
           fill = "Distribution",
           caption = caption) +
      guides(colour = "none",
             fill = guide_legend(override.aes = list(size = 3,
                                                     alpha = 1,
                                                     colour = c("#0072B2",
                                                                "grey40")))) +
      theme_classic() +
      theme(axis.text = element_text(size = 13),
            axis.title = element_text(size = 13, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 13, face = "bold"))

  } else if (distr == "halfnormal") {


    ## Give a distinct name to each function
    # For heter_prior1
    name1 <- paste0("HN(", heter_prior1[1], ", ", heter_prior1[2], ")")

    # For heter_prior2
    name2 <- paste0("HN(", heter_prior2[1], ", ", heter_prior2[2], ")")


    ## Bring into a data-frame
    dataset_dist <-
      data.frame(distr = c(name1, name2),
                 mean = c(heter_prior1[1], heter_prior2[1]),
                 sd = c(heter_prior1[2], heter_prior2[2]))


    ## Define scale parameter
    scale1 <- ifelse (heter_prior1[2] == 1, sqrt(pi / 2),
                      sqrt(pi / 2) * (1 / heter_prior1[2]))
    scale2 <- ifelse (heter_prior2[2] == 1, sqrt(pi / 2),
                      sqrt(pi / 2) * (1 / heter_prior2[2]))


    ## Add to 'dataset_dist'
    dataset_dist$scale <- c(scale1, scale2)


    ## Define tau range
    values <- seq(0, 3, 0.01)


    ## Obtain the pdf
    prob_dens <-
      lapply(1:dim(dataset_dist)[1],
             function(x) dhalfnorm(values, dataset_dist[x, 4]))


    ## Bring together in a data-frame
    dataset_pdf <-
      data.frame(value_x = rep(values, 2),
                 prob_dens_y = unlist(prob_dens),
                 distr = rep(dataset_dist$distr, length(values)))
    dataset_pdf$distr <- factor(dataset_pdf$distr, levels = c(name1, name2))


    ## Define breaks (for ggplot2)
    breaks_tau <-
      unique(sort(unlist(lapply(1:dim(dataset_dist)[1],
                                function(x) qhalfnorm(c(0, 0.50,  0.99),
                                                      dataset_dist[x, 4])))))


    ## Get table with quartiles
    tab0 <- lapply(1:dim(dataset_dist)[1],
                   function(x) qhalfnorm(c(0.025, 0.25, 0.50, 0.75, 0.975),
                                         dataset_dist[x, 4]))
    tab <- matrix(label_value(unlist(tab0)),
                  nrow = 5, ncol = 2, byrow = FALSE)
    colnames(tab) <- c(name1, name2)
    rownames(tab) <- c("2.5%", "25%", "50%", "75%", "97.5%")


    ## Get density plots
    plot <-
      ggplot(dataset_pdf,
             aes(x = value_x,
                 y = prob_dens_y)) +
      stat_function(fun = function(z) {dhalfnorm(z, dataset_dist[1, 4])},
                    xlim = c(min(breaks_tau),
                             qhalfnorm(0.5, dataset_dist[1, 4])),
                    geom = "area",
                    fill = "#0072B2",
                    alpha = 0.2) +
      stat_function(fun = function(z) {dhalfnorm(z, dataset_dist[1, 4])},
                    col = "#0072B2",
                    linewidth = 1.3) +
      stat_function(fun = function(z) {dhalfnorm(z, dataset_dist[2, 4])},
                    xlim = c(min(breaks_tau),
                             qhalfnorm(0.5, dataset_dist[2, 4])),
                    geom = "area",
                    fill = "grey20",
                    alpha = 0.2) +
      stat_function(fun = function(z) {dhalfnorm(z, dataset_dist[2, 4])},
                    col = "grey40",
                    linewidth = 1.3) +
      geom_point(data = dataset_pdf,
                 aes(x = value_x,
                     y = prob_dens_y,
                     fill = distr),
                 alpha = 0) +
      scale_x_continuous(breaks = breaks_tau,
                         labels = sprintf("%.2f", breaks_tau),
                         #limits = c(min(breaks_tau), max(breaks_tau)),
                         guide = guide_axis(check.overlap = TRUE)) +
      labs(x = "Between-study standard deviation",
           y = "Density",
           fill = "Distribution",
           captiotn = "") +
      guides(colour = "none",
             fill = guide_legend(override.aes = list(size = 3,
                                                     alpha = 1,
                                                     colour = c("#0072B2",
                                                                "grey40")))) +
      theme_classic() +
      theme(axis.text = element_text(size = 13),
            axis.title = element_text(size = 13, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 13),
            legend.title = element_text(size = 13, face = "bold"))

  }

  return(list(Density_plots = plot,
              tabulated_percentiles =
                knitr::kable(tab,
                             align = "cc",
                             col.names = c("Percentiles", colnames(tab)),
                             caption = "Percentiles of each distribution")))
}
