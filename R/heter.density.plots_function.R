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
#' @param heter1 Character string indicating the heterogeneity parameter
#'   for \code{heter_prior1}. Set \code{heter1} equal to one of the following:
#'   \code{"tau"}, or \code{"tau_omega"}, which refers to a between-study
#'   heterogeneity or between-design heterogeneity (inconsistency),
#'   respectively. This argument is relevant only when
#'   \code{distr = "lognormal"} or \code{distr = "logt"}. The default is
#'   \code{"tau"}.
#' @param heter2 Character string indicating the heterogeneity parameter
#'   for \code{heter_prior2}. Set \code{heter2} equal to one of the following:
#'   \code{"tau"}, or \code{"tau_omega"}, which refers to a between-study
#'   heterogeneity or between-design heterogeneity (inconsistency),
#'   respectively. This argument is relevant only when
#'   \code{distr = "lognormal"} or \code{distr = "logt"}. The default is
#'   \code{"tau"}.
#' @param caption Logical to indicate whether to report a caption at the bottom
#'   right of the plot. It is relevant only when \code{distr = "lognormal"} and
#'   \code{distr = "logt"}. The default is \code{FALSE} (do not report).
#' @param x_axis_name Logical to indicate whether to present the title of x-axis
#'   ('Between-study standard deviation'). The default is \code{TRUE} (report).
#' @param y_axis_name Logical to indicate whether to present the title of y-axis
#'   ('Density'). The default is \code{TRUE} (report).
#' @param title_name Text for the title of the plot. \code{title_name}
#'   determines the labs argument of the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param axis_title_size A positive integer for the font size of axis title.
#'   \code{axis_title_size} determines the axis.title argument found in the
#'   theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'   The default option is 13.
#' @param axis_text_size A positive integer for the font size of axis text.
#'   \code{axis_text_size} determines the axis.text argument found in the
#'   theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'   The default option is 13.
#' @param legend_title_size A positive integer for the font size of legend
#'   title. \code{legend_text_size} determines the legend.text argument found in
#'   the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'   The default option is 13.
#' @param legend_text_size A positive integer for the font size of legend text.
#'   \code{legend_text_size} determines the legend.text argument found in the
#'   theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'   The default option is 13.
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
#' \dontrun{
#' ## Two empirical priors for between-study variance of log odds ratio.
#' heter_density_plot(distr = "lognormal",
#'                    heter_prior1 = c(-2.56, 1.74),  # General healthcare setting
#'                    heter_prior2 = c(-1.83, 1.52))  # Pain and pharma vs. placebo/ctrl
#'
#' ## Two empirical priors for between-study variance of standardised mean
#' ## difference.
#' heter_density_plot(distr = "logt",
#'                    heter_prior1 = c(-3.44, 2.59),  # General healthcare setting
#'                    heter_prior2 = c(-0.60, 2.61))  # Pain and pharma vs. placebo/ctrl for cancer
#'
#' ## Two half-normal prior distributions for between-study standard deviation
#' heter_density_plot(distr = "halfnormal",
#'                    heter_prior1 = c(0, 1),
#'                    heter_prior2 = c(0, 0.5))
#' }
#'
#' @export
heter_density_plot <- function (distr,
                                heter_prior1,
                                heter_prior2,
                                heter1 = "tau",
                                heter2 = "tau",
                                caption = FALSE,
                                x_axis_name = TRUE,
                                y_axis_name = TRUE,
                                title_name = NULL,
                                axis_title_size = 13,
                                axis_text_size = 13,
                                legend_title_size = 13,
                                legend_text_size = 13) {


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
  het10 <- if (missing(heter1)) {
    "tau"
  } else if (!is.element(heter1, c("tau", "tau_omega"))) {
    stop("Insert 'tau', or 'tau_omega'", call. = FALSE)
  } else if (is.element(heter1, c("tau", "tau_omega"))) {
    heter1
  }
  het1 <- if (distr == "lognormal" & het10 == "tau") {
    "\u03C4\u00b2"
  } else if (distr == "lognormal" & het10 == "tau_omega") {
    "\u03C9\u00b2"
  } else if (distr == "logt" & het10 == "tau") {
    "ln(\u03C4\u00b2)"
  } else if (distr == "logt" & het10 == "tau_omega") {
    "ln(\u03C9\u00b2)"
  }
  het20 <- if (missing(heter2)) {
    "tau"
  } else if (!is.element(heter2, c("tau", "tau_omega"))) {
    stop("Insert 'tau', or 'tau_omega'", call. = FALSE)
  } else if (is.element(heter2, c("tau", "tau_omega"))) {
    heter2
  }
  het2 <- if (distr == "lognormal" & het20 == "tau") {
    "\u03C4\u00b2"
  } else if (distr == "lognormal" & het20 == "tau_omega") {
    "\u03C9\u00b2"
  } else if (distr == "logt" & het20 == "tau") {
    "ln(\u03C4\u00b2)"
  } else if (distr == "logt" & het20 == "tau_omega") {
    "ln(\u03C9\u00b2)"
  }
  caption <- if (caption == TRUE) {
    "*Note: Variance values are plotted on the logarithmic scale."
  } else {
    ""
  }
  x_axis_name <- if (x_axis_name == TRUE) {
    "Between-study standard deviation"
  } else {
    ""
  }
  y_axis_name <- if (y_axis_name == TRUE) {
    "Density"
  } else {
    ""
  }


  ## Function to power values close to zero
  label_value <- function (x) {
    unlist(lapply(x, function(x) if (x < 0.01)
      format(x, scientific = TRUE, digits = 2) else sprintf("%.2f", x)))
  }


  ## Self-written function for density probability of half-normal distribution
  dhnorm2 <- function (x, sigma) {

    theta <- sqrt(pi / 2) * (1 / sigma)
    res <- ((2 * theta) / pi) * exp(-(x * x * theta * theta)/pi)
    return(res)
  }


  value_x <- prob_dens_y <- NULL
  if (distr == "lognormal") {

    ## Give a distinct name to each function
    # For heter_prior1
    name1 <-
      paste0(het1, "~LN(", sprintf("%.2f", heter_prior1[1]), ", ",
             sprintf("%.2f", heter_prior1[2]), "\u00b2)")

    # For heter_prior2
    name2 <-
      paste0(het2, "~LN(", sprintf("%.2f", heter_prior2[1]), ", ",
             sprintf("%.2f", heter_prior2[2]), "\u00b2)")


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
                   function(x) sqrt(qlnorm(c(0.025, 0.25, 0.50, 0.75, 0.975),
                                           dataset_dist[x, 2],
                                           dataset_dist[x, 3])))
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
        fill = "#B4AF46",
        alpha = 0.2) +
      stat_function(fun = function(z) {
        dnorm(z, dataset_dist[1, 2], dataset_dist[1, 3])},
        xlim = c(min(log(breaks_tau2)),
                 log(qlnorm(0.999, dataset_dist[1, 2], dataset_dist[1, 3]))),
        col = "#B4AF46",
        linewidth = 1.3) +
      stat_function(fun = function(z) {
        dnorm(z, dataset_dist[2, 2], dataset_dist[2, 3])},
        xlim = c(min(log(breaks_tau2)),
                 log(qlnorm(0.5, dataset_dist[2, 2], dataset_dist[2, 3]))),
        geom = "area",
        fill = "#B4464B",
        alpha = 0.2) +
      stat_function(fun = function(z) {
        dnorm(z, dataset_dist[2, 2], dataset_dist[2, 3])},
        xlim = c(min(log(breaks_tau2)),
                 log(qlnorm(0.999, dataset_dist[2, 2], dataset_dist[2, 3]))),
        col = "#B4464B",
        linewidth = 1.3) +
      geom_point(data = dataset_pdf,
                 aes(x = value_x,
                     y = prob_dens_y,
                     fill = distr),
                 alpha = 0) +
      scale_x_continuous(breaks = log(breaks_tau2),
                         labels = label_value(sqrt(breaks_tau2)),
                         #limits = c(min(log(breaks_tau2)),
                         #           max(log(breaks_tau2))),
                         guide = guide_axis(check.overlap = TRUE)) +
      labs(x = x_axis_name,
           y = y_axis_name,
           fill = "Distribution",
           title_name = title_name,
           caption = caption) +
      guides(colour = "none",
             fill = guide_legend(override.aes = list(size = 3,
                                                     alpha = 1,
                                                     colour = c("#B4AF46",
                                                                "#B4464B")))) +
      theme_classic() +
      theme(plot.title = element_text(size = axis_title_size, face = "bold"),
            axis.text = element_text(size = axis_text_size),
            axis.title = element_text(size = axis_title_size, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = legend_text_size),
            legend.title = element_text(size = legend_title_size,
                                        face = "bold"))

  } else if (distr == "logt") {

    ## Give a distinct name to each function
    # For heter_prior1
    name1 <-
      paste0(het1, "~t(", sprintf("%.2f", heter_prior1[1]), ", ",
             sprintf("%.2f", heter_prior1[2]), "\u00b2, 5)")

    # For heter_prior2
    name2 <-
      paste0(het2, "~t(", sprintf("%.2f", heter_prior2[1]), ", ",
             sprintf("%.2f", heter_prior2[2]), "\u00b2, 5)")


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
                     sqrt(exp((qt(c(0.025, 0.25, 0.50, 0.75, 0.975), 5) *
                                 dataset_dist[x, 3]) + dataset_dist[x, 2])))
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
        fill = "#B4AF46",
        alpha = 0.2) +
      stat_function(fun = function(z) {(1 / dataset_dist[1, 3]) *
          dt((z - dataset_dist[1, 2]) / dataset_dist[1, 3], 5)},
          xlim = c(min(log(breaks_tau2)),
                   (qt(0.999, 5) * dataset_dist[1, 3]) + dataset_dist[1, 2]),
        col = "#B4AF46",
        linewidth = 1.3) +
      stat_function(fun = function(z) {(1 / dataset_dist[2, 3]) *
          dt((z - dataset_dist[2, 2]) / dataset_dist[2, 3], 5)},
        xlim = c(min(log(breaks_tau2)),
                 (qt(0.5, 5) * dataset_dist[2, 3]) + dataset_dist[2, 2]),
        geom = "area",
        fill = "#B4464B",
        alpha = 0.2) +
      stat_function(fun = function(z) {(1 / dataset_dist[2, 3]) *
          dt((z - dataset_dist[2, 2]) / dataset_dist[2, 3], 5)},
          xlim = c(min(log(breaks_tau2)),
                   (qt(0.999, 5) * dataset_dist[2, 3]) + dataset_dist[2, 2]),
        col = "#B4464B",
        linewidth = 1.3) +
      geom_point(data = dataset_pdf,
                 aes(x = value_x,
                     y = prob_dens_y,
                     fill = distr),
                 alpha = 0) +
      scale_x_continuous(breaks = log(breaks_tau2),
                         labels = label_value(sqrt(breaks_tau2)),
                         #limits = c(min(log(breaks_tau2)),
                         #           max(log(breaks_tau2))),
                         guide = guide_axis(check.overlap = TRUE)) +
      labs(x = x_axis_name,
           y = y_axis_name,
           fill = "Distribution",
           title_name = title_name,
           caption = caption) +
      guides(colour = "none",
             fill = guide_legend(override.aes = list(size = 3,
                                                     alpha = 1,
                                                     colour = c("#B4AF46",
                                                                "#B4464B")))) +
      theme_classic() +
      theme(plot.title = element_text(size = axis_title_size, face = "bold"),
            axis.text = element_text(size = axis_text_size),
            axis.title = element_text(size = axis_title_size, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = legend_text_size),
            legend.title = element_text(size = legend_title_size,
                                        face = "bold"))

  } else if (distr == "halfnormal") {


    ## Give a distinct name to each function
    # For heter_prior1
    name1 <- paste0("\u03C4~HN(", sprintf("%.2f", heter_prior1[1]), ", ",
                    sprintf("%.2f", heter_prior1[2]), ")")

    # For heter_prior2
    name2 <- paste0("\u03C4~HN(", sprintf("%.2f", heter_prior2[1]), ", ",
                    sprintf("%.2f", heter_prior2[2]), ")")


    ## Bring into a data-frame
    dataset_dist <-
      data.frame(distr = c(name1, name2),
                 mean = c(heter_prior1[1], heter_prior2[1]),
                 sd = c(heter_prior1[2], heter_prior2[2]))


    ## Define scale parameter
    #scale1 <- ifelse (heter_prior1[2] == 1, sqrt(pi / 2),
    #                  sqrt(pi / 2) * (1 / heter_prior1[2]))
    #scale2 <- ifelse (heter_prior2[2] == 1, sqrt(pi / 2),
    #                  sqrt(pi / 2) * (1 / heter_prior2[2]))


    ## Add to 'dataset_dist'
    #dataset_dist$scale <- c(scale1, scale2)


    ## Define tau range
    values <- seq(0, 2.0, 0.01)


    ## Obtain the pdf
    prob_dens <-
      lapply(1:dim(dataset_dist)[1],
             function(x) dhnorm2(values, dataset_dist[x, 3]))


    ## Bring together in a data-frame
    dataset_pdf <-
      data.frame(value_x = rep(values, 2),
                 prob_dens_y = unlist(prob_dens),
                 distr = rep(dataset_dist$distr, length(values)))
    dataset_pdf$distr <- factor(dataset_pdf$distr, levels = c(name1, name2))


    ## Define breaks (for ggplot2)
    breaks_tau <-
      unique(sort(unlist(lapply(1:dim(dataset_dist)[1],
                                function(x) qnorm(c(0.50, 0.75,  0.995),
                                                  0, dataset_dist[x, 3])))))


    ## Get table with quartiles
    tab0 <- lapply(1:dim(dataset_dist)[1],
                   function(x) qnorm(c(0.5125, 0.625, 0.75, 0.875, 0.9875),
                                     0, dataset_dist[x, 3]))
    tab <- matrix(label_value(unlist(tab0)),
                  nrow = 5, ncol = 2, byrow = FALSE)
    colnames(tab) <- c(name1, name2)
    rownames(tab) <- c("2.5%", "25%", "50%", "75%", "97.5%")


    ## Get density plots
    plot <-
      ggplot(dataset_pdf,
             aes(x = value_x,
                 y = prob_dens_y)) +
      stat_function(fun = function(z) {dhnorm2(z, dataset_dist[1, 3])},
                    xlim = c(min(breaks_tau),
                             qnorm(0.75, 0, dataset_dist[1, 3])),
                    geom = "area",
                    fill = "#B4AF46",
                    alpha = 0.2) +
      stat_function(fun = function(z) {dhnorm2(z, dataset_dist[1, 3])},
                    col = "#B4AF46",
                    linewidth = 1.3) +
      stat_function(fun = function(z) {dhnorm2(z, dataset_dist[2, 3])},
                    xlim = c(min(breaks_tau),
                             qnorm(0.75, 0, dataset_dist[2, 3])),
                    geom = "area",
                    fill = "#B4464B",
                    alpha = 0.2) +
      stat_function(fun = function(z) {dhnorm2(z, dataset_dist[2, 3])},
                    col = "#B4464B",
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
      labs(x = x_axis_name,
           y = y_axis_name,
           fill = "Distribution",
           title_name = title_name,
           captiotn = "") +
      guides(colour = "none",
             fill = guide_legend(override.aes = list(size = 3,
                                                     alpha = 1,
                                                     colour = c("#B4AF46",
                                                                "#B4464B")))) +
      theme_classic() +
      theme(plot.title = element_text(size = axis_title_size, face = "bold"),
            axis.text = element_text(size = axis_text_size),
            axis.title = element_text(size = axis_title_size, face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = legend_text_size),
            legend.title = element_text(size = legend_title_size,
                                        face = "bold"))

  }

  return(list(Density_plots = suppressWarnings({plot}),
              tabulated_percentiles =
                knitr::kable(tab,
                             align = "cc",
                             #col.names = c("Percentiles", colnames(tab)),
                             caption = "Percentiles in standard deviation")))
}
