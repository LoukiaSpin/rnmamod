#' Predictive distributions for the between-study variance in a future
#' meta-analysis on odds ratio or standardised mean difference
#'
#' @description
#' A table with the hyperparameters of the predictive distributions for the
#' between-study variance developed by Turner et al. (2015) and Rhodes et al.
#' (2015): log-normal distribution and t-distribution (with 5 degrees of
#' freedom) when the outcome data are analysed in the odds ratio or
#' standardised mean difference scale, respectively.
#'
#' @param measure Character string indicating the effect measure with possible
#'   values \code{"OR"} for odds ratio and \code{"SMD"} for standardised mean
#'   difference.
#' @param area Character string indicating the medical area relating to the
#'   predictive distributions for standardised mean difference with possible
#'   values \code{"cancer"} for medical areas of cancer, \code{"respiratory"}
#'   for medical areas of respiratory diseases, and \code{"other"} for medical
#'   areas other than cancer or respiratory diseases. The argument is not
#'   relevant for odds ratio.
#'
#' @return A cross-sectional table as a heatmap showing the hyperparameters
#' (mean and standard deviation) of the corresponding predictive distribution
#' for all combinations between the outcome types and treatment-comparison types
#' and according to the selected medical area (only relevant with standardised
#' mean difference) as defined by Turner et al. (2015) and Rhodes et al. (2015).
#' The tiles are coloured with different shades according to the corresponding
#' median value: the larger the median, the darker the colour.
#'
#' @details
#' This table aids in selecting the hyperparameters for the function
#' \code{\link{heterogeneity_param_prior}} when considering an informative prior
#' distribution for the between-study variance parameter based on the two
#' publications mentioned above (relevant for the function
#' \code{\link{run_model}} to conduct random-effects network meta-analysis).
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{heterogeneity_param_prior}}, \code{\link{run_model}}
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
#' @export
table_tau2_prior <- function (measure, area) {


  ## Default arguments
  measure <- if (missing(measure)) {
    stop("The argument 'measure' must be defined.", call. = FALSE)
  } else if (!is.element(measure, c("OR", "SMD"))) {
    stop("Insert 'OR', or 'SMD'.", call. = FALSE)
  } else if (is.element(measure, c("OR", "SMD"))) {
    measure
  }
  area <- if (measure == "OR") {
    NULL
  } else if (measure == "SMD" & missing(area)) {
    stop("The argument 'area' must be defined.", call. = FALSE)
  } else if (measure == "SMD" & !is.element(area, c("cancer", "respiratory", "other"))) {
    stop("Insert 'cancer', 'respiratory', or 'other'.", call. = FALSE)
  } else if (measure == "SMD" & is.element(area, c("cancer", "respiratory", "other"))) {
    area
  }


  if (measure == "OR") {
    ## Outcome types as suggested by Turner et al.
    outcome_type <- c("All-cause mortality",
                      "Obstetric outcomes",
                      "Cause-specific mortality/major morbidity event/composite",
                      "Resource use/hospital stay/process",
                      "Surgical/device related success/failure",
                      "Withdrawals/drop-outs",
                      "Internal/external structure-related outcomes",
                      "General physical health indicators",
                      "Adverse events",
                      "Infection/onset of new disease",
                      "Signs/symptoms reflecting continuation/end of condition",
                      "Pain",
                      "Quality of life/functioning (dichotomised)",
                      "Mental health indicators",
                      "Biological markers (dichotomised)",
                      "Subjective outcomes (various)",
                      "General healthcare setting")


    ## Treatment-comparison types as suggested by Turner et al.
    interv_compat_type <- c("Pharma vs. Placebo/control",
                            "Pharma vs. Pharma",
                            "Non-pharma vs. Placebo/control",
                            "Non-pharma vs. Pharma",
                            "Non-pharma vs. Non-pharma")


    ## Mean values from the prior log-normal distribution for the between-study variance (Turner et al.)
    mean_value <- c(-3.95, -3.52, -3.71, -2.34, -2.14, -2.99, -2.71, -2.29, -1.87, -2.49, -2.06, -1.83, -2.54, -2.12, -1.77, -2.70, -2.56,
                    -4.18, -3.75, -3.95, -2.58, -2.37, -3.23, -2.94, -2.53, -2.10, -2.73, -2.29, -2.06, -2.78, -2.35, -2.00, -2.93, -2.56,
                    -4.17, -3.74, -3.93, -2.56, -2.36, -3.21, -2.93, -2.51, -2.10, -2.71, -2.28, -2.05, -2.77, -2.34, -1.99, -2.92, -2.56,
                    -2.92, -2.49, -2.68, -1.31, -1.11, -1.96, -1.67, -1.26, -0.84, -1.46, -1.03, -0.80, -1.51, -1.09, -0.74, -1.67, -2.56,
                    -3.50, -3.08, -3.27, -1.90, -1.69, -2.55, -2.26, -1.85, -1.43, -2.05, -1.61, -1.38, -2.10, -1.67, -1.33, -2.26, -2.56)


    ## Standard deviation values from the prior log-normal distribution for the between-study variance (Turner et al.)
    sd_value <- c(1.34, 1.74, 1.74, 1.74, 1.74, 1.74, 1.74, 1.53, 1.52, 1.52, 1.51, 1.52, 1.54, 1.53, 1.52, 1.52, 1.74,
                  1.41, 1.79, 1.79, 1.79, 1.79, 1.79, 1.79, 1.58, 1.58, 1.58, 1.58, 1.58, 1.60, 1.60, 1.58, 1.58, 1.74,
                  1.55, 1.91, 1.91, 1.91, 1.91, 1.91, 1.92, 1.72, 1.71, 1.71, 1.71, 1.71, 1.73, 1.72, 1.71, 1.71, 1.74,
                  1.02, 1.50, 1.51, 1.50, 1.50, 1.51, 1.51, 1.25, 1.24, 1.24, 1.24, 1.25, 1.27, 1.27, 1.24, 1.25, 1.74,
                  1.26, 1.68, 1.68, 1.68, 1.68, 1.68, 1.68, 1.46, 1.45, 1.45, 1.45, 1.45, 1.47, 1.47, 1.45, 1.45, 1.74)


    ## Median value (in *standard deviation* scale) based on its prior (log-normal) distribution
    median_value <- sqrt(qlnorm(0.5, mean_value, sqrt(1 / sd_value)))


    ## Bring hyperparameters wtith distribution type into a cross-sectional table
    table_priors_text <-
      matrix(paste0("LN(", sprintf("%.2f", mean_value), ", ", sprintf("%.2f", sd_value), "\u00b2)"),
             nrow = length(outcome_type),
             ncol = length(interv_compat_type),
             byrow = FALSE)
    rownames(table_priors_text) <- outcome_type
    colnames(table_priors_text) <- interv_compat_type


    ## Bring median values into a cross-sectional table
    table_median <-
      matrix(median_value,
             nrow = length(outcome_type),
             ncol = length(interv_compat_type),
             byrow = FALSE)
    rownames(table_median) <- outcome_type
    colnames(table_median) <- interv_compat_type


    ## The empirical priors do *not* refer to a specific medical area
    caption <- "(doi: 10.1002/sim.6381)"


    ## Define the colour of the upper bound of median values
    upper_level <- "#F0E442"


    ## Publication source
    author <- "Turner et al. (2015)"

  } else if (measure == "SMD" & area == "cancer") {
    outcome_type <- c("Obstetric outcomes",
                      "Resource use/hospital stay/process",
                      "Internal/external structure-related outcomes",
                      "General physical health indicators",
                      "Adverse events",
                      "Pain",
                      "Quality of life/functioning (dichotomised)",
                      "Infection/onset of new disease",
                      "Signs/symptoms reflecting continuation/end of condition",
                      "Mental health indicators",
                      "Biological markers (dichotomised)",
                      "Subjective outcomes (various)",
                      "General healthcare setting")


    ## Treatment-comparison types as suggested by Rhodes et al.
    interv_compat_type <- c("Pharma vs. Placebo/control",
                            "Pharma vs. Pharma",
                            "Non-pharma vs. Placebo/control",
                            "Non-pharma vs. Pharma",
                            "Non-pharma vs. Non-pharma")


    ## Mean values from the prior t-distribution for the *log* between-study variance (Rhodes et al.)
    mean_value <- c(-1.57, -0.01, -0.13, -0.60, -0.60, -0.60, -0.60, -0.44, -0.44, -0.43, -0.85, -0.20, -3.44,
                    -1.85, -0.27, -0.14, -0.88, -0.88, -0.88, -0.88, -0.71, -0.71, -0.71, -1.13, -0.48, -3.44,
                    -1.43,  0.15,  0.27, -0.46, -0.46, -0.46, -0.46, -0.30, -0.30, -0.29, -0.71, -0.06, -3.44,
                    -1.43,  0.15,  0.27, -0.46, -0.46, -0.46, -0.46, -0.30, -0.30, -0.29, -0.71, -0.06, -3.44,
                    -1.43,  0.15,  0.27, -0.46, -0.46, -0.46, -0.46, -0.30, -0.30, -0.29, -0.71, -0.06, -3.44)


    ## Standard deviation values from the prior t-distribution for the *log* between-study variance (Rhodes et al.)
    sd_value <- c(2.45, 2.83, 2.61, 2.61, 2.61, 2.61, 2.61, 2.60, 2.60, 2.28, 2.93, 2.68, 2.59,
                  2.41, 2.79, 2.56, 2.55, 2.55, 2.55, 2.55, 2.57, 2.57, 2.25, 2.87, 2.68, 2.59,
                  2.24, 2.68, 2.45, 2.40, 2.40, 2.40, 2.40, 2.46, 2.46, 2.08, 2.78, 2.53, 2.59,
                  2.24, 2.68, 2.45, 2.40, 2.40, 2.40, 2.40, 2.46, 2.46, 2.08, 2.78, 2.53, 2.59,
                  2.24, 2.68, 2.45, 2.40, 2.40, 2.40, 2.40, 2.46, 2.46, 2.08, 2.78, 2.53, 2.59)


    ## Median value (in *standard deviation* scale) based on its prior (t) distribution
    median_value <- sqrt(exp(qt(0.5, df = 5) * sqrt(1 / sd_value) + mean_value))


    ## Bring hyperparameters wtith distribution type into a cross-sectional table
    table_priors_text <-
      matrix(paste0("t(", sprintf("%.2f", mean_value), ", ", sprintf("%.2f", sd_value), "\u00b2, 5)"),
             nrow = length(outcome_type),
             ncol = length(interv_compat_type),
             byrow = FALSE)
    rownames(table_priors_text) <- outcome_type
    colnames(table_priors_text) <- interv_compat_type


    ## Bring median values into a cross-sectional table
    table_median <-
      matrix(median_value,
             nrow = length(outcome_type),
             ncol = length(interv_compat_type),
             byrow = FALSE)
    rownames(table_median) <- outcome_type
    colnames(table_median) <- interv_compat_type


    ## The empirical priors refer to a specif medical area (cancer)
    caption <- "Medical areas of cancer (doi: 10.1016/j.jclinepi.2014.08.012)"


    ## Define the colour of the upper bound of median values
    upper_level <- "#009E73"


    ## Publication source
    author <- "Rhodes et al. (2015)"

  } else if (measure == "SMD" & area == "respiratory") {
    outcome_type <- c("Obstetric outcomes",
                      "Resource use/hospital stay/process",
                      "Internal/external structure-related outcomes",
                      "General physical health indicators",
                      "Adverse events",
                      "Pain",
                      "Quality of life/functioning (dichotomised)",
                      "Infection/onset of new disease",
                      "Signs/symptoms reflecting continuation/end of condition",
                      "Mental health indicators",
                      "Biological markers (dichotomised)",
                      "Subjective outcomes (various)",
                      "General healthcare setting")


    ## Treatment-comparison types as suggested by Rhodes et al.
    interv_compat_type <- c("Pharma vs. Placebo/control",
                            "Pharma vs. Pharma",
                            "Non-pharma vs. Placebo/control",
                            "Non-pharma vs. Pharma",
                            "Non-pharma vs. Non-pharma")


    ## Mean values from the prior t-distribution for the *log* between-study variance (Rhodes et al.)
    mean_value <- c(-6.03, -4.46, -4.33, -5.07, -5.07, -5.07, -5.07, -4.90, -4.90, -4.90, -5.31, -4.66, -3.44,
                    -6.31, -4.73, -4.61, -5.34, -5.34, -5.34, -5.34, -5.18, -5.18, -5.17, -5.59, -4.94, -3.44,
                    -5.89, -4.32, -4.19, -4.93, -4.93, -4.93, -4.93, -4.76, -4.76, -4.76, -5.17, -4.52, -3.44,
                    -5.89, -4.32, -4.19, -4.93, -4.93, -4.93, -4.93, -4.76, -4.76, -4.76, -5.17, -4.52, -3.44,
                    -5.89, -4.32, -4.19, -4.93, -4.93, -4.93, -4.93, -4.76, -4.76, -4.76, -5.17, -4.52, -3.44)


    ## Standard deviation values from the prior t-distribution for the *log* between-study variance (Rhodes et al.)
    sd_value <- c(2.36, 2.74, 2.51, 2.51, 2.51, 2.51, 2.51, 2.50, 2.50, 2.17, 2.83, 2.59, 2.59,
                  2.31, 2.70, 2.46, 2.45, 2.45, 2.45, 2.45, 2.47, 2.47, 2.14, 2.78, 2.59, 2.59,
                  2.21, 2.57, 2.33, 2.28, 2.28, 2.28, 2.28, 2.33, 2.33, 1.94, 2.66, 2.41, 2.59,
                  2.21, 2.57, 2.33, 2.28, 2.28, 2.28, 2.28, 2.33, 2.33, 1.94, 2.66, 2.41, 2.59,
                  2.21, 2.57, 2.33, 2.28, 2.28, 2.28, 2.28, 2.33, 2.33, 1.94, 2.66, 2.41, 2.59)


    ## Median value (in *standard deviation* scale) based on its prior (t) distribution
    median_value <- sqrt(exp(qt(0.5, df = 5) * sqrt(1 / sd_value) + mean_value))


    ## Bring hyperparameters wtith distribution type into a cross-sectional table
    table_priors_text <-
      matrix(paste0("t(", sprintf("%.2f", mean_value), ", ", sprintf("%.2f", sd_value), "\u00b2, 5)"),
             nrow = length(outcome_type),
             ncol = length(interv_compat_type),
             byrow = FALSE)
    rownames(table_priors_text) <- outcome_type
    colnames(table_priors_text) <- interv_compat_type


    ## Bring median values into a cross-sectional table
    table_median <-
      matrix(median_value,
             nrow = length(outcome_type),
             ncol = length(interv_compat_type),
             byrow = FALSE)
    rownames(table_median) <- outcome_type
    colnames(table_median) <- interv_compat_type


    ## The empirical priors refer to a specif medical area (respiratory)
    caption <- "Medical areas of respiratory diseases (doi: 10.1016/j.jclinepi.2014.08.012)"


    ## Define the colour of the upper bound of median values
    upper_level <- "#56B4E9"


    ## Publication source
    author <- "Rhodes et al. (2015)"

  } else if (measure == "SMD" & area == "other") {
    outcome_type <- c("Obstetric outcomes",
                      "Resource use/hospital stay/process",
                      "Internal/external structure-related outcomes",
                      "General physical health indicators",
                      "Adverse events",
                      "Pain",
                      "Quality of life/functioning (dichotomised)",
                      "Infection/onset of new disease",
                      "Signs/symptoms reflecting continuation/end of condition",
                      "Mental health indicators",
                      "Biological markers (dichotomised)",
                      "Subjective outcomes (various)",
                      "General healthcare setting")


    ## Treatment-comparison types as suggested by Rhodes et al.
    interv_compat_type <- c("Pharma vs. Placebo/control",
                            "Pharma vs. Pharma",
                            "Non-pharma vs. Placebo/control",
                            "Non-pharma vs. Pharma",
                            "Non-pharma vs. Non-pharma")


    ## Mean values from the prior t-distribution for the *log* between-study variance (Rhodes et al.)
    mean_value <- c(-4.13, -2.55, -2.43, -3.16, -3.16, -3.16, -3.16, -3.00, -3.00, -2.99, -3.41, -2.76, -3.44,
                    -4.40, -2.83, -2.70, -3.44, -3.44, -3.44, -3.44, -3.27, -3.27, -3.27, -3.68, -3.03, -3.44,
                    -3.99, -2.41, -2.29, -3.02, -3.02, -3.02, -3.02, -2.86, -2.86, -3.85, -3.27, -2.62, -3.44,
                    -3.99, -2.41, -2.29, -3.02, -3.02, -3.02, -3.02, -2.86, -2.86, -3.85, -3.27, -2.62, -3.44,
                    -3.99, -2.41, -2.29, -3.02, -3.02, -3.02, -3.02, -2.86, -2.86, -3.85, -3.27, -2.62, -3.44)


    ## Standard deviation values from the prior t-distribution for the *log* between-study variance (Rhodes et al.)
    sd_value <- c(2.34, 2.73, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.16, 2.83, 2.58, 2.59,
                  2.31, 2.70, 2.46, 2.44, 2.44, 2.44, 2.44, 2.47, 2.47, 2.14, 2.78, 2.59, 2.59,
                  2.11, 2.57, 2.32, 2.27, 2.27, 2.27, 2.27, 2.33, 2.33, 1.93, 2.66, 2.41, 2.59,
                  2.11, 2.57, 2.32, 2.27, 2.27, 2.27, 2.27, 2.33, 2.33, 1.93, 2.66, 2.41, 2.59,
                  2.11, 2.57, 2.32, 2.27, 2.27, 2.27, 2.27, 2.33, 2.33, 1.93, 2.66, 2.41, 2.59)


    ## Median value (in *standard deviation* scale) based on its prior (t) distribution
    median_value <- sqrt(exp(qt(0.5, df = 5) * sqrt(1 / sd_value) + mean_value))


    ## Bring hyperparameters wtith distribution type into a cross-sectional table
    table_priors_text <-
      matrix(paste0("t(", sprintf("%.2f", mean_value), ", ", sprintf("%.2f", sd_value), "\u00b2, 5)"),
             nrow = length(outcome_type),
             ncol = length(interv_compat_type),
             byrow = FALSE)
    rownames(table_priors_text) <- outcome_type
    colnames(table_priors_text) <- interv_compat_type


    ## Bring median values into a cross-sectional table
    table_median <-
      matrix(median_value,
             nrow = length(outcome_type),
             ncol = length(interv_compat_type),
             byrow = FALSE)
    rownames(table_median) <- outcome_type
    colnames(table_median) <- interv_compat_type


    ## The empirical priors refer to a specif medical area (other)
    caption <- "Medical areas other than cancer and respiratory diseases (doi: 10.1016/j.jclinepi.2014.08.012)"


    ## Define the colour of the upper bound of median values
    upper_level <- "#D55E00"


    ## Publication source
    author <- "Rhodes et al. (2015)"

  }


  ## Prepare dataset for ggplot2
  # Median values
  comparisons <- NULL
  dataset_median <- melt(table_median)
  colnames(dataset_median) <- c("outcome", "comparisons", "value")

  # Distributions
  distribution <- NULL
  dataset_text <- melt(table_priors_text)
  colnames(dataset_text) <- c("outcome", "comparisons", "distribution")


  ## Create the plot
  plot <-
    ggplot(dataset_median,
           aes(x = comparisons,
               y = outcome,
               fill = value)) +
    geom_tile(colour = "white") +
    geom_text(data = dataset_text,
              aes(x = comparisons,
                  y = outcome,
                  label = distribution),
              size = rel(4.5),
              inherit.aes = FALSE) +
    scale_fill_gradientn(colours =  c("white", upper_level),
                         guide = "none",
                         limits = c(min(dataset_median$value),
                                    max(dataset_median$value))) +
    scale_x_discrete(position = "top") +
    labs(x = "Treatment-comparison type",
         y = "Outcome type",
         caption = caption) +
    ggtitle(paste("Predictive distributions for between-study variance based on",
                  author)) +
    theme_bw() +
    theme(legend.position = "none",
          plot.title = element_text(size = 13, face = "bold", colour = "black"),
          axis.title = element_text(size = 12, face = "bold", colour = "black"),
          axis.text = element_text(size = 12),
          plot.caption = element_text(size = 12, face = "italic",hjust = 0.001))

  return(plot)
}
