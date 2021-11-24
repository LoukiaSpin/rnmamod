#' Robustness index
#'
#' @description
#'   Calculates the robustness index, a novel index that quantifies the overall
#'   divergence of the sensitivity analysis results from the primary analysis
#'   results, and considers objective decision rules to infer the presence or
#'   lack of robustness of the primary analysis results when conducting a
#'   sensitivity analysis (Spineli et al., 2021).
#'
#' @param sens An object of S3 class \code{\link{run_sensitivity}} when
#'   sensitivity analysis refers to different scenarios about the average
#'   missingness parameter. See 'Value' in \code{\link{run_sensitivity}}. For a
#'   general sensitivity analysis, insert a list of three elements with the
#'   following order: 1) a row-bind matrix with the estimated summary effect
#'   measure, \code{EM}, obtained from the \code{link{run_model}} function for
#'   each re-analysis (the first \code{EM} should refer to the estimated summary
#'   effect measure under the primary analysis); 2) the effect measure with
#'   values \code{"OR"}, \code{"MD"}, \code{"SMD"}, or \code{"ROM"} for the odds
#'   ratio, mean difference, standardised mean difference and ratio of means,
#'   respectively; and 3) a character vector that indicates the name of the
#'   re-analyses employed.
#' @param threshold A number indicating the threshold of robustness, that is,
#'   the minimally allowed deviation between the primary analysis and
#'   re-analysis results. See 'Details' below.
#'
#' @return \code{robustness_index} prints on the R console a message in green
#'   text on the threshold of robustness determined by the user.
#'   Then, the function returns the following list of elements:
#'   \tabular{ll}{
#'    \code{robust_index} \tab A a numeric scalar or vector on the robustness
#'    index values. In the case of a pairwise meta-analysis (PMA),
#'    \code{robust_index} is scalar as only one summary effect size is obtained.
#'    In the case of network meta-analysis (NMA), \code{robust_index} is a
#'    vector with length equal to the number of possible pairwise comparisons;
#'    one robustness index per possible comparison.\cr
#'    \tab \cr
#'    \code{robust} \tab A character or character vector (of same length with
#'    \code{robust_index}) on whether the primary analysis results are
#'    \emph{robust} or \emph{frail} to the different re-analyses.\cr
#'    \tab \cr
#'    \code{kld} \tab A vector or matrix on the Kullback-Leibler divergence
#'    (KLD) measure in the summary effect size from a subsequent re-analysis to
#'    the primary analysis. In the case of a PMA, \code{kld} is a vector with
#'    length equal to the number of total analyses. The latter equals the square
#'    of the number of scenarios indicated in the argument \code{mean_scenarios}
#'    of \code{\link{run_sensitivity}}, in the case of missing participant
#'    outcome data, or the length of the character vector in argument
#'    \code{sens}. Therefore, one KLD value per analysis. In the case of NMA,
#'    \code{robust_index} is a matrix with number of rows equal to the number of
#'    total analyses and number of columns equal to the number of possible
#'    pairwise comparisons; one KLD value per analysis and possible
#'    comparison.\cr
#'   }
#'
#' @details The user may consider the values 0.28 and 0.17 in the argument
#'   \code{threshold} for binary and continuous outcome data (the default
#'   values), respectively, or consider other plausible values.
#'   Spineli et al., (2021) offers a discussion on specifying the
#'   \code{threshold} of robustness.
#'
#'   In the case of missing participant outcome data, the primary analysis is
#'   considered to be the middle of the numbers in the argument
#'   \code{mean_scenarios} of \code{\link{run_sensitivity}} (see 'Arguments'
#'   and 'Details' in \code{link{run_sensitivity}}).
#'
#'   In \code{robust}, the value \code{"robust"} appears when
#'   \code{robust_index} is less than \code{threshold}); otherwise, the value
#'   \code{"frail"} appears.
#'
#'   In the case of missing participant outcome data, \code{robustness_index}
#'   can be used only when missing participant outcome data have been
#'   extracted for at least one trial. Otherwise, the execution of the function
#'   will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_model}}, \code{\link{run_sensitivity}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of
#' primary analysis results: A case study on missing outcome data in pairwise
#' and network meta-analysis.
#' \emph{Res Synth Methods} 2021;\bold{12}(4):475--490.
#' \doi{10.1002/jrsm.1478}
#'
#' Kullback S, Leibler RA. On information and sufficiency.
#' \emph{Ann Math Stat} 1951;\bold{22}(1):79--86.
#'
#' @examples
#' data("nma.baker2009")
#'
#' # Read results from 'run_sensitivity' (using the default arguments)
#' res_sens <- readRDS(system.file('extdata/res_sens_baker.rds',
#'                     package = 'rnmamod'))
#'
#' # Calculate the robustness index
#' robustness_index(sens = res_sens,
#'                  threshold = 0.28)
#'
#' @export
robustness_index <- function(sens, threshold) {


  if (is.null(sens$EM)) {
    es_mat <- as.matrix(sens[[1]])
    measure <- sens[[2]]
    scenarios <- sens[[3]]
    n_scenar <- length(scenarios)
    primary_scenar <- 1
  } else {
    es_mat <- sens$EM
    measure <- sens$measure
    scenarios <- sens$scenarios
    n_scenar <- length(scenarios)^2
    primary_scenar <- median(seq_len(n_scenar))

    if (any(is.na(sens))) {
      stop("Missing participant outcome data have *not* been collected.
           This function cannot be used.", call. = FALSE)
    }
  }

  # The quadratic formula for the roots of the general quadratic equation
  nt <- (1 + sqrt(1 + 8 * (length(es_mat[, 1]) / sqrt(n_scenar)^2))) / 2
  poss_comp <- (nt * (nt - 1)) / 2

  if (missing(threshold) & is.element(measure, "OR")) {
    threshold <- 0.28
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt =
                       "The value 0.28 was assigned on 'threshold' by default",
                       "\033[0m", "\n")))
  } else if (missing(threshold) & is.element(measure, c("MD", "SMD", "ROM"))) {
    threshold <- 0.17
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt =
                       "The value 0.17 was assigned on 'threshold' by default",
                       "\033[0m", "\n")))
  } else {
    threshold <- threshold
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt = paste("The value", threshold,
                       "was assigned on 'threshold' for",
                                   effect_measure_name(measure)),
                       "\033[0m", "\n")))
  }

  # Function for the Kullback-Leibler Divergence (two normal distributions)
  kld_measure_univ <- function(mean_y, sd_y, mean_x, sd_x) {
    # x is the 'truth' (e.g. the MAR assumption)
    kld_xy <- 0.5 * (((sd_x / sd_y)^2) + ((mean_y - mean_x)^2)
                     / (sd_y^2) - 1 + 2 * log(sd_y / sd_x))

    return(kld_xy)
  }

  # A matrix of estimates for all possible comparisons under each scenario
  mean_mat <- matrix(rep(NA, length(es_mat[, 1])), nrow = n_scenar)
  sd_mat <- mean_mat
  for (i in 1:n_scenar) {
    for (j in 1:poss_comp) {
      mean_mat[i, j] <- es_mat[j + poss_comp * (i - 1), 1]
      sd_mat[i, j] <- es_mat[j + poss_comp * (i - 1), 2]
    }
  }

  kldxy <- list()
  robust_index <- rep(NA, poss_comp)
  for (i in 1:poss_comp) {
    kldxy[[i]] <- rep(NA, n_scenar)
    for (j in 1:n_scenar) {
      # Returns the KLD of informative scenario j compared with primary analysis
      # for comparison i
      kldxy[[i]][j] <- kld_measure_univ(mean_mat[j, i],
                                        sd_mat[j, i],
                                        mean_mat[primary_scenar, i],
                                        sd_mat[primary_scenar, i])
    }
    # This refers to the primary analysis
    kldxy[[i]][primary_scenar] <- 0

    # Returns the Robustness Index of comparison i across all
    # informative scenarios
    robust_index[i] <- sqrt(round(t(kldxy[[i]]) %*% kldxy[[i]], 3))
  }
  kld <- matrix(unlist(kldxy), ncol = n_scenar, byrow = TRUE)
  robust <- ifelse(robust_index < threshold, "robust", "frail")

  return(list(robust_index = robust_index,
              robust = robust,
              kld = kld,
              measure = measure,
              threshold = threshold,
              scenarios = scenarios))
}
