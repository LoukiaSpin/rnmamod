#' Function for the Kullback-Leibler Divergence of two normally distributed
#' treatment effects for the same pairwise comparison
#'
#' @description
#' The user specify the (posterior) mean and standard error (or posterior
#' standard deviation) of two estimated treatment effects, X and Y, that refer
#' to the \strong{same} pairwise comparison and are assumed to follow a normal
#' distribution. The function returns the Kullback-Leibler Divergence (KLD)
#' measure of 1) approximating X with Y, 2) approximating Y with X, and 3) their
#' average.
#'
#' @param mean_y A real number that refers to the mean of the estimated
#'   treatment effect Y on the scale of the selected effect measure
#'   (in logarithmic scale for relative effect measures).
#' @param sd_y A positive integer that refers to the posterior standard
#'   deviation or the standard error of the estimated treatment effect Y on the
#'   scale of the selected effect measure (in logarithmic scale for relative
#'   effect measures).
#' @param mean_x A real number that refers to the mean of the estimated
#'   treatment effect X on the scale of the selected effect measure
#'   (in logarithmic scale for relative effect measures).
#' @param sd_x A positive integer that refers to the posterior standard
#'   deviation or the standard error of the estimated treatment effect X on the
#'   scale of the selected effect measure (in logarithmic scale for relative
#'   effect measures).
#'
#' @return
#' The function return the following numeric results:
#' \tabular{ll}{
#'    \strong{kld_sym} \tab The symmetric KLD value as the average of two KLD
#'    values .\cr
#'    \tab \cr
#'    \strong{kld_x_true} \tab The KLD value when approximating X by Y
#'    (X is the 'truth').\cr
#'    \tab \cr
#'    \strong{kld_y_true} \tab The KLD value when approximating Y by X
#'    (Y is the 'truth').\cr
#'   }
#'
#' @seealso \code{\link{kld_inconsistency}},
#'   \code{\link{kld_inconsistency_user}}, \code{\link{robustness_index}},
#'   \code{\link{robustness_index_user}}
#'
#' @references
#' Kullback S, Leibler RA. On information and sufficiency.
#' \emph{Ann Math Stat} 1951;\bold{22}(1):79--86. doi: 10.1214/aoms/1177729694
#'
#' @export
kld_measure <- function(mean_y, sd_y, mean_x, sd_x) {
  # x is the 'truth' (e.g. the direct estimate)
  kld_xy <- 0.5 * (((sd_x / sd_y)^2) + ((mean_y - mean_x)^2)
                   / (sd_y^2) - 1 + 2 * log(sd_y / sd_x))

  # y is the 'truth' (e.g. the indirect estimate)
  kld_yx <- 0.5 * (((sd_y / sd_x)^2) + ((mean_x - mean_y)^2)
                   / (sd_x^2) - 1 + 2 * log(sd_x / sd_y))

  # Symmetric KLD, also known as J-divergence
  sym_kld <- (kld_xy + kld_yx) / 2

  return(list(kld_sym = sym_kld,
              kld_x_true = kld_xy,
              kld_y_true = kld_yx))
}
