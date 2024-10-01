#' Calculate study percentage contributions to summary treatment effects or
#' regression coefficients
#'
#' @description
#' A data-frame on the percentage contributions of each study to every possible
#' pairwise comparison in the investigated network. Study percentage
#' contributions are based on the proposed methodology of Donegan and colleagues
#' (2018).
#'
#' @param study_name A vector of labels with the names of the studies included
#'   in the investigated network. For multi-arm studies, the study name should
#'   appear as many times as the number of possible comparisons among the
#'   compared treatments.
#' @param base_t A vector of numbers referring to the treatment identifier for
#'   the baseline arm (comparator) of each study.
#' @param exp_t A vector of numbers referring to the treatment identifier for
#'   the experimental arm of each study.
#' @param ref_t A scalar for the selected reference treatment in the network.
#' @param obs_se A vector of numbers referring to the estimated standard error
#'   of the treatment effect of each study. For multi-arm studies, the standard
#'   error of the treatment effect of each possible  comparison among the
#'   compared treatments should be included.
#' @param obs_cov A vector of numbers referring to the covariance in the block
#'    variance-covariance matrix of the estimated treatments effects for the
#'    multi-arm studies only. This argument should be left unspecified if there
#'    are no multi-arm studies in the network.
#' @param covar A vector of numbers referring to a continuous covariate that
#'   indicates a study characteristic or summary patient characteristic.
#' @param covar_assum Character string indicating the structure of the
#'   treatment-by-covariate interaction, as described in Cooper et al. (2009) if
#'   interest also lies on the study percentage contributions to the estimated
#'   regression coefficients. Set \code{covar_assumption} equal to \code{"no"},
#'   \code{"exchangeable"}, \code{"independent"}, or \code{"common"}. When
#'   \code{covar_assum = "no"}, only the study percentage contributions to the
#'   summary treatment effects will be calculated. There is no default argument.
#' @param model Character string indicating the analysis model with values
#'   \code{"RE"}, or \code{"FE"} for the random-effects and fixed-effect model,
#'   respectively. There is no default argument.
#' @param tau A scalar referring to the estimated between-study standard
#'   deviation obtained from network meta-analysis,
#'   if \code{covar_assum = "no"}, or network meta-regression for a specific
#'   treatment-by-covariate interaction assumption.  This argument should be
#'   left unspecified when \code{model = "FE"}.
#' @param tau_beta A scalar referring to the estimated standard
#'   deviation of the exchangeable regression coefficients obtained from
#'   network meta-regression with exchangeable treatment-by-covariate
#'   interaction. This argument should be left unspecified when
#'   \code{covar_assum} is not "exchangeable". There is no default argument.
#'
#' @return A list of the following two elements:
#'   \item{perc_contribute}{A data-frame with four columns referring to the
#'   study name, baseline and experimental treatment arm, and the covariate and
#'   as many columns as the number of possible comparisons with the study
#'   percentage contributions to summary treatment effects referring to the
#'   basic parameters, and followed by the functional parameters. If interest
#'   lies also on the regression coefficients, extra columns appear referring to
#'   the study percentage contributions to the regression coefficients specified
#'   from the argument \code{covar_assum}.}
#'   \item{covar_assumption}{The estimated summary odd ratio in the logarithmic scale when
#'   \code{measure = "RR"} or \code{measure = "RD"}.}
#'
#' @details
#' Note that the columns referring to the study percentage contributions to
#' summary treatment effects are indicated by the letter 'd' with two numbers in
#' decreasing order for the comparison: the first number refers to the
#' comparator and the second number refers to the experimental treatment of the
#' comparison. If interest lies also on the regression coefficients, the
#' correspoding columns are indicated by 'beta'.
#'
#' The function centers the covariate to the mean but presents the original
#' version of the covariate.
#'
#' @author {Loukia M. Spineli}
#'
#' @references
#' Cooper NJ, Sutton AJ, Morris D, Ades AE, Welton NJ. Addressing between-study
#' heterogeneity and inconsistency in mixed treatment comparisons: Application
#' to stroke prevention treatments in individuals with non-rheumatic atrial
#' fibrillation. \emph{Stat Med} 2009;\bold{28}(14):1861--81.
#' doi: 10.1002/sim.3594
#'
#' Donegan S, Dias S, Tudur-Smith C, Marinho V, Welton NJ. Graphs of study
#' contributions and covariate distributions for network meta-regression.
#' \emph{Res Synth Methods} 2018;\bold{9}(2):243--60. doi: 10.1002/jrsm.1292
#'
#' @examples
#' data("nma.malaria.donegan2018")
#'
#' # Get study contributions to fixed-effect network meta-regression
#' # results under the assumption of independent treatment-by-covariate
#' # interaction
#' study_perc_contrib(study_name = nma.malaria.donegan2018$s,
#'                    base_t = nma.malaria.donegan2018$t1,
#'                    exp_t = nma.malaria.donegan2018$t2,
#'                    ref_t = 1,
#'                    obs_se = nma.malaria.donegan2018$se,
#'                    covar = nma.malaria.donegan2018$x,
#'                    covar_assum = "independent",
#'                    model = "FE")
#'
#' @export
study_perc_contrib <- function (study_name,
                                base_t,
                                exp_t,
                                ref_t,
                                obs_se,
                                obs_cov = NULL,
                                covar,
                                covar_assum,
                                model,
                                tau = NULL,
                                tau_beta = NULL) {


  ## Default arguments
  # The covariate variable
  covar <- if(missing(covar)) {
    stop("'covar' must be defined", call. = FALSE)
  } else {
    covar
  }
  # Create a centered variable if 'covar' is metric
  covar_val <-
    if (length(unique(covar)) > 2) { # Center metric variable
    covar - mean(covar)
  } else if (length(unique(covar)) == 2) { # Binary
    covar
  }
  # Treatment-by-covariate Interaction assumption
  covar_assum <- if (missing(covar_assum)) {
    stop("'covar_assum' must be defined", call. = FALSE)
  } else if (!is.element(covar_assum,
                         c("no", "common", "exchangeable", "independent"))) {
    aa <- "'common', 'exchangeable', or 'independent'."
    stop(paste("'covar_assum' must be any of the following: 'no',", aa),
         call. = FALSE)
  } else {
    covar_assum
  }
  # The meta-analysis model
  model <- if (missing(model)) {
    stop("The argument 'model' must be specified.", call. = FALSE)
  } else if (!is.element(model, c("FE", "RE"))) {
    stop("Insert 'FE' or 'RE' for the argument model'.", call. = FALSE)
  } else {
    model
  }
  # The between-study standard deviation value
  tau <- if (model == "RE" & is.null(tau)) {
    stop("The argument 'tau' must be specified", call. = FALSE)
  } else {
    tau
  }
  # The between-study standard deviation value for beta (exchangeable interaction)
  tau_beta <- if (covar_assum == "exchangeable" & is.null(tau_beta)) {
    stop("The argument 'tau_beta' must be specified", call. = FALSE)
  } else {
    tau_beta
  }


  ## Design matrices X and Z
  desing_mat <-
    design_matrices(study_name, base_t, exp_t, ref_t, covar = covar_val, covar_assum)


  ## Find the unique treatments and sort in ascending order
  unique_treats <- sort(unique(c(cbind(base_t, exp_t))))

  if (model == "FE") {

    # Calculate the contribution matrix
    if (is.element(covar_assum, c("no", "common", "independent"))) {

      # *Observed* Variance-covariance matrix
      V_mat <- matrix(diag(obs_se^2), ncol = length(base_t))

      # Required contribution matrix
      contribute <-
        desing_mat$design_Z %*% (solve((t(desing_mat$design_X) %*% solve(V_mat)) %*% desing_mat$design_X)
                                 %*% t(desing_mat$design_X) %*% solve(V_mat))
    } else if (covar_assum == "exchangeable") {
      # Design matrix star
      design_X_star <- as.matrix(bdiag(desing_mat$design_X, rep(c(0, -1), each = length(unique_treats) - 1), 1))
      diag(design_X_star[(length(base_t) + 1):(length(base_t) + 2 * (length(unique_treats) - 1)),
                         1:(dim(desing_mat$design_X)[2] + 1)]) <- 1

      # *Observed* Variance-covariance matrix
      var_cov_obs <- matrix(diag(obs_se^2), ncol = length(base_t))

      # *Beta-effects* Variance-covariance matrix
      var_cov_beta <- rep(c(10^4, tau_beta^2), each = length(unique_treats) - 1) *
        diag(x = 1, nrow = 2 * (length(unique_treats) - 1))

      # Variance-covariance matrix star
      V_mat_star0 <- as.matrix(bdiag(var_cov_obs, var_cov_beta))
      V_mat_star <- cbind(rbind(V_mat_star0, rep(0, dim(V_mat_star0)[2])),
                          c(rep(0, dim(V_mat_star0)[1]), 10^4))
      # Contribution matrix A
      A_mat <- solve((t(design_X_star) %*% solve(V_mat_star)) %*% design_X_star) %*%
        t(design_X_star) %*% solve(V_mat_star)

      # Required contribution matrix (part of the matrix A)
      contribute <- desing_mat$design_Z %*% A_mat[1:(2 * (length(unique_treats) - 1)), 1:length(base_t)]
    }

  } else if (model == "RE") {

    # *Observed* Variance-covariance matrix (multi-arm trials addressed)
    if (!is.null(obs_cov)) {
      split_cov_obs <- lapply(split(obs_cov, study_name), function(x) unique(x))
      var_cov_obs <- lapply(split(cbind(obs_se^2, obs_cov), study_name),
                            function(x) diag(x[1:(length(x) / 2)], nrow = length(x) / 2))
      var_cov_obs_fin0 <-
        lapply(1:length(var_cov_obs),
               function(x) ifelse(row(var_cov_obs[[x]]) == col(var_cov_obs[[x]]), var_cov_obs[[x]], split_cov_obs[[x]]))
      var_cov_obs_fin <- as.matrix(bdiag(var_cov_obs_fin0))

      # *Random-effects* Variance-covariance matrix (multi-arm trials addressed)
      split_cov_tau <- lapply(split(ifelse(!is.na(obs_cov), (tau^2) / 2, NA), study_name), function(x) unique(x))
      var_cov_tau <- lapply(split(cbind(rep(tau^2, length(base_t)),
                                        ifelse(!is.na(obs_cov), (tau^2) / 2, NA)),
                                  study_name),
                            function(x) diag(x[1:(length(x) / 2)], nrow = length(x) / 2))
      var_cov_tau_fin0 <-
        lapply(1:length(var_cov_tau),
               function(x) ifelse(row(var_cov_tau[[x]]) == col(var_cov_tau[[x]]), var_cov_tau[[x]], split_cov_tau[[x]]))
      var_cov_tau_fin <- as.matrix(bdiag(var_cov_tau_fin0))

    } else {
      var_cov_obs_fin <- matrix(diag(obs_se^2), ncol = length(base_t))
      var_cov_tau_fin <- (tau^2) * diag(x = 1, nrow = length(base_t))
    }

    # Calculate the contribution matrix
    if (is.element(covar_assum, c("no", "common", "independent"))) {

      # Design matrix star
      design_X_star <- as.matrix(bdiag(diag(x = 1, nrow = length(base_t)), (-1) * desing_mat$design_X))
      diag(design_X_star[(length(base_t) + 1):(length(base_t) * 2), 1:length(base_t)]) <- 1

      # Variance-covariance matrix star
      V_mat_star <- as.matrix(bdiag(var_cov_obs_fin, var_cov_tau_fin))

    } else if (covar_assum == "exchangeable") {

      # Design matrix star
      design_X_star1 <- as.matrix(bdiag(diag(x = 1, nrow = length(base_t)), (-1) * desing_mat$design_X))
      diag(design_X_star1[(length(base_t) + 1):(length(base_t) * 2), 1:length(base_t)]) <- 1

      design_X_star2 <- rbind(cbind(matrix(0, nrow = 2 *(length(unique_treats) - 1), ncol = length(base_t)),
                              diag(x = 1, nrow = 2 *(length(unique_treats) - 1))),
                              c(rep(0, length(base_t)), rep(0, 2 *(length(unique_treats) - 1))))
      design_X_star3 <- c(rep(0, 2 * length(base_t)), rep(c(0, -1), each = length(unique_treats) - 1), 1)
      design_X_star <- cbind(rbind(design_X_star1, design_X_star2), design_X_star3) # It is correct! :-)

      # *Beta-effects* Variance-covariance matrix
      var_cov_beta <- rep(c(10^4, tau_beta^2), each = length(unique_treats) - 1) *
        diag(x = 1, nrow = 2 * (length(unique_treats) - 1))

      # Variance-covariance matrix star
      V_mat_star <- as.matrix(bdiag(var_cov_obs_fin, var_cov_tau_fin, var_cov_beta, 10^4)) # It is correct! :-)
    }

    # Contribution matrix A
    A_mat <- solve((t(design_X_star) %*% solve(V_mat_star)) %*% design_X_star) %*% t(design_X_star) %*% solve(V_mat_star)

    # Required contribution matrix (part of the matrix A)
    contribute <- if (is.element(covar_assum, c("independent", "exchangeable"))) {
      desing_mat$design_Z %*%
        A_mat[(length(base_t) + 1):(length(base_t) + (2 * (length(unique_treats) - 1))), 1:length(base_t)]
    } else if (covar_assum == "common") {
      desing_mat$design_Z %*%
        A_mat[(length(base_t) + 1):(length(base_t) + length(unique_treats)), 1:length(base_t)]
    } else if (covar_assum == "no") {
      desing_mat$design_Z %*%
        A_mat[(length(base_t) + 1):(length(base_t) + (length(unique_treats) - 1)), 1:length(base_t)]
    }

  }

  # Sum the absolute contributions by row (namely, across the studies)
  denom <- apply(abs(contribute), 1, sum)

  # Calculate the percentage contribution matrix :-)
  perc_contribute0 <- round((abs(contribute) / denom) * 100, 3)
  colnames(perc_contribute0) <- paste("study", study_name)

  # Transpose for readability
  perc_contribute <-
    data.frame(study_name = gsub(" ", "",
                                 sub("study", "", colnames(perc_contribute0))),
               comparator_arm = base_t,
               experimental_arm = exp_t,
               t(perc_contribute0),
               covariate = covar)

  results <- list(perc_contribute = perc_contribute,
                  covar_assumption = covar_assum)

  class(results) <- "study_perc_contrib"

  return(results)
}
