#' Scatterplot of SUCRA values
#'
#' @description Creates a scatterplot of the SUCRA values from the
#'   network meta-analysis and the network meta-regression for a specified level
#'   or value of the investigated covariate.
#'
#' @param full An object of S3 class \code{\link{run_model}}.
#'   See 'Value' in \code{\link{run_model}}.
#' @param reg An object of S3 class \code{\link{run_metareg}}. See 'Value' in
#'   \code{\link{run_metareg}}.
#' @param cov_value A list of two elements in the following order: 1) a number
#'   for the covariate value of interest (see 'Arguments' in
#'   \code{\link{run_metareg}}), and 2) a character to indicate the name of
#'   the covariate. See also 'Details'.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#'
#' @return A scatterplot of the SUCRA values under the network meta-analysis
#'   (y-axis) against the SUCRA values under the network meta-regression
#'   (x-axis) for a specified level or value of the investigated covariate.
#'
#' @details The names of the interventions appear above each point in the plot.
#'   Three coloured rectangles are drawn in the scatterplot: a red rectangle for
#'   SUCRA values up to 50\%, a yellow rectangular for SUCRA values between
#'   50\% and 80\%, and a green rectangle for SUCRA values over 80\%.
#'   Interventions falling at the green area are considered as the highest
#'   ranked interventions, whilst interventions falling at the red area are
#'   considered as the lowest ranked interventions.
#'
#'   When the covariate is binary, specify in the second element of
#'   \code{cov_value} the name of the level for which the scatterplot will be
#'   created.
#'
#'   \code{scatterplot_sucra} is integrated in \code{\link{metareg_plot}}.
#'
#'   \code{scatterplot_sucra} can be used only for a network of interventions.
#'   Otherwise, the execution of the function will be stopped and an error
#'   message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{metareg_plot}}, \code{\link{run_metareg}},
#'   \code{\link{run_model}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries
#' for presenting results from multiple-treatment meta-analysis: an overview and
#' tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71.
#' doi: 10.1016/j.jclinepi.2010.03.016
#'
#' @export
scatterplot_sucra <- function(full, reg, cov_value, drug_names) {

  if (length(unique(reg$covariate)) < 3 &
      !is.element(cov_value[[1]], reg$covariate)) {
    stop("The first element of the argument 'cov_value' is out of the value
         range of the analysed covariate", call. = FALSE)
  } else if (length(unique(reg$covariate)) > 2 &
             (cov_value[[1]] < min(reg$covariate) |
              cov_value[[1]] > max(reg$covariate))) {
    stop("The first element of the argument 'cov_value' is out of the value
         range of the analysed covariate", call. = FALSE)
  }

  drug_names <- if (missing(drug_names)) {
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt = "The argument 'drug_names' has not been defined.
                       The intervention ID, as specified in 'data' is used as
                       intervention names", "\033[0m", "\n")))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug_names
  }

  if (length(drug_names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis",
         call. = FALSE)
  }

  cov_value <- if (!is.null(reg$beta_all) & missing(cov_value)) {
    stop("The argument 'cov_value' has not been defined", call. = FALSE)
  } else if (!is.null(reg$beta_all) & length(cov_value) < 2) {
    stop("The argument 'cov_value' must be a list with elements a number and a
         character", call. = FALSE)
  } else if (!is.null(reg$beta_all) & length(cov_value) == 2) {
    cov_value
  }

  covar <- if (length(unique(reg$covariate)) < 3) {
    cov_value[[1]]
  } else {
    cov_value[[1]] - mean(reg$covariate)
  }

  sucra_full <- full$SUCRA

  #Source: https://rdrr.io/github/nfultz/stackoverflow/man/reflect_triangle.html
  reflect_triangle <- function(m, from = c("lower", "upper")) {
    ix <- switch(match.arg(from),
                 lower = upper.tri,
                 upper = lower.tri)(m, diag = FALSE)
    m[ix] <- t(-m)[ix]
    m
  }

  sucra_nma <- round(sucra_full[, 1] * 100, 0)
  par_mean <- reg$EM[, 1] + reg$beta_all[, 1] * covar
  par_sd <- sqrt(((reg$EM[, 2])^2) + ((reg$beta_all[, 2] * covar)^2))
  z_test <- par_mean / par_sd
  z_test_mat <- matrix(NA, nrow = length(drug_names), ncol = length(drug_names))
  z_test_mat[lower.tri(z_test_mat, diag = FALSE)] <- z_test * (-1)
  z_test_mat <- reflect_triangle(z_test_mat, from = "lower")
  prob_diff <- if (full$D == 0) {
    pnorm(z_test_mat)
  } else {
    1 - pnorm(z_test_mat)
  }
  # The p-scores per intervention
  sucra_nmr <- round(apply(prob_diff, 1, sum, na.rm = TRUE) /
                       (length(drug_names) - 1) * 100, 0)
  dataset <- data.frame(sucra_nma, sucra_nmr, drug_names)

  ggplot(dataset, aes(x = sucra_nmr, y = sucra_nma)) +
    geom_rect(aes(xmin = 80, xmax = 100, ymin = 80, ymax = 100),
              fill = "#009E73", alpha = 0.1) +
    geom_rect(aes(xmin = 50, xmax = 80, ymin = 50, ymax = 80),
              fill = "orange", alpha = 0.1) +
    geom_rect(aes(xmin = 0, xmax = 50, ymin = 0, ymax = 50),
              fill = "#D55E00", alpha = 0.1) +
    geom_point(size = 2, color = "blue") +
    geom_abline(lty = 1, size = 1, col = "grey60") +
    geom_text_repel(data = dataset,
                    aes(x = sucra_nmr, y = sucra_nma, label = drug_names),
                    color = "black",
                    fontface = "bold",
                    hjust = "right",
                    size = 3.8,
                    max.overlaps = Inf,
                    nudge_x = -0.1,
                    direction = "both") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    labs(x = "Network meta-regression SUCRA values (%)",
         y = "Network meta-analysis SUCRA values (%)",
         caption = paste("Results for",
                         cov_value[[2]],
                         ifelse(length(unique(reg$covariate)) < 3, " ",
                                cov_value[[1]]))) +
    theme_classic() +
    theme(axis.title = element_text(color = "black", face = "bold", size = 12),
          axis.text = element_text(color = "black", size = 12),
          plot.caption = element_text(color = "black", size = 11, hjust = 0.01))
}
