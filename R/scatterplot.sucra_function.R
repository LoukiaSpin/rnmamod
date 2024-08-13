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
#' @param cov_name A character or text to indicate the name of the covariate.
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
scatterplot_sucra <- function(full,
                              reg,
                              cov_name = "covariate value",
                              drug_names) {

  if (!inherits(full, "run_model") || is.null(full)) {
    stop("'full' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  if (!inherits(reg, "run_metareg") || is.null(reg)) {
    stop("'reg' must be an object of S3 class 'run_metareg'.",
         call. = FALSE)
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

  sucra_full <- full$SUCRA
  sucra_reg <- reg$SUCRA
  sucra_nma <- round(sucra_full[, 1] * 100, 0)
  sucra_nmr <- round(sucra_reg[, 1] * 100, 0)

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
         caption = paste("Results for", cov_name)) +
         #caption = paste("Results for",
         #                cov_name,
         #                ifelse(length(unique(reg$covariate)) < 3, " ",
         #                       round(reg$cov_value, 2)))) +
    theme_classic() +
    theme(axis.title = element_text(color = "black", face = "bold", size = 12),
          axis.text = element_text(color = "black", size = 12),
          plot.caption = element_text(color = "black", size = 11, hjust = 0.01))
}
