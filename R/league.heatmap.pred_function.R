#' League heatmap for prediction
#'
#' @description
#'   Creates a heatmap with the predicted effec size of all possible
#'   comparisons of interventions in the network. \code{league_heatmap_pred} can
#'    be used only for a random-effects network meta-analysis and network
#'    meta-regression. \code{league_heatmap_pred} is applied for one outcome
#'    only.
#'
#' @param full An object of S3 class \code{\link{run_model}} for network
#'   meta-analysis or \code{\link{run_metareg}} for network meta-regression.
#'   See 'Value' in \code{\link{run_model}} and \code{\link{run_metareg}}.
#' @param cov_value A list of two elements in the following order: a number
#'   for the covariate value of interest and a character for the name of the
#'   covariate. See also 'Details'.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#'
#' @return A league heatmap of the posterior mean and 95\% predictive interval
#'   of the effect measure (according to the argument \code{measure} defined in
#'   \code{\link{run_model}}) for all possible comparisons in the off-diagonals,
#'   and the posterior mean of the SUCRA values in the diagonal.
#'
#' @details The rows and columns of the heatmap display the names of
#'   interventions which are sorted by decreasing order from the best to the
#'   worst based on their SUCRA value (Salanti et al., 2011). The main diagonal
#'   contains the SUCRA values of the corresponding interventions when the
#'   argument \code{full} refers to the \code{\link{run_model}} function.
#'   When the argument \code{full} refers to the \code{\link{run_metareg}}
#'   function, the p-score (Ruecker and Schwarzer, 2015) is calculated for each
#'   intervention while taking into account the covariate value in
#'   the argument \code{cov_value}. P-score is the 'frequentist analogue to
#'   SUCRA' (Ruecker and Schwarzer, 2015).
#'
#'   Results in the lower triangle refer to comparisons in the opposite
#'   direction after converting negative values into positive values
#'   (in absolute or logarithmic scale), and vice versa.
#'   Darker shades of red and green correspond to larger treatment effects in
#'   the upper and lower triangle, respectively, for a beneficial outcome, and
#'   vice versa for a harmful outcome.
#'   Odds ratio and ratio of means are reported in the original scale after
#'   exponentiation of the logarithmic scale.
#'
#'   Comparisons between interventions should be read from left to right.
#'   Therefore, each cell refers to the corresponding row-defining intervention
#'   against the column-defining intervention. Results that indicate strong
#'   evidence in favour of the row-defining intervention (i.e. the respective
#'   95\% predictive interval does not include the null value) are indicated in
#'   bold.
#'
#'   In the case of network meta-regression, when the covariate is binary,
#'   specify in the second element of \code{cov_value} the name of the level for
#'   which the heatmap will be created.
#'
#'   \code{league_heatmap_pred} can be used only for a network of interventions.
#'   In the case of two interventions, the execution of the function will be
#'   stopped and an error message will be printed on the R console. Similarly,
#'   when the function is executed for a fixed-effect network meta-analysis or
#'   network meta-regression.
#'
#' @author {Loukia M. Spineli}, {Chrysostomos Kalyvas},
#'   {Katerina Papadimitropoulou}
#'
#' @seealso \code{\link{run_metareg}}, \code{\link{run_model}}
#'
#' @references
#' Ruecker G, Schwarzer G. Ranking treatments in frequentist network
#' meta-analysis works without resampling methods.
#' \emph{BMC Med Res Methodol} 2015;\bold{15}:58.
#' \doi{10.1186/s12874-015-0060-8}
#'
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries
#' for presenting results from multiple-treatment meta-analysis: an overview and
#' tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71.
#' \doi{10.1016/j.jclinepi.2010.03.016}
#'
#' @examples
#' data("nma.liu2013")
#'
#' # Read results from 'run_model' (using the default arguments)
#' res <- readRDS(system.file('extdata/res_liu.rds', package = 'rnmamod'))
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "pramipexole", "serotonin-norepinephrine
#'                   reuptake inhibitor", "serotonin reuptake inhibitor",
#'                   "tricyclic antidepressant", "pergolide")
#'
#' # Create the league heatmap
#' league_heatmap_pred(full = res,
#'                     drug_names = interv_names)
#'
#' @export
league_heatmap_pred <- function(full, cov_value = NULL, drug_names) {


  if (full$model == "FE") {
    stop("This function cannot be used for a fixed-effect model",
         call. = FALSE)
  }

  drug_names <- if (missing(drug_names)) {
    aa <- "The argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as specified in 'data' is used as"
    cc <- "intervention names"
    message(cat(paste0("\033[0;", col = 32, "m", aa, " ", bb, " ", cc,
                       "\033[0m", "\n")))
    as.character(seq_len(length(full$SUCRA[, 1])))
  } else {
    drug_names
  }

  if (length(drug_names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis",
         call. = FALSE)
  }

  measure <- full$measure

  #Source: https://rdrr.io/github/nfultz/stackoverflow/man/reflect_triangle.html
  reflect_triangle <- function(m, from = c("lower", "upper")) {
    ix <- switch(match.arg(from),
                 lower = upper.tri,
                 upper = lower.tri)(m, diag = FALSE)
    m[ix] <- t(-m)[ix]
    m
  }

  if (is.null(full$beta_all)) {
    par <- full$EM_pred
    sucra <- full$SUCRA[, 1]
  } else {

    cov_value <- if (missing(cov_value)) {
      stop("The argument 'cov_value' has not been defined", call. = FALSE)
    } else if (length(cov_value) != 2) {
      aa <- "The argument 'cov_value' must be a list with elements a number and"
      stop(paste(aa, "a character"), call. = FALSE)
    } else if (length(cov_value) == 2) {
      cov_value
    }

    if (length(unique(full$covariate)) < 3 &
        !is.element(cov_value[[1]], full$covariate)) {
      aa <- "The first element of the argument 'cov_value' is out of the value"
      stop(paste(aa, "range of the analysed covariate"), call. = FALSE)
    } else if (length(unique(full$covariate)) > 2 &
               (cov_value[[1]] < min(full$covariate) |
                cov_value[[1]] > max(full$covariate))) {
      aa <- "The first element of the argument 'cov_value' is out of the value"
      stop(paste(aa, "range of the analysed covariate"), call. = FALSE)
    }

    covar <- if (length(unique(full$covariate)) < 3) {
      cov_value[[1]]
    } else {
      cov_value[[1]] - mean(full$covariate)
    }

    par <- full$EM_pred
    par_mean <- full$EM_pred[, 1] + full$beta_all[, 1] * covar
    par_sd <- sqrt(((full$EM_pred[, 2])^2) + ((full$beta_all[, 2] * covar)^2))
    par_lower <- par_mean - 1.96 * par_sd
    par_upper <- par_mean + 1.96 * par_sd
    par <- data.frame(par_mean,
                      par_sd,
                      par_lower,
                      full$EM_pred[, 4:6],
                      par_upper)
    par_mean_em <- full$EM[, 1] + full$beta_all[, 1] * covar
    par_sd_em <- sqrt(((full$EM[, 2])^2) + ((full$beta_all[, 2] * covar)^2))
    z_test <- par_mean_em / par_sd_em
    z_test_mat <- matrix(NA,
                         nrow = length(drug_names),
                         ncol = length(drug_names))
    z_test_mat[lower.tri(z_test_mat, diag = FALSE)] <- z_test * (-1)
    z_test_mat <- reflect_triangle(z_test_mat, from = "lower")
    prob_diff <- if (full$D == 0) {
      pnorm(z_test_mat)
    } else {
      1 - pnorm(z_test_mat)
    }
    # The p-scores per intervention
    sucra <- apply(prob_diff, 1, sum, na.rm = TRUE) / (length(drug_names) - 1)
  }

  # Matrix of effect measure for all possible comparisons
  # Lower triangle
  point0 <- matrix(NA, nrow = length(drug_names), ncol = length(drug_names))
  lower0 <- upper0 <- point0
  point0[lower.tri(point0, diag = FALSE)] <- round(par[, 1], 2)
  # Incorporate upper triangle
  point1 <- reflect_triangle(point0, from = "lower")

  # Matrix of lower and upper bound of effect measure (all possible comparisons)
  # Lower triangle
  lower0[lower.tri(lower0, diag = FALSE)] <- round(par[, 3], 2)
  upper0[lower.tri(upper0, diag = FALSE)] <- round(par[, 7], 2)
  # Incorporate upper triangle
  lower1 <- reflect_triangle(upper0, from = "lower")
  lower1[lower.tri(lower1, diag = FALSE)] <- round(par[, 3], 2)
  upper1 <- reflect_triangle(lower0, from = "lower")
  upper1[lower.tri(upper1, diag = FALSE)] <- round(par[, 7], 2)

  # Interventions order based on their SUCRA value (from best to worst)
  drug_order_col <- drug_order_row <- order(-sucra)

  # Order interventions based on their SUCRA value (from best to worst)
  order_drug <- drug_names[order(-sucra)]
  len_drug <- length(order_drug)

  # Symmetric matrix for effect measure and its bounds after ordering rows and
  # columns from the best to the worst intervention
  if (measure != "OR" & measure != "ROM") {
    point <- point1[drug_order_col, drug_order_row]
    lower <- lower1[drug_order_col, drug_order_row]
    upper <- upper1[drug_order_col, drug_order_row]

    # Spot the statistically significant comparisons (i.e. the 95% CrI does not
    # include the value of no difference)
    (signif_status <- melt(ifelse(upper < 0 | lower > 0,
                                  "significant", "non-significant"),
                           na.rm = FALSE)[3])
    signif_status[is.na(signif_status)] <- 0
  } else {
    point <- round(exp(point1[drug_order_col, drug_order_row]), 2)
    lower <- round(exp(lower1[drug_order_col, drug_order_row]), 2)
    upper <- round(exp(upper1[drug_order_col, drug_order_row]), 2)

    # Spot the statistically significant comparisons (i.e. the 95% CrI does not
    # include the value of no difference)
    signif_status <- melt(ifelse(upper < 1 | lower > 1,
                                 "significant", "non-significant"),
                          na.rm = FALSE)[3]
    signif_status[is.na(signif_status)] <- 1
  }

  # Merge point estimate with 95% credible interval in a new symmetric matric
  final <- matrix(
    paste0(sprintf("%.2f", point),  "\n", "(", sprintf("%.2f", lower), ",", " ",
           sprintf("%.2f", upper), ")"),
    nrow = length(drug_names),
    ncol = length(drug_names))
  colnames(final) <- rownames(final) <- order_drug

  # Include SUCRA values in the diagonal of the new matrix
  diag(final) <- paste0(round(sort(sucra * 100, decreasing = TRUE), 1), "%")

  # Preparing the dataset for the ggplot2
  mat_new1 <- melt(final, na.rm = FALSE)

  # Merge both datasets to be used for ggplot2
  mat <- point
  mat_new <- cbind(mat_new1, melt(mat, na.rm = FALSE)[, 3])
  colnames(mat_new) <- c("Var1", "Var2", "value", "value2")

  # The final dataset for ggplot2
  diag(mat) <- NA
  final_col <- melt(mat)
  mat_new$value_sucra <- final_col$value

  caption <- if (!is.null(full$beta_all) & length(unique(full$covariate)) > 2) {
    paste("Posterior mean (95% predictive interval) for", cov_value[[2]],
          cov_value[[1]])
  } else if (!is.null(full$beta_all) & length(unique(full$covariate)) < 3) {
    paste("Posterior mean (95% predictive interval) for", cov_value[[2]])
  } else if (is.null(full$beta_all)) {
    paste("Posterior mean (95% predictive interval)")
  }

  # The league table as a heatmap
  p <- ggplot(mat_new,
              aes(factor(Var2, levels = order_drug[seq_len(len_drug)]),
                  factor(Var1, levels = order_drug[rev(seq_len(len_drug))]),
                  fill = value2)) +
    geom_tile(aes(fill = value_sucra)) +
    geom_fit_text(aes(factor(Var2, levels = order_drug[seq_len(len_drug)]),
                      factor(Var1, levels = order_drug[rev(seq_len(len_drug))]),
                      label = value), reflow = TRUE) +
    geom_fit_text(aes(factor(Var2, levels = order_drug[seq_len(len_drug)]),
                      factor(Var1, levels = order_drug[rev(seq_len(len_drug))]),
                      label = value,
                      fontface = ifelse(signif_status == "significant",
                                        "bold", "plain")),
                  reflow = TRUE) +
    scale_fill_gradientn(colours = c("#009E73", "white", "#D55E00"),
                         values = rescale(
                           c(min(mat_new$value2, na.rm = TRUE),
                             ifelse(measure != "OR" & measure != "ROM", 0, 1),
                             max(mat_new$value2, na.rm = TRUE))),
                         limits = c(min(mat_new$value2, na.rm = TRUE),
                                    max(mat_new$value2, na.rm = TRUE))) +
    scale_x_discrete(position = "top") +
    labs(x = "", y = "", caption = caption) +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.x = element_text(size = 12, angle = 50, hjust = 0.0),
          axis.text.y = element_text(size = 12),
          plot.caption = element_text(hjust = 0.01))

  return(p)
}
