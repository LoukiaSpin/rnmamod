#' League heatmap for estimation
#'
#' @description
#'   For one outcome, it creates a heatmap of the estimated effect measure for
#'   all possible comparisons of interventions in the network.
#'   For two outcomes, the heatmap illustrates these two outcomes for the same
#'   effect measure in the upper and lower off-diagonals for all
#'   possible comparisons of interventions in the network.
#'   The function can also be used to illustrate the results of two different
#'   models on the same outcome and effect measure.
#'   \code{league_heatmap} can be used for a random-effects or fixed-effect
#'   network meta-analysis, network meta-regression, and series of pairwise
#'   meta-analyses.
#'
#' @param full1 An object of S3 class \code{\link{run_model}} for network
#'   meta-analysis, or \code{\link{run_metareg}} for network meta-regression.
#'   See 'Value' in \code{\link{run_model}} and \code{\link{run_metareg}}.
#' @param full2 An object of S3 class \code{\link{run_model}} for network
#'   meta-analysis, \code{\link{run_metareg}} for network meta-regression, or
#'   \code{\link{run_series_meta}} for a series of pairwise meta-analyses.
#'   See 'Value' in \code{\link{run_model}}, \code{\link{run_metareg}}, and
#'   \code{\link{run_series_meta}}.
#' @param cov_value A list of two elements in the following order: a number
#'   for the covariate value of interest and a character for the name of the
#'   covariate. See also 'Details'.
#' @param drug_names1 A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}} for \code{full1}.
#' @param drug_names2 A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}} for \code{full2}. The elements must be a subset of
#'   \code{drug_names1}.
#' @param name1 The text for the title of the results that refer to
#'   the outcome or model under \code{full1}.
#' @param name2 The text for the title of the results that refer to
#'   the outcome or model under \code{full2}.
#' @param show A vector of at least three character strings that refer to the
#'   names of the interventions \emph{exactly} as defined in \code{drug_names1}.
#'   Then, the league table will be created for these interventions only.
#'   If \code{show} is not defined, the league table will present all
#'   interventions as defined in \code{drug_names1}.
#'
#' @return A heatmap of the league table showing the posterior mean and 95\%
#'   credible interval of the comparisons in the off-diagonals, and
#'   the posterior mean of the SUCRA values in the diagonal.
#'
#' @details \code{heatmap_league} offers the following options to display
#'   \bold{one} estimated effect measure for all (or some) pairwise comparisons:
#'   \itemize{
#'    \item one outcome, with results in the lower triangle referring to
#'    comparisons in the opposite direction after converting negative values
#'    into positive values (in absolute or logarithmic scale), and vice versa.
#'    Comparisons between interventions should be read from left to right.
#'    Therefore, each cell refers to the corresponding row-defining intervention
#'    against the column-defining intervention.
#'    Results that indicate strong evidence in favour of the
#'    row-defining intervention (i.e. the respective 95\% credible interval does
#'    not include the null value) are indicated in bold. A message is printed on
#'    the R console on how to read the heatmap;
#'    \item two outcomes for the same model, namely, network meta-analysis (via
#'    \code{\link{run_model}}) or network meta-regression (via
#'    \code{\link{run_metareg}}).
#'    When one of the outcomes includes more interventions, the argument
#'    \code{full1} should be considered for that outcome.
#'    Comparisons between interventions should be read as follows: for the upper
#'    diagonal, each cell refers to the corresponding row-defining intervention
#'    against the column-defining intervention, and for the lower diagonal, each
#'    cell refers to the corresponding column-defining intervention against the
#'    row-defining intervention. Results that indicate strong evidence (i.e. the
#'    respective 95\% credible interval does not include the null value) are
#'    indicated in bold. A message is printed on the R console on how to read
#'    the heatmap;
#'    \item two models for the same outcome, namely, network meta-analysis
#'    versus network meta-regression, or network meta-analysis versus series of
#'    pairwise meta-analyses.
#'    The instructions to read the heatmap are in line with the previous point.
#'    A message is printed on the R console on how to read the heatmap.
#'   }
#'
#'   For a beneficial outcome, red favours the first intervention of the
#'   comparison, and blue favours the second intervention. For a harmful
#'   outcome, blue favours the first intervention of the comparison, and red
#'   favours the second intervention. The larger the treatment effect, the
#'   darker the colour shade.
#'
#'   The function displays the effect measure as inherited by the argument
#'   \code{full1}. For binary outcome, it can display the odds ratio,
#'   relative risk, and risk difference. See 'Details' in
#'   \code{\link{run_model}} for the relative risk, and risk difference.
#'   For continuous outcome, it can display the mean difference, standardised
#'   mean difference, and ratio of means. Odds ratios, relative risk and ratio
#'   of means are reported in the original scale after exponentiation of the
#'   logarithmic scale.
#'
#'   The rows and columns of the heatmap display the names of
#'   interventions  sorted by decreasing order from the best to the worst
#'   based on their SUCRA value (Salanti et al., 2011) for the outcome or model
#'   under the argument \code{full1}. The off-diagonals contain the posterior
#'   mean and 95\% credible interval of the effect measure (according to the
#'   argument \code{measure} as inherited in the argument \code{full1}) of the
#'   corresponding comparisons.
#'
#'   The main diagonal contains the posterior mean of SUCRA of the corresponding
#'   interventions when the arguments \code{full1} refers to the
#'   \code{\link{run_model}} function. When the arguments \code{full1}
#'   refers to the \code{\link{run_metareg}} function, the p-score
#'   (Ruecker and Schwarzer, 2015) is calculated for each intervention while
#'   taking into account the covariate value in the argument \code{cov_value}.
#'   P-score is the 'frequentist analogue to SUCRA'
#'   (Ruecker and Schwarzer, 2015).
#'
#'   In the case of network meta-regression, when the covariate is binary,
#'   specify in the second element of \code{cov_value} the name of the level
#'   for which the heatmap will be created.
#'
#'   \code{league_heatmap} can be used only for a network of interventions.
#'   In the case of two interventions, the execution of the function will be
#'   stopped and an error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}, {Chrysostomos Kalyvas},
#'   {Katerina Papadimitropoulou}
#'
#' @seealso \code{\link{run_metareg}}, \code{\link{run_model}},
#'   \code{\link{run_series_meta}}
#'
#' @references
#' Ruecker G, Schwarzer G. Ranking treatments in frequentist network
#' meta-analysis works without resampling methods.
#' \emph{BMC Med Res Methodol} 2015;\bold{15}:58.
#' doi: 10.1186/s12874-015-0060-8
#'
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries
#' for presenting results from multiple-treatment meta-analysis: an overview and
#' tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71.
#' doi: 10.1016/j.jclinepi.2010.03.016
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
#' league_heatmap(full1 = res,
#'                drug_names1 = interv_names)
#'
#' @export
league_heatmap <- function(full1,
                           full2 = NULL,
                           cov_value = NULL,
                           drug_names1,
                           drug_names2 = NULL,
                           name1 = NULL,
                           name2 = NULL,
                           show = NULL) {

  if (!is.element(full1$type, c("nma", "nmr")) || is.null(full1$type)) {
    stop("'full1' must be an object of S3 class 'run_model', or 'run_metareg'.",
         call. = FALSE)
  }

  if (!is.null(full2) & (!is.element(full2$type, c("nma", "nmr", "series")) ||
                                     is.null(full2$type))) {
    aa <- "'run_metareg', or 'run_series_meta'."
    stop(paste("'full2' must be an object of S3 class 'run_model',", aa),
         call. = FALSE)
  }

  # Both objects must refer to the same effect measure
  measure <- if (is.null(full2) || (!is.null(full2) &
                                full1$measure == full2$measure)) {
    full1$measure
  } else if (!is.null(full2) & full1$measure != full2$measure) {
    stop("'full1' and 'full2' must have the same effect measure.", call. = FALSE)
  }

  # Forcing to define 'drug_names1' & 'drug_names2' so that 'show' can be used
  drug_names1 <- if (missing(drug_names1)) {
    stop("The argument 'drug_names1' has not been defined.", call. = FALSE)
  } else {
    drug_names1
  }

  if (length(drug_names1) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis.",
         call. = FALSE)
  }

  drug_names2 <- if (!is.null(full2) & is.null(drug_names2)) {
    stop("The argument 'drug_names2' has not been defined.", call. = FALSE)
  } else if (!is.null(full2) & !is.null(drug_names2)) {
    drug_names2
  }

  if (length(unique(is.element(drug_names2, drug_names1))) > 1) {
    stop("The argument 'drug_names2' must be a subset of 'drug_names1'.",
         call. = FALSE)
  }

  drug_names0 <- if (is.null(full2) ||
                     (!is.null(full2) &
                      length(drug_names1) >= length(drug_names2))) {
    drug_names1
  } else if (!is.null(full2) & length(drug_names1) < length(drug_names2)) {
    stop("'drug_names1' must have greater length than 'drug_names2'.",
         call. = FALSE)
  }

  name1 <- if (!is.null(full2) & is.null(name1)) {
    stop("The argument 'name1' has not been defined.", call. = FALSE)
  } else if (!is.null(full2) & !is.null(name1)) {
    name1
  }

  name2 <- if (!is.null(full2) & is.null(name2)) {
    stop("The argument 'name2' has not been defined.", call. = FALSE)
  } else if (!is.null(full2) & !is.null(name2)) {
    name2
  }

  show0 <- if (length(unique(!is.element(show, drug_names0))) > 1) {
    aa <- "All elements of the argument 'show' must be found in 'drug_names1'"
    bb <- "or 'drug_names2'."
    stop(paste(aa, bb), call. = FALSE)
  } else if (length(unique(!is.element(show, drug_names0))) == 1 &
             length(show) < 3) {
    stop("The argument 'show' must have length greater than 2.", call. = FALSE)
  } else if (length(unique(!is.element(show, drug_names0))) == 1 &
             length(show) > 2) {
    cbind(combn(show, 2)[2,], combn(show, 2)[1,])
  }

  drug_names <- if (is.null(show0)) {
    drug_names0
  } else {
    subset(drug_names0, is.element(drug_names0, show))
  }

  if (is.null(full2)) {
    message("Tips to read the table: row versus column.")
  } else {
    aa <- "Tips to read the table: upper triangle, row versus column;"
    bb <- "lower triangle, column versus row."
    message(paste(aa, bb))
  }

  #Source: https://rdrr.io/github/nfultz/stackoverflow/man/reflect_triangle.html
  reflect_triangle <- function(m, from = c("lower", "upper")) {
    ix <- switch(match.arg(from),
                 lower = upper.tri,
                 upper = lower.tri)(m, diag = FALSE)
    m[ix] <- t(-m)[ix]
    m
  }

  select <- cbind(combn(drug_names0, 2)[2, ], combn(drug_names0, 2)[1, ])
  ## Prepare the first outcome or model
  if (is.null(full1$beta_all)) {
    par <- if (is.null(show0)) {
      full1$EM
    } else {
      na.omit(subset(data.frame(full1$EM, select),
                     is.element(select[, 1], show) &
                       is.element(select[, 2], show)))
    }
    sucra <- if (is.null(show0)) {
      full1$SUCRA[, 1]
    } else {
      na.omit(subset(data.frame(full1$SUCRA[, 1], drug_names0),
                     is.element(drug_names0, show)))[, 1]
    }

  } else {
    cov_value <- if (missing(cov_value)) {
      stop("The argument 'cov_value' has not been defined.", call. = FALSE)
    } else if (length(cov_value) != 2) {
      aa <- "The argument 'cov_value' must be a list with elements a number"
      stop(paste(aa, "and a character."), call. = FALSE)
    } else if (length(cov_value) == 2) {
      cov_value
    }

    if (length(unique(full1$covariate)) < 3 &
        !is.element(cov_value[[1]], full1$covariate)) {
      aa <- "The first element of the argument 'cov_value' is out of the value"
      stop(paste(aa, "range of the analysed covariate."), call. = FALSE)
    } else if (length(unique(full1$covariate)) > 2 &
               (cov_value[[1]] < min(full1$covariate) |
                cov_value[[1]] > max(full1$covariate))) {
      aa <- "The first element of the argument 'cov_value' is out of the value"
      stop(paste(aa, "range of the analysed covariate."), call. = FALSE)
    }

    covar <- if (length(unique(full1$covariate)) < 3) {
      cov_value[[1]]
    } else {
      cov_value[[1]] - mean(full1$covariate)
    }

    par_mean <- full1$EM[, 1] + full1$beta_all[, 1] * covar
    par_sd <- sqrt(((full1$EM[, 2])^2) + ((full1$beta_all[, 2] * covar)^2))
    par_lower <- par_mean - 1.96 * par_sd
    par_upper <- par_mean + 1.96 * par_sd
    par0 <- data.frame(par_mean, par_sd, par_lower, full1$EM[, 4:6], par_upper)
    par <- if (is.null(show0)) {
      par0
    } else {
      na.omit(subset(data.frame(par0, select),
                     is.element(select[, 1], show) &
                       is.element(select[, 2], show)))
    }
    z_test <- par0[, 1] / par0[, 2]
    z_test_mat <- matrix(NA,
                         nrow = length(drug_names0),
                         ncol = length(drug_names0))
    z_test_mat[lower.tri(z_test_mat, diag = FALSE)] <- z_test * (-1)
    z_test_mat <- reflect_triangle(z_test_mat, from = "lower")
    prob_diff <- if (full1$D == 0) {
      pnorm(z_test_mat)
    } else {
      1 - pnorm(z_test_mat)
    }
    # The p-scores per intervention
    sucra0 <- apply(prob_diff, 1, sum, na.rm = TRUE) / (length(drug_names0) - 1)
    sucra <- if (is.null(show0)) {
      sucra0
    } else {
      na.omit(subset(data.frame(sucra0, drug_names0),
                     is.element(drug_names0, show)))[, 1]
    }
  }

  # Interventions order based on their SUCRA value (from best to worst)
  drug_order <- order(-sucra)

  # Order interventions based on their SUCRA value (from best to worst)
  order_drug <- drug_names[order(-sucra)]
  len_drug <- length(order_drug)

  # First: Matrix of effect measure for all possible comparisons
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

  # First: Symmetric matrix for effect measure and its bounds after ordering
  # rows and columns from the best to the worst intervention
  if (!is.element(measure, c("OR", "RR", "ROM"))) {
    point <- point1[drug_order, drug_order]
    lower <- lower1[drug_order, drug_order]
    upper <- upper1[drug_order, drug_order]

    # Spot the statistically significant comparisons (i.e. the 95% CrI does not
    # include the value of no difference)
    signif <- ifelse(upper < 0 | lower > 0, 1, 0)
    signif[is.na(signif)] <- 0
    signif_status <- melt(signif, na.rm = FALSE)[3]
  } else {
    point <- round(exp(point1[drug_order, drug_order]), 2)
    lower <- round(exp(lower1[drug_order, drug_order]), 2)
    upper <- round(exp(upper1[drug_order, drug_order]), 2)

    # Spot the statistically significant comparisons (i.e. the 95% CrI does not
    # include the value of no difference)
    signif <- ifelse(upper < 1 | lower > 1, 1, 0)
    signif[is.na(signif)] <- 1
    signif_status <- melt(signif, na.rm = FALSE)[3]
  }

  if (!is.null(full2) & length(full2$EM[1, ]) == 11) {
    names <- full2$EM[, 1:2]
    for (i in 1:length(names[, 1])) {
      names[i, 1] <- drug_names1[full2$EM[i, 2]]
      names[i, 2] <- drug_names1[full2$EM[i, 1]]
    }
  }

  # Prepare the second outcome or model
  if (!is.null(full2) & is.null(full2$beta_all)) {
    par2 <- if (is.null(show0) & length(full2$EM[1, ]) < 11) {
     full2$EM
    } else if (!is.null(show0) & length(full2$EM[1, ]) < 11) {
      na.omit(subset(data.frame(full2$EM, select),
                     is.element(select[, 1], show) &
                       is.element(select[, 2], show)))
    } else if (is.null(show0) & length(full2$EM[1, ]) == 11) {
      data.frame(full2$EM[, c(-1, -2)], names)
    } else if (!is.null(show0) & length(full2$EM[1, ]) == 11) {
      na.omit(subset(data.frame(full2$EM[, c(-1, -2)], names),
                     is.element(names[, 1], show) &
                       is.element(names[, 2], show)))
    }
  } else if (!is.null(full2) & !is.null(full2$beta_all)) {
    cov_value <- if (missing(cov_value)) {
      stop("The argument 'cov_value' has not been defined", call. = FALSE)
    } else if (length(cov_value) != 2) {
      aa <- "The argument 'cov_value' must be a list with elements a number"
      stop(paste(aa, "and a character"), call. = FALSE)
    } else if (length(cov_value) == 2) {
      cov_value
    }

    # Covariate value shall be the same with that for the first outcome or model
    if (length(unique(full2$covariate)) < 3 &
        !is.element(cov_value[[1]], full2$covariate)) {
      aa <- "The first element of the argument 'cov_value' is out of the value"
      stop(paste(aa, "range of the analysed covariate"), call. = FALSE)
    } else if (length(unique(full2$covariate)) > 2 &
               (cov_value[[1]] < min(full2$covariate) |
                cov_value[[1]] > max(full2$covariate))) {
      aa <- "The first element of the argument 'cov_value' is out of the value"
      stop(paste(aa, "range of the analysed covariate"), call. = FALSE)
    }

    covar <- if (length(unique(full2$covariate)) < 3) {
      cov_value[[1]]
    } else {
      cov_value[[1]] - mean(full2$covariate)
    }

    par_mean2 <- full2$EM[, 1] + full2$beta_all[, 1] * covar
    par_sd2 <- sqrt(((full2$EM[, 2])^2) + ((full2$beta_all[, 2] * covar)^2))
    par_lower2 <- par_mean2 - 1.96 * par_sd2
    par_upper2 <- par_mean2 + 1.96 * par_sd2
    par20 <- data.frame(par_mean2, par_sd2, par_lower2, full2$EM[, 4:6],
                        par_upper2)
    par2 <- if (is.null(show0)) {
      par20
    } else {
      na.omit(subset(data.frame(par20, select),
                     is.element(select[, 1], show) &
                       is.element(select[, 2], show)))
    }
  }

  # Second: Matrix of effect measure for all possible comparisons
  if (!is.null(full2) & length(full2$EM[1, ]) < 11) {
    comp0 <- t(combn(drug_names, 2))
    colnames(comp0) <- c("t1", "t2")
    comp <- cbind(comp0[, 2], comp0[, 1])
  } else if (!is.null(full2) & length(full2$EM[1, ]) == 11) {
    comp0 <- par2[, c("t1", "t2")]
    comp <- subset(comp0,
                   is.element(comp0[, 1], drug_names) &
                     is.element(comp0[, 2], drug_names))
  }

  reflect_triangle2 <- function(m, from = c("lower", "upper")) {
    ix <- switch(match.arg(from),
                 lower = upper.tri,
                 upper = lower.tri)(m, diag = FALSE)
    m[ix] <- t(m)[ix]
    m
  }

  if (!is.null(full2)) {
    point20 <- matrix(NA,
                      nrow = length(drug_names),
                      ncol = length(drug_names))
    lower20 <- upper20 <- point20
    rownames(point20) <- colnames(point20) <- drug_names
    rownames(lower20) <- colnames(lower20) <- drug_names
    rownames(upper20) <- colnames(upper20) <- drug_names
    for (i in 1:length(comp[, 1])) {
      point20[comp[i, 1], comp[i, 2]] <- round(par2[i, 1], 2)
      # Lower triangle
      lower20[comp[i, 1], comp[i, 2]] <- round(par2[i, 3], 2)
      upper20[comp[i, 1], comp[i, 2]] <- round(par2[i, 7], 2)
    }

    # Incorporate upper triangle
    point_21 <- reflect_triangle(point20, from = "lower")

    # Matrix of lower and upper bound of effect measure (all possible comparisons)
    # Incorporate upper triangle
    lower_210 <- reflect_triangle(upper20, from = "lower")
    upper_210 <- reflect_triangle(lower20, from = "lower")
    #lower_21[lower.tri(lower_21, diag = FALSE)] <- round(par2[, 3], 2)
    lower_21 <- lower_210
    upper_21 <- upper_210
    lower_21[lower.tri(lower_21, diag = FALSE)] <- upper_210[lower.tri(upper_210, diag = FALSE)]
    #upper_21[lower.tri(upper_21, diag = FALSE)] <- round(par2[, 7], 2)
    upper_21[lower.tri(upper_21, diag = FALSE)] <- lower_210[lower.tri(lower_210, diag = FALSE)]

    # Second: Symmetric matrix for effect measure and its bounds after ordering
    # rows and columns from the best to the worst intervention
    point02 <- point_21[drug_order, drug_order]
    lower02 <- lower_21[drug_order, drug_order]
    upper02 <- upper_21[drug_order, drug_order]
    if (!is.element(measure, c("OR", "RR", "ROM"))) {
      point2 <- reflect_triangle2(point02, from = "upper")
      lower2 <- reflect_triangle2(lower02, from = "upper")
      upper2 <- reflect_triangle2(upper02, from = "upper")

      # Spot the statistically significant comparisons (i.e. the 95% CrI does
      # not include the value of no difference)
      signif2 <- ifelse(upper2 < 0 | lower2 > 0, 1, 0)
      signif2[is.na(signif2)] <- 0
    } else {
      point2 <- round(exp(reflect_triangle2(point02, from = "upper")), 2)
      lower2 <- round(exp(reflect_triangle2(lower02, from = "upper")), 2)
      upper2 <- round(exp(reflect_triangle2(upper02, from = "upper")), 2)

      # Spot the statistically significant comparisons (i.e. the 95% CrI does not
      # include the value of no difference)
      signif2 <- ifelse(upper2 < 1 | lower2 > 1, 1, 0)
      signif2[is.na(signif2)] <- 1
    }
  }

  if (is.null(full2)) {
    point_f <- point
    lower_f <- lower
    upper_f <- upper
  } else {
    ## First outcome/model in upper, second outcome/model in lower diagonal
    point_f <- point2
    lower_f <- lower2
    upper_f <- upper2
    point_f[upper.tri(point_f, diag = FALSE)] <- point[upper.tri(point,
                                                                 diag = FALSE)]
    lower_f[upper.tri(lower_f, diag = FALSE)] <- lower[upper.tri(lower,
                                                                 diag = FALSE)]
    upper_f[upper.tri(upper_f, diag = FALSE)] <- upper[upper.tri(upper,
                                                                 diag = FALSE)]
    # Bring both into a matrix (Second: lower triangle; First: upper triangle)
    signif2[upper.tri(signif2, diag = FALSE)] <- signif[upper.tri(signif,
                                                                  diag = FALSE)]
    signif_status <- melt(signif2, na.rm = FALSE)[3]
  }

  # Merge point estimate with 95% credible interval in a new symmetric matrix
  final <- matrix(
    paste0(sprintf("%.2f", point_f),  "\n", "(",
           sprintf("%.2f", lower_f), ",", " ",
           sprintf("%.2f", upper_f), ")"),
    nrow = length(drug_names),
    ncol = length(drug_names))
  colnames(final) <- rownames(final) <- order_drug

  # Include SUCRA values in the diagonal of the new matrix
  diag(final) <- paste0(round(sort(sucra * 100, decreasing = TRUE), 1), "%")

  # Preparing the dataset for the ggplot2
  mat_new1 <- melt(final, na.rm = FALSE)

  # When is.null(show) replace the corresponding cells with "".
  mat_new1[, 3] <- ifelse(mat_new1[, 3] == "NA\n(NA, NA)", " ", mat_new1[, 3])

  # Merge both datasets to be used for ggplot2
  mat <- point_f
  diag(mat) <- ifelse(!is.element(measure, c("OR", "RR", "ROM")), 0, 1)
  mat_new <- cbind(mat_new1, melt(mat, na.rm = FALSE)[, 3])
  colnames(mat_new) <- c("Var1", "Var2", "value", "value2")

  caption0 <- if (!is.null(full1$beta_all) &
                 length(unique(full1$covariate)) > 2) {
    paste("Posterior mean of", effect_measure_name(measure, lower = TRUE),
          "(95% credible interval) for", cov_value[[2]], cov_value[[1]])
  } else if (!is.null(full1$beta_all) & length(unique(full1$covariate)) < 3) {
    paste("Posterior mean of", effect_measure_name(measure, lower = TRUE),
          "(95% credible interval) for", cov_value[[2]])
  } else if (is.null(full1$beta_all)) {
    paste("Posterior mean of", effect_measure_name(measure, lower = TRUE),
          "(95% credible interval)")
  }

  ## To create the orders of the lower diagonal
  xmin1 <- rep(seq(0.5, len_drug - 0.5, 1), each = len_drug)
  xmax1 <- xmin1 + 1
  ymin1 <- rep(seq(len_drug - 0.5, 0.5, -1), each = len_drug)
  ymax1 <- ymin1 + 1

  # Argument in scale_fill_gradientn
  min_value <- min(mat_new$value2, na.rm = TRUE)
  max_value <- max(mat_new$value2, na.rm = TRUE)
  values <- rescale(c(min_value,
                      ifelse(!is.element(measure,
                                         c("OR", "RR", "ROM")),
                             0.0001, 1.0001),
                      max_value))

  val <- if (!is.element(measure, c("OR", "RR", "ROM"))) {
    0
  } else {
    1
  }

  colours <- if (is.null(full2) ||
                 !is.null(full2) & min_value < val & max_value > val) {
   c("blue", "white", "#D55E00")
  } else if (!is.null(full2) & min_value < val & max_value == val) {
   c("blue", "white", "white")
  } else if (!is.null(full2) & min_value == val & max_value > val) {
    c("white", "white", "#D55E00")
  }

  # The league table as a heatmap
  p <- ggplot(mat_new,
              aes(factor(Var2, levels = order_drug[seq_len(len_drug)]),
                  factor(Var1, levels = order_drug[rev(seq_len(len_drug))]),
                  fill = value2)) +
    geom_tile() +
    geom_fit_text(aes(factor(Var2, levels = order_drug[seq_len(len_drug)]),
                      factor(Var1, levels =
                               order_drug[rev(seq_len(len_drug))]),
                      label = value),
                  reflow = TRUE) +
    geom_fit_text(aes(factor(Var2, levels = order_drug[seq_len(len_drug)]),
                      factor(Var1, levels =
                               order_drug[rev(seq_len(len_drug))]),
                      label = value,
                      fontface = ifelse(signif_status == 1, "bold", "plain")),
                  reflow = TRUE) +
    scale_fill_gradientn(colours = colours,
                         na.value = "grey70",
                         guide = "none",
                         values = values,
                         limits = c(min(mat_new$value2, na.rm = TRUE),
                                    max(mat_new$value2, na.rm = TRUE))) +
    geom_rect(aes(xmin = xmin1, xmax = xmax1, ymin = ymin1, ymax = ymax1),
              fill = "transparent", color = "black", size = 0.55) +
    scale_x_discrete(position = "top") +
    labs(x = name1, y = name2, caption = caption0) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 12, face = "bold",
                                      colour = "black"),
          axis.title.y = element_text(size = 12, face = "bold",
                                      colour = "black"),
          axis.text.x = element_text(size = 12, angle = 50, hjust = 0.0), #0.5
          axis.text.y = element_text(size = 12),
          plot.caption = element_text(hjust = 0.01))
  return(p)
}
