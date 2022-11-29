#' League table for relative and absolute effects (user defined)
#'
#' @description
#'   In line with \code{\link{league_table_absolute}}, provides a league table
#'   of the estimated odds ratio, and risk difference per 1000 participants for
#'   all possible comparisons of interventions in the network.
#'   The main diagonal of the table presents the absolute risk for each
#'   intervention in the network. \code{league_table_absolute_user} requires
#'   users to input the summary effect and 95\% credible or confidence interval
#'   of the basic parameters in the reported effect measure. This function
#'   should be used when the user has access to the results of a published
#'   systematic review rather than the raw trial-level data. In the latter case,
#'   the user should consider the function \code{\link{league_table_absolute}}.
#'   \code{league_table_absolute_user} is applied for one binary outcome only.
#'
#' @param data A data-frame with the summary effects of comparisons with the
#'   reference intervention of the network, known as basic parameters. The
#'   data-frame has \code{T} rows (\code{T} is the number of interventions in
#'   the network) and four columns that contain the point estimate, the lower
#'   and upper bound of the 95\% (confidence or credible) interval of the
#'   corresponding basic parameters, and a ranking measure to indicate the order
#'   of the interventions in the hierarchy from the best to the worst with
#'   possible choices a non-zero positive integer for the rank, the SUCRA value
#'   (Salanti et al., 2011) or p-score value (Ruecker and Schwarzer, 2015).
#'   The first row of the data-frame refers to the selected reference
#'   intervention and should include (1) the null value three times at the
#'   investigated effect measure (i.e. 1 for odds ratio and relative risk, and 0
#'   for risk difference), and (2) the value of the ranking measure.
#' @param measure Character string indicating the effect measure of \code{data}.
#'   For a binary outcome, the following can be considered: \code{"OR"},
#'   \code{"RR"} or \code{"RD"} for the odds ratio, relative risk, and risk
#'   difference, respectively.
#' @param base_risk A number in the interval (0, 1) that indicates the baseline
#'   risk for the selected reference intervention.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data}.
#'   The first intervention should be the selected reference intervention.
#' @param show A vector of at least three character strings that refer to the
#'   names of the interventions \emph{exactly} as defined in \code{drug_names}.
#'   Then, the league table will be created for these interventions only.
#'   If \code{show} is not defined, the league table will present all
#'   interventions as defined in \code{drug_names}.
#' @param save_xls Logical to indicate whether to export the tabulated results
#'   to an 'xlsx' file (via the \code{\link[writexl:write_xlsx]{write_xlsx}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}) to the working
#'   directory of the user. The default is \code{FALSE} (do not export).
#'
#' @return A league table showing the estimate and 95\% confidence interval of
#'   the odds ratio (upper off-diagonals), risk difference per 1000
#'   participants (lower off-diagonals), and absolute risks per 1000
#'   participants (main diagonal).
#'
#' @details
#'   When the published results are reported in the relative risk scale
#'   (i.e., \code{measure = "RR"}), the function calculates odds ratios and risk
#'   differences (point estimate and 95\% confidence interval) for all possible
#'   pairwise comparisons in the network based on the obtained absolute risks
#'   and the selected baseline risk. Likewise, when the published results are in
#'   the odds ratio or risk difference scale (i.e., \code{measure = "OR"} or
#'   \code{measure = "RD"}, respectively), the function calculates risk
#'   differences or odds ratios (point estimate and 95\% confidence interval),
#'   respectively, for all possible pairwise comparisons in the network based on
#'   the obtained absolute risks and the selected baseline risk.
#'
#'   The rows and columns of the league table display the names of the
#'   interventions  sorted by decreasing order from the best to the worst
#'   based on the ranking measure in the fourth column of the argument
#'   \code{data}. The upper off-diagonals contain the estimate and 95\%
#'   confidence interval of the odds ratio, the lower off-diagonals contain the
#'   estimate and 95\% confidence interval of the risk difference (per 1000
#'   participants), and the main diagonal comprises the absolute risks and their
#'   95\% confidence interval (per 1000 participants) of the corresponding
#'   non-reference interventions. The reference intervention of the network
#'   (which the baseline risk has been selected for) is indicated in the main
#'   diagonal with a black, thick frame.
#'
#'   Comparisons between interventions should be read from left to right.
#'   Results that indicate strong evidence in favour of the row-defining
#'   intervention (i.e. the respective 95\% confidence interval does not include
#'   the null value) are indicated in bold.
#'
#'   Furthermore, \code{league_table_absolute_user} exports
#'   \code{table_relative_absolute_effect}, a table with the relative and
#'   absolute effects of the basic parameters, as an 'xlsx' file (via the
#'   \code{\link[writexl:write_xlsx]{write_xlsx}} function) to the working
#'   directory of the user.
#'
#'   To obtain unique absolute risks for each intervention, we have considered
#'   the transitive risks framework, namely, an intervention has the same
#'   absolute risk regardless of the comparator intervention(s) in a trial
#'   (Spineli et al., 2017).
#'
#'   \code{league_table_absolute_user} can be used only for a network of
#'   interventions. In the case of two interventions, the execution of the
#'   function will be stopped and an error message will be printed in the R
#'   console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso code{\link{league_table_absolute}},
#'   \code{\link[writexl:write_xlsx]{write_xlsx}}
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
#' Spineli LM, Brignardello-Petersen R, Heen AF, Achille F, Brandt L,
#' Guyatt GH, et al. Obtaining absolute effect estimates to facilitate shared
#' decision making in the context of multiple-treatment comparisons.
#' Abstracts of the Global Evidence Summit, Cape Town, South Africa.
#' \emph{Cochrane Database of Systematic Reviews} 2017;\bold{9}(Suppl 1):1891.
#'
#' @export
league_table_absolute_user <- function(data,
                                       measure,
                                       base_risk,
                                       drug_names,
                                       show = NULL,
                                       save_xls) {

  data <- if (missing(data)) {
    stop("The argument 'data' needs to be defined", call. = FALSE)
  } else if (dim(data)[2] != 4) {
    stop("The argument 'data' must have four columns", call. = FALSE)
  } else {
    data
  }

  measure <- if (missing(measure)) {
    aa <- "Insert 'OR', 'RR', or 'RD'."
    stop(paste("The argument 'measure' needs to be defined.", aa),
         call. = FALSE)
  } else if (!is.element(measure, c("OR", "RR", "RD"))) {
    stop("Insert 'OR', 'RR', or 'RD'", call. = FALSE)
  } else {
    measure
  }

  base_risk <- if (missing(base_risk)) {
    stop("The argument 'base_risk' needs to be defined", call. = FALSE)
  } else if (base_risk <= 0 || base_risk >= 1) {
    stop("The argument 'base_risk' must be defined in (0, 1).", call. = FALSE)
  } else {
    base_risk
  }

  drug_names0 <- if (missing(drug_names)) {
    stop("The argument 'drug_names' has not been defined.", call. = FALSE)
  } else {
    drug_names
  }

  show0 <- if (length(unique(!is.element(show, drug_names0))) > 1) {
    stop("All elements of the argument 'show' must be found in 'drug_names'",
         call. = FALSE)
  } else if (length(unique(!is.element(show, drug_names0))) == 1 &
             length(show) < 3) {
    stop("The argument 'show' must have length greater than 2.", call. = FALSE)
  } else if (length(unique(!is.element(show, drug_names0))) == 1 &
             length(show) > 2) {
    cbind(combn(show, 2)[2,], combn(show, 2)[1,])
  }

  save_xls <- if (missing(save_xls)) {
    FALSE
  } else {
    save_xls
  }

  drug_names <- if (is.null(show0)) {
    drug_names0
  } else {
    subset(drug_names0, is.element(drug_names0, show))
  }

  if (length(drug_names0) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis",
         call. = F)
  } else {
    message("Tips to read the table: row versus column.")
  }

  # Function to obtain the absolute risks
  absol_risk_fun <- function (estimate, baseline, effect_measure) {

    if (measure == "OR") {
      round((estimate * baseline) / (1 + baseline * (estimate - 1)), 3)
    } else if (measure == "RR") {
      round(estimate * baseline, 3)
    } else {
      round((estimate + baseline), 3)
    }
  }

  # Obtain the absolute risks specific to the effect measure
  absol_risk <- absol_risk_fun(estimate = data[, -4],
                               baseline = base_risk,
                               effect_measure = measure)

  # Use basic parameters to calculate the functional parameters
  # via the consistency equation (point estimate and standard error)
  comb0 <- t(combn(1:length(drug_names0), 2))
  comb <- cbind(comb0[, 2], comb0[, 1])  # Intervention versus control
  #data_new <- cbind(data[, 1], (data[, 3] - data[, 2])/3.92) # Add the Std Error
  full_point <- full_se <- full_lower <- full_upper <- rep(NA, dim(comb)[1])
  for(i in 1:dim(comb)[1]) {
    if (is.element(measure, c("OR", "RR"))) {
      # Add the Std Error
      data_new <- cbind(data[, 1], (log(data[, 3]) - log(data[, 2]))/3.92)
      full_point[i] <- data_new[comb[i, 1], 1]/data_new[comb[i, 2], 1]
      full_se[i] <- sqrt((data_new[comb[i, 1], 2]^2) +
                           (data_new[comb[i, 2], 2]^2))
      full_lower[i] <- full_point[i] * exp(-1.96 * full_se[i])
      full_upper[i] <- full_point[i] * exp(1.96 * full_se[i])
    } else {
      # Add the Std Error
      data_new <- cbind(data[, 1], (data[, 3] - data[, 2])/3.92)
      full_point[i] <- data_new[comb[i, 1], 1] - data_new[comb[i, 2], 1]
      full_se[i] <- sqrt((data_new[comb[i, 1], 2]^2) +
                           (data_new[comb[i, 2], 2]^2))
      full_lower[i] <- full_point[i] - 1.96 * full_se[i]
      full_upper[i] <- full_point[i] + 1.96 * full_se[i]
    }
  }
  full <- cbind(full_point, full_lower, full_upper)
  full[1:(length(drug_names0) - 1), ] <- data[2:length(drug_names0), -4]

  # Obtain OR and RD based on the effect measure and absolute risks
  if (measure == "OR") {
    full_rd_point <- absol_risk_fun(full[, 1],
                                    absol_risk[comb[, 2], 1],
                                    measure) - absol_risk[comb[, 2], 1]
    full_rd_lower <- absol_risk_fun(full[, 2],
                                    absol_risk[comb[, 2], 2],
                                    measure) - absol_risk[comb[, 2], 2]
    full_rd_upper <- absol_risk_fun(full[, 3],
                                    absol_risk[comb[, 2], 3],
                                    measure) - absol_risk[comb[, 2], 3]
    full_rd <- cbind(full_rd_point, full_rd_lower, full_rd_upper)
    # Odds ratios (log scale; from published results)
    full_lor <- log(full)
  } else if (measure == "RR") {
    # Risk differences (as a function of odds ratios (full) and absolute risks)
    full_rd_point <- absol_risk_fun(full[, 1],
                                    absol_risk[comb[, 2], 1],
                                    measure) - absol_risk[comb[, 2], 1]
    full_rd_lower <- absol_risk_fun(full[, 2],
                                    absol_risk[comb[, 2], 2],
                                    measure) - absol_risk[comb[, 2], 2]
    full_rd_upper <- absol_risk_fun(full[, 3],
                                    absol_risk[comb[, 2], 3],
                                    measure) - absol_risk[comb[, 2], 3]
    full_rd <- cbind(full_rd_point, full_rd_lower, full_rd_upper)
    # Odds ratios (log scale; as a function of relative risks (full) and
    # the absolute risks)
    full_lor_point <- log((exp(full[, 1]) * absol_risk[comb[, 2], 1]) /
                            (1 + absol_risk[comb[, 2], 1] *
                               (exp(full[, 1]) - 1)))
    full_lor_lower <- log((exp(full[, 2]) * absol_risk[comb[, 2], 2]) /
                            (1 + absol_risk[comb[, 2], 2] *
                               (exp(full[, 2]) - 1)))
    full_lor_upper <- log((exp(full[, 3]) * absol_risk[comb[, 2], 3]) /
                            (1 + absol_risk[comb[, 2], 3] *
                               (exp(full[, 3]) - 1)))
    full_lor <- cbind(full_lor_point, full_lor_lower, full_lor_upper)
  } else if (measure == "RD") {
    # Risk differences (from published results)
    full_rd <- full
    # Relative risks (log scale; as a function of risk differences (full) and
    # the absolute risks)
    full_lrr_point <- (full[, 1] + absol_risk[comb[, 2], 1]) /
      absol_risk[comb[, 2], 1]
    full_lrr_lower <- (full[, 2] + absol_risk[comb[, 2], 2]) /
      absol_risk[comb[, 2], 2]
    full_lrr_upper <- (full[, 3] + absol_risk[comb[, 2], 3]) /
      absol_risk[comb[, 2], 3]
    full_lrr <- cbind(full_lrr_point, full_lrr_lower, full_lrr_upper)
    # Odds ratios (log scale; as a function of relative risks (full_lrr_X) and
    # the absolute risks)
    full_lor_point <- log((exp(full_lrr[, 1]) * absol_risk[comb[, 2], 1]) /
                            (1 + absol_risk[comb[, 2], 1] *
                               (exp(full_lrr[, 1]) - 1)))
    full_lor_lower <- log((exp(full_lrr[, 2]) * absol_risk[comb[, 2], 2]) /
                            (1 + absol_risk[comb[, 2], 2] *
                               (exp(full_lrr[, 2]) - 1)))
    full_lor_upper <- log((exp(full_lrr[, 3]) * absol_risk[comb[, 2], 3]) /
                            (1 + absol_risk[comb[, 2], 3] *
                               (exp(full_lrr[, 3]) - 1)))
    full_lor <- cbind(full_lor_point, full_lor_lower, full_lor_upper)
  }

  select <- cbind(combn(drug_names0, 2)[2, ], combn(drug_names0, 2)[1, ])
  par_or <- if (is.null(show0)) {
    full_lor
  } else {
    na.omit(subset(data.frame(full_lor, select),
                   is.element(select[, 1], show) &
                     is.element(select[, 2], show)))
  }

  par_rd <- if (is.null(show0)) {
    full_rd
  } else {
    na.omit(subset(data.frame(full_rd, select),
                   is.element(select[, 1], show) &
                     is.element(select[, 2], show)))
  }

  par_absol <- if (is.null(show0)) {
    absol_risk
  } else {
    na.omit(subset(data.frame(absol_risk, drug_names0),
                   is.element(drug_names0, show)))
  }

  hiera <- if (is.null(show0)) {
    data[, 4]
  } else {
    na.omit(subset(data.frame(data[, 4], drug_names0),
                   is.element(drug_names0, show)))[, 1]
  }

  nt <- length(hiera)

  # Source:https://rdrr.io/github/nfultz/stackoverflow/man/reflect_triangle.html
  reflect_triangle <- function(m, from = c("lower", "upper")) {
    ix <- switch(match.arg(from),
                 lower = upper.tri, upper = lower.tri)(m,  diag = FALSE)
    m[ix] <- t(-m)[ix]
    m
  }

  # Interventions order according to their order value (from best to worst)
  if (all.equal(hiera, as.integer(hiera)) == FALSE) {
    drug_order_col <- drug_order_row <- order(-hiera)
    order_drug <- drug_names[order(-hiera)]
  } else {
    drug_order_col <- drug_order_row <- order(hiera)
    order_drug <- drug_names[order(hiera)]
  }

  # Working on log odds ratios
  point_or0 <- matrix(NA, nrow = length(drug_names), ncol = length(drug_names))
  lower_or0 <- upper_or0 <- point_or0
  point_or0[lower.tri(point_or0, diag = FALSE)] <- par_or[, 1]
  # Incorporate upper triangle
  point_or1 <- reflect_triangle(point_or0, from = "lower")

  # Matrix of lower and upper bound of effect measure (all possible comparisons)
  # Lower triangle
  lower_or0[lower.tri(lower_or0, diag = FALSE)] <- par_or[, 2]
  upper_or0[lower.tri(upper_or0, diag = FALSE)] <- par_or[, 3]
  # Incorporate upper triangle
  lower_or1 <- reflect_triangle(upper_or0, from = "lower")
  lower_or1[lower.tri(lower_or1, diag = FALSE)] <- par_or[, 2]
  upper_or1 <- reflect_triangle(lower_or0, from = "lower")
  upper_or1[lower.tri(upper_or1, diag = FALSE)] <- par_or[, 3]

  # Working on risk differences
  point_rd0 <- matrix(NA, nrow = length(drug_names), ncol = length(drug_names))
  lower_rd0 <- upper_rd0 <- point_rd0
  point_rd0[lower.tri(point_rd0, diag = FALSE)] <- par_rd[, 1]
  # Incorporate upper triangle
  point_rd1 <- reflect_triangle(point_rd0, from = "lower")

  # Matrix of lower and upper bound of effect measure (all possible comparisons)
  # Lower triangle
  lower_rd0[lower.tri(lower_rd0, diag = FALSE)] <- par_rd[, 2]
  upper_rd0[lower.tri(upper_rd0, diag = FALSE)] <- par_rd[, 3]
  # Incorporate upper triangle
  lower_rd1 <- reflect_triangle(upper_rd0, from = "lower")
  lower_rd1[lower.tri(lower_rd1, diag = FALSE)] <- par_rd[, 2]
  upper_rd1 <- reflect_triangle(lower_rd0, from = "lower")
  upper_rd1[lower.tri(upper_rd1, diag = FALSE)] <- par_rd[, 3]

  # Symmetric matrix for effect measure and its bounds
  # after ordering rows and columns from the best to the worst intervention
  # ODDS RATIO
  point_or <- round(exp(point_or1[drug_order_col, drug_order_row]), 2)
  lower_or <- round(exp(lower_or1[drug_order_col, drug_order_row]), 2)
  upper_or <- round(exp(upper_or1[drug_order_col, drug_order_row]), 2)
  # RISK DIFFERENCE
  point_rd <- point_rd1[drug_order_col, drug_order_row]
  lower_rd <- lower_rd1[drug_order_col, drug_order_row]
  upper_rd <- upper_rd1[drug_order_col, drug_order_row]

  # Odds ratio in upper diagonal, risk difference in lower diagonal
  point <- round(point_rd*1000, 0)
  lower <- round(lower_rd*1000, 0)
  upper <- round(upper_rd*1000, 0)
  point[upper.tri(point, diag = FALSE)] <- point_or[upper.tri(point_or,
                                                              diag = FALSE)]
  lower[upper.tri(lower, diag = FALSE)] <- lower_or[upper.tri(lower_or,
                                                              diag = FALSE)]
  upper[upper.tri(upper, diag = FALSE)] <- upper_or[upper.tri(upper_or,
                                                              diag = FALSE)]

  # Spot the statistically significant comparisons
  # ODDS RATIO
  sign_or <- ifelse(upper_or < 1 | lower_or > 1, 1, 0)
  sign_or[is.na(sign_or)] <- 1
  # RISK DIFFERENCE
  sign_rd <- ifelse(upper_rd < 0 | lower_rd > 0, 1, 0)
  sign_rd[is.na(sign_rd)] <- 0

  # Bring both into a matrix (RD: lower triangle; OR: upper triangle)
  sign_rd[upper.tri(sign_rd, diag = FALSE)] <- sign_or[upper.tri(sign_or,
                                                                 diag = FALSE)]
  signif_status <- melt(sign_rd, na.rm = FALSE)[3]

  # Absolute risks presented in 1000 and rounded up to 0 decimals
  # Ordered according to their order value (from the best to the worst)
  if (all.equal(hiera, as.integer(hiera)) == FALSE) {
    point_abs_risk <- round(par_absol[order(-hiera), 1]*1000, 0)
    lower_abs_risk <- round(par_absol[order(-hiera), 2]*1000, 0)
    upper_abs_risk <- round(par_absol[order(-hiera), 3]*1000, 0)
  } else {
    point_abs_risk <- round(par_absol[order(hiera), 1]*1000, 0)
    lower_abs_risk <- round(par_absol[order(hiera), 2]*1000, 0)
    upper_abs_risk <- round(par_absol[order(hiera), 3]*1000, 0)
  }

  # Merge point estimate with 95% credible interval in a new symmetric matric
  final <- matrix(paste0(point,  "\n", "(",
                         lower, ",", " ",
                         upper, ")"),
                  nrow = length(drug_names), ncol = length(drug_names))
  colnames(final) <- order_drug; rownames(final) <- order_drug

  # Include order values in the diagonal of the new matrix
  diag(final) <- ifelse(lower_abs_risk == upper_abs_risk,
                        paste0(point_abs_risk, "\n", ""), #(reference)
                        paste0(point_abs_risk,  "\n", "(",
                               lower_abs_risk, ",", " ", upper_abs_risk, ")"))

  # Preparing the dataset for the ggplot2
  mat_new1 <- melt(final, na.rm = FALSE)

  # Merge both datasets to be used for ggplot2
  mat <- point
  mat_new <- cbind(mat_new1, melt(mat, na.rm = FALSE)[, 3])
  colnames(mat_new) <- c("Var1", "Var2", "value", "value2")

  # Spot the reference intervention
  ref <- drug_names0[1]

  # To create the orders of the lower diagonal
  xmin1 <- rep(seq(0.5, nt - 0.5, 1), each = nt)
  xmax1 <- rep(seq(1.5, nt + 0.5, 1), each = nt)
  ymin1 <- rep(seq(nt - 0.5, 0.5, -1), each = nt)
  ymax1 <- ymin1

  # The league table as a heatmap
  fig <- if (is.element(drug_names0[1], drug_names)) {
    ggplot(mat_new,
           aes(factor(Var2, levels = order_drug[1:nt]),
               factor(Var1, levels = order_drug[nt:1]))) +
      geom_tile(aes(fill = value2)) +
      geom_tile(aes(x = drug_names0[1],
                    y = drug_names0[1]),
                colour = "black",
                fill = "grey90",
                size = 2) +
      geom_fit_text(aes(factor(Var2, levels = order_drug[1:nt]),
                        factor(Var1, levels = order_drug[nt:1]),
                        label = value),
                    fontface = ifelse(signif_status == 1, "bold", "plain"),
                    reflow = TRUE) +
      scale_fill_gradient(low = "white",
                          high = "white",
                          na.value = "grey90") +
      geom_rect(aes(xmin = xmin1,
                    xmax = xmax1,
                    ymin = ymin1,
                    ymax = ymax1),
                color = "black", size = 1) +
      geom_rect(aes(xmin = ymin1,
                    xmax = ymax1,
                    ymin = xmin1,
                    ymax = xmax1),
                color = "black", size = 1) +
      scale_x_discrete(position = "top") +
      labs(x = "", y = "") +
      theme_classic() +
      theme(legend.position = "none",
            axis.title.x = element_text(size = 12,
                                        face = "bold",
                                        colour = "black"),
            axis.title.y = element_text(size = 12,
                                        face = "bold",
                                        colour = "black"),
            axis.text.x = element_text(size = 12,
                                       angle = 50,
                                       hjust = 0.0), #0.5
            axis.text.y = element_text(size = 12),
            plot.caption = element_text(hjust = 0.01))
  } else {
    ggplot(mat_new,
           aes(factor(Var2, levels = order_drug[1:nt]),
               factor(Var1, levels = order_drug[nt:1]))) +
      geom_tile(aes(fill = value2)) +
      geom_fit_text(aes(factor(Var2, levels = order_drug[1:nt]),
                        factor(Var1, levels = order_drug[nt:1]),
                        label = value),
                    fontface = ifelse(signif_status == 1, "bold", "plain"),
                    reflow = TRUE) +
      scale_fill_gradient(low = "white",
                          high = "white",
                          na.value = "grey90") +
      geom_rect(aes(xmin = xmin1,
                    xmax = xmax1,
                    ymin = ymin1,
                    ymax = ymax1),
                color = "black", size = 1) +
      geom_rect(aes(xmin = ymin1,
                    xmax = ymax1,
                    ymin = xmin1,
                    ymax = xmax1),
                color = "black", size = 1) +
      scale_x_discrete(position = "top") +
      labs(x = "", y = "") +
      theme_classic() +
      theme(legend.position = "none",
            axis.title.x = element_text(size = 12,
                                        face = "bold",
                                        colour = "black"),
            axis.title.y = element_text(size = 12,
                                        face = "bold",
                                        colour = "black"),
            axis.text.x = element_text(size = 12,
                                       angle = 50,
                                       hjust = 0.0), #0.5
            axis.text.y = element_text(size = 12),
            plot.caption = element_text(hjust = 0.01))
  }

  # Tabulate relative and absolute effects for the basic parameters
  n_t <- length(drug_names0)
  tab0 <- data.frame(drug_names0,
                     round(rbind(rep(1, 3), exp(full_lor[1:(n_t - 1), ])), 2),
                     absol_risk * 1000,
                     rbind(rep(0, 3), round(full_rd[1:(n_t - 1), ] * 1000, 0)))
  colnames(tab0) <- c("Interventions",
                      "OR", "lower", "upper",
                      "AR", "lower", "upper",
                      "RD", "lower", "upper")
  tab <- if (all.equal(hiera, as.integer(hiera)) == FALSE) {
    subset(tab0, is.element(tab0$Interventions, drug_names))[order(-hiera), ]
  } else {
    subset(tab0, is.element(tab0$Interventions, drug_names))[order(hiera), ]
  }
  rownames(tab) <- NULL

  # Write the table as .xlsx
  if (save_xls == TRUE) {
    write_xlsx(tab, paste0("Table relative $ absolute", ".xlsx"))
  }

  # Collect results
  results <- list(table_relative_absolute_effects =
                    knitr::kable(tab,
                                 align = "lcccccc",
                                 caption = "Relative and absolute effects"),
                  league_table = fig)

  return(results)
}
