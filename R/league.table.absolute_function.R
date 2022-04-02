#' League table for relative and absolute effects
#'
#' @description
#'   Provides a league table of the estimated odds ratio, and risk difference
#'   per 1000 participants for all possible comparisons of interventions in the
#'   network. The main diagonal of the table presents the absolute risk for each
#'   intervention in the network. \code{league_table_absolute} can be used for a
#'   random-effects or fixed-effect network meta-analysis.
#'   It is applied for one binary outcome only.
#'
#' @param full An object of S3 class \code{\link{run_model}}.
#'   See 'Value' in \code{\link{run_model}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}.
#' @param show A vector of at least three character strings that refer to the
#'   names of the interventions \emph{exactly} as defined in \code{drug_names}.
#'   Then, the league table will be created for these interventions only.
#'   If \code{show} is not defined, the league table will present all
#'   interventions as defined in \code{drug_names}.
#'
#' @return A league table showing the posterior estimate and 95\% credible
#'   interval of the odds ratio (upper off-diagonals), risk difference per 1000
#'   participants (lower off-diagonals), and absolute risks per 1000
#'   participants (main diagonal).
#'
#' @details The user must define the argument \code{measure = "RD"} in
#'   \code{\link{run_model}}; otherwise, the function will be stopped and an
#'   error message will be printed in the R console.
#'
#'   The rows and columns of the league table display the names of the
#'   interventions  sorted by decreasing order from the best to the worst
#'   based on their SUCRA value (Salanti et al., 2011) for the odds ratio. The
#'   upper off-diagonals contain the posterior median and 95\% credible interval
#'   of the odds ratio, the lower off-diagonals contain the posterior median and
#'   95\% credible interval of the risk difference (per 1000 participants), and
#'   the main diagonal comprises the posterior median and 95\% credible interval
#'   of the absolute risks (per 1000 participants) of the  corresponding
#'   interventions. The reference intervention of the network (which the
#'   baseline risk has been selected for) is indicated in the main diagonal with
#'   a homonymous label.
#'
#'   Comparisons between interventions should be read from left to right.
#'   Results that indicate strong evidence in favor of the row-defining
#'   intervention (i.e. the respective 95\% credible interval does not include
#'   the null value) are indicated in bold and dark blue.
#'
#'   To obtain unique absolute risks for each intervention, the network
#'   meta-analysis model has been extended to incorporate the transitive risks
#'   framework, namely, an intervention has the same absolute risk regardless of
#'   the comparator intervention(s) in a trial (Spineli et al., 2017).
#'   The absolute risks are a function of the odds ratio (the \strong{base-case}
#'   effect measure for a binary outcome) and the selected baseline risk for the
#'   reference intervention (Appendix in Dias et al., 2013). See 'Arguments' in
#'   \code{\link{run_model}}. We advocate using the odds ratio as an effect
#'   measure for its desired mathematical properties. Then, the risk difference
#'   can be obtained as a function of the absolute risks of the corresponding
#'   interventions in the comparison of interest.
#'
#'   \code{league_table_absolute} can be used only for a network of
#'   interventions. In the case of two interventions, the execution of the
#'   function will be stopped and an error message will be printed in the R
#'   console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_model}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries
#' for presenting results from multiple-treatment meta-analysis: an overview and
#' tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71.
#' \doi{10.1016/j.jclinepi.2010.03.016}
#'
#' Spineli LM, Brignardello-Petersen R, Heen AF, Achille F, Brandt L,
#' Guyatt GH, et al. Obtaining absolute effect estimates to facilitate shared
#' decision making in the context of multiple-treatment comparisons.
#' Abstracts of the Global Evidence Summit, Cape Town, South Africa.
#' \emph{Cochrane Database of Systematic Reviews} 2017;\bold{9}(Suppl 1):1891.
#'
#' @export
league_table_absolute <- function(full, drug_names, show = NULL) {


  if ((full$type != "nma") || is.null(full$type)) {
    stop("'full' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  if (full$measure != "RD") {
    stop("The argument 'measure' in 'run_model' must be 'RD'", call. = FALSE)
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

  select <- cbind(combn(drug_names0, 2)[2, ], combn(drug_names0, 2)[1, ])
  par_or <- if (is.null(show0)) {
    full$EM_LOR
  } else {
    na.omit(subset(data.frame(full$EM_LOR, select),
                   is.element(select[, 1], show) &
                     is.element(select[, 2], show)))
  }
  par_rd <- if (is.null(show0)) {
    full$EM
  } else {
    na.omit(subset(data.frame(full$EM, select),
                   is.element(select[, 1], show) &
                     is.element(select[, 2], show)))
  }
  sucra <- if (is.null(show0)) {
    full$SUCRA[, 1]
  } else {
    na.omit(subset(data.frame(full$SUCRA[, 1], drug_names0),
                   is.element(drug_names0, show)))[, 1]
  }
  abs_risk <- if (is.null(show0)) {
    full$abs_risk
  } else {
    na.omit(subset(data.frame(full$abs_risk, drug_names0),
                   is.element(drug_names0, show)))
  }
  nt <- length(sucra)

  # Source:https://rdrr.io/github/nfultz/stackoverflow/man/reflect_triangle.html
  reflect_triangle <- function(m, from = c("lower", "upper")) {
    ix <- switch(match.arg(from),
                 lower = upper.tri, upper = lower.tri)(m,  diag = FALSE)
    m[ix] <- t(-m)[ix]
    m
  }

  # Interventions order according to their SUCRA value (from best to worst)
  drug_order_col <- drug_order_row <- order(-sucra)

  # Order interventions according to their SUCRA value (from best to worst)
  order_drug <- drug_names[order(-sucra)]

  # Working on log odds ratios
  point_or0 <- matrix(NA, nrow = length(drug_names), ncol = length(drug_names))
  lower_or0 <- upper_or0 <- point_or0
  point_or0[lower.tri(point_or0, diag = FALSE)] <- par_or[, 5]
  # Incorporate upper triangle
  point_or1 <- reflect_triangle(point_or0, from = "lower")

  # Matrix of lower and upper bound of effect measure (all possible comparisons)
  # Lower triangle
  lower_or0[lower.tri(lower_or0, diag = FALSE)] <- par_or[, 3]
  upper_or0[lower.tri(upper_or0, diag = FALSE)] <- par_or[, 7]
  # Incorporate upper triangle
  lower_or1 <- reflect_triangle(upper_or0, from = "lower")
  lower_or1[lower.tri(lower_or1, diag = FALSE)] <- par_or[, 3]
  upper_or1 <- reflect_triangle(lower_or0, from = "lower")
  upper_or1[lower.tri(upper_or1, diag = FALSE)] <- par_or[, 7]

  # Working on risk differences
  point_rd0 <- matrix(NA, nrow = length(drug_names), ncol = length(drug_names))
  lower_rd0 <- upper_rd0 <- point_rd0
  point_rd0[lower.tri(point_rd0, diag = FALSE)] <- par_rd[, 5]
  # Incorporate upper triangle
  point_rd1 <- reflect_triangle(point_rd0, from = "lower")

  # Matrix of lower and upper bound of effect measure (all possible comparisons)
  # Lower triangle
  lower_rd0[lower.tri(lower_rd0, diag = FALSE)] <- par_rd[, 3]
  upper_rd0[lower.tri(upper_rd0, diag = FALSE)] <- par_rd[, 7]
  # Incorporate upper triangle
  lower_rd1 <- reflect_triangle(upper_rd0, from = "lower")
  lower_rd1[lower.tri(lower_rd1, diag = FALSE)] <- par_rd[, 3]
  upper_rd1 <- reflect_triangle(lower_rd0, from = "lower")
  upper_rd1[lower.tri(upper_rd1, diag = FALSE)] <- par_rd[, 7]

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
  # Ordered according to their SUCRA value (from the best to the worst)
  point_abs_risk <- round(abs_risk[order(-sucra), 5]*1000, 0)
  lower_abs_risk <- round(abs_risk[order(-sucra), 3]*1000, 0)
  upper_abs_risk <- round(abs_risk[order(-sucra), 7]*1000, 0)

  # Merge point estimate with 95% credible interval in a new symmetric matric
  final <- matrix(paste0(point,  "\n", "(",
                         lower, ",", " ",
                         upper, ")"),
                  nrow = length(drug_names), ncol = length(drug_names))
  colnames(final) <- order_drug; rownames(final) <- order_drug

  # Include SUCRA values in the diagonal of the new matrix
  diag(final) <- ifelse(lower_abs_risk == upper_abs_risk,
                        paste0(point_abs_risk, "\n", "(reference)"),
                        paste0(point_abs_risk,  "\n", "(",
                               lower_abs_risk, ",", " ", upper_abs_risk, ")"))

  # Preparing the dataset for the ggplot2
  mat_new1 <- melt(final, na.rm = FALSE)

  # Merge both datasets to be used for ggplot2
  mat <- point
  mat_new <- cbind(mat_new1, melt(mat, na.rm = FALSE)[, 3])
  colnames(mat_new) <- c("Var1", "Var2", "value", "value2")

  # To create the orders of the lower diagonal
  xmin1 <- rep(seq(0.5, nt - 0.5, 1), each = nt)
  xmax1 <- rep(seq(1.5, nt + 0.5, 1), each = nt)
  ymin1 <- rep(seq(nt - 0.5, 0.5, -1), each = nt)
  ymax1 <- ymin1

  # The league table as a heatmap
  ggplot(mat_new,
         aes(factor(Var2, levels = order_drug[1:nt]),
             factor(Var1, levels = order_drug[nt:1]))) +
    geom_tile(aes(fill = value2)) +
    geom_fit_text(aes(factor(Var2, levels = order_drug[1:nt]),
                      factor(Var1, levels = order_drug[nt:1]),
                      label = value),
                  fontface = ifelse(signif_status == 1, "bold", "plain"),
                  colour = ifelse(signif_status == 1, "blue", "black"),
                  reflow = TRUE) +
    scale_fill_gradient(low = "white", high = "white", na.value = "grey90") +
    geom_rect(aes(xmin = xmin1, xmax = xmax1, ymin = ymin1, ymax = ymax1),
              color = "black", size = 1) +
    geom_rect(aes(xmin = ymin1, xmax = ymax1, ymin = xmin1, ymax = xmax1),
              color = "black", size = 1) +
    scale_x_discrete(position = "top") +
    labs(x = "", y = "") +
    theme_classic() +
    theme(legend.position = "none",
          axis.title.x = element_text(size = 12, face = "bold",
                                      colour = "black"),
          axis.title.y = element_text(size = 12, face = "bold",
                                      colour = "black"),
          axis.text.x = element_text(size = 12, hjust = 0.5), #angle = 50,
          axis.text.y = element_text(size = 12),
          plot.caption = element_text(hjust = 0.01))

}
