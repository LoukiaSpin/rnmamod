#' Barplot for the Kullback-Leibler divergence measure
#'
#' @description Produces a barplot with the Kullback-Leibler divergence measure
#'   from each re-analysis to the primary analysis for a pairwise
#'   comparison. Currently, \code{kld_barplot} is used concerning the impact of
#'   missing participant outcome data.
#'
#' @param robust An object of S3 class \code{\link{robustness_index}}.
#'   See 'Value' in \code{\link{robustness_index}}.
#' @param compar A character vector with two elements that indicates the
#'   pairwise comparison of interest. The first element refers to the
#'   'experimental' intervention and the second element refers to the
#'   'control' intervention of the comparison.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#'
#' @return \code{kld_barplot} returns a panel of barplots on the
#'   Kullback-Leibler divergence measure for each re-analysis.
#'
#' @details \code{kld_barplot} uses the scenarios inherited by
#'   \code{\link{robustness_index}} via the \code{\link{run_sensitivity}}
#'   function. The scenarios for the missingness parameter (see 'Details' in
#'   \code{\link{run_sensitivity}}) in the compared interventions are split to
#'   \emph{Extreme}, \emph{Sceptical}, and \emph{Optimistic} following the
#'   classification of Spineli et al. (2021). In each class, bars will green,
#'   orange, and red colour refer to scenarios without distance, less distant,
#'   and more distant from the primary analysis
#'   (the missing-at-random assumption).
#'
#'   \code{kld_barplot} can be used only when missing participant outcome
#'   data have been extracted for at least one trial. Otherwise, the execution
#'   of the function will be stopped and an error message will be printed on
#'   the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{robustness_index}}, \code{\link{run_model}},
#'   \code{\link{run_sensitivity}}
#
#' @references
#' Kullback S, Leibler RA. On information and sufficiency.
#' \emph{Ann Math Stat} 1951;\bold{22}(1):79--86.
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of
#' primary analysis results: A case study on missing outcome data in pairwise
#' and network meta-analysis.
#' \emph{Res Synth Methods} 2021;\bold{12}(4):475--490. \doi{10.1002/jrsm.1478}
#'
#' @examples
#' data("pma.taylor2004")
#'
#' # Read results from 'run_sensitivity' (using the default arguments)
#' res_sens <- readRDS(system.file('extdata/res_sens_taylor.rds',
#'                     package = 'rnmamod'))
#'
#' # Calculate the robustness index
#' robust <- robustness_index(sens = res_sens,
#'                            threshold = 0.28)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "inositol")
#'
#' # Create the barplot for the comparison 'inositol versus placebo'
#' kld_barplot(robust = robust,
#'             compar = c("inositol", "placebo"),
#'             drug_names = interv_names)
#'
#' @export
kld_barplot <- function(robust, compar, drug_names) {

  if (any(is.na(robust))) {
    aa <- "Missing participant outcome data have *not* been collected."
    stop(paste(aa, "This function cannot be used."), call. = FALSE)
  }

  if (dim(robust$kld)[2] < 25) {
    stop("Use *only* for different scenarios about the missingness parameter.",
         call. = FALSE)
  }

  drug_names <- if (missing(drug_names)) {
    stop("The argument 'drug_names' has not been defined.", call. = FALSE)
  } else {
    drug_names
  }

  compar <- if (missing(compar)) {
    stop("The argument 'compar' needs to be defined.", call. = FALSE)
  } else if (length(drug_names) < 3 & missing(compar)) {
    c(drug_names[2], drug_names[1])
  } else if (!is.element(compar[1], drug_names) |
             !is.element(compar[2], drug_names)) {
    stop("The value of 'compar' is not found in the argument 'drug_names'.",
         call. = FALSE)
  } else if (is.element(compar[1], drug_names) &
             is.element(compar[2], drug_names) &
             match(compar[1], drug_names) < match(compar[2], drug_names)) {
    stop("Re-arrange the order of the element in the argument 'compar'.",
         call. = FALSE)
  } else if (is.element(compar[1], drug_names) &
             is.element(compar[2], drug_names) &
             match(compar[1], drug_names) > match(compar[2], drug_names)) {
    compar
  }

  # Indicate all possible comparisons (necessary for NMA)
  comparison <- matrix(combn(drug_names, 2),
                       nrow = length(combn(drug_names, 2)) / 2,
                       ncol = 2,
                       byrow = TRUE)
  compar_id <- which(comparison[, 1] == compar[2] &
                       comparison[, 2] == compar[1])
  experim <- comparison[compar_id, 2]
  control <- comparison[compar_id, 1]


  # Define the scenarios
  measure <- robust$measure
  scenarios <- if (is.element(measure, c("OR", "ROM"))) {
    cbind(rep(as.character(fractions(exp(robust$scenarios))),
              each = length(robust$scenarios)),
          rep(as.character(fractions(exp(robust$scenarios))),
              times = length(robust$scenarios)))
  } else if (!is.element(measure, c("OR", "ROM"))) {
    cbind(rep(robust$scenarios, each = length(robust$scenarios)),
          rep(robust$scenarios, times = length(robust$scenarios)))
  }
  colnames(scenarios) <- c("active", "ctrl")
  kld <- robust$kld[compar_id, ]
  len.scen <- length(unique(scenarios[, 1]))

  # Rank the scenarios to calculate their distance in the compared arms
  ranked_scenarios <- cbind(rep(rank(1:len.scen), each = len.scen),
                            rep(rank(1:len.scen), times = len.scen))
  distance <-
    ifelse(abs(ranked_scenarios[, 1] - ranked_scenarios[, 2]) > 1,
           "more distant",
           ifelse(abs(ranked_scenarios[, 1] - ranked_scenarios[, 2]) == 1,
                  "less distant", "no distance"))

  # Characterise the scenarios to extreme, sceptical, and optimistic with
  # respect to their position from MAR
  plausibility <- factor(
    c("Extreme", rep("Sceptical", len.scen - 2), "Extreme",
      rep(c("Sceptical", rep("Optimistic", len.scen - 2), "Sceptical"),
          len.scen - 2),
      "Extreme", rep("Sceptical", len.scen - 2), "Extreme"),
    levels = c("Extreme", "Sceptical", "Optimistic"))

  # Dataset for the barplot
  dataset_new <- data.frame(kld[-13],
                            paste0(scenarios[-13, 1], ",", scenarios[-13, 2]),
                            plausibility[-13],
                            distance[-13])
  colnames(dataset_new) <- c("kld",
                             "scenarios",
                             "plausibility",
                             "distance")

  # In each facet, x-axis is sorted by KLD in descending order
  barplot <- ggplot(dataset_new,
                    aes(x = reorder(scenarios, -kld),
                        y = kld,
                        fill = distance)) +
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_manual(breaks = c("more distant", "less distant", "no distance"),
                      values = c("#D55E00", "orange", "#009E73")) +
    facet_grid(. ~  plausibility, scales = "free_x", space = "free") +
    labs(x = "Scenarios (active vs control)",
         y = "Kullback-Leibler divergence measure",
         fill = "Distance between the scenarios") +
    ylim(0, ifelse(max(dataset_new$kld) > 0.3, max(dataset_new$kld), 0.3)) +
    ggtitle(paste(experim, "versus", control)) +
    theme_classic() +
    theme(axis.title = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 10.5),
          axis.text.x = element_text(size = 10.5, angle = 45, hjust = 1),
          legend.position = "bottom",
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 12),
          strip.text = element_text(size = 12),
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

  return(barplot)
}
