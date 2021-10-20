#' Barplot for the Kullback-Leibler divergence measure
#'
#' @description Produces a barplot with the Kullback-Leibler divergence (KLD)
#'   measure from each re-analysis to primary analysis for a pairwise
#'   comparison. Currently, \code{kld_barplot} is used concerning the impact of
#'   missing participant outcome data.
#'
#' @param robust An object of S3 class \code{\link{robustness_index}}.
#'   See 'Value' in \code{\link{robustness_index}}.
#' @param compar A character vector with two characters that indicates the
#'   pairwise comparison of interest. The first character refers to the
#'   'experimental' intervention and the second character refers to the
#'   'control' intervention of the comparison.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If the argument \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#'
#' @return \code{kld_barplot} returns a panel of barplots on the KLD measure for
#'   each analysis.
#'
#' @details The scenarios for the missingness parameter in compared
#'   interventions are split to \emph{Extreme}, \emph{Sceptical}, and
#'   \emph{Optimistic} following the classification of Spineli et al. (2021).
#'   In each class, bars will green, orange, and red colour refer to scenarios
#'   without distance, less distant, and more distant scenarios from the primary
#'    analysis (the missing-at-random assumption).
#'
#'   \code{kld_barplot} can be used only for when missing participant outcome
#'   data have been extracted for at least one trial. Otherwise, the execution
#'   of the function will be stopped and an error message will be printed in
#'   the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{robustness_index}}, \code{\link{run_model}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of
#' primary analysis results: A case study on missing outcome data in pairwise
#' and network meta-analysis.
#' \emph{Res Synth Methods} 2021;\bold{12}(4):475--490.
#' [\doi{10.1002/jrsm.1478}]
#'
#' Kullback S, Leibler RA. On information and sufficiency.
#' \emph{Ann Math Stat} 1951;\bold{22}(1):79--86.
#' [\doi{10.1214/aoms/1177729694}]
#'
#' @examples
#' data("nma.baker2009")
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' res <- run_model(data = nma.baker2009,
#'                  measure = "OR",
#'                  model = "RE",
#'                  assumption = "IDE-ARM",
#'                  heter_prior = list("halfnormal", 0, 1),
#'                  mean_misspar = c(0, 0),
#'                  var_.misspar = 1,
#'                  D = 1,
#'                  n_chains = 3,
#'                  n_iter = 10000,
#'                  n_burnin = 1000,
#'                  n_thin = 1)
#'
#' # Perform the sensitivity analysis (missing-at-random assumption)
#' res_sens <- run_sensitivity(full = res,
#'                             var_misspar = 1,
#'                             n_chains = 3,
#'                             n_iter = 10000,
#'                             n_burnin = 1000,
#'                             n_thin = 1)
#'
#' # Calculate the robustness index
#' robust <- robustness_index(sens = res_sens,
#'                            threshold = 0.28)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "budesonide", "budesonide plus formoterol",
#'                   "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium")
#'
#' # Crate the barplot for the comparison 'tiotropium versus salmeterol'
#' kld_barplot(robust = robust,
#'             compar = c("tiotropium", "salmeterol"),
#'             drug_names = interv_names)
#' }
#' @export
kld_barplot <- function(robust, compar, drug_names) {

  if (any(is.na(robust))) {
    aa <- "Missing participant outcome data have *not* been collected."
    stop(paste(aa, "This function cannot be used."), call. = F)
  }

  if (dim(robust$kld)[2] < 25) {
    stop("Use *only* for different scenarios about the missingness parameter.",
         call. = F)
  }

  drug_names <- if (missing(drug_names)) {
    aa <- "The argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as specified in 'data' is used as"
    cc <- "intervention names"
    message(cat(paste0("\033[0;", col = 32, "m", aa, " ", bb, " ", cc,
                       "\033[0m", "\n")))
    nt <- (1 + sqrt(1 + 8 * length(robust$robust))) / 2
    as.character(1:nt)
  } else {
    drug_names
  }

  compar <- if (missing(compar)) {
    stop("The argument 'compar' needs to be defined", call. = F)
  } else if (length(drug_names) < 3 & missing(compar)) {
    c(drug_names[2], drug_names[1])
  } else if (!is.element(compar, drug_names)) {
    stop("The value of 'compar' is not found in the argument 'drug_names'",
         call. = F)
  } else if (is.element(compar, drug_names)) {
    compar
  }

  # Indicate all possible comparisons (necessary for NMA)
  comparison <- matrix(combn(drug_names, 2),
                       nrow = length(combn(drug_names, 2)) / 2,
                       ncol = 2,
                       byrow = T)
  compar_id <- which(comparison[, 1] == compar[2] &
                       comparison[, 2] == compar[1])
  experim <- comparison[compar_id, 2]
  control <- comparison[compar_id, 1]


  # Define the scenarios
  scenarios <- if (is.element(robust$measure, c("OR", "ROM"))) {
    cbind(rep(as.character(fractions(exp(robust$scenarios))),
              each = length(robust$scenarios)),
          rep(as.character(fractions(exp(robust$scenarios))),
              times = length(robust$scenarios)))
  } else if (is.element(robust$measure, c("MD", "SMD"))) {
    cbind(rep(robust$scenarios, each = length(robust$scenarios)),
          rep(robust$scenarios, times = length(robust$scenarios)))
  }
  colnames(scenarios) <- c("active", "ctrl")
  kld <- robust$kld[compar_id, ]

  # Rank the scenarios to calculate their distance in the compared arms
  ranked_scenarios <- cbind(rep(rank(1:5), each = 5), rep(rank(1:5), times = 5))
  distance <-
    ifelse(abs(ranked_scenarios[, 1] - ranked_scenarios[, 2]) > 1,
           "more distant",
           ifelse(abs(ranked_scenarios[, 1] - ranked_scenarios[, 2]) == 1,
                  "less distant", "no distance"))

  # Characterise the scenarios to extreme, sceptical, and optimistic with
  # respect to their position from MAR
  plausibility <- factor(c("Extreme", rep("Sceptical", 3), "Extreme",
                           "Sceptical", rep("Optimistic", 3), "Sceptical",
                           "Sceptical", rep("Optimistic", 3), "Sceptical",
                           "Sceptical", rep("Optimistic", 3), "Sceptical",
                           "Extreme", rep("Sceptical", 3), "Extreme"),
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
