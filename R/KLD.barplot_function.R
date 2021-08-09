#' Barplot for the Kullback-Leibler divergence measure: Investigate the impact of missing participant outcome data
#'
#' @description This function produces a barplot with the Kullback-Leibler divergence (KLD) measure from each re-analysis to primary analysis for a pairwise comparison.
#'   Currently, \code{KLD.barplot} is used concerning the impact of missing participant outcome data.
#'
#' @param robust An object of S3 class \code{\link{robustness.index}}. See 'Value' in \code{\link{robustness.index}}.
#' @param compar A character vector with two characters that indicates the pairwise comparison of interest. The first character refers to the 'experimental' intervention
#'   and the second character refers to the 'control' intervention of the comparison.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @return \code{KLD.barplot} returns a panel of barplots on the KLD measure for each analysis.
#'
#' @details The scenarios for the missingness parameter in compared interventions are split to \emph{Extreme}, \emph{Sceptical}, and \emph{Optimistic} following the classification of Spineli et al. (2021).
#'   In each class, bars will green, orange, and red colour refer to scenarios without distance, less distant, and more distant scenarios from the primary analysis (the missing-at-random assumption).
#'
#'   \code{KLD.barplot} can be used only for when missing participant outcome data have been extracted for at least one trial. Otherwise, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{robustness.index}}, \code{\link{run.model}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of primary analysis results: A case study on missing outcome data in pairwise and network meta-analysis. \emph{Res Synth Methods} 2021;\bold{12}(4):475--490. [\doi{10.1002/jrsm.1478}]
#'
#' Kullback S, Leibler RA. On information and sufficiency. \emph{Ann Math Stat} 1951;\bold{22}(1):79--86. [\doi{10.1214/aoms/1177729694}]
#'
#' @examples
#' data("nma.baker2009")
#'
#' # Perform a random-effects network meta-analysis
#' res1 <- run.model(data = nma.baker2009,
#'                   measure = "OR",
#'                   model = "RE",
#'                   assumption = "IDE-ARM",
#'                   heter.prior = list("halfnormal", 0, 1),
#'                   mean.misspar = 0,
#'                   var.misspar = 1,
#'                   D = 1,
#'                   n.chains = 3,
#'                   n.iter = 10000,
#'                   n.burnin = 1000,
#'                   n.thin = 1)
#'
#' # Perform the sensitivity analysis (using the 'default' of the argument 'mean.scenarios')
#' res.sens <- run.sensitivity(full = res1,
#'                             assumption = "IDE-ARM",
#'                             var.misspar = 1,
#'                             n.chains = 3,
#'                             n.iter = 10000,
#'                             n.burnin = 1000,
#'                             n.thin = 1)
#'
#' # Calculate the robustness index
#' robust <- robustness.index(sens = res.sens, threshold = 0.28)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("budesodine", "budesodine plus formoterol", "fluticasone", "fluticasone plus
#'                   salmeterol", "formoterol", "salmeterol", "tiotropium", "placebo")
#'
#' # Crate the barplot for the comparison 'tiotropium versus salmeterol'
#' KLD.barplot(robust = robust, compar = c("tiotropium", "salmeterol"), drug.names = interv.names)
#'
#' \dontshow{
#' closeAllConnections()
#' }
#'
#' @export
KLD.barplot <- function(robust, compar, drug.names){


  if (any(is.na(robust))) {
    stop("Missing participant outcome data have *not* been collected. This function cannot be used.", call. = F)
  }


  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- (1 + sqrt(1 + 8*length(robust$robust)))/2
    as.character(1:nt)
  } else {
    drug.names
  }


  ## Indicate all possible comparisons (necessary for NMA)
  comparison <- matrix(combn(drug.names, 2), nrow = length(combn(drug.names, 2))/2, ncol = 2, byrow = T)
  compar.id <- which(comparison[, 1] == compar[2] & comparison[, 2] == compar[1])
  experim <- comparison[compar.id, 2]
  control <- comparison[compar.id, 1]

  compar <- if (length(drug.names) > 2 & missing(compar)) {
    stop("The argument 'compar' needs to be defined", call. = F)
  } else if (length(drug.names) < 3 & missing(compar)) {
    c(comparison[2], comparison[2])
  } else {
    compar
  }


  ## Define the scenarios
  scenarios <- if (is.element(robust$measure, c("OR", "ROM"))) {
    cbind(rep(as.character(fractions(exp(robust$scenarios))), each = length(robust$scenarios)),
          rep(as.character(fractions(exp(robust$scenarios))), times = length(robust$scenarios)))
  } else if (is.element(robust$measure, c("MD", "SMD"))) {
    cbind(rep(robust$scenarios, each = length(robust$scenarios)),
          rep(robust$scenarios, times = length(robust$scenarios)))
  }


  colnames(scenarios) <- c("active", "ctrl")


  KLD <- robust$KLD[compar.id, ]


  ## Rank the scenarios to calculate their distance in the compared arms
  (ranked.scenarios <- cbind(rep(rank(1:5), each = 5), rep(rank(1:5), times = 5)))
  (distance <- ifelse(abs(ranked.scenarios[, 1] - ranked.scenarios[, 2]) > 1, "more distant",
                      ifelse(abs(ranked.scenarios[, 1] - ranked.scenarios[, 2]) == 1, "less distant", "no distance")))


  ## Characterise the scenarios to extreme, sceptical, and optimistic with respect to their position from MAR
  plausibility <- factor(c("Extreme", rep("Sceptical", 3), "Extreme", "Sceptical", rep("Optimistic", 3), "Sceptical", "Sceptical", rep("Optimistic", 3),
                           "Sceptical", "Sceptical", rep("Optimistic", 3), "Sceptical", "Extreme", rep("Sceptical", 3), "Extreme"),
                         levels = c("Extreme", "Sceptical", "Optimistic"))


  ## Dataset for the barplot
  dataset.new <- data.frame(KLD[-13], paste0(scenarios[-13, 1], ",", scenarios[-13, 2]), plausibility[-13], distance[-13])
  colnames(dataset.new) <- c("KLD", "scenarios", "plausibility", "distance")


  ## In each facet, x-axis is sorted by KLD in descending order
  barplot <- ggplot(dataset.new, aes(x = reorder(scenarios, -KLD), y = KLD, fill = distance)) +
    geom_bar(stat = "identity", width = 0.5) +
    scale_fill_manual(breaks = c("more distant", "less distant", "no distance"), values = c("#D55E00", "orange", "#009E73")) +
    facet_grid(. ~  plausibility, scales = "free_x", space = "free") +
    labs(x = "Scenarios (active vs control)", y = "Kullback-Leibler divergence measure", fill = "Distance between the scenarios") +
    ylim(0, ifelse(max(dataset.new$KLD) > 0.3, max(dataset.new$KLD), 0.3)) +
    ggtitle(paste(experim, "versus", control)) +
    theme_classic() +
    theme(axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10.5),
                   axis.text.x = element_text(size = 10.5, angle = 45, hjust = 1),
          legend.position = "bottom", legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12),
          strip.text = element_text(size = 12), plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

  return(barplot)
}








