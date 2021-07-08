#' Barplot for the Kullback-Leibler divergence measure: Investigate the impact of missing participant outcome data
#'
#' @description This function produces a barplot with the Kullback-Leibler divergence measure for all combinations of scenarios for the interventions in a pairwise comparison.
#'   Currently, \code{KLD.barplot} is used concerning the impact of missing participant outcome data.
#'
#' @param robust An object of S3 class \code{\link{robustness.index}}. See 'Value' in \code{\link{robustness.index}}.
#' @param compar A character vector with two characters that indicates the pairwise comparison of interest. The first character refers to the 'experimental' intervention
#'   and the second character refers to the 'control' intervention of the comparison.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @return \code{KLD.barplot} prints on the R console a message on the threshold of robustness determined by the user in green text.
#' Then, the function returns the following list of elements:
#' \tabular{ll}{
#'  \code{RI} \tab A a numeric scalar or vector on the robustnessindex values. In the case of a pairwise meta-analysis (PMA), \code{RI} is scalar as only one summary effect size is obtained.
#'    In the case of network meta-analysis (NMA), \code{RI} is a vector with length equal to the number of possible pairwise comparisons; one robustness index per possible comparison.\cr
#'  \tab \cr
#'  \code{robust} \tab A character vector on whether the primary analysis results are \emph{robust} or \emph{frail} to the different re-analyses.\cr
#'  \tab \cr
#'  \code{KLD} \tab A vector or matrix on the Kullback-Leibler divergence (KLD) measure in the summary effect size from a subsequent re-analysis to the primary analysis.
#'    In the case of a PMA, \code{KLD} is a vector with length equal to the number of total analyses.
#'    The latter equals the square of the number of scenarios indicated in argument \code{mean.scenarios} of \code{run.sensitivity}. Therefore, one KLD value per analysis.
#'    In the case of NMA, \code{RI} is a matrix with number of rows equal to the number of total analyses and number of columns equal to the number of possible pairwise comparisons;
#'    one KLD value per analysis and possible comparison.\cr
#' }
#'
#' @details The user may consider the values 0.28 and 0.17 in the argument \code{threshold} for binary and continuous outcome data (the default values), respectively, or consider other plausible values.
#'   Spineli et al. (2021) offers a discussion on specifying the \code{threshold} of robustness.
#'
#'   In \code{robust}, the value \code{"robust"} appears when \code{RI} \eqn{<} \code{threshold}); otherwise, the value \code{"frail"} appears.
#'
#'   \code{robustness.index} can be used only for when missing participant outcome data have been extracted for at least one trial. Otherwise, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.sensitivity}}, \code{\link{robustness.index}}, \code{\link{run.model}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of primary analysis results: A case study on missing outcome data in pairwise and network meta-analysis. \emph{Res Synth Methods} 2021;\bold{12}(4):475--490. [\doi{10.1002/jrsm.1478}]
#'
#' Kullback S, Leibler RA. On information and sufficiency. \emph{Ann Math Stat} 1951;\bold{22}(1):79--86. [\doi{10.1214/aoms/1177729694}]
#'
#' @examples
#' data("nma.Baker2009.RData")
#'
#' # Perform the sensitivity analysis (using the 'default' of the argument 'mean.scenarios')
#' res.sens <- run.sensitivity(full = res1, assumption = "IDE-ARM", var.misspar = 1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' # Calculate the robustness index
#' robuat <- robustness.index(sens = res.sens, primary.scenar = 13, threshold = 0.28)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("budesodine", "budesodine plus formoterol", "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium", "placebo")
#'
#' # Crate the barplot for the comparison 'tiotropium versus salmeterol'
#' KLD.barplot(robust = robust, compar = c("tiotropium", "salmeterol"), drug.names = interv.names)
#'
#' @export
KLD.barplot <- function(robust, compar, drug.names){


  if (is.na(robust)) {
    stop("Missing participant outcome data have *not* been collected. This function cannot be used.", call. = F)
  }


  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- (1 + sqrt(1 + 8*length(robust$robust)))/2
    as.character(1:nt)
  } else {
    drug.names
  }


  compar <- if (length(drug.names) > 2 & missing(compar)) {
    stop("The argument 'compar' needs to be defined", call. = F)
  } else if (length(drug.names) < 3 & missing(compar)) {
    1
  } else {
    compar
  }


  outcome <- if (is.element(robust$measure, c("MD", "SMD", "ROM"))) {
    "continuous"
  } else {
    "binary"
  }


  comparisons <- t(combn(drug.names, 2))


  if(outcome == "binary"){

    ## Define the scenarios for IMOR
    scenarios <- cbind(rep(c(0.3, 0.50, 1, 2, 3), each = 5), rep(c(0.3, 0.50, 1, 2, 3), times = 5))

  } else {

    ## Define the scenarios for IMDoM
    scenarios <- cbind(rep(c(-2, -1, 0, 1, 2), each = 5), rep(c(-2, -1, 0, 1, 2), times = 5))

  }
  colnames(scenarios) <- c("active", "ctrl")

  KLD <- robust$KLD[compar, ]


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
    #geom_hline(yintercept = 0.28, linetype = 2) +
    ylim(0, ifelse(max(dataset.new$KLD) > 0.3, max(dataset.new$KLD), 0.3)) +
    ggtitle(paste(comparisons[compar, 2], "versus", comparisons[compar, 1])) +
    theme_classic() +
    theme(axis.title = element_text(size = 12, face = "bold"), axis.text = element_text(size = 10.5), axis.text.x = element_text(size = 10.5, angle = 45, hjust = 1),
          legend.position = "bottom", legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 11),
          strip.text = element_text(size = 12), plot.title = element_text(size = 14, face = "bold", hjust = 0.5))

  return(barplot)
}








