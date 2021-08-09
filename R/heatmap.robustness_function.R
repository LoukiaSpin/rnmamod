#' Heatmap of robsutness index: Investigate the impact of missing participant outcome data
#'
#' @description This functions facilitates the detection of comparisons that are sensible to different assumptions about the informative missingness parameter
#'   for the compared interventions. The heatmap is based on the robustness index of each possible pairwise comparison in the investigated network (see \code{robustness.index}).
#'   Currently, \code{heatmap.robustness} is used concerning the impact of missing participant outcome data.
#'
#' @param robust An object of S3 class \code{\link{robustness.index}}. See 'Value' in \code{\link{robustness.index}}.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @details The heatmap illustrates the robustness index for each possible pairwise comparison in the network.
#'   The pairwise comparisons are read from left to right and are highlighted either with a green or red colour.
#'   Comparisons highlighted with green or red colour imply robust or frail conclusions for the primary analysis.
#'   This corresponds to robustness index below or at least the selected threshold of robustness (see 'Details'in \code{robustness.index}).
#'   The robustness index value of each pairwise comparison also appears in the corresponding cell.
#'   When there is at least one comparison with frail conclusions, the primary analysis results may be questionable for the whole network (Spineli et al., 2021).
#'
#'   \code{heatmap.robustness} uses the threshold of robustness selected in the \code{robustness.index} function.
#'
#'   \code{heatmap.robustness} can be used only for a network of interventions and when missing participant outcome data have been extracted for at least one trial.
#'   Otherwise, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @return \code{heatmap.robustness} first prints on the R console a message on the threshold of robustness determined by the user in the \code{robustness.index} function.
#'   Then, it returns a lower triangular heatmap matrix with the robustness index value of all possible pairwise comparisons.
#'
#' @import reshape2 ggplot2
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{robustness.index}}, \code{\link{run.model}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of primary analysis results: A case study on missing outcome data in pairwise and network meta-analysis. \emph{Res Synth Methods} 2021;\bold{12}(4):475--490. [\doi{10.1002/jrsm.1478}]
#'
#' @examples
#' data("nma.liu2013")
#'
#' # Perform a random-effects network meta-analysis
#' res1 <- run.model(data = nma.liu2013,
#'                   measure = "OR",
#'                   model = "RE",
#'                   assumption = "IDE-ARM",
#'                   heter.prior = list("halfnormal", 0, 1),
#'                   mean.misspar = 0,
#'                   var.misspar = 1,
#'                   D = 1,
#'                   n.chains = 2,
#'                   n.iter = 1000,
#'                   n.burnin = 100,
#'                   n.thin = 1)
#'
#' # Perform the sensitivity analysis (using the 'default' of the argument 'mean.scenarios')
#' res.sens <- run.sensitivity(full = res1,
#'                             assumption = "IDE-ARM",
#'                             var.misspar = 1,
#'                             n.chains = 2,
#'                             n.iter = 1000,
#'                             n.burnin = 100,
#'                             n.thin = 1)
#'
#' # Calculate the robustness index
#' robust <- robustness.index(sens = res.sens, threshold = 0.28)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("budesodine", "budesodine plus formoterol", "fluticasone", "fluticasone plus
#'                   salmeterol", "formoterol", "salmeterol", "tiotropium", "placebo")
#'
#' # Create the heatmap of robustness
#' heatmap.robustness(robust = robust, drug.names = interv.names)
#'
#' \dontshow{
#' closeAllConnections()
#' }
#'
#' @export
heatmap.robustness <- function(robust, drug.names){


  if (any(is.na(robust))) {
    stop("Missing participant outcome data have *not* been collected. This function cannot be used.", call. = F)
  }


  RI <- robust$RI; threshold <- robust$threshold

  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The  argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- (1 + sqrt(1 + 8*length(robust$robust)))/2
    as.character(1:nt)
  } else {
    drug.names
  }


  if(length(drug.names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }


  if (missing(threshold) & is.element(robust$measure, "OR")) {
    threshold <- 0.28
    #message("The value 0.28 was assigned on 'threshold' by default")
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The value 0.28 was assigned on 'threshold' by default", "\033[0m", "\n")))
  } else if (missing(threshold) & is.element(robust$measure, c("MD", "SMD", "ROM"))) {
    threshold <- 0.17
    #message("The value 0.17 was assigned on 'threshold' by default")
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The value 0.17 was assigned on 'threshold' by default", "\033[0m", "\n")))
  } else {
    threshold <- threshold
    #message(paste("The value", threshold, "was assigned on 'threshold' for", effect.measure.name(full$measure)))
    message(cat(paste0("\033[0;", col = 32, "m", txt = paste("The value", threshold, "was assigned on 'threshold' for", effect.measure.name(robust$measure)), "\033[0m", "\n")))
  }


  ## Lower triangular heatmap matrix - Comparisons are read from the left to the right
  ## CAREFUL: The interventions in the drug.names should follow the order you considered to run NMA pattern-mixture model!
  mat <- matrix(NA, nrow = length(drug.names) - 1, ncol = length(drug.names) - 1)
  mat[lower.tri(mat, diag = T)] <- sprintf("%.2f", RI)
  colnames(mat) <- drug.names[1:(length(drug.names) - 1)]; rownames(mat) <- drug.names[2:length(drug.names)]
  mat.new <- melt(mat, na.rm = T)


  ## Create the heatmap for one network of interventions
  p <- ggplot(mat.new, aes(factor(Var2, levels = drug.names[1:(length(drug.names) - 1)]), factor(Var1, levels = drug.names[length(drug.names):2]), fill = ifelse(value < threshold, "high", "poor"))) +
         geom_tile(colour = "white") +
         geom_text(aes(factor(Var2, levels = drug.names[1:(length(drug.names) - 1)]), factor(Var1, levels = drug.names[length(drug.names):2]), label = value, fontface = "bold"), size = rel(4.5)) +
         scale_fill_manual(breaks = c("high", "poor"), values = c("#009E73", "#D55E00")) +
         scale_x_discrete(position = "top") +
         labs(x = "", y = "") +
         theme(legend.position = "none", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

  return(p)
}

