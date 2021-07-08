#' Enhanced balloon plot: Investigate the impact of missing participant outcome data
#'
#' @description This functions creates the enhanced balloon plot with the summary results of a pairwise comparison under different scenarios about the missingness parameter
#'   for a pair of interventions. Currently, \code{balloon.plot.mod} supports the illustration of the summary effect size and its standard error.
#'
#' @param sens An object of S3 class \code{\link{run.sensitivity}}. See 'Value' in \code{\link{run.sensitivity}}.
#' @param compar A positive integer that indicates the pairwise comparison of interest. For a pairwise meta-analysis, \code{compar} is inherently equal to 1.
#'   For a network meta-analysis, the user must first use the \code{possible.observed.comparisons} function
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
#'   Then, it returns a lower triangular heatmap matrix with the robustness index valueof all possible pairwise comparisons.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.sensitivity}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of primary analysis results: A case study on missing outcome data in pairwise and network meta-analysis. \emph{Res Synth Methods} 2021;\bold{12}(4):475--490. [\doi{10.1002/jrsm.1478}]
#'
#' @examples
#' data("nma.liu2013.RData")
#'
#' # Perform the sensitivity analysis (using the 'default' of the argument 'mean.scenarios')
#' res.sens <- run.sensitivity(full = res1, assumption = "IDE-ARM", var.misspar = 1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("budesodine", "budesodine plus formoterol", "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium", "placebo")
#'
#' # Create the enhanced balloon plot for the comparison 'X'
#' balloon.plot.mod(sens = res.sens, compar = 1, drug.names = interv.names)
#'
#' @export
balloon.plot.mod <- function(sens, compar, drug.names){


  if (is.na(sens)) {
    stop("Missing participant outcome data have *not* been collected. This function cannot be used.", call. = F)
  }

  ES.all <- sens$EM; D <- sens$D; scenarios <- sens$scenarios


  ## Define the position and number of the scenarios
  nt <- (1 + sqrt(1 + 8*(length(ES.all[, 1])/length(scenarios)^2)))/2  # The quadratic formula for the roots of the general quadratic equation


  outcome <- if (is.element(sens$measure, c("MD", "SMD", "ROM"))) {
    "continuous"
  } else {
    "binary"
  }


  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
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


  ## Each parameter is a matrix with rows referring to the scenarios and columns referring to the possible comparisons
  signif <- upper.ES <- lower.ES <- ES.stand <- sd.ES <- ES <- matrix(NA, nrow = length(scenarios)^2, ncol = (nt*(nt - 1))/2)
  for(i in 1:(length(scenarios)^2)){
    for(j in 1:(nt*(nt - 1))/2){

      # Effect estimate (e.g. posterior mean of an effect measure)
      ES[i, j] <- round(ES.all[j + (nt*(nt - 1)/2)*(i - 1), 1], 2)

      # Uncertainty arounf the effect estimate (e.g. posterior standard deviation of an effect measure)
      sd.ES[i, j] <- round(ES.all[j + (nt*(nt - 1)/2)*(i - 1), 2], 2)

      # Standardise the effect estimate according to the outcome type (positive or negative)
      ES.stand[i, j] <- ifelse(D == 1, ES[i, j]/sd.ES[i, j], -ES[i, j]/sd.ES[i, j])

      # Lower bound of the 95% (credible or confidence) interval of the effect estimate
      lower.ES[i, j] <- ES.all[j + (nt*(nt - 1)/2)*(i - 1), 3]

      # Upper bound of the 95% (credible or confidence) interval of the effect estimate
      upper.ES[i, j] <- ES.all[j + (nt*(nt - 1)/2)*(i - 1), 4]

      # Dummy variable to indicate presence or absence of statistical significance
      signif[i, j] <- ifelse(lower.ES[i, j] < 0 & upper.ES[i, j] > 0, "yes", "no")

    }
  }


  ## Normalise the standardised effect estimate (ES.stand). We need this to weight the bubbles in the balloon plot (see, geom_point below).
  ES.normalised <- (ES.stand - min(ES.stand))/(max(ES.stand) - min(ES.stand))


  ## Indicate all possible comparisons (necessary for NMA)
  ## CAREFUL: The interventions in the drug.names should follow the order you considered to run NMA pattern-mixture model!
  comparison <- matrix(combn(drug.names, 2), nrow = length(combn(drug.names, 2))/2, ncol = 2, byrow = T)
  comparison.numeric <- matrix(combn(1:nt, 2), nrow = length(combn(1:nt, 2))/2, ncol = 2, byrow = T)


  ## Indicate all combinations of scenarios for the 'active vs control' comparison
  missp <- data.frame(rep(scenarios, each = 5), rep(scenarios, 5));colnames(missp) <- c("active", "control")


  ## Now, bring all necessary input data in a dataframe to proceed with the creation of the balloon plot (via ggplot2)
  ## Do not forget to check the heatmap, first, to decide on the 'compar' for NMA.
  mat <- data.frame(missp, ES[, compar], sd.ES[, compar], signif[, compar]); colnames(mat) <- c("active", "control", "value", "sd.value", "significance")


  ## Create the proposed balloon plot (separately, for binary and continuous outcomes)
  if(outcome == "binary"){

    bubble <- ggplot(mat, aes(x = active, y = control, color = sd.value, label = sprintf("%.2f", round(exp(value), 2)))) +
                geom_rect(mapping = aes(NULL, NULL, xmin = 1, xmax = 5), ymin = 1, ymax = 5, color = "grey93", fill = "grey93", alpha = 0.1) +
                geom_rect(mapping = aes(NULL, NULL, xmin = 2, xmax = 4), ymin = 2, ymax = 4, color = "grey100", fill = "grey100", alpha = 0.1) +
                geom_point(aes(size = ES.normalised[, compar]), stroke = 2, shape = ifelse(signif[, compar] == "yes", "circle", "circle plus")) +
                scale_size(range = c(0, 30)) +
                geom_text(colour = "black", fontface = "bold", size = 4.5) +
                geom_label(aes(3, 3, label = round(exp(mat[13, 3]), 2)), colour = "black", fontface = "bold",  size = 4.5) +
                scale_color_gradient(low = "deepskyblue", high = "brown1") +
                scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("1/3", "1/2", "1", "2", "3"), position = "bottom", expand = c(0.2, 0)) +
                scale_y_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("1/3", "1/2", "1", "2", "3"), expand = c(0.2, 0)) +
                coord_cartesian(ylim = c(1, 5), clip = 'off') +
                labs(x = paste("IMOR scenario in", comparison[compar, 2]), y = paste("IMOR scenario in", comparison[compar, 1]), color = "") +
                guides(shape = F, size = F) +
                theme_bw() +
                theme(axis.text.x = element_text(size = 12, angle = 360, vjust = 0.8, hjust = 0.5), axis.text.y = element_text(size = 12, vjust = 0.5, hjust = 1),
                      axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, angle = 90, face = "bold"),
                      legend.position = "bottom", legend.text = element_text(size = 12), legend.key.width = unit(1.5, "cm"),
                      legend.title = element_text(size = 12, face = "bold"), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "grey86"))

  } else {

    bubble <- ggplot(mat, aes(x = active, y = control, color = sd.value, label = sprintf("%.2f", round(value, 2)))) +
                geom_rect(mapping = aes(NULL, NULL, xmin = 1, xmax = 5), ymin = 1, ymax = 5, color = "grey93", fill = "grey93", alpha = 0.1) +
                geom_rect(mapping = aes(NULL, NULL, xmin = 2, xmax = 4), ymin = 2, ymax = 4, color = "grey100", fill = "grey100", alpha = 0.1) +
                geom_point(aes(size = ES.normalised[, compar]), stroke = 2, shape = ifelse(signif[, compar] == "yes", "circle", "circle plus")) +
                scale_size(range = c(0, 30)) +
                geom_text(colour = "black", fontface = "bold", size = 4.5) +
                geom_label(aes(3, 3, label = mat[13, 3]), colour = "black", fontface = "bold",  size = 4.5) +
                scale_color_gradient(low = "deepskyblue", high = "brown1") +
                scale_x_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("-2", "-1", "0", "1", "2"), position = "bottom", expand = c(0.2, 0)) +
                scale_y_continuous(breaks = c(1, 2, 3, 4, 5), labels = c("-2", "-1", "0", "1", "2"), expand = c(0.2, 0)) +
                coord_cartesian(ylim = c(1, 5), clip = 'off') +
                labs(x = paste("IMDoM scenario in", comparison[compar, 2]), y = paste("IMDoM scenario in", comparison[compar, 1]), color = "") +
                guides(shape = F, size = F) +
                theme_bw() +
                theme(axis.text.x = element_text(size = 12, angle = 360, vjust = 0.8, hjust = 0.5), axis.text.y = element_text(size = 12, vjust = 0.5, hjust = 1),
                      axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, angle = 90, face = "bold"),
                      legend.position = "bottom", legend.text = element_text(size = 12), legend.key.width = unit(1.5, "cm"),
                      legend.title = element_text(size = 12, face = "bold"), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "grey86"))

  }

  return(bubble)
}
