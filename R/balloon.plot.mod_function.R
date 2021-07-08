#' Enhanced balloon plot: Investigate the impact of missing participant outcome data
#'
#' @description This functions creates the enhanced balloon plot with the summary results of a pairwise comparison under different scenarios about the missingness parameter
#'   for a pair of interventions. Currently, \code{balloon.plot.mod} supports the illustration of the summary effect size and its standard error.
#'
#' @param sens An object of S3 class \code{\link{run.sensitivity}}. See 'Value' in \code{\link{run.sensitivity}}.
#' @param compar A character vector with two characters that indicates the pairwise comparison of interest. The first character refers to the 'experimental' intervention
#'   and the second character refers to the 'control' intervention of the comparison.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @details The different colours and bubble sizes reflect the posterior standard deviation of the summary effect size and the magnitude of the posterior mean of the summary effect size for one pairwise comparison.
#'   A colour key appears below the plot that matches the colours with the corresponding posterior standard deviation value.
#'   The size of the bubble is proportional to the magnitude of the corresponding effect size.
#'   Crossed bubbles indicate scenrios with conclusive evidence (i.e., the 95% credible interval excludes the null value), and filled bubbles indicate scenarios with inconclusive evidence
#'   (i.e., the 95% credible interval includes the null value).
#'   The missing-at-random assumption (the primary analysis) is labelled in a white frame.
#'   Both axes illustrate the scenarios as specified in the argument \code{mean.scenarios} of the \code{run.sensitivity} function:
#'   the x-axis refers to the 'experimental' intervention, and the y-axis refers to the 'control' intervention.
#'
#'   \code{heatmap.robustness} can be used only for a network of interventions and when missing participant outcome data have been extracted for at least one trial.
#'   Otherwise, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @return \code{balloon.plot.mod} a balloon plot for one comparison which is enhanced with colours and bubbles of different size and shape. See 'Details'.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.sensitivity}}, \code{\link{run.model}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of primary analysis results: A case study on missing outcome data in pairwise and network meta-analysis. \emph{Res Synth Methods} 2021;\bold{12}(4):475--490. [\doi{10.1002/jrsm.1478}]
#'
#' @examples
#' data("nma.baker2009.RData")
#'
#' #' # Perform a random-effects network meta-analysis
#' res1 <- run.model(data = nma.baker2009, measure = "OR", model = "RE", assumption = "IDE-ARM", heter.prior = list("halfnormal", 0, 1), mean.misspar = 0, var.misspar = 1, D = 1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' # Perform the sensitivity analysis (using the 'default' of the argument 'mean.scenarios')
#' res.sens <- run.sensitivity(full = res1, assumption = "IDE-ARM", var.misspar = 1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("budesodine", "budesodine plus formoterol", "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium", "placebo")
#'
#' # Create the enhanced balloon plot for the comparison 'tiotropium versus salmeterol'
#' balloon.plot.mod(sens = res.sens, compar = c("tiotropium", "salmeterol"), drug.names = interv.names)
#'
#' @export
balloon.plot.mod <- function(sens, compar, drug.names){


  if (is.na(sens)) {
    stop("Missing participant outcome data have *not* been collected. This function cannot be used.", call. = F)
  }

  ES.all <- sens$EM; D <- sens$D; scenarios <- sens$scenarios; measure <- sens$measure


  ## Define the position and number of the scenarios
  nt <- (1 + sqrt(1 + 8*(length(ES.all[, 1])/length(scenarios)^2)))/2  # The quadratic formula for the roots of the general quadratic equation


  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
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


  ## Indicate all combinations of scenarios for the 'active vs control' comparison
  missp <- data.frame(rep(1:length(scenarios), each = length(scenarios)), rep(1:length(scenarios), length(scenarios)));colnames(missp) <- c("active", "control")


  ## Now, bring all necessary input data in a dataframe to proceed with the creation of the balloon plot (via ggplot2)
  ## Do not forget to check the heatmap, first, to decide on the 'compar' for NMA.
  mat <- data.frame(missp, ES[, compar.id], sd.ES[, compar.id], signif[, compar.id]); colnames(mat) <- c("active", "control", "value", "sd.value", "significance")


  ## Create the proposed balloon plot (separately, for binary and continuous outcomes)
  bubble <- if (is.element(measure, c("OR", "ROM"))) {
   ggplot(mat, aes(x = active, y = control, color = sd.value, label = sprintf("%.2f", round(exp(value), 2)))) +
     geom_rect(mapping = aes(NULL, NULL, xmin = 1, xmax = length(scenarios)), ymin = 1, ymax = length(scenarios), color = "grey93", fill = "grey93", alpha = 0.1) +
     geom_rect(mapping = aes(NULL, NULL, xmin = 2, xmax = length(scenarios) - 1), ymin = 2, ymax = length(scenarios) - 1, color = "grey100", fill = "grey100", alpha = 0.1) +
     geom_point(aes(size = ES.normalised[, compar.id]), stroke = 2, shape = ifelse(signif[, compar.id] == "yes", "circle", "circle plus")) +
     scale_size(range = c(0, 30)) +
     geom_text(colour = "black", fontface = "bold", size = 4.5) +
     geom_label(aes(median(order(scenarios)), median(order(scenarios)), label = round(exp(mat[median(1:(length(scenarios)^2)), 3]), 2)), colour = "black", fontface = "bold",  size = 4.5) +
     scale_color_gradient(low = "deepskyblue", high = "brown1") +
     scale_x_continuous(breaks = 1:length(scenarios), labels = as.character(fractions(exp(scenarios))), position = "bottom", expand = c(0.2, 0)) +
     scale_y_continuous(breaks = 1:length(scenarios), labels = as.character(fractions(exp(scenarios))), expand = c(0.2, 0)) +
     coord_cartesian(ylim = c(1, length(scenarios)), clip = 'off') +
     labs(x = paste("IMOR scenario in", experim), y = paste("IMOR scenario in", control), color = "") +
     guides(shape = F, size = F) +
     theme_bw() +
     theme(axis.text.x = element_text(size = 12, angle = 360, vjust = 0.8, hjust = 0.5), axis.text.y = element_text(size = 12, vjust = 0.5, hjust = 1),
           axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, angle = 90, face = "bold"),
           legend.position = "bottom", legend.text = element_text(size = 12), legend.key.width = unit(1.5, "cm"),
           legend.title = element_text(size = 12, face = "bold"), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "grey86"))
  } else if (is.element(measure, c("MD", "SMD"))) {
    ggplot(mat, aes(x = active, y = control, color = sd.value, label = sprintf("%.2f", round(value, 2)))) +
      geom_rect(mapping = aes(NULL, NULL, xmin = 1, xmax = length(scenarios)), ymin = 1, ymax = length(scenarios), color = "grey93", fill = "grey93", alpha = 0.1) +
      geom_rect(mapping = aes(NULL, NULL, xmin = 2, xmax = length(scenarios) - 1), ymin = 2, ymax = length(scenarios) - 1, color = "grey100", fill = "grey100", alpha = 0.1) +
      geom_point(aes(size = ES.normalised[, compar.id]), stroke = 2, shape = ifelse(signif[, compar.id] == "yes", "circle", "circle plus")) +
      scale_size(range = c(0, 30)) +
      geom_text(colour = "black", fontface = "bold", size = 4.5) +
      geom_label(aes(median(order(scenarios)), median(order(scenarios)), label = mat[median(1:(length(scenarios)^2)), 3]), colour = "black", fontface = "bold",  size = 4.5) +
      scale_color_gradient(low = "deepskyblue", high = "brown1") +
      scale_x_continuous(breaks = 1:length(scenarios), labels = as.character(scenarios), position = "bottom", expand = c(0.2, 0)) +
      scale_y_continuous(breaks = 1:length(scenarios), labels = as.character(scenarios), expand = c(0.2, 0)) +
      coord_cartesian(ylim = c(1, length(scenarios)), clip = 'off') +
      labs(x = paste("IMDoM scenario in", experim), y = paste("IMDoM scenario in", control), color = "") +
      guides(shape = F, size = F) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 12, angle = 360, vjust = 0.8, hjust = 0.5), axis.text.y = element_text(size = 12, vjust = 0.5, hjust = 1),
            axis.title.x = element_text(size = 12, face = "bold"), axis.title.y = element_text(size = 12, angle = 90, face = "bold"),
            legend.position = "bottom", legend.text = element_text(size = 12), legend.key.width = unit(1.5, "cm"),
            legend.title = element_text(size = 12, face = "bold"), panel.grid.minor = element_blank(), panel.background = element_rect(fill = "grey86"))
  }

  return(bubble)}
