#' Enhanced balloon plot (missing participant outcome data)
#'
#' @description
#'   Creates the enhanced balloon plot for the summary effect size and
#'   between-trial standard deviation, \eqn{\tau}, under different scenarios
#'   about the missingness parameter for a pair of interventions.
#'   \code{balloon_plot} uses the scenarios considered in
#'   \code{\link{run_sensitivity}}.
#'
#' @param sens An object of S3 class \code{\link{run_sensitivity}}.
#'   See 'Value' in \code{\link{run_sensitivity}}.
#' @param compar A character vector with two elements indicating the pairwise
#'   comparison of interest. The first element refers to the 'experimental'
#'   and the second element to the 'control' intervention of the comparison.
#' @param drug_names A vector of labels with the name of the interventions
#'   in the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If the argument \code{drug_names} is not defined,
#'   the order of the interventions as they appear in
#'   \code{data} is used, instead.
#'
#' @details
#'   For the summary effect size of the selected pairwise comparison,
#'   the different colour and size of the bubbles reflect the posterior
#'   standard deviation and the magnitude of the posterior mean, respectively.
#'   A colour key appears below the plot.
#'   The size of the bubble is proportional to the corresponding posterior mean.
#'   Crossed bubbles indicate scenarios with conclusive evidence (the 95\%
#'   credible interval excludes the null value), and filled bubbles indicate
#'   scenarios with inconclusive evidence (the 95\% credible interval includes
#'   the null value). The missing-at-random assumption (primary analysis) is
#'   labeled in a white frame.
#'   Both axes illustrate the scenarios as specified in the argument
#'   \code{mean_scenarios} of the \code{run_sensitivity} function:
#'   the x-axis refers to the 'experimental' intervention, and the y-axis refers
#'   to the 'control' intervention.
#'
#'   The same enhanced balloon plot is created for \eqn{\tau}.
#'   However, filled bubbles indicate low statistical heterogeneity
#'   (the posterior median of \eqn{\tau} is lower than the median of the
#'   prior distribution for the heterogeneity parameter),
#'   and crossed bubbles indicate considerable statistical heterogeneity
#'   (the posterior median of \eqn{\tau} exceeds the median of the prior
#'   distribution).
#'
#'   \code{balloon_plot_mod} can be used only for a network of interventions
#'   and when missing outcome data have been extracted for at least one trial.
#'   Otherwise, the execution of the function will be stopped and an error
#'   message will be printed in the R console.
#'
#' @return \code{balloon_plot} A list of two enhanced balloon plots for one
#'   comparison (see 'Details' for the description of the enhanced
#'   balloon plot):
#'   \tabular{ll}{
#'    \code{Plot_effect_size} \tab The enhanced balloon plot for the
#'    summary effect size of one pairwise comparison.\cr
#'    \tab \cr
#'    \code{Plot_tau} \tab The enhanced balloon plot for \eqn{\tau}. When the
#'    fixed-effect model has been performed in \code{run_sensitivity},
#'    the function will not return the \code{Plot_tau}.\cr
#' }
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_sensitivity}}, \code{\link{run_model}}
#'
#' @references
#'   Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of
#'   primary analysis results: A case study on missing outcome data in pairwise
#'   and network meta-analysis. \emph{Res Synth Methods}
#'   2021;\bold{12}(4):475--490. [\doi{10.1002/jrsm.1478}]
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
#'                  var_misspar = 1,
#'                  D = 1,
#'                  n_chains = 3,
#'                  n_iter = 10000,
#'                  n_burnin = 1000,
#'                  n_thin = 1)
#'
#' # Perform the sensitivity analysis (missing-at-random assumption)
#' res.sens <- run_sensitivity(full = res,
#'                             var_misspar = 1,
#'                             n_chains = 3,
#'                             n_iter = 10000,
#'                             n_burnin = 1000,
#'                             n_thin = 1)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("placebo", "budesonide", "budesonide plus formoterol",
#'                   "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium")
#'
#' # Create the enhanced balloon plot for 'tiotropium versus salmeterol'
#' balloon_plot(sens = res.sens,
#'              compar = c("tiotropium", "salmeterol"),
#'              drug_names = interv.names)
#' }
#' @export
balloon_plot <- function(sens, compar, drug_names) {

  if (any(is.na(sens))) {
    stop("No missing outcome data collected; function cannot be used.",
         call. = F)
  }

  es_all <- sens$EM
  D <- sens$D
  scenarios <- sens$scenarios
  measure <- sens$measure
  tau_all <- sens$tau

  # Define the position and number of the scenarios
  nt <- (1 + sqrt(1 + 8 * (length(es_all[, 1]) / length(scenarios)^2))) / 2

  drug_names <- if (missing(drug_names)) {
    message(cat(paste0("\033[0;", col = 32, "m",
                       txt = "The argument 'drug_names' has not been defined.
                       The intervention ID, as specified in 'data' is used as
                       intervention names",
                       "\033[0m", "\n")))
    as.character(1:nt)
  } else {
    drug_names
  }

  compar <- if (length(drug_names) > 2 & missing(compar)) {
    stop("The argument 'compar' needs to be defined", call. = F)
  } else if (length(drug_names) < 3 & missing(compar)) {
    c(comparison[2], comparison[2])
  } else {
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

  # A matrix: rows for the scenarios and columns for the possible comparisons
  es <- matrix(NA, nrow = length(scenarios)^2, ncol = (nt * (nt - 1)) / 2)
  signif <- upper_es <- lower_es <- es_stand <- sd_es <- es
  for (i in 1:(length(scenarios)^2)) {
    for (j in 1:(nt * (nt - 1)) / 2) {
      # Posterior mean of an effect measure
      es[i, j] <- round(es_all[j + (nt * (nt - 1) / 2) * (i - 1), 1], 3)

      # Posterior standard deviation of an effect measure
      sd_es[i, j] <- round(es_all[j + (nt * (nt - 1) / 2) * (i - 1), 2], 3)

      # Standardise the effect estimate based on the outcome direction
      es_stand[i, j] <- ifelse(D == 1, es[i, j] / sd_es[i, j],
                               -es[i, j] / sd_es[i, j])

      # Lower bound of the 95% credible interval of the effect estimate
      lower_es[i, j] <- es_all[j + (nt * (nt - 1) / 2) * (i - 1), 3]

      # Upper bound of the 95% credible interval of the effect estimate
      upper_es[i, j] <- es_all[j + (nt * (nt - 1) / 2) * (i - 1), 4]

      # Dummy variable to indicate present or absent statistical significance
      signif[i, j] <- ifelse(lower_es[i, j] < 0 & upper_es[i, j] > 0,
                             "yes", "no")
    }
  }

  # Normalise the standardised effect estimate (es_stand)
  # We need this to weight the bubbles in the balloon plot (via geom_point)
  es_normalised <- round((es_stand - min(es_stand)) /
                           (max(es_stand) - min(es_stand)), 2)

  # All combinations of scenarios for the (active vs control) comparison
  missp <- data.frame(rep(seq_len(length(scenarios)), each = length(scenarios)),
                      rep(seq_len(length(scenarios)), length(scenarios)))
  colnames(missp) <- c("active", "control")

  # Prepare dataframe to create the enhanced balloon plot (via ggplot2)
  mat <- data.frame(missp,
                    es[, compar_id],
                    sd_es[, compar_id],
                    signif[, compar_id])
  colnames(mat) <- c("active", "control", "value", "sd.value", "significance")

  # Enhanced balloon plot for summary effect size of the selected comparison
  bubble_es <- if (is.element(measure, c("OR", "ROM"))) {
   ggplot(mat, aes(x = active,
                   y = control,
                   color = sd.value,
                   label = sprintf("%.2f", round(exp(value), 2)))) +
      geom_rect(mapping = aes(x = NULL,
                              y = NULL,
                              xmin = 1,
                              xmax = length(scenarios),
                              ymin = 1,
                              ymax = length(scenarios)),
                color = "grey93", fill = "grey93", alpha = 0.1) +
      geom_rect(mapping = aes(x = NULL,
                              y = NULL,
                              xmin = 2,
                              xmax = length(scenarios) - 1,
                              ymin = 2,
                              ymax = length(scenarios) - 1),
                color = "grey100", fill = "grey100", alpha = 0.1) +
      geom_point(aes(size = es_normalised[, compar_id]),
                 stroke = 2,
                 shape = ifelse(signif[, compar_id] == "yes", "circle",
                                "circle plus")) +
      scale_size(range = c(0, 30)) +
      geom_text(colour = "black", fontface = "bold", size = 4.5) +
      geom_label(aes(median(order(scenarios)),
                     median(order(scenarios)),
                     label = round(exp(mat[median(1:(length(scenarios)^2)), 3]),
                                   2)),
                 colour = "black", fontface = "bold",  size = 4.5) +
      scale_color_gradient(low = "deepskyblue", high = "#D55E00") +
      scale_x_continuous(breaks = seq_len(length(scenarios)),
                         labels = as.character(fractions(exp(scenarios))),
                         position = "bottom",
                         expand = c(0.2, 0)) +
      scale_y_continuous(breaks = seq_len(length(scenarios)),
                         labels = as.character(fractions(exp(scenarios))),
                         expand = c(0.2, 0)) +
      coord_cartesian(ylim = c(1, length(scenarios)), clip = "off") +
      labs(x = paste("IMOR scenario in", experim),
           y = paste("IMOR scenario in", control),
           color = "") +
      guides(shape = F, size = F) +
      ggtitle(paste("Summary", effect.measure.name(measure))) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 12,
                                       angle = 360,
                                       vjust = 0.8,
                                       hjust = 0.5),
            axis.text.y = element_text(size = 12,
                                       vjust = 0.5,
                                       hjust = 1),
            axis.title.x = element_text(size = 12,
                                        face = "bold"),
            axis.title.y = element_text(size = 12,
                                        angle = 90,
                                        face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.key.width = unit(1.5, "cm"),
            legend.title = element_text(size = 12,
                                        face = "bold"),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "grey86"))
  } else if (is.element(measure, c("MD", "SMD"))) {
    ggplot(mat, aes(x = active,
                    y = control,
                    color = sd.value,
                    label = sprintf("%.2f", round(value, 2)))) +
      geom_rect(mapping = aes(x = NULL,
                              y = NULL,
                              xmin = 1,
                              xmax = length(scenarios),
                              ymin = 1,
                              ymax = length(scenarios)),
                color = "grey93", fill = "grey93", alpha = 0.1) +
      geom_rect(mapping = aes(x = NULL,
                              y = NULL,
                              xmin = 2,
                              xmax = length(scenarios) - 1,
                              ymin = 2,
                              ymax = length(scenarios) - 1),
                color = "grey100", fill = "grey100", alpha = 0.1) +
      geom_point(aes(size = es_normalised[, compar_id]),
                 stroke = 2,
                 shape = ifelse(signif[, compar_id] == "yes", "circle",
                                "circle plus")) +
      scale_size(range = c(0, 30)) +
      geom_text(colour = "black", fontface = "bold", size = 4.5) +
      geom_label(aes(median(order(scenarios)),
                     median(order(scenarios)),
                     label = mat[median(1:(length(scenarios)^2)), 3]),
                 colour = "black", fontface = "bold",  size = 4.5) +
      scale_color_gradient(low = "deepskyblue", high = "#D55E00") +
      scale_x_continuous(breaks = seq_len(length(scenarios)),
                         labels = as.character(scenarios),
                         position = "bottom",
                         expand = c(0.2, 0)) +
      scale_y_continuous(breaks = seq_len(length(scenarios)),
                         labels = as.character(scenarios),
                         expand = c(0.2, 0)) +
      coord_cartesian(ylim = c(1, length(scenarios)), clip = "off") +
      labs(x = paste("IMDoM scenario in", experim),
           y = paste("IMDoM scenario in", control),
           color = "") +
      guides(shape = F, size = F) +
      ggtitle(paste("Summary", effect.measure.name(measure))) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 12,
                                       angle = 360,
                                       vjust = 0.8,
                                       hjust = 0.5),
            axis.text.y = element_text(size = 12,
                                       vjust = 0.5,
                                       hjust = 1),
            axis.title.x = element_text(size = 12,
                                        face = "bold"),
            axis.title.y = element_text(size = 12,
                                        angle = 90,
                                        face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.key.width = unit(1.5, "cm"),
            legend.title = element_text(size = 12,
                                        face = "bold"),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "grey86"))
  }

  if (!is.null(sens$heter)) {
    # Posterior mean of tau
    tau <- tau_all[, 1]

    # Posterior standard deviation of tau
    sd_tau <- tau_all[, 2]

    # Dummy variable to indicate the extent of tau
    median_extent <- if (sens$heter[3] == 1) {
      fdrtool::qhalfnorm(0.5, theta = (1 / sens$heter[2]) * sqrt(pi / 2))
    } else if (sens$heter[3] == 2) {
      qunif(10000, 0, sens$heter[2])
    } else if (sens$heter[3] == 3) {
      sqrt(qlnorm(0.5, sens$heter[1], sqrt(1 / sens$heter[2])))
    } else if (sens$heter[3] == 4) {
      sqrt(exp(qt(0.5, df = 5) * sqrt(1 / sens$heter[2]) + sens$heter[1]))
    }
    extent_tau <- ifelse(tau < median_extent, "low", "considerable")

    # Normalise tau2. We need this to weight the bubbles in the balloon plot
    tau_normalised <- round((tau - min(tau)) / (max(tau) - min(tau)), 2)

    # Prepare dataset for the ggplot2
    mat_tau <- data.frame(missp, tau, sd_tau, extent_tau)
    colnames(mat_tau) <- c("active", "control", "value", "sd.value", "extent")
  }

  # Enhanced balloon plot for the between-trial standard deviation
  bubble_tau <- if (!is.null(sens$heter)) {
    ggplot(mat_tau, aes(x = active,
                        y = control,
                        color = sd.value,
                        label = sprintf("%.2f", value))) +
      geom_rect(mapping = aes(x = NULL,
                              y = NULL,
                              xmin = 1,
                              xmax = length(scenarios),
                              ymin = 1,
                              ymax = length(scenarios)),
                color = "grey93", fill = "grey93", alpha = 0.1) +
      geom_rect(mapping = aes(x = NULL,
                              y = NULL,
                              xmin = 2,
                              xmax = length(scenarios) - 1,
                              ymin = 2,
                              ymax = length(scenarios) - 1),
                color = "grey100", fill = "grey100", alpha = 0.1) +
      geom_point(aes(size = tau_normalised),
                 stroke = 2,
                 shape = ifelse(extent_tau == "low", "circle", "circle plus")) +
      scale_size(range = c(0, 30)) +
      geom_text(colour = "black", fontface = "bold", size = 4.5) +
      geom_label(aes(median(order(scenarios)),
                     median(order(scenarios)),
                     label = sprintf("%.2f",
                                     mat_tau[median(1:(length(scenarios)^2)),
                                             3])),
                 colour = "black", fontface = "bold",  size = 4.5) +
      scale_color_gradient(low = "deepskyblue", high = "#D55E00") +
      scale_x_continuous(breaks = seq_len(length(scenarios)),
                         labels = as.character(fractions(exp(scenarios))),
                         position = "bottom",
                         expand = c(0.2, 0)) +
      scale_y_continuous(breaks = seq_len(length(scenarios)),
                         labels = as.character(fractions(exp(scenarios))),
                         expand = c(0.2, 0)) +
      coord_cartesian(ylim = c(1, length(scenarios)), clip = "off") +
      labs(x = paste("IMOR scenario in experimental"),
           y = paste("IMOR scenario in control"),
           color = "") +
      ggtitle("Between-trial standard deviation") +
      guides(shape = F, size = F) +
      theme_bw() +
      theme(axis.text.x = element_text(size = 12,
                                       angle = 360,
                                       vjust = 0.8,
                                       hjust = 0.5),
            axis.text.y = element_text(size = 12,
                                       vjust = 0.5,
                                       hjust = 1),
            axis.title.x = element_text(size = 12,
                                        face = "bold"),
            axis.title.y = element_text(size = 12,
                                        angle = 90,
                                        face = "bold"),
            legend.position = "bottom",
            legend.text = element_text(size = 12),
            legend.key.width = unit(1.5, "cm"),
            legend.title = element_text(size = 12,
                                        face = "bold"),
            panel.grid.minor = element_blank(),
            panel.background = element_rect(fill = "grey86"))
  } else {
    NA
  }

  results <- if (!is.null(sens$heter)) {
    list(Plot_effect_size = bubble_ES,
         Plot_tau = bubble_tau)
    cat("Plot_effect_size: Enhance balloon plot for
        the effect measure\n")
    invisible(bubble_ES)
    cat("\n")

    cat("Plot_tau: Enhance balloon plot for
        the between-trial standard deviation\n")
    invisible(bubble_tau)
  } else {
    bubble_ES
    cat("Plot_effect_size: Enhance balloon plot for
        the effect measure \n")
    invisible(bubble_ES)
  }

  return(results)
}
