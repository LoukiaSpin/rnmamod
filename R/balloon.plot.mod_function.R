#' Enhanced balloon plot
#'
#' @description
#'   Creates the enhanced balloon plot for the summary effect size and
#'   between-trial standard deviation, \emph{tau}, under different scenarios
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
#' @return \code{balloon_plot} returns two enhanced balloon plots for one
#'   comparison (see 'Details'):
#'   \tabular{ll}{
#'    \code{Plot_effect_size} \tab The enhanced balloon plot for the
#'    summary effect size (according to the argument \code{measure} inherited
#'    by \code{\link{run_sensitivity}}) of one pairwise comparison.\cr
#'    \tab \cr
#'    \code{Plot_tau} \tab The enhanced balloon plot for tau. When the
#'    fixed-effect model has been performed in \code{\link{run_sensitivity}},
#'    the function will not return the \code{Plot_tau}.\cr
#' }
#'
#' @details
#'   For the \code{Plot_effect_size} of the selected pairwise
#'   comparison, the different colours and sizes of the bubbles reflect the
#'   posterior standard deviation and the posterior mean, respectively.
#'   A colour key appears below the plot.
#'   The size of the bubble is proportional to the corresponding posterior mean.
#'   Crossed bubbles indicate scenarios with conclusive evidence (the
#'   95\% credible interval excludes the null value), and filled bubbles
#'   indicate scenarios with inconclusive evidence (the 95\% credible interval
#'   includes the null value). The missing-at-random assumption (primary
#'   analysis) is labeled in a white frame.
#'   Both axes illustrate the scenarios as specified in the argument
#'   \code{mean_scenarios} of the \code{\link{run_sensitivity}}:
#'   the x-axis refers to the 'experimental' intervention, and the y-axis refers
#'   to the 'control' intervention.
#'
#'   The same enhanced balloon plot is created for tau (\code{Plot_tau}).
#'   However, filled bubbles indicate low statistical heterogeneity
#'   (the posterior median of tau is lower than the median of the
#'   prior distribution for the heterogeneity parameter),
#'   and crossed bubbles indicate considerable statistical heterogeneity
#'   (the posterior median of tau exceeds the median of the prior
#'   distribution for the heterogeneity parameter).
#'
#'   \code{balloon_plot_mod} can be used only when missing participant
#'   outcome data have been extracted for at least one trial.
#'   Otherwise, the execution of the function will be stopped and an error
#'   message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_model}}, \code{\link{run_sensitivity}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of
#' primary analysis results: A case study on missing outcome data in pairwise
#' and network meta-analysis. \emph{Res Synth Methods}
#' 2021;\bold{12}(4):475--490. \doi{10.1002/jrsm.1478}
#'
#' @examples
#' data("pma.taylor2004")
#'
#' # Read results from 'run_sensitivity' (using the default arguments)
#' res_sens <- readRDS(system.file('extdata/res_sens_taylor.rds',
#'                     package = 'rnmamod'))
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "inositol")
#'
#' # Create the enhanced balloon plot for 'inositol versus placebo'
#' balloon_plot(sens = res_sens,
#'              compar = c("inositol", "placebo"),
#'              drug_names = interv_names)
#'
#' @export
balloon_plot <- function(sens, compar, drug_names) {

  if (any(is.na(sens))) {
    stop("No missing outcome data collected; function cannot be used.",
         call. = FALSE)
  }

  es_all <- sens$EM
  D <- sens$D
  scenarios <- sens$scenarios
  measure <- sens$measure
  tau_all <- sens$tau

  # Define the position and number of the scenarios
  nt <- (1 + sqrt(1 + 8 * (length(es_all[, 1]) / length(scenarios)^2))) / 2

  drug_names <- if (missing(drug_names)) {
    stop("The argument 'drug_names' has not been defined", call. = FALSE)
  } else {
    drug_names
  }

  compar <- if (missing(compar)) {
    stop("The argument 'compar' needs to be defined", call. = FALSE)
  } else if (length(drug_names) < 3 & missing(compar)) {
    c(drug_names[2], drug_names[1])
  } else if (!is.element(compar[1], drug_names) |
             !is.element(compar[2], drug_names)) {
    stop("The value of 'compar' is not found in the argument 'drug_names'",
         call. = FALSE)
  } else if (is.element(compar[1], drug_names) &
             is.element(compar[2], drug_names) &
             match(compar[1], drug_names) < match(compar[2], drug_names)) {
    stop("Re-arrange the order of the element in the argument 'compar'",
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

  # A matrix: rows for the scenarios and columns for the possible comparisons
  es <- matrix(NA, nrow = length(scenarios)^2, ncol = (nt * (nt - 1)) / 2)
  signif <- upper_es <- lower_es <- es_stand <- sd_es <- es
  for (i in 1:(length(scenarios)^2)) {
    for (j in 1:(nt * (nt - 1)) / 2) {
      # Posterior mean of an effect measure
      es[i, j] <- round(es_all[j + (nt * (nt - 1) / 2) * (i - 1), 1], 2)

      # Posterior standard deviation of an effect measure
      sd_es[i, j] <- round(es_all[j + (nt * (nt - 1) / 2) * (i - 1), 2], 2)

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
  colnames(mat) <- c("active", "control", "value", "sd_value", "significance")

  # Enhanced balloon plot for summary effect size of the selected comparison
  labels <- if (is.element(measure, c("OR", "ROM"))) {
    as.character(fractions(exp(scenarios)))
  } else {
    as.character(scenarios)
  }

  imp <- if (is.element(measure, "OR")) {
    "IMOR"
  } else if (is.element(measure, c("MD", "SMD"))) {
    "IMDoM"
  } else if (is.element(measure, "ROM")) {
    "IMRoM"
  }

  bubble_es <- if (is.element(measure, c("OR", "ROM"))) {
   ggplot(mat, aes(x = active,
                   y = control,
                   color = sd_value,
                   label = sprintf("%.2f", exp(value)))) +
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
                         labels = labels,
                         position = "bottom",
                         expand = c(0.2, 0)) +
      scale_y_continuous(breaks = seq_len(length(scenarios)),
                         labels = labels,
                         expand = c(0.2, 0)) +
      coord_cartesian(ylim = c(1, length(scenarios)), clip = "off") +
      labs(x = paste(imp, "scenario in", experim),
           y = paste(imp, "scenario in", control),
           color = "") +
      guides(shape = "none", size = "none") +
      ggtitle(paste("Summary", effect_measure_name(measure))) +
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
                    color = sd_value,
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
                         labels = labels,
                         position = "bottom",
                         expand = c(0.2, 0)) +
      scale_y_continuous(breaks = seq_len(length(scenarios)),
                         labels = labels,
                         expand = c(0.2, 0)) +
      coord_cartesian(ylim = c(1, length(scenarios)), clip = "off") +
      labs(x = paste(imp, "scenario in", experim),
           y = paste(imp, "scenario in", control),
           color = "") +
      guides(shape = "none", size = "none") +
      ggtitle(paste("Summary", effect_measure_name(measure))) +
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
    tau <- round(tau_all[, 1], 2)

    # Posterior standard deviation of tau
    sd_tau <- round(tau_all[, 2], 2)

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
    colnames(mat_tau) <- c("active", "control", "value", "sd_value", "extent")
  }

  # Enhanced balloon plot for the between-trial standard deviation
  axis.name.x <- if (is.element(measure, "OR")) {
    "IMOR scenario in experimental"
  } else if (is.element(measure, c("MD", "SMD"))) {
    "IMDoM scenario in experimental"
  } else if (is.element(measure, "ROM")) {
    "IMRoM scenario in experimental"
  }
  axis.name.y <- if (is.element(measure, "OR")) {
    "IMOR scenario in control"
  } else if (is.element(measure, c("MD", "SMD"))) {
    "IMDoM scenario in control"
  } else if (is.element(measure, "ROM")) {
    "IMRoM scenario in control"
  }

  bubble_tau <- if (!is.null(sens$heter)) {
    ggplot(mat_tau, aes(x = active,
                        y = control,
                        color = sd_value,
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
                         labels = labels,
                         position = "bottom",
                         expand = c(0.2, 0)) +
      scale_y_continuous(breaks = seq_len(length(scenarios)),
                         labels = labels,
                         expand = c(0.2, 0)) +
      coord_cartesian(ylim = c(1, length(scenarios)), clip = "off") +
      labs(x = axis.name.x,
           y = axis.name.y,
           color = "") +
      ggtitle("Between-trial standard deviation") +
      guides(shape = "none", size = "none") +
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
    list(plot_effect_size = bubble_es,
         plot_tau = bubble_tau)
  } else {
    bubble_es
  }

  return(results)
}
