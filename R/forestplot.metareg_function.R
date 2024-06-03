#' Comparator-specific forest plot for network meta-regression
#'
#' @description
#'   Provides a forest plot with the posterior median and 95\% credible
#'   and prediction intervals for comparisons with the selected intervention
#'   (comparator) in the network under the network meta-analysis \emph{and}
#'   network meta-regression for a specified level or value of the investigated
#'   covariate, and a forest plot with the corresponding SUCRA values.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param reg An object of S3 class \code{\link{run_metareg}}. See 'Value' in
#'   \code{\link{run_metareg}}.
#' @param compar A character to indicate the comparator intervention. It must
#'   be any name found in \code{drug_names}.
#' @param cov_value A list of two elements in the following order: a number
#'   for the covariate value of interest (see 'Arguments' in
#'   \code{\link{run_metareg}}), and a character to indicate the name of
#'   the covariate. See also 'Details'.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#'
#' @return A panel of two forest plots: (1) a forest plot on the effect
#'   estimates and predictions of comparisons with the selected intervention in
#'   the network under the network meta-analysis and network meta-regression for
#'   a specified level or value of the investigated covariate, and (2) a forest
#'   plot on the posterior mean and 95\% credible interval of SUCRA values of
#'   the interventions (Salanti et al., 2011).
#'
#' @details The y-axis of the forest plot on \bold{effect sizes} displays the
#'   labels of the interventions in the network; the selected intervention that
#'   comprises the \code{compar} argument is annotated in the plot with the
#'   label 'Comparator intervention'.
#'   For each comparison with the selected intervention, the 95\% credible and
#'   prediction intervals are displayed as overlapping lines. Black lines refer
#'   to estimation under both analyses. Green and red lines refer to prediction
#'   under network meta-analysis and network meta-regression, respectively. The
#'   corresponding numerical results are displayed above each line:
#'   95\% credible intervals are found in parentheses, and 95\% predictive
#'   intervals are found in brackets. Odds ratios, relative risks, and ratio of
#'   means are reported in the original scale after exponentiation of the
#'   logarithmic scale.
#'
#'   The y-axis for the forest plot on \bold{SUCRA} values displays the
#'   labels of the interventions in the network.
#'   The corresponding numerical results are displayed above each line.
#'   Three coloured rectangles appear in the forest plot: a red rectangle for
#'   SUCRA values up to 50\%, a yellow rectangular for SUCRA values between
#'   50\% and 80\%, and a green rectangle for SUCRA values over 80\%.
#'   Interventions falling at the green area are considered as the highest
#'   ranked interventions, whilst interventions falling at the red area are
#'   considered as the lowest ranked interventions.
#'
#'   In both plots, the interventions are sorted in the descending order of
#'   their SUCRA values based on the network meta-analysis.
#'
#'   To obtain the posterior distribution of SUCRAs under the network
#'   meta-regression for a specified level or value of the investigated
#'   covariate, a two-step procedured is followed. First, the posterior median and
#'   standard deviation of the basic parameters and corresponding beta
#'   coefficients are considered to calculate the meadin and standard deviation of
#'   the basic parameters at the selected level or value of the investigated
#'   covariate: d + beta * cov_value. Then, these calculated values are fed into
#'   the BUGS code to get the posterior distribution of SUCRA. The progress of
#'   the simulation appears on the R console.
#'
#'   \code{forestplot_metareg} is integrated in \code{\link{metareg_plot}}.
#'
#'   \code{forestplot_metareg} can be used only for a network of interventions.
#'   In the case of two interventions, the execution of the function will be
#'   stopped and an error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{metareg_plot}}, \code{\link{run_metareg}},
#'   \code{\link{run_model}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries
#' for presenting results from multiple-treatment meta-analysis: an overview and
#' tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71.
#' \doi{10.1016/j.jclinepi.2010.03.016}
#'
#' @export
forestplot_metareg <- function(full, reg, compar, cov_value, drug_names) {

  if (!inherits(full, "run_model") || is.null(full)) {
    stop("The argument 'full' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  if (!inherits(reg, "run_metareg") || is.null(reg)) {
    stop("The argument 'reg' must be an object of S3 class 'run_metareg'.",
         call. = FALSE)
  }

  if (length(unique(reg$covariate)) < 3 &
      !is.element(cov_value[[1]], reg$covariate)) {
    aa <- "The first element of the argument 'cov_value' is out of the value"
    bb <- "range of the analysed covariate."
    stop(paste(aa, bb), call. = FALSE)
  } else if (length(unique(reg$covariate)) > 2 &
             (cov_value[[1]] < min(reg$covariate) |
              cov_value[[1]] > max(reg$covariate))) {
    aa <- "The first element of the argument 'cov_value' is out of the value"
    bb <- "range of the analysed covariate."
    stop(paste(aa, bb), call. = FALSE)
  }

  drug_names <- if (missing(drug_names)) {
    aa <- "The argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as specified in 'data' is used, instead."
    message(paste(aa, bb))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug_names
  }
  len_drug <- length(drug_names)

  if (length(drug_names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis.",
         call. = FALSE)
  }

  # Sort the drugs by their SUCRA in decreasing order and remove the reference
  # intervention (number 1)
  drug_names_sorted <- drug_names[order(full$SUCRA[, 1], decreasing = TRUE)]

  compar <- if (missing(compar)) {
    stop("The argument 'compar' has not been defined.", call. = FALSE)
  } else if (!is.element(compar, drug_names)) {
    stop("The value of the argument 'compar' is not found in the 'drug_names'.",
         call. = FALSE)
  } else if (is.element(compar, drug_names)) {
    compar
  }

  cov_value <- if (!is.null(reg$beta_all) & missing(cov_value)) {
    stop("The argument 'cov_value' has not been defined.", call. = FALSE)
  } else if (!is.null(reg$beta_all) & length(cov_value) < 2) {
    aa <- "The argument 'cov_value' must be a list with elements a number and"
    stop(paste(aa, "a character."), call. = FALSE)
  } else if (!is.null(reg$beta_all) & length(cov_value) == 2) {
    cov_value
  }

  model <- full$model
  measure <- full$measure
  em_full <- full$EM
  pred_full <- full$EM_pred
  sucra_full <- full$SUCRA

  cov_val <- ifelse(length(unique(reg$covariate)) < 3,
                    cov_value[[1]],
                    cov_value[[1]] - mean(reg$covariate))

  # A matrix with all possible comparisons in the network
  poss_pair_comp1 <- data.frame(exp = t(combn(drug_names, 2))[, 2],
                                comp = t(combn(drug_names, 2))[, 1])
  poss_pair_comp2 <- data.frame(exp = t(combn(drug_names, 2))[, 1],
                                comp = t(combn(drug_names, 2))[, 2])
  poss_pair_comp <- rbind(poss_pair_comp1, poss_pair_comp2)

  # Effect size of all possible pairwise comparisons (NMA)
  em_ref00_nma <- cbind(rbind(data.frame(median = em_full[, 5],
                                         lower = em_full[, 3],
                                         upper = em_full[, 7]),
                              data.frame(median = em_full[, 5] * (-1),
                                         lower = em_full[, 7] * (-1),
                                         upper = em_full[, 3] * (-1))),
                        poss_pair_comp)
  em_subset_nma <- subset(em_ref00_nma, em_ref00_nma[5] == compar)
  em_ref0_nma <- rbind(em_subset_nma[, c(1:3)], c(rep(NA, 3)))
  sucra_new <- data.frame(sucra_full[, 1],
                          drug_names)[order(match(data.frame(sucra_full[, 1],
                                                             drug_names)[, 2],
                                                  em_subset_nma[, 4])), 1]
  em_ref_nma <- em_ref0_nma[order(sucra_new, decreasing = TRUE), ]

  # Effect size of all possible pairwise comparisons (NMR)
  par_median <- as.vector(c(reg$EM[, 5] + reg$beta_all[, 5] * cov_val,
                          (reg$EM[, 5] * (-1)) +
                            (reg$beta_all[, 5] * (-1) * cov_val)))
  par_sd <- as.vector(c(sqrt(((reg$EM[, 2])^2) +
                               ((reg$beta_all[, 2] * cov_val)^2)),
                        sqrt(((reg$EM[, 2])^2) +
                               ((reg$beta_all[, 2] * cov_val)^2))))

  em_ref00_nmr <- cbind(median = par_median,
                        sd = par_sd,
                        lower = par_median - 1.96 * par_sd,
                        upper = par_median + 1.96 * par_sd,
                        poss_pair_comp)
  em_subset_nmr <- subset(em_ref00_nmr, em_ref00_nmr[6] == compar)
  em_ref0_nmr <- rbind(em_subset_nmr[, c(1, 3, 4)], c(rep(NA, 3)))
  em_ref_nmr <- em_ref0_nmr[order(sucra_new, decreasing = TRUE), ]
  rownames(em_ref_nma) <- rownames(em_ref_nmr) <- NULL

  # Posterior results on the predicted estimates of comparisons with the
  # selected comparator as reference
  if (model == "RE") {
    pred_ref00_nma <- cbind(rbind(data.frame(median = pred_full[, 5],
                                             lower = pred_full[, 3],
                                             upper = pred_full[, 7]),
                                  data.frame(median = pred_full[, 5] * (-1),
                                             lower = pred_full[, 7] * (-1),
                                             upper = pred_full[, 3] * (-1))),
                            poss_pair_comp)
    pred_subset_nma <- subset(pred_ref00_nma, pred_ref00_nma[5] == compar)
    pred_ref0_nma <- rbind(pred_subset_nma[, 1:3], c(rep(NA, 3)))
    par_median <- as.vector(c(reg$EM_pred[, 5] + reg$beta_all[, 5] * cov_val,
                            (reg$EM_pred[, 5] * (-1)) +
                              (reg$beta_all[, 5] * (-1) * cov_val)))
    par_sd <- as.vector(c(sqrt(((reg$EM_pred[, 2])^2) +
                                 ((reg$beta_all[, 2] * cov_val)^2)),
                          sqrt(((reg$EM_pred[, 2])^2) +
                                 ((reg$beta_all[, 2] * cov_val)^2))))

    pred_ref00_nmr <-  cbind(data.frame(median = par_median,
                                        lower = par_median - 1.96 * par_sd,
                                        upper = par_median + 1.96 * par_sd),
                           poss_pair_comp)
    pred_subset_nmr <- subset(pred_ref00_nmr, pred_ref00_nmr[5] == compar)
    pred_ref0_nmr <- rbind(pred_subset_nmr[, 1:3], c(rep(NA, 3)))

    # Sort by SUCRA in decreasing order and remove the reference intervention
    pred_ref_nma <- pred_ref0_nma[order(sucra_new, decreasing = TRUE), ]
    pred_ref_nmr <- pred_ref0_nmr[order(sucra_new, decreasing = TRUE), ]
    rownames(pred_ref_nma) <- rownames(pred_ref_nmr) <- NULL
  }

  # Create a data-frame with credible and prediction intervals of comparisons
  # with the reference intervention
  if (!is.element(measure, c("OR", "RR", "ROM")) & model == "RE") {
    prepare_em_nma <- data.frame(as.factor(rep(rev(seq_len(len_drug)), 2)),
                                 rep(drug_names_sorted, 2),
                                 round(rbind(em_ref_nma, pred_ref_nma), 2),
                                 rep(c("Estimation", "Prediction"),
                                     each = length(drug_names)))
    prepare_em_nmr <- data.frame(as.factor(rep(rev(seq_len(len_drug)), 2)),
                                 rep(drug_names_sorted, 2),
                                 round(rbind(em_ref_nmr, pred_ref_nmr), 2),
                                 rep(c("Estimation", "Prediction"),
                                     each = length(drug_names)))
    colnames(prepare_em_nma) <- c("order",
                                  "comparison",
                                  "median", "lower", "upper",
                                  "interval")
    colnames(prepare_em_nmr) <- colnames(prepare_em_nma)
  } else if (is.element(measure, c("OR", "RR", "ROM")) & model == "RE") {
    prepare_em_nma <- data.frame(as.factor(rep(rev(seq_len(len_drug)), 2)),
                                 rep(drug_names_sorted, 2),
                                 round(rbind(exp(em_ref_nma),
                                             exp(pred_ref_nma)), 2),
                                 rep(c("Estimation", "Prediction"),
                                     each = length(drug_names)))
    prepare_em_nmr <- data.frame(as.factor(rep(rev(seq_len(len_drug)), 2)),
                                 rep(drug_names_sorted, 2),
                                 round(rbind(exp(em_ref_nmr),
                                             exp(pred_ref_nmr)), 2),
                                 rep(c("Estimation", "Prediction"),
                                     each = length(drug_names)))
    colnames(prepare_em_nma) <- c("order",
                                  "comparison",
                                  "median", "lower", "upper",
                                  "interval")
    colnames(prepare_em_nmr) <- colnames(prepare_em_nma)
  } else if (!is.element(measure, c("OR", "RR", "ROM")) & model == "FE") {
    prepare_em_nma <- data.frame(as.factor(rev(seq_len(len_drug))),
                                 drug_names_sorted,
                                 round(em_ref_nma, 2))
    prepare_em_nmr <- data.frame(as.factor(rev(seq_len(len_drug))),
                                 drug_names_sorted,
                                 round(em_ref_nmr, 2))
    colnames(prepare_em_nma) <- c("order",
                                  "comparison",
                                  "median", "lower", "upper")
    colnames(prepare_em_nmr) <- colnames(prepare_em_nma)
  } else if (is.element(measure, c("OR", "RR", "ROM")) & model == "FE") {
    prepare_em_nma <- data.frame(as.factor(rev(seq_len(len_drug))),
                                 drug_names_sorted,
                                 round(exp(em_ref_nma), 2))
    prepare_em_nmr <- data.frame(as.factor(rev(seq_len(len_drug))),
                                 drug_names_sorted,
                                 round(exp(em_ref_nmr), 2))
    colnames(prepare_em_nma) <- c("order",
                                  "comparison",
                                  "median", "lower", "upper")
    colnames(prepare_em_nmr) <- colnames(prepare_em_nma)
  }
  prepare_em <- cbind(rbind(prepare_em_nma, prepare_em_nmr),
                      analysis = rep(c("Network meta-analysis",
                                       "Network meta-regression"),
                                     each = length(prepare_em_nma[, 1])))

  # Forest plots on credible/prediction intervals of comparisons with the
  # selected comparator
  subtitle <- paste("Results for", cov_value[[2]],
                   ifelse(length(unique(reg$covariate)) < 3, " ",
                          cov_value[[1]]))

  measure2 <- effect_measure_name(measure, lower = FALSE)
  caption <- if (full$D == 0 & is.element(measure, c("OR", "RR", "ROM"))) {
    paste(measure2, "< 1, favours the first arm.",
          measure2, "> 1, favours", compar)
  } else if (full$D == 1 & is.element(measure, c("OR", "RR", "ROM"))) {
    paste(measure2, "< 1, favours", compar,
          ".", measure2, "> 1, favours the first arm")
  } else if (full$D == 0 & !is.element(measure, c("OR", "RR", "ROM"))) {
    paste(measure2, "< 0, favours the first arm.",
          measure2, "> 0, favours", compar)
  } else if (full$D == 1 & !is.element(measure, c("OR","RR",  "ROM"))) {
    paste(measure2, "< 0, favours", compar,
          ".", measure2, "> 0, favours the first arm")
  }


  p1 <- if (model == "RE") {
    ggplot(data = prepare_em[prepare_em$interval == "Prediction", ],
           aes(x = order,
               y = median,
               ymin = lower,
               ymax = upper,
               group = analysis)) +
      geom_hline(yintercept = ifelse(!is.element(
        measure, c("OR", "RR", "ROM")), 0, 1),
        lty = 1,
        linewidth = 1,
        col = "grey60") +
      geom_linerange(aes(color = analysis),
                     linewidth = 2,
                     position = position_dodge(width = 0.5)) +
      geom_errorbar(data = prepare_em[prepare_em$interval == "Estimation", ],
                    aes(x = order,
                        y = median,
                        ymin = lower,
                        ymax = upper,
                        group = analysis),
                    size = 2,
                    position = position_dodge(width = 0.5),
                    width = 0.0) +
      geom_point(size = 1.5,
                 colour = "white",
                 stroke = 0.3,
                 position = position_dodge(width = 0.5)) +
      geom_text(data = prepare_em[prepare_em$interval == "Estimation", ],
                aes(x = order,
                    y = median,
                    group = analysis,
                    colour = analysis,
                    label = paste0(sprintf("%.2f", median), " ", "(",
                                   sprintf("%.2f", lower),
                                   ",",
                                   " ",
                                   sprintf("%.2f", upper),
                                   ")",
                                   " ",
                                   "[",
                                   prepare_em[
                                     (length(drug_names_sorted) + 1):
                                       (length(drug_names_sorted) * 2)
                                     & prepare_em$interval == "Prediction",
                                              4],
                                   ",",
                                   " ",
                                   prepare_em[
                                     (length(drug_names_sorted) + 1):
                                       (length(drug_names_sorted) * 2)
                                     & prepare_em$interval == "Prediction",
                                              5],
                                   "]"),
                    hjust = 0,
                    vjust = -0.5),
                #colour = "black",
                size = 3.8,
                check_overlap = FALSE,
                parse = FALSE,
                position = position_dodge(width = 0.5),
                inherit.aes = TRUE,
                na.rm = TRUE) +
      labs(x = "",
           y = measure2,
           colour = "Analysis",
           subtitle = subtitle,
           caption = caption) +
      scale_x_discrete(breaks = as.factor(seq_len(len_drug)),
                       labels = drug_names_sorted[rev(seq_len(len_drug))]) +
      scale_colour_manual(name = "Analysis",
                          breaks = c("Network meta-analysis",
                                     "Network meta-regression"),
                         values = c("#009E73", "#D55E00")) +
      geom_label(aes(x = unique(order[is.na(median)]),
                     y = ifelse(!is.element(
                       measure, c("OR", "RR", "ROM")), 0, 1), # -0.2, 0.65
                     hjust = 0,
                     vjust = 1,
                     label = "Comparator intervention"),
                 fill = "beige",
                 colour = "black",
                 fontface = "plain",
                 size = 4) +
      scale_y_continuous(trans = ifelse(!is.element(
        measure, c("OR", "RR", "ROM")), "identity", "log10")) +
      guides(colour = guide_legend(nrow = 1)) +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold",
                                        size = 12),
            legend.position = "bottom",
            legend.justification = c(0.43, 0),
            legend.text =  element_text(color = "black", size = 12),
            legend.title = element_text(color = "black", face = "bold",
                                        size = 12),
            plot.caption = element_text(hjust = 0.01))
  } else {
    ggplot(data = prepare_em,
           aes(x = order,
               y = median,
               ymin = lower,
               ymax = upper,
               group = analysis,
               colour = analysis)) +
      geom_hline(yintercept = ifelse(!is.element(
        measure, c("OR", "RR", "ROM")), 0, 1),
        lty = 1,
        linewidth = 1,
        col = "grey60") +
      geom_linerange(linewidth = 2,
                     position = position_dodge(width = 0.5)) +
      geom_point(size = 1.5,
                 colour = "white",
                 stroke = 0.3,
                 position = position_dodge(width = 0.5)) +
      geom_text(aes(x = order,
                    y = median,
                    label = paste0(sprintf("%.2f", median), " ", "(",
                                   sprintf("%.2f", lower),
                                   ",",
                                   " ",
                                   sprintf("%.2f", upper),
                                   ")"),
                    hjust = 0,
                    vjust = -0.5),
                color = "black",
                size = 3.8,
                check_overlap = TRUE,
                parse = FALSE,
                position = position_dodge(width = 0.5),
                inherit.aes = TRUE,
                na.rm = TRUE) +
      labs(x = "",
           y = measure2,
           colour = "Analysis",
           subtitle = subtitle,
           caption = caption) +
      scale_x_discrete(breaks = as.factor(seq_len(len_drug)),
                       labels = drug_names_sorted[rev(seq_len(len_drug))]) +
      scale_color_manual(breaks = c("Network meta-analysis",
                                    "Network meta-regression"),
                         values = c("#009E73", "#D55E00")) +
      geom_label(aes(x = unique(order[is.na(median)]),
                     y = ifelse(!is.element(
                       measure, c("OR", "RR", "ROM")), 0, 1), #-0.2, 0.65
                     hjust = 0,
                     vjust = 1,
                     label = "Comparator intervention"),
                 fill = "beige",
                 colour = "black",
                 fontface = "plain",
                 size = 4) +
      scale_y_continuous(trans = ifelse(!is.element(
        measure, c("OR", "RR", "ROM")), "identity", "log10")) +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12),
            strip.text = element_text(color = "black", face = "bold",
                                      size = 12),
            axis.title.x = element_text(color = "black", face = "bold",
                                        size = 12),
            legend.position = "bottom",
            legend.justification = c(0.23, 0),
            legend.text =  element_text(color = "black", size = 12),
            legend.title = element_text(color = "black", face = "bold",
                                        size = 12),
            plot.caption = element_text(hjust = 0.01))
  }

  ## Get SUCRA for the selected covariate value
  # Include the 'compar' in the dataset
  basic_par <- subset(em_ref00_nmr, em_ref00_nmr[6] == drug_names[full$ref])
  indexing <- match(c(drug_names[full$ref], basic_par[, 5]), drug_names)
  em_subset_nmr_new <- matrix(NA, nrow = len_drug, ncol = 2)
  em_subset_nmr_new[indexing, 1] <- c(NA, basic_par[, 1])
  em_subset_nmr_new[indexing, 2] <- c(NA, basic_par[, 2])

  # Data for BUGS
  data_jag <- list("d_cov_mean" = em_subset_nmr_new[, 1],
                   "d_cov_prec" = 1/(em_subset_nmr_new[, 2]^2),
                   "nt" = len_drug,
                   "D" = full$D,
                   "ref" = full$ref)

  # Parameters to monitor
  param_jags <- c("d.new", "SUCRA")

  # Run the BUGS code for SUCRA
  suppressWarnings({
  jagsfit <- jags(data = data_jag,
                  parameters.to.save = param_jags,
                  DIC = FALSE,
                  model.file = textConnection('
                    model {
                       d.new[ref] <- 0
                       for (t in 1:(ref - 1)) {
                         d.new[t] ~ dnorm(d_cov_mean[t], d_cov_prec[t])
                       }
                       for (t in (ref + 1):nt) {
                         d.new[t] ~ dnorm(d_cov_mean[t], d_cov_prec[t])
                       }
                       sorted <- rank(d.new[])
                       for (t in 1:nt) {
                         order[t] <- (nt + 1 - sorted[t])*equals(D, 1) + sorted[t]*(1 - equals(D, 1))
                         most.effective[t] <- equals(order[t], 1)
                         for (l in 1:nt) {
                           effectiveness[t, l] <- equals(order[t], l)
                           cumeffectiveness[t, l] <- sum(effectiveness[t, 1:l])
                         }
                        SUCRA[t] <- sum(cumeffectiveness[t, 1:(nt - 1)])/(nt - 1)
                       }
                    }
                                              '),
                  n.chains = full$n_chains,
                  n.iter = full$n_iter,
                  n.burnin = full$n_burnin,
                  n.thin = full$n_thin)
  })

  # Turn R2jags object into a data-frame
  get_results <- as.data.frame(t(jagsfit$BUGSoutput$summary))

  # Obtain posterior distribution from parameters of interest
  sucra_res <-
    t(get_results)[startsWith(rownames(t(get_results)), "SUCRA["), c(1, 3, 7)]

  # Order by SUCRA of NMA model
  sucra_nmr <- sucra_res[order(sucra_full[, 1], decreasing = TRUE), ]

  # SUCRA of NMA model ordered
  sucra_nma <- sucra_full[order(sucra_full[, 1], decreasing = TRUE), c(1, 3, 7)]
  colnames(sucra_nmr) <- colnames(sucra_nma) <- c("mean", "lower", "upper")

  # Prepare dataset for SUCRA forest plot
  prepare_sucra <- data.frame(as.factor(rev(seq_len(len_drug))),
                              rep(drug_names_sorted, 2),
                              rbind(sucra_nma, sucra_nmr),
                              rep(c("Network meta-analysis",
                                    "Network meta-regression"),
                                  each = len_drug))
  colnames(prepare_sucra) <- c("order",
                               "intervention",
                               "mean", "lower", "upper",
                               "analysis")

  # Forest plots of SUCRA per intervention and analysis
  p2 <- ggplot(data = prepare_sucra,
               aes(x = order,
                   y = mean,
                   ymin = lower,
                   ymax = upper,
                   group = analysis)) +
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 0.5,
                  fill = "lowest"),
              alpha = 0.01) +
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0.5, ymax = 0.8,
                  fill = "intermediate"),
              alpha = 0.01) +
    geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0.8, ymax = 1.0,
                  fill = "highest"),
              alpha = 0.01) +
    geom_linerange(aes(colour = analysis),
                   linewidth = 2,
                   position = position_dodge(width = 0.5)) +
    geom_point(size = 1.5,
               colour = "white",
               stroke = 0.3,
               position = position_dodge(width = 0.5)) +
    geom_text(aes(x = order, # as.factor(order)
                  y = mean,
                  label = paste0(round(mean * 100, 0),
                                 " ",
                                 "(",
                                 round(lower * 100, 0),
                                 ",",
                                 " ",
                                 round(upper * 100, 0), ")"),
                  hjust = 0,
                  vjust = -0.5),
              color = "black",
              size = 3.8,
              check_overlap = FALSE,
              parse = FALSE,
              position = position_dodge(width = 0.5),
              inherit.aes = TRUE) +
    scale_color_manual(breaks = c("Network meta-analysis",
                                  "Network meta-regression"),
                       values = c("#009E73", "#D55E00"),
                       guide = "none") +
    scale_fill_manual(name = "Ranked",
                      values = c("lowest" = "#D55E00",
                                 "intermediate" = "orange",
                                 "highest" = "#009E73")) +
    labs(x = "",
         y = "Surface under the cumulative ranking curve value",
         caption = " ") +
    scale_x_discrete(breaks = as.factor(seq_len(len_drug)),
                     labels = prepare_sucra$intervention[rev(
                       seq_len(len_drug))]) +
    scale_y_continuous(labels = percent) +
    coord_flip() +
    theme_classic() +
    theme(axis.text.x = element_text(color = "black", size = 12),
          axis.text.y = element_text(color = "black", size = 12),
          axis.title.x = element_text(color = "black", face = "bold",
                                      size = 12),
          legend.position = "bottom",
          legend.box="vertical",
          legend.text =  element_text(color = "black", size = 12),
          legend.title =  element_text(color = "black", face = "bold",
                                       size = 12))

  # Bring together both forest-plots
  forest_plots <- suppressWarnings(
    ggarrange(p1, p2 + guides(fill = guide_legend(override.aes =
                                                    list(alpha = 0.4))),
              nrow = 1, ncol = 2, labels = c("A)", "B)"),
              common.legend = FALSE, legend = "bottom")
  )

 return(forest_plots)
}
