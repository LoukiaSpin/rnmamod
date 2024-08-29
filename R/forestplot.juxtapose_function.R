#' Forest plot of juxtaposing several network meta-analysis models
#'
#' @description
#'   Provides a forest plot with the posterior median and 95\% credible
#'   and prediction intervals for comparisons with the selected intervention
#'   (comparator) in the network under several network meta-analyses models, as
#'   well as a forest plot with the corresponding SUCRA values.
#'
#' @param results A list of at least two objects of S3 class
#'   \code{\link{run_model}} or \code{\link{run_metareg}}. See 'Value' in
#'   \code{\link{run_model}} and \code{\link{run_metareg}}.
#' @param compar A character to indicate the comparator intervention. It must
#'   be any name found in \code{drug_names}.
#' @param name A vector of characters referring to the juxtaposed models. If the
#'   argument is left unspecified, the names of models appear as 'Model X' with
#'   'X' being the order/position of each model in the argument \code{results}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#' @param axis_title_size A positive integer for the font size of x axis title.
#'   \code{axis_title_size} determines the axis.title argument
#'   found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param axis_text_size A positive integer for the font size of axis text (both
#'   axes). \code{axis_text_size} determines the axis.text argument found in the
#'   theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param caption_text_size A positive integer for the font size of caption
#'   text. \code{caption_text_size} determines the plot.caption argument found
#'   in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param label_size A positive integer for the font size of labels appearing on
#'   each interval. \code{label_size} determines the size argument found in the
#'   geom's aesthetic properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param position_width A positive integer specifying the vertical position of
#'   the intervals. \code{position_width} is found in the geom's aesthetic
#'   properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'
#' @return A list of the following two figures:
#'   \item{forest_plots}{A panel of two forest plots: (1) a forest plot on the
#'   posterior median and 95\% credible and prediction intervals for comparisons
#'   with the selected comparator treatment (specified with \code{compar}), and
#'   (2) a forest plot on the posterior mean and 95\% credible interval of SUCRA
#'   values of the treatments (Salanti et al., 2011).}
#'   \item{tau_plot}{A forest plot on the posterior median and 95\% credible
#'   interval of the between-study standard deviation.}
#'
#' @details The y-axis of the forest plot on \bold{forest_plots} displays the
#'   labels of the treatments in the network; the selected treatment that
#'   comprises the \code{compar} argument is annotated in the plot with the
#'   label 'Comparator intervention'.
#'   For each comparison with the selected treatment, the 95\% credible and
#'   prediction intervals are displayed as overlapping lines. Black lines refer
#'   to estimation under both analyses. Coloured lines refer to prediction
#'   under each model, respectively. The corresponding numerical results are
#'   displayed above each line: 95\% credible intervals are found in
#'   parentheses, and 95\% predictive intervals are found in brackets.
#'   Odds ratios, relative risks, and ratio of means are reported in the
#'   original scale after exponentiation of the logarithmic scale.
#'
#'   If one of the models refer to network meta-regression
#'   (\code{\link{run_metareg}}) the results on treatment effects (estimation
#'   and prediction) and SUCRA values refer to the covariate value selected
#'   when employing \code{\link{run_metareg}}.
#'
#'   The y-axis for the forest plot on \bold{SUCRA} values displays the
#'   labels of the treatments in the network. The corresponding numerical
#'   results are displayed above each line.
#'
#'   In \bold{forest_plots} and \bold{tau_plot}, the treatments are sorted in
#'   the descending order of their SUCRA values based on the first model
#'   specified in \code{results}.
#'
#'   \bold{Important note:} \code{forestplot_juxtapose} should be used to
#'   compare the results from several network meta-analysis models that contain
#'   the same treatments, have the same meta-analysis model (fixed-effect or
#'   random-effects) and the same effect measure; otherwise the execution of the
#'   function will be stopped and an error message will be printed on the R
#'   console.
#'
#'   \code{forestplot_juxtapose} is used only for a network of treatments.
#'   In the case of two treatments, the execution of the function will be
#'   stopped and an error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_metareg}}, \code{\link{run_model}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries
#' for presenting results from multiple-treatment meta-analysis: an overview and
#' tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71.
#' \doi{10.1016/j.jclinepi.2010.03.016}
#'
#' @export
forestplot_juxtapose <- function(results,
                                 compar,
                                 name,
                                 drug_names,
                                 axis_title_size = 12,
                                 axis_text_size = 12,
                                 caption_text_size = 9,
                                 label_size = 3.5,
                                 position_width = 0.8) {


  ## Check default
  type <- unlist(lapply(results, function(x) class(x)))
  if (any(!is.element(type, c("run_model", "run_metareg"))) == TRUE) {
    stop("A list of objects of S3 class 'run_model' or 'run_metareg'",
         call. = FALSE)
  } else if (length(unique(lapply(results,
                                  function(x) length(x$SUCRA[, 1])))) != 1) {
    stop("All elements of the list 'results' must contain the same treatments",
         call. = FALSE)
  } else if (length(unique(lapply(results,
                                  function(x) length(x$model)))) != 1) {
    stop("The argument must refer to the same meta-analysis model (RE or FE)",
         call. = FALSE)
  } else if (length(unique(lapply(results,
                                  function(x) length(x$measure)))) != 1) {
    stop("The argument must refer to the same effect measure",
         call. = FALSE)
  }
  name <- if (missing(name)) {
    paste("Model", 1:length(results))
  } else if (length(name) != length(results)) {
    stop("The argument must have the same length with 'results'",
         call. = FALSE)
  } else {
    name
  }
  drug_names <- if (missing(drug_names)) {
    aa <- "The argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as specified in 'data' is used, instead."
    message(paste(aa, bb))
    nt <- length(results[[1]]$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug_names
  }
  len_drug <- length(drug_names)
  compar <- if (missing(compar)) {
    stop("The argument 'compar' has not been defined.", call. = FALSE)
  } else if (!is.element(compar, drug_names)) {
    stop("The value of the argument 'compar' is not found in the 'drug_names'.",
         call. = FALSE)
  } else if (is.element(compar, drug_names)) {
    compar
  }


  ## The function is suitable not for meta-analysis
  if (length(drug_names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis.",
         call. = FALSE)
  }


  ## Sort the drugs by their SUCRA in decreasing order of the first model
  drug_names_sorted <-
    drug_names[order(results[[1]]$SUCRA[, 1], decreasing = TRUE)]


  ## A matrix with all possible comparisons in the network
  poss_pair_comp1 <- data.frame(exp = t(combn(drug_names, 2))[, 2],
                                comp = t(combn(drug_names, 2))[, 1])
  poss_pair_comp2 <- data.frame(exp = t(combn(drug_names, 2))[, 1],
                                comp = t(combn(drug_names, 2))[, 2])
  poss_pair_comp <- rbind(poss_pair_comp1, poss_pair_comp2)


  ## Prepare dataset with comparisons with the selected comparator
  # Effect size of all possible pairwise comparisons
  em_ref00_nma <-
    lapply(results,
           function(x) cbind(rbind(data.frame(median = x$EM[, 5],
                                              lower = x$EM[, 3],
                                              upper = x$EM[, 7]),
                                   data.frame(median = x$EM[, 5] * (-1),
                                              lower = x$EM[, 7] * (-1),
                                              upper = x$EM[, 3] * (-1))),
                             poss_pair_comp))

  # Restrict to comparisons with the selected comparator
  em_subset_nma <- lapply(em_ref00_nma, function(x) subset(x, x[5] == compar))

  # Replace with NA the results on the comparator versus comparator
  em_ref0_nma <- lapply(em_subset_nma,
                        function(x) rbind(x[, c(1:3)], c(rep(NA, 3))))

  # Move the SUCRA of the comparator at the end
  sucra_new <-
    lapply(1:length(results),
           function(x)
             data.frame(results[[x]]$SUCRA[, 1], drug_names)
           [order(match(data.frame(results[[x]]$SUCRA[, 1], drug_names)[, 2],
                        em_subset_nma[[x]][, 4])), 1])

  # Sort the summary effects by the SUCRA of the irst model in decreasing order
  em_ref_nma <-
    lapply(1:length(results),
           function(x) em_ref0_nma[[x]][order(sucra_new[[1]],
                                              decreasing = TRUE), ])


  ## Posterior results on the predicted estimates of comparisons with the
  # selected comparator as reference
  if (results[[1]]$model == "RE") {

    # Predicted effect size of all possible pairwise comparisons (NMA)
    pred_ref00_nma <-
      lapply(results,
             function(x) cbind(rbind(data.frame(median = x$EM_pred[, 5],
                                                lower = x$EM_pred[, 3],
                                                upper = x$EM_pred[, 7]),
                                     data.frame(median = x$EM_pred[, 5] * (-1),
                                                lower = x$EM_pred[, 7] * (-1),
                                                upper = x$EM_pred[, 3] * (-1))),
                               poss_pair_comp))
    # Restrict to comparisons with the selected comparator
    pred_subset_nma <-
      lapply(pred_ref00_nma, function(x) subset(x, x[5] == compar))

    # Replace with NA the results on the comparator versus comparator
    pred_ref0_nma <- lapply(pred_subset_nma,
                            function(x) rbind(x[, c(1:3)], c(rep(NA, 3))))

    # Sort by SUCRA in decreasing order and remove the reference intervention
    pred_ref_nma <-
      lapply(1:length(results),
             function(x) pred_ref0_nma[[x]][order(sucra_new[[1]],
                                                  decreasing = TRUE), ])
    rownames(pred_ref_nma) <- NULL
  }


  ## Create a data-frame with credible and prediction intervals of comparisons
  # with the reference intervention
  if (!is.element(results[[1]]$measure, c("OR", "RR", "ROM")) &
      results[[1]]$model == "RE") {
    prepare_em_nma <-
      lapply(1:length(results),
             function(x)
               data.frame(as.factor(rep(rev(seq_len(len_drug)), 2)),
                          rep(drug_names_sorted, 2),
                          round(rbind(em_ref_nma[[x]], pred_ref_nma[[x]]), 2),
                          rep(c("Estimation", "Prediction"),
                              each = length(drug_names))) )
  } else if (is.element(results[[1]]$measure, c("OR", "RR", "ROM")) &
             results[[1]]$model == "RE") {
    prepare_em_nma <-
      lapply(1:length(results),
             function(x) data.frame(as.factor(rep(rev(seq_len(len_drug)), 2)),
                                    rep(drug_names_sorted, 2),
                                    round(rbind(exp(em_ref_nma[[x]]),
                                                exp(pred_ref_nma[[x]])), 2),
                                    rep(c("Estimation", "Prediction"),
                                        each = length(drug_names))))
  } else if (!is.element(results[[1]]$measure, c("OR", "RR", "ROM")) &
             results[[1]]$model == "FE") {
    prepare_em_nma <-
      lapply(1:length(results),
             function(x) data.frame(as.factor(rev(seq_len(len_drug))),
                                    drug_names_sorted,
                                    round(em_ref_nma[[x]], 2)))
  } else if (is.element(results[[1]]$measure, c("OR", "RR", "ROM")) &
             results[[1]]$model == "FE") {
    prepare_em_nma <-
      lapply(1:length(results),
             function(x) data.frame(as.factor(rev(seq_len(len_drug))),
                                    drug_names_sorted,
                                    round(exp(em_ref_nma[[x]]), 2)))
  }


  ## Bring all model results into one data-frame
  # Bind the 'prepare_em_nma' data-frames by row
  prepare_em <- do.call("rbind", prepare_em_nma)
  colnames(prepare_em) <-
    if (results[[1]]$model == "RE") {
      c("order",  "comparison", "median", "lower", "upper", "interval")
    } else {
      c("order",  "comparison", "median", "lower", "upper")
    }

  # Add the model name indicator
  prepare_em$analysis <- rep(name, each = length(prepare_em_nma[[1]][, 1]))


  ## Define the
  measure2 <- effect_measure_name(results[[1]]$measure, lower = FALSE)
  caption <- if (results[[1]]$D == 0 &
                 is.element(results[[1]]$measure, c("OR", "RR", "ROM"))) {
    paste0(measure2, " < 1, favours the first arm. ",
           measure2, " > 1, favours ", compar, ".")
  } else if (results[[1]]$D == 1 &
             is.element(results[[1]]$measure, c("OR", "RR", "ROM"))) {
    paste0(measure2, " < 1, favours ", compar, ". ",
           measure2, " > 1, favours the first arm.")
  } else if (results[[1]]$D == 0 &
             !is.element(results[[1]]$measure, c("OR", "RR", "ROM"))) {
    paste0(measure2, " < 0, favours the first arm. ",
           measure2, " > 0, favours ", compar, ".")
  } else if (results[[1]]$D == 1 &
             !is.element(results[[1]]$measure, c("OR","RR",  "ROM"))) {
    paste0(measure2, " < 0, favours ", compar, ". ",
           measure2, " > 0, favours the first arm.")
  }

  p1 <-
    if (results[[1]]$model == "RE") {
      ggplot(data = prepare_em[prepare_em$interval == "Prediction", ],
             aes(x = order,
                 y = median,
                 ymin = lower,
                 ymax = upper,
                 group = factor(analysis, levels = rev(name)))) +
        geom_hline(yintercept = ifelse(!is.element(
          results[[1]]$measure, c("OR", "RR", "ROM")), 0, 1),
          lty = 1,
          linewidth = 1,
          col = "grey60") +
        geom_linerange(aes(color = factor(analysis, levels = rev(name))),
                       linewidth = 2,
                       position = position_dodge(width = position_width)) +
        geom_errorbar(data = prepare_em[prepare_em$interval == "Estimation", ],
                      aes(x = order,
                          y = median,
                          ymin = lower,
                          ymax = upper,
                          group = factor(analysis, levels = rev(name))),
                      linewidth = 2,
                      position = position_dodge(width = position_width),
                      width = 0.0) +
        geom_point(size = 1.5,
                   colour = "white",
                   stroke = 0.3,
                   position = position_dodge(width = position_width)) +
        geom_text(data = prepare_em[prepare_em$interval == "Estimation", ],
                  aes(x = order,
                      y = median,
                      group = factor(analysis, levels = rev(name)),
                      #colour = analysis,
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
                                         (length(drug_names_sorted) * 2) &
                                         prepare_em$interval == "Prediction",
                                       4],
                                     ",",
                                     " ",
                                     prepare_em[
                                       (length(drug_names_sorted) + 1):
                                         (length(drug_names_sorted) * 2) &
                                         prepare_em$interval == "Prediction",
                                       5],
                                     "]"),
                      hjust = 0,
                      vjust = -0.5),
                  colour = "black",
                  size = label_size,
                  check_overlap = FALSE,
                  parse = FALSE,
                  position = position_dodge(width = position_width),
                  inherit.aes = TRUE,
                  na.rm = TRUE) +
        labs(x = "",
             y = measure2,
             colour = "Analysis",
             caption = caption) +
        scale_x_discrete(breaks = as.factor(seq_len(len_drug)),
                         labels = drug_names_sorted[rev(seq_len(len_drug))]) +
        geom_label(aes(x = unique(order[is.na(median)]),
                       y = ifelse(!is.element(
                         results[[1]]$measure, c("OR", "RR", "ROM")), 0, 1), # -0.2, 0.65
                       hjust = 0,
                       vjust = 1,
                       label = "Comparator intervention"),
                   fill = "beige",
                   colour = "black",
                   fontface = "plain",
                   size = label_size) +
        scale_y_continuous(trans = ifelse(!is.element(
          results[[1]]$measure, c("OR", "RR", "ROM")), "identity", "log10")) +
        scale_colour_manual(breaks = as.factor(name),
                            values = hue_pal()(length(name))) +
        guides(colour = guide_legend(nrow = 1)) +
        coord_flip() +
        theme_classic() +
        theme(axis.text = element_text(color = "black",
                                       size = axis_text_size),
              axis.title = element_text(color = "black",
                                        face = "bold",
                                        size = axis_title_size),
              legend.position = "bottom",
              legend.text =  element_text(color = "black",
                                          size = axis_text_size),
              legend.title = element_text(color = "black",
                                          face = "bold",
                                          size = axis_title_size),
              plot.caption = element_text(size = caption_text_size,
                                          hjust = 0.01))
    } else {
      ggplot(data = prepare_em,
             aes(x = order,
                 y = median,
                 ymin = lower,
                 ymax = upper,
                 group = factor(analysis, levels = rev(name)),
                 colour = factor(analysis, levels = rev(name)))) +
        geom_hline(yintercept = ifelse(!is.element(
          results[[1]]$measure, c("OR", "RR", "ROM")), 0, 1),
          lty = 1,
          linewidth = 1,
          col = "grey60") +
        geom_linerange(linewidth = 2,
                       position = position_dodge(width = position_width)) +
        geom_point(size = 1.5,
                   colour = "white",
                   stroke = 0.3,
                   position = position_dodge(width = position_width)) +
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
                  size = label_size,
                  check_overlap = TRUE,
                  parse = FALSE,
                  position = position_dodge(width = position_width),
                  inherit.aes = TRUE,
                  na.rm = TRUE) +
        labs(x = "",
             y = measure2,
             colour = "Analysis",
             caption = caption) +
        scale_x_discrete(breaks = as.factor(seq_len(len_drug)),
                         labels = drug_names_sorted[rev(seq_len(len_drug))]) +
        geom_label(aes(x = unique(order[is.na(median)]),
                       y = ifelse(!is.element(
                         results[[1]]$measure, c("OR", "RR", "ROM")), 0, 1), #-0.2, 0.65
                       hjust = 0,
                       vjust = 1,
                       label = "Comparator intervention"),
                   fill = "beige",
                   colour = "black",
                   fontface = "plain",
                   size = label_size) +
        scale_y_continuous(trans = ifelse(!is.element(
          results[[1]]$measure, c("OR", "RR", "ROM")), "identity", "log10")) +
        scale_colour_manual(breaks = as.factor(name),
                            values = hue_pal()(length(name))) +
        coord_flip() +
        theme_classic() +
        theme(axis.text = element_text(color = "black",
                                       size = axis_text_size),
              axis.title = element_text(color = "black",
                                        face = "bold",
                                        size = axis_title_size),
              legend.position = "bottom",
              legend.text =  element_text(color = "black",
                                          size = axis_text_size),
              legend.title = element_text(color = "black",
                                          face = "bold",
                                          size = axis_title_size),
              plot.caption = element_text(size = caption_text_size,
                                          hjust = 0.01))
    }


  # SUCRA of NMA model ordered
  sucra_ordered <-
    lapply(results,
           function(x) x$SUCRA[order(results[[1]]$SUCRA[, 1],
                                     decreasing = TRUE), c(1, 3, 7)])

  # Bring all in one dataset (bind by row)
  sucra_ordered_all <- do.call("rbind", sucra_ordered)
  colnames(sucra_ordered_all) <- c("mean", "lower", "upper")

  # Prepare dataset for SUCRA forest plot
  prepare_sucra <- data.frame(as.factor(rev(seq_len(len_drug))),
                              rep(drug_names_sorted, length(results)),
                              sucra_ordered_all,
                              rep(name, each = length(drug_names)))
  colnames(prepare_sucra) <- c("order",
                               "intervention",
                               "mean", "lower", "upper",
                               "analysis")

  # Forest plot of SUCRA per intervention and analysis
  p2 <-
    ggplot(data = prepare_sucra,
           aes(x = order,
               y = mean,
               ymin = lower,
               ymax = upper,
               group = factor(analysis, levels = rev(name)))) +
    geom_linerange(aes(colour = factor(analysis, levels = rev(name))),
                   linewidth = 2,
                   position = position_dodge(width = position_width)) +
    geom_point(size = 1.5,
               colour = "white",
               stroke = 0.3,
               position = position_dodge(width = position_width)) +
    geom_text(aes(x = order, # as.factor(order)
                  y = mean,
                  label = paste0(round(mean * 100, 0),
                                 " ",
                                 "(",
                                 round(lower * 100, 0),
                                 ",",
                                 " ",
                                 round(upper * 100, 0), ")"),
                  hjust = ifelse(mean < 0.80, 0, 1),
                  vjust = -0.5),
              color = "black",
              size = label_size,
              check_overlap = FALSE,
              parse = FALSE,
              position = position_dodge(width = position_width),
              inherit.aes = TRUE) +
    scale_colour_manual(breaks = as.factor(name),
                        values = hue_pal()(length(name))) +
    labs(x = "",
         y = "Surface under the cumulative ranking curve value",
         colour = "Analysis",
         caption = " ") +
    scale_x_discrete(breaks = as.factor(seq_len(len_drug)),
                     labels = prepare_sucra$intervention[rev(
                       seq_len(len_drug))]) +
    scale_y_continuous(labels = percent) +
    coord_flip() +
    theme_classic() +
    theme(axis.text = element_text(color = "black",
                                   size = axis_text_size),
          axis.title = element_text(color = "black",
                                    face = "bold",
                                    size = axis_title_size),
          legend.position = "bottom",
          legend.text =  element_text(color = "black",
                                      size = axis_text_size),
          legend.title =  element_text(color = "black",
                                       face = "bold",
                                       size = axis_title_size))


  ## Results on between-study standard deviation (tau)
  # Prepare dataset
  tau_res <- data.frame(do.call("rbind",
                                lapply(results, function(x) x$tau[c(5, 3, 7)])),
                        order = length(name):1,
                        analysis = name)
  colnames(tau_res)[1:3] <- c("median", "lower", "upper")

  # Get forestplot on tau
  tau_plot <-
    ggplot(data = tau_res,
           aes(x = order,
               y = median,
               ymin = lower,
               ymax = upper,
               group = factor(analysis, levels = rev(name)),
               colour = factor(analysis, levels = rev(name)))) +
    geom_linerange(linewidth = 2,
                   position = position_dodge(width = position_width)) +
    geom_point(size = 1.5,
               colour = "white",
               stroke = 0.3,
               position = position_dodge(width = position_width)) +
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
              size = label_size + 1,
              check_overlap = TRUE,
              parse = FALSE,
              position = position_dodge(width = position_width),
              inherit.aes = TRUE,
              na.rm = TRUE) +
    labs(x = "",
         y = "Between-study standard deviation",
         colour = "Analysis") +
    scale_x_continuous(breaks = seq_len(length(name)),
                      labels = rev(name)) +
    scale_colour_manual(breaks = as.factor(name),
                        values = hue_pal()(length(name))) +
    guides(colour = guide_legend(nrow = 1, byrow = TRUE)) +
    coord_flip() +
    theme_classic() +
    theme(axis.text = element_text(color = "black",
                                     size = axis_text_size),
          axis.title = element_text(color = "black",
                                    face = "bold",
                                      size = axis_title_size),
          legend.position = "none")


  # Bring together both forest-plots
  forest_plots <- suppressWarnings(
    ggarrange(p1, p2,
              nrow = 1, ncol = 2, labels = c("A)", "B)"),
              common.legend = TRUE, legend = "bottom"))


  ## Collect results
  collect_results <- list(forest_plots = forest_plots,
                          tau_plot = tau_plot)

  return(collect_results)
}
