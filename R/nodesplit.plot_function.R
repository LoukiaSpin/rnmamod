#' End-user-ready results for the node-splitting approach
#'
#' @description \code{nodesplit_plot} hosts a toolkit of functions that
#'   facilitates the comparison of the consistency model
#'   (via \code{\link{run_model}}) with the node-splitting approach
#'   (via \code{\link{run_nodesplit}}) regarding the posterior summaries of the
#'   direct and indirect effects and inconsistency factor of the split
#'   nodes, the between-trial standard deviation and model assessment
#'   parameters (Spiegelhalter et al., 2002) after each split node in the
#'   network.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param node An object of S3 class \code{\link{run_nodesplit}}. See 'Value'
#'   in \code{\link{run_nodesplit}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#' @param save_xls Logical to indicate whether to export the tabulated results
#'   to an 'xlsx' file (via the \code{\link[writexl:write_xlsx]{write_xlsx}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}) at the working
#'   directory of the user. The default is \code{FALSE} (do not export).
#'
#' @return \code{nodesplit_plot} returns the following list of elements:
#'   \item{table_effect_size}{A data-frame with the posterior median,
#'   posterior standard deviation and 95\% credible interval of the direct and
#'   indirect effect and the inconsistency factor of each split node.}
#'   \item{table_model_assessment}{A data-frame with the model assessment
#'   parameters (DIC, posterior mean of total residual deviance, and number of
#'   effective parameters), the posterior median, posterior standard deviation
#'   and 95\% credible interval of \emph{tau} under the consistency model and
#'   after each split node. See 'Details'.}
#'   \item{intervalplot_inconsistency_factor}{A panel of interval plots on
#'   the direct and indirect effect of the split nodes and the corresponding
#'   inconsistency factor. See 'Details'.}
#'   \item{intervalplot_tau}{An interval plot on \emph{tau} after each
#'   split node. See 'Details'.}
#'
#' @details \code{intervalplot_inconsistency_factor} includes as many interval
#'   plots as the number of split nodes in the network. Each interval plot
#'   illustrates the posterior median and 95\% credible interval of the direct and
#'   indirect effect of the split nodes and the corresponding inconsistency
#'   factor.
#'   The line that corresponds to the inconsistency factor is highlighted with
#'   green, when it does not cross the vertical line of no difference (between
#'   the direct and indirect effect), and red otherwise. If there are more than
#'   30 split nodes, the function presents the interval plots on split nodes
#'   with conclusive inconsistency factor (green intervals) or those with
#'   an opposite sign in the direct and indirect effects.
#'
#'   \code{intervalplot_tau} is an interval plot on the median and 95\% credible
#'   interval of \emph{tau} after each split node. The lines that correspond to
#'   the split nodes are sorted in ascending order of the deviance information
#'   criterion (DIC) which appears at the top of each line.
#'   The estimated median and 95\% credible intervals of \emph{tau} under the
#'   consistency model appear in the interval plot as a solid and two dotted
#'   parallel blue lines, respectively. The different levels of heterogeneity
#'   appear as green, yellow, orange, and red rectangulars to indicate a low,
#'   reasonable, fairly high, and fairly extreme heterogeneity, respectively,
#'   following the classification of Spiegelhalter et al. (2004).
#'   When a fixed-effect model has been performed, \code{nodesplit_plot} does
#'   not return the \code{intervalplot_tau}.
#'
#'   \code{table_model_assessment} also includes the column
#'   \emph{DIC-based better fit} that indicates the preferred model in terms of
#'   parsimony for each split node. Therefore, the DIC of the model after each
#'   split node is compared with the DIC of the consistency model
#'   (Dias et al., 2010). If the difference in DIC exceeds 5, the consistency
#'   model is preferred; if the difference in DIC is less than -5, the model
#'   after the split node is preferred; otherwise, there is little to choose
#'   between the compared models.
#'
#'   For a binary outcome, when \code{measure} is "RR" (relative risk) or "RD"
#'   (risk difference) in \code{\link{run_model}}, \code{nodesplit_plot}
#'   currently presents the results in the odds ratio scale. This is because,
#'   the odds ratio is used as the 'best-case' effect measure in
#'   \code{\link{run_model}}. Then, relative risk, and risk difference are
#'   obtained as a function of the odds ratio and the selected baseline risk
#'   (See 'Details' in \code{\link{run_model}}).
#'
#'   The split nodes have been automatically selected via the
#'   \code{\link[gemtc:mtc.nodesplit.comparisons]{mtc.nodesplit.comparisons}}
#'   function of the R-package
#'   \href{https://CRAN.R-project.org/package=gemtc}{gemtc}.
#'   See 'Details' in \code{\link{run_nodesplit}}.
#'
#'   Furthermore, \code{nodesplit_plot} exports both data-frames to separate
#'   'xlsx' files (via the \code{\link[writexl:write_xlsx]{write_xlsx}} function
#'   of the R-package
#'   \href{https://CRAN.R-project.org/package=writexl}{writexl}) to the working
#'   directory of the user.
#'
#'   \code{nodesplit_plot} can be used only for a network of interventions and
#'   when there is at least one split node. Otherwise, the execution of the
#'   function will be stopped and an error message will be printed on the R
#'   console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso
#'   \code{\link[gemtc:mtc.nodesplit.comparisons]{mtc.nodesplit.comparisons}},
#'   \code{\link{run_model}}, \code{\link{run_nodesplit}},
#'   \code{\link[writexl:write_xlsx]{write_xlsx}}
#'
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
#' treatment comparison meta-analysis.
#' \emph{Stat Med} 2010;\bold{29}(7-8):932--44. doi: 10.1002/sim.3767
#'
#' Spiegelhalter DJ, Abrams KR, Myles JP. Bayesian approaches to clinical trials
#' and health-care evaluation. John Wiley and Sons, Chichester, 2004.
#'
#' Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A. Bayesian measures of
#' model complexity and fit. \emph{J R Stat Soc B} 2002;\bold{64}(4):583--616.
#' doi: 10.1111/1467-9868.00353
#'
#' @examples
#' data("nma.baker2009")
#'
#' # Read results from 'run_model' (using the default arguments)
#' res <- readRDS(system.file('extdata/res_baker.rds', package = 'rnmamod'))
#'
#' # Read results from 'run_nodesplit' (using the default arguments)
#' node <- readRDS(system.file('extdata/node_baker.rds', package = 'rnmamod'))
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "budesonide", "budesonide plus formoterol",
#'                   "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium")
#'
#' # Plot the results from both models
#' nodesplit_plot(full = res,
#'                node = node,
#'                drug_names = interv_names)
#'
#' @export
nodesplit_plot <- function(full, node, drug_names, save_xls) {

  if (!inherits(full, "run_model") || is.null(full)) {
    stop("'full' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  if (!inherits(node, "run_nodesplit") || is.null(node)) {
    stop("'node' must be an object of S3 class 'run_nodesplit'.",
         call. = FALSE)
  }

  if (is.null(node)) {
    stop("There is no split node.", call. = FALSE)
  }

  save_xls <- if (missing(save_xls)) {
    FALSE
  } else {
    save_xls
  }

  data <- full$data
  measure <- if (is.element(full$measure, c("RR", "RD"))) {
    "OR"
  } else {
    full$measure
  }
  item <- data_preparation(data, measure)
  if (item$nt < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis.",
         call. = FALSE)
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

  # Keep tau and model assessment measures from NMA model
  tau_values <- full$tau[c(5, 2, 3, 7)]
  model_assess_nma <- full$model_assessment
  data_points <- model_assess_nma[4]

  # Effect measure
  measure2 <- effect_measure_name(measure, lower = FALSE)

  # Analysis model
  model <- full$model

  # Direct effects (node split)
  direct0 <- node$direct

  # Inirect effects (node split)
  indirect0 <- node$indirect

  # Inconsistency factor
  incons_factor0 <- node$diff

  # Model assessment (node split)
  model_assess <- node$model_assessment

  # Sort by DIC in ascending order
  direct <- direct0[order(model_assess$DIC), ]
  indirect <- indirect0[order(model_assess$DIC), ]
  incons_factor <- incons_factor0[order(model_assess$DIC), ]
  if (model == "RE") {
    tau <- node$tau[order(model_assess$DIC), ]
  } else {
    tau <- NA
  }
  model_assess_sort <- model_assess[order(model_assess$DIC), ]

  # Interventions' name: Replace code with original names
  # For treat1 (non-baseline arm)
  for (i in sort(unique(unlist(direct[, 1])))) {
    direct[direct$treat1 == i, 1] <- drug_names[i]
  }
  # For treat2 (baseline arm)
  for (i in sort(unique(unlist(direct[, 2])))) {
    direct[direct$treat2 == i, 2] <- drug_names[i]
  }

  # Prepare the dataset to create the panel of interval plots
  # on the 'direct evidence', 'indirect evidence', and 'inconsistency factor'
  # for each split node
  comp <- paste(direct[, 1], "vs", direct[, 2])
  prepare <- data.frame(rep(comp, 3),
                        rbind(direct[, c(3, 5:6)], indirect[, c(3, 5:6)],
                              incons_factor[, c(3, 5:6)]),
                        rep(c("direct", "indirect", "IF"),
                            each = length(direct[, 1])))
  colnames(prepare) <- c("node", "median", "lower", "upper", "evidence")
  prepare$stat_sign <- ifelse(prepare$lower > 0 | prepare$upper < 0,
                              "strong evidence",
                              "weak evidence")
  prepare$stat_sign <- ifelse(prepare$evidence != "IF", NA, prepare$stat_sign)
  prepare$DIC <- sort(model_assess$DIC)

  # Create the panel of interval plots
  if (length(unique(comp)) <= 30) {
    p1 <- ggplot(data = prepare,
                 aes(x = factor(evidence,
                                levels = c("IF", "indirect", "direct")),
                     y = median,
                     ymin = lower,
                     ymax = upper,
                     colour = stat_sign)) +
            geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
            geom_hline(yintercept = 0, lty = 1, size = 1, col = "grey") +
            geom_point(size = 1.5,
                       colour = "white",
                       stroke = 0.3,
                       position = position_dodge(width = 0.5)) +
            geom_text(aes(x = as.factor(evidence),
                          y = round(median, 2),
                          label = paste0(sprintf("%.2f", median), " ", "(",
                                         sprintf("%.2f", lower), ",", " ",
                                         sprintf("%.2f", upper), ")"),
                          hjust = 0,
                          vjust = -0.5),
                      color = "black",
                      size = 4.0,
                      check_overlap = FALSE,
                      parse = FALSE,
                      position = position_dodge(width = 0.5),
                      inherit.aes = TRUE) +
            geom_label(aes(x = 3.5,
                           y = -Inf,
                           hjust = 0,
                           vjust = 1,
                           label = sprintf("%.2f", DIC)),
                       fill = "beige",
                       colour = "black",
                       fontface = "plain",
                       size = 3.1) +
            scale_y_continuous(trans = "identity") +
            facet_wrap(vars(factor(node, levels = unique(prepare$node))),
                       scales = "fixed") +
            labs(x = "",
                 y = ifelse(
                   is.element(measure, c("OR", "ROM")),
                   paste(measure2, "(in logarithmic scale)"), measure2),
                 colour = "Evidence on inconsistency") +
            coord_flip() +
            scale_color_manual(breaks = c("strong evidence",
                                          "weak evidence"),
                               values = c("#009E73", "#D55E00"),
                               na.value = "black") +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12),
                  axis.text.y = element_text(color = "black", size = 12),
                  axis.title.x = element_text(color = "black", size = 12,
                                              face = "bold"),
                  legend.position = "bottom",
                  legend.title = element_text(color = "black", size = 12,
                                              face = "bold"),
                  legend.text = element_text(color = "black", size = 12),
                  strip.text = element_text(size = 11))
  } else {
    # Keep nodes with strong evidence inconsistency (selection1) or with
    # inconsistent sign in the direct and indirect estimate (selection2)
    selection1 <- subset(prepare, stat_sign == "strong evidence")
    selection2 <- subset(prepare, (median[evidence == "direct"] < 0 &
                             median[evidence == "indirect"] > 0) |
                          (median[evidence == "direct"] > 0 &
                             median[evidence == "indirect"] < 0))
    selection  <- rbind(subset(prepare, is.element(node, selection1[, 1])),
                        selection2)
    p1 <- ggplot(data = selection,
                 aes(x = factor(evidence,
                                levels = c("IF", "indirect", "direct")),
                     y = median,
                     ymin = lower,
                     ymax = upper,
                     colour = stat_sign)) +
            geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
            geom_hline(yintercept = 0, lty = 1, size = 1, col = "grey") +
            geom_point(size = 1.5,
                       colour = "white",
                       stroke = 0.3,
                       position = position_dodge(width = 0.5)) +
            geom_text(aes(x = as.factor(evidence),
                          y = round(median, 2),
                          label = paste0(sprintf("%.2f", median), " ", "(",
                                         sprintf("%.2f", lower), ",", " ",
                                         sprintf("%.2f", upper), ")"),
                          hjust = 0,
                          vjust = -0.5),
                      color = "black",
                      size = 4.0,
                      check_overlap = FALSE,
                      parse = FALSE,
                      position = position_dodge(width = 0.5),
                      inherit.aes = TRUE) +
            geom_label(aes(x = 3.5,
                           y = -Inf,
                           hjust = 0,
                           vjust = 1,
                           label = round(DIC, 0)),
                       fill = "beige",
                       colour = "black",
                       fontface = "plain",
                       size = 3.1) +
            facet_wrap(vars(factor(node, levels = unique(prepare$node))),
                       scales = "fixed") +
            scale_y_continuous(trans = "identity") +
            labs(x = "",
                 y = ifelse(
                   is.element(measure, c("OR", "ROM")),
                   paste(measure2, "(in logarithmic scale)"), measure2),
                 colour = "Evidence on inconsistency") +
            coord_flip() +
            scale_color_manual(breaks = c("strong evidence",
                                          "weak evidence"),
                               values = c("#009E73", "#D55E00"),
                               na.value = "black") +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12),
                  axis.title.x = element_text(color = "black", size = 12,
                                              face = "bold"),
                  axis.text.y = element_text(color = "black", size = 12),
                  legend.position = "bottom",
                  legend.title = element_text(color = "black", size = 12,
                                              face = "bold"),
                  legend.text = element_text(color = "black", size = 12),
                  strip.text = element_text(size = 11))
  }

  # Create a table on the direct, indirect and IF per split node
  cri_direct <- paste0("(", round(direct[, 5], 2), ",",
                       " ", round(direct[, 6], 2), ")",
                       ifelse(direct[, 5] > 0 | direct[, 6] < 0, "*", " "))
  cri_indirect <- paste0("(", round(indirect[, 5], 2), ",",
                         " ", round(indirect[, 6], 2), ")",
                         ifelse(indirect[, 5] > 0 |
                                  indirect[, 6] < 0, "*", " "))
  cri_if <- paste0("(", round(incons_factor[, 5], 2), ",", " ",
                   round(incons_factor[, 6], 2), ")",
                   ifelse(incons_factor[, 5] > 0 |
                            incons_factor[, 6] < 0, "*", " "))
  table_em <- data.frame(comp, round(direct[, 3:4], 2), cri_direct,
                         round(indirect[, 3:4], 2), cri_indirect,
                         round(incons_factor[, 3:4], 2), cri_if)
  colnames(table_em) <- c("Node",
                          "Median direct",
                          "SD direct",
                          "95% CrI direct",
                          "Median indirect",
                          "SD indirect",
                          "95% CrI indirect",
                          "Median IF",
                          "SD IF",
                          "95% CrI IF")
  rownames(table_em) <- NULL

  # Find whether at least one split node improve the fit of the model
  model_selection <- data.frame(comp,
                                model_assess_sort[, 4] -
                                  rep(model_assess_nma$DIC,
                                      length(model_assess_sort[, 4])))
  colnames(model_selection) <- c("Comparison", "dic_diff")
  better_fit <- c(ifelse(model_selection$dic_diff > 5, "Consistency model",
                         ifelse(model_selection$dic_diff < -5,
                                "After split node", "Little to choose")))

  ## Table on the model assessment measures and between-trial standard deviation
  # per split node
  if (model == "RE") {
    cri_tau <- paste0("(", round(tau[, 5], 2), ",", " ",
                      round(tau[, 6], 2), ")")
    table_assess0 <- data.frame(comp,
                                round(model_assess_sort[, -c(1:2)], 2),
                                better_fit,
                                round(tau[, 3:4], 2),
                                cri_tau)
    colnames(table_assess0) <- c("Approach",
                                 "Residual deviance", "DIC", "pD",
                                 "DIC-based better fit",
                                 "Median tau", "SD tau",
                                 "95% CrI tau")
    add <- data.frame("NMA", round(model_assess_nma[c(3, 1, 2)], 2), "-",
                      round(tau_values[1], 2), round(tau_values[2], 2),
                      paste0("(", round(tau_values[3], 2), ",", " ",
                             round(tau_values[4], 2), ")"))
    colnames(add) <- colnames(table_assess0)
    table_assess <- rbind(add, table_assess0)
  } else {
    table_assess0 <- data.frame(comp,
                                round(model_assess_sort[, -c(1:2)], 2),
                                better_fit)
    colnames(table_assess0) <- c("Approach",
                                 "Residual deviance", "DIC", "pD",
                                 "DIC-based better fit")
    add <- data.frame("NMA", round(model_assess_nma[c(3, 1, 2)], 2), "-")
    colnames(add) <- colnames(table_assess0)
    table_assess <- rbind(add, table_assess0)
  }
  rownames(table_assess) <- NULL

  ## Dataset to create the panel of interval plot on the 'between-trial standard
  # deviation' for each split node
  if (model == "RE") {
    prepare_tau <- data.frame(comp, tau[, c(3, 5:6)], sort(model_assess$DIC))
    colnames(prepare_tau) <- c("comp", "median", "lower", "upper", "DIC")
  } else {
    prepare_tau <- NA
  }

  # Create the interval plot for 'between-trial standard deviation'
  p2 <- if (model == "RE") {
    ggplot(data = prepare_tau,
           aes(x = as.factor(seq_len(length(comp))),
               y = median, ymin = lower, ymax = upper)) +
      #geom_rect(aes(xmin = 0,
      #              xmax = Inf,
      #              ymin = tau_values[3],
      #              ymax = tau_values[4]),
      #          fill = "#D55E00",
      #          alpha = 0.1) +
      geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0, ymax = 0.099,
                    fill = "low"),
                alpha = 0.02) +
      geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0.1, ymax = 0.5,
                    fill = "reasonable"),
                alpha = 0.02) +
      geom_rect(aes(xmin = 0, xmax = Inf, ymin = 0.5, ymax = 1.0,
                    fill = "fairly high"),
                alpha = 0.02) +
      geom_rect(aes(xmin = 0, xmax = Inf, ymin = 1.0, ymax = Inf,
                    fill = "fairly extreme"),
                alpha = 0.02) +
      geom_hline(yintercept = tau_values[1],
                 lty = 1,
                 size = 1,
                 col = "#006CD1") +
      geom_hline(yintercept = tau_values[3],
                 lty = 3,
                 size = 1,
                 col = "#006CD1") +
      geom_hline(yintercept = tau_values[4],
                 lty = 3,
                 size = 1,
                 col = "#006CD1") +
      geom_linerange(size = 2,
                     position = position_dodge(width = 0.5)) +
      #geom_hline(yintercept = tau_values[1],
      #           lty = 1,
      #           size = 1,
      #           col = "#D55E00") +
      geom_point(size = 1.5,
                 colour = "white",
                 stroke = 0.3,
                 position = position_dodge(width = 0.5)) +
      geom_text(aes(x = as.factor(seq_len(length(comp))),
                    y = round(median, 2),
                    label = sprintf("%.2f", median),
                    hjust = -0.2,
                    vjust = 0.3),
                color = "black",
                size = 4.0,
                check_overlap = FALSE,
                parse = FALSE,
                position = position_dodge(width = 0.5),
                inherit.aes = TRUE) +
      geom_label(aes(x = as.factor(seq_len(length(comp))),
                     y = upper,
                     label = sprintf("%.2f", DIC)),
                 fill = "beige",
                 colour = "black",
                 fontface = "plain",
                 size = 3.1) +
      scale_x_discrete(breaks = as.factor(seq_len(length(comp))),
                       labels = comp[seq_len(length(comp))]) +
      scale_fill_manual(name = "Heterogeneity",
                        values = c("low" = "#009E73",
                                   "reasonable" = "orange",
                                   "fairly high" = "#D55E00",
                                   "fairly extreme" = "red"),
                        limits = c("low", "reasonable", "fairly high",
                                   "fairly extreme")) +
      labs(x = "Split nodes (sorted by DIC in ascending order)",
           y = "Between-trial standard deviation") +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12, angle = 45,
                                       hjust = 1),
            axis.text.y = element_text(color = "black", size = 12),
            legend.position = "bottom",
            legend.text =  element_text(color = "black", size = 12),
            legend.title =  element_text(color = "black", face = "bold",
                                         size = 12))
  } else {
    NA
  }

  # Write the table with the EMs from both models as .xlsx
  if (save_xls == TRUE) {
    write_xlsx(table_em, paste0("Table NMA vs Node-Split", ".xlsx"))
    write_xlsx(table_assess, paste0("Table assesssment Node-Split", ".xlsx"))
  }

  # Return results
  results <- if (model == "RE") {
    list(table_effect_size =
           knitr::kable(table_em,
                        align = "lccccccccc",
                        caption = "Estimates for the split nodes"),
         table_model_assessment =
           knitr::kable(table_assess,
                        align = "lccclccc",
                        caption = paste0("Model assessment parameters (",
                                         data_points, " ",
                                         "unconstrained data points)")),
         intervalplot_inconsistency_factor = p1,
         intervalplot_tau = p2 +
           guides(fill = guide_legend(override.aes = list(alpha = 0.4))))
  } else {
    list(table_effect_size =
           knitr::kable(table_em,
                        align = "lccccccccc",
                        caption = "Estimates fot split nodes"),
         table_model_assessment =
           knitr::kable(table_assess,
                        align = "lcccl",
                        caption =paste0("Model assessment parameters (",
                                        data_points, " ",
                                        "unconstrained data points)")),
         intervalplot_inconsistency_factor = p1)
  }

  return(results)
}
