#' End-user-ready results: consistency model versus node-splitting approach
#'
#' @description This function illustrates the direct and indirect effects and inconsistency factor of the split nodes in a panel of interval plots and also exports these results in an Excel format.
#'   Furthermore, \code{nodesplit.plot} facilitates the comparison of the consistency model (via \code{run.model}) with the node-splitting approach (via \code{run.nodesplit}) regarding between-trial standard deviation (\eqn{\tau}) and model assessment parameters
#'   after each split node of the network.
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param node An object of S3 class \code{\link{run.nodesplit}}. See 'Value' in \code{\link{run.nodesplit}}.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @return \code{nodesplit.plot} returns a panel of as many interval plots as the number of split nodes in the network. Each interval plot illustrates the posterior mean and 95\% credible interval of the direct and indirect effect of the split nodes
#'   and the corresponding inconsistency factor. The line that corresponds to the inconsistency factor is highlighted with green, when it does not cross the vertical line of no difference (between the direct and indirect effect), and red otherwise.
#'   If there are more than 30 split nodes, the function presents the interval plots on split nodes with statistically significant inconsistency factor or those with inconsistent sign in the direct and indirect effect.
#'
#'   Furthermore, the function returns a interval plot on the median and 95\% credible interval of \eqn{\tau} after each split node. The lines that correspond to the split nodes are sorted in ascending order of the
#'   deviance information criterion (DIC) which appears at the top of each line. The 95\% credible interval of \eqn{\tau} under the consistency model appears as a rectangle in the interval plot.
#'   When a fixed-effect model has been performed, \code{nodesplit.plot} does not return the interval plot of \eqn{\tau}.
#'
#'   The R console prints the data-frame with the posterior mean, posterior standard deviation and 95\% credible interval of the direct and indirect effect and the inconsistency factor of each split.
#'   The console also prints the data-frame with the model assessment parameters (DIC, posterior mean of total residual deviance, and number of effective parameters), the posterior median, posterior standard deviation and 95\% credible interval of \eqn{\tau}
#'   under the consistency model and after each split node. The DIC of the model after each split node is compared with the DIC of the consistency model (Spiegelhalter et al. (2002), Dias et al. (2010)). If the difference in DIC exceeds 5, the consistency model is preferred; if the difference in DIC is less than -5,
#'   the model after split node is preferred; otherwise, there is little to choose between the compared models.
#'
#'   Furthermore, \code{nodesplit.plot} exports both data-frames to an Excel 'xlsx' format (via the \code{\link[writexl]{write_xlsx}} function) to the working directory of the user.
#'
#' @details \code{nodesplit.plot} can be used only for a network of interventions. In the case of two interventions, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{run.nodepslit}}
#'
#' @references
#' Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed treatment comparison meta-analysis. \emph{Stat Med} 2010;\bold{29}(7-8):932--44. [\doi{10.1002/sim.3767}]
#'
#' Spiegelhalter DJ, Best NG, Carlin BP, van der Linde A. Bayesian measures of model complexity and fit. \emph{J R Stat Soc B} 2002;\bold{64}:583--616. [\doi{10.1111/1467-9868.00353}]
#'
#' @examples
#' data("nma.baker2009.RData")
#'
#' # Perform a random-effects network meta-analysis
#' res1 <- run.model(data = nma.baker2009, measure = "OR", model = "RE", assumption = "IDE-ARM", heter.prior = list("halfnormal", 0, 1), mean.misspar = 0, var.misspar = 1, D = 1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' # Run random-effects network meta-analysis with node-splitting approachs
#' node1 <- run.nodesplit(data = nma.baker2009, full = res1, n.chains = 3, n.iter = 10000, n.burnin = 1000, n.thin = 1)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("budesodine", "budesodine plus formoterol", "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium", "placebo")
#'
#' # Plot the results from the consistency model and the node-splitting approach
#' nodesplit.plot(full = res1, node = node1, drug.names = interv.names)
#'
#' @export
nodesplit.plot <- function(full, node, drug.names) {


  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug.names
  }


  if(length(drug.names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }

  ## Keep tau and model assessment measures from NMA model
  tau.values <- full$tau[c(5, 2, 3, 7)]
  model.assess.NMA <- full$model.assessment

  # Effect measure
  measure <- effect.measure.name(full$measure)



  # Analysis model
  model <- full$model



  ## Keep results on 'direct evidence', 'indirect evidence', 'inconsistency factor', 'between-trial standard deviation',
  ## and model assessment measures (i.e., DIC, posterior mean of refisual deviance, and pD)
  direct0 <- node$direct; indirect0 <- node$indirect; IF0 <- node$diff; model.assess <- node$model.assessment


  ## Sort 'direct evidence', 'indirect evidence', 'inconsistency factor', and 'between-trial standard deviation' by DIC in ascending order
  direct <- direct0[order(model.assess$DIC), ]
  indirect <- indirect0[order(model.assess$DIC), ]
  IF <- IF0[order(model.assess$DIC), ]
  if (model == "RE") {
    tau <- node$tau[order(model.assess$DIC), ]
  } else {
    tau <- NA
  }
  model.assess.sort <- model.assess[order(model.assess$DIC), ]


  ## Interventions' name: Replace code with original names
  # For treat1 (non-baseline arm)
  for(i in sort(unique(unlist(direct[, 1])))) {
    direct[direct$treat1 == i, 1] <- drug.names[i]
  }

  # For treat2 (baseline arm)
  for(i in sort(unique(unlist(direct[, 2])))) {
    direct[direct$treat2 == i, 2] <- drug.names[i]
  }


  ## Prepare the dataset to create the panel of interval plots on the 'direct evidence', 'indirect evidence', and 'inconsistency factor' for each split node
  comp <- paste(direct[, 1], "vs", direct[, 2])
  prepare <- data.frame(rep(comp, 3), rbind(direct[, c(3, 5:6)], indirect[, c(3, 5:6)], IF[, c(3, 5:6)]), rep(c("direct", "indirect", "IF"), each = length(direct[, 1])))
  colnames(prepare) <- c("node", "mean", "lower", "upper", "evidence")
  prepare$stat.signif <- ifelse(prepare$lower > 0 | prepare$upper < 0  , "statistically significant", "statistically non-significant")
  prepare$stat.signif <- ifelse(prepare$evidence != "IF", NA, prepare$stat.signif)
  prepare$DIC <- sort(model.assess$DIC)


  ## Create the panel of interval plots
  if(length(unique(comp)) <= 30) {

    p1 <- ggplot(data = prepare, aes(x = factor(evidence, levels = c("IF", "indirect", "direct")), y = mean, ymin = lower, ymax = upper, colour = stat.signif) ) +
            geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
            geom_hline(yintercept = 0, lty = 2, size = 1.5, col = "grey") +
            geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
            geom_text(aes(x = as.factor(evidence), y = round(mean, 2), label = round(mean, 2), hjust = 0, vjust = -0.4), color = "black", size = 4.0,
                      check_overlap = F, parse = F, position = position_dodge(width = 0.5),  inherit.aes = T) +
            geom_label(aes(x = 3.5, y = -Inf, hjust = 0, vjust = 1, label = round(DIC, 0)), fill = "beige", colour = "black", fontface = "plain", size = 3.1) +
            scale_y_continuous(trans = "identity") +
            facet_wrap(vars(factor(node, levels = unique(prepare$node))), scales = "free_x") +
            labs(x = "", y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), paste(measure, "(in logarithmic scale)"), measure), colour = "") +
            coord_flip() +
            scale_color_manual(breaks = c("statistically significant", "statistically non-significant"), values = c("#009E73", "#D55E00"), na.value = "black") +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), axis.title.x = element_text(color = "black", size = 12, face = "bold"),
                  legend.position = "none", legend.title = element_text(color = "black", size = 12, face = "bold"), legend.text = element_text(color = "black", size = 12),
                  strip.text = element_text(size = 11))

  } else {

    # Keep nodes with statistically significant inconsistency OR with inconsistent sign in the direct and indirect estimate
    selection <- subset(prepare, stat.signif == "statistically significant" |
               (mean[evidence == "direct"] < 0 & mean[evidence == "indirect"] > 0) |
               (mean[evidence == "direct"] > 0 & mean[evidence == "indirect"] < 0))

    p1 <- ggplot(data = selection, aes(x = factor(evidence, levels = c("IF", "indirect", "direct")), y = mean, ymin = lower, ymax = upper, colour = stat.signif) ) +
            geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
            geom_hline(yintercept = 0, lty = 2, size = 1.5, col = "grey") +
            geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
            geom_text(aes(x = as.factor(evidence), y = round(mean, 2), label = round(mean, 2), hjust = 0, vjust = -0.4), color = "black", size = 4.0,
                      check_overlap = F, parse = F, position = position_dodge(width = 0.5),  inherit.aes = T) +
            geom_label(aes(x = 3.5, y = -Inf, hjust = 0, vjust = 1, label = round(DIC, 0)), fill = "beige", colour = "black", fontface = "plain", size = 3.1) +
            facet_wrap(vars(factor(node, levels = unique(prepare$node))), scales = "free_x") +
            scale_y_continuous(trans = "identity") +
            labs(x = "", y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), paste(measure, "(in logarithmic scale)"), measure), colour = "Evidence on inconsistency") +
            coord_flip() +
            scale_color_manual(breaks = c("statistically significant", "statistically non-significant"), values = c("#009E73", "#D55E00"), na.value = "black") +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.title.x = element_text(color = "black", size = 12, face = "bold"), axis.text.y = element_text(color = "black", size = 12),
                  legend.position = "bottom", legend.title = element_text(color = "black", size = 12, face = "bold"), legend.text = element_text(color = "black", size = 12),
                  strip.text = element_text(size = 11))
  }



  ## Create a table on the direct, indirect and IF per split node
  CrI.direct <- paste0("(", round(direct[, 5], 2), ",", " ", round(direct[, 6], 2), ")", ifelse(direct[, 5] > 0 | direct[, 6] < 0, "*", " "))
  CrI.indirect <- paste0("(", round(indirect[, 5], 2), ",", " ", round(indirect[, 6], 2), ")", ifelse(indirect[, 5] > 0 | indirect[, 6] < 0, "*", " "))
  CrI.IF <- paste0("(", round(IF[, 5], 2), ",", " ", round(IF[, 6], 2), ")", ifelse(IF[, 5] > 0 | IF[, 6] < 0, "*", " "))
  table.EM <- data.frame(comp, round(direct[, 3:4], 2), CrI.direct, round(indirect[, 3:4], 2), CrI.indirect, round(IF[, 3:4], 2), CrI.IF)
  colnames(table.EM) <- c("Split node", "Post. mean dir.", "Post. SD dir.", "95% CrI dir.", "Post. mean indir", "Post. SD indir.", "95% CrI indir.",
                         "Post. mean IF", "Post. SD IF", "95% CrI IF")



  ## Find whether at least one split node improve the fit of the model
  model.selection <- data.frame(comp, model.assess.sort[, 3] - rep(model.assess.NMA$DIC, length(model.assess.sort[, 3])))
  colnames(model.selection) <- c("Comparison", "DIC.diff")
  Better.fit <- c(ifelse(model.selection$DIC.diff > 5, "Consistency model", ifelse(model.selection$DIC.diff < -5, "After split node", "Little to choose")))



  ## Create a table on the model assessment measures and between-trial standard deviation per split node
  if (model == "RE") {
    CrI.tau <- paste0("(", round(tau[, 5], 2), ",", " ", round(tau[, 6], 2), ")")
    table.assess0 <- data.frame(comp, round(model.assess.sort[, -c(1:2)], 2), Better.fit, round(tau[, 3:4], 2), CrI.tau)
    colnames(table.assess0) <- c("Approach", "DIC", "Post. mean dev.", "pD", "DIC-based better fit", "Post. median tau", "Post. SD tau", "95% CrI tau")
    add <- data.frame("NMA", round(model.assess.NMA[c(1, 3, 2)], 2), "-", round(tau.values[1], 2), round(tau.values[2], 2), paste0("(", round(tau.values[3], 2), ",", " ", round(tau.values[4], 2), ")"))
    colnames(add) <- colnames(table.assess0)
    table.assess <- rbind(add, table.assess0)
  } else {
    table.assess0 <- data.frame(comp, round(model.assess.sort[, -c(1:2)], 2), Better.fit)
    colnames(table.assess0) <- c("Approach", "DIC", "Post. mean dev.", "pD", "DIC-based better fit")
    add <- data.frame("NMA", round(model.assess.NMA[c(1, 3, 2)], 2), "-")
    colnames(add) <- colnames(table.assess0)
    table.assess <- rbind(add, table.assess0)
  }



  ## Prepare the dataset to create the panel of interval plot on the 'between-trial standard deviation' for each split node
  if (model == "RE") {
    prepare.tau <- data.frame(comp, tau[, c(3, 5:6)], sort(model.assess$DIC))
    colnames(prepare.tau) <- c("node", "median", "lower", "upper", "DIC")
  } else {
    prepare.tau <- NA
  }


  ## Create the interval plot for 'between-trial standard deviation'
  p2 <- if (model == "RE") {
    ggplot(data = prepare.tau, aes(x = as.factor(1:length(prepare.tau$node)), y = median, ymin = lower, ymax = upper) ) +
      geom_rect(aes(xmin = 0, xmax = Inf, ymin = tau.values[3], ymax = tau.values[4]), fill = "#D55E00", alpha = 0.1) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = tau.values[1], lty = 1, size = 1, col = "#D55E00") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = as.factor(1:length(prepare.tau$node)), y = round(median, 2), label = round(median, 2), hjust = -0.2, vjust = 0.3), color = "black", size = 4.0,
                check_overlap = F, parse = F, position = position_dodge(width = 0.5),  inherit.aes = T) +
      geom_label(aes(x = as.factor(1:length(prepare.tau$node)), y = upper, label = round(DIC, 0)), fill = "beige", colour = "black", fontface = "plain",  size = 3.1) +
      scale_x_discrete(breaks = as.factor(1:length(prepare.tau$node)), labels = prepare.tau$node[1:length(prepare.tau$node)]) +
      labs(x = "Split nodes (sorted by DIC in ascending order)", y = "Between-trial standard deviation") +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1), axis.text.y = element_text(color = "black", size = 12), legend.position = "none")
  } else {
    NA
  }



  ## Write the table with the EMs from both models as .xlsx
  write_xlsx(table.EM, paste0(getwd(),"Table NMA vs Node-Split.xlsx"))
  write_xlsx(table.assess, paste0(getwd(),"Table assesssment Node-Split.xlsx"))


  ## Return results
  results <- if (model == "RE") {
    list(table.EM = table.EM,
         table.assess = table.assess,
         node.split.intervalplot = p1,
         tau.intervalplot = p2)
  } else {
    list(table.EM = table.EM,
         table.assess = table.assess,
         node.split.intervalplot = p1)
  }

  return(results)
}
