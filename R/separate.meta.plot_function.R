#' Plot the results from the unrelated mean effects model
#'
#' @export
separate.meta.plot <- function(full, meta, drug.names) {


  if(length(drug.names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }

  ## The results on the following parameters will be used:
  # Posterior results on the effect estimates under NMA
  EM.full <-  full$EM

  # Posterior results on the effect estimates under separate random-effect pairwise meta-analysis (RE-MAs)
  EM.meta <- meta$EM

  # Posterior results on between-trial standard deviation under NMA
  tau.full <- full$tau

  # Posterior results on between-trial standard deviation under RE-MAs
  tau.meta <- meta$tau

  # Effect measure
  measure <- if (full$measure != meta$measure) {
    stop("The argument 'measure' differs in 'run.model' and 'run.separate.meta'. Specify the same 'measure' and run the analysis again")
  } else {
    effect.measure.name(full$measure)
  }

  # Analysis model
  model <- if (full$model != meta$model) {
    stop("The argument 'model' differs in 'run.model' and 'run.separate.meta'. Specify the same 'model' and run the analysis again")
  } else {
    full$model
  }


  # Possible and observed comparisons
  possible.comp <- possible.observed.frail.comparisons(drug.names, obs.comp = paste0(meta$EM[, "t2"], "vs", meta$EM[, "t1"]))

  # Observed comparisons
  obs.comp <- possible.comp$obs.comp[, 3]


  ## Keep only the effect estimates according to the 'poss.pair.comp.clean' - Consistency model
  EM.full.clean <- format(round(EM.full[is.element(possible.comp$poss.comp[, 4], obs.comp), c(1:3, 7)], 2), nsmall = 2)



  ## Effect estimate of separate RE-MAs
  EM.meta.clean <- format(round(EM.meta[, c(3:5, 9)], 2), nsmall = 2)



  ## Between-trial standard deviation of separate RE-MAs
  if (model == "RE") {
    tau.meta.clean <- format(round(tau.meta[, c(3:5, 9)], 2), nsmall = 2)
  } else {
    tau.meta.clean <- NA
  }



  ## Keep only the 95% credible intervals (CrI) according to the 'poss.pair.comp.clean' - Consistency model
  CrI.full.clean <- paste0("(", EM.full.clean[, 3], ",", " ", EM.full.clean[, 4], ")", ifelse(as.numeric(EM.full.clean[, 3]) > 0 | as.numeric(EM.full.clean[, 4]) < 0, "*", " "))



  ## The 95% CrIs of the effect estimate of separate RE-MAs
  CrI.meta.clean <- paste0("(", EM.meta.clean[, 3], ",", " ", EM.meta.clean[, 4], ")", ifelse(as.numeric(EM.meta.clean[, 3]) > 0 | as.numeric(EM.meta.clean[, 4]) < 0, "*", " "))



  ## The 95% CrIs of the between-trial standard deviation of separate RE-MAs
  CrI.tau.meta <- if (model == "RE") {
    paste0("(", tau.meta.clean[, 3], ",", " ", tau.meta.clean[, 4], ")")
  } else {
    NA
  }



  ## Create a data-frame with effect estimates on both models
  if (model == "RE") {
    EM.both.models <- data.frame(possible.comp$obs.comp[, 4], EM.full.clean[, 1:2], CrI.full.clean, EM.meta.clean[, 1:2], CrI.meta.clean, tau.meta.clean[, 1:2], CrI.tau.meta)
    colnames(EM.both.models) <- c("Comparison", "Posterior mean NMA", "Posterior SD NMA", "95% CrI NMA", "Posterior mean MA",
                                  "Posterior SD MA", "95% CrI MA", "Posterior median tau", "Posterior SD tau", "95% CrI tau")
  } else {
    EM.both.models <- data.frame(possible.comp$obs.comp[, 4], EM.full.clean[, 1:2], CrI.full.clean, EM.meta.clean[, 1:2], CrI.meta.clean)
    colnames(EM.both.models) <- c("Comparison", "Posterior mean NMA", "Posterior SD NMA", "95% CrI NMA", "Posterior mean MA",
                                  "Posterior SD MA", "95% CrI MA")
  }
  rownames(EM.both.models) <- NULL



  ## Prepare dataset for ggplot2
  # Effect estimate
  #if (is.element(measure, c("Odds ratio", "Ratio of means"))) {
  #  prepare <- data.frame(rep(1:length(obs.comp), 2), rep(possible.comp$obs.comp[, 4], 2), rbind(apply(apply(EM.full.clean[, -2], 2, as.numeric), 2, exp), apply(apply(EM.meta.clean[, -2], 2, as.numeric), 2, exp)), rep(c("Network meta-analysis", "Paiwise meta-analysis"), each = length(obs.comp)))
  #} else  {
    prepare <- data.frame(rep(1:length(obs.comp), 2), rep(possible.comp$obs.comp[, 4], 2), rbind(apply(EM.full.clean[, -2], 2, as.numeric), apply(EM.meta.clean[, -2], 2, as.numeric)), rep(c("Network meta-analysis", "Paiwise meta-analysis"), each = length(obs.comp)))
  #}
  colnames(prepare) <- c("order", "comparison", "mean", "lower", "upper", "analysis")
  rownames(prepare) <- NULL

  # Between-trial standard deviation
  if (model == "RE") {
    prepare.tau <- data.frame(possible.comp$obs.comp[, 4], tau.meta.clean[, -2])
    colnames(prepare.tau) <- c("comparison", "median", "lower", "upper")
  } else {
    prepare.tau <- NA
  }




  ## Forest plots of comparisons on effect estimate
  p1 <- ggplot(data = prepare, aes(x = as.factor(order), y = mean, ymin = lower, ymax = upper, colour = analysis, group = analysis)) +
          geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
          #geom_hline(yintercept = ifelse(measure != "Odds ratio" & measure != "Ratio of means", 0, 1), lty = 2, size = 1.3, col = "grey53") +
          geom_hline(yintercept = 0, lty = 2, size = 1.3, col = "grey53") +
          geom_point(size = 1.5,  colour = "black", stroke = 0.3, position = position_dodge(width = 0.5)) +
          geom_text(aes(x = as.factor(order), y = round(as.numeric(mean), 2), label = round(as.numeric(mean), 2)),
                    color = "black", hjust = -0.1, vjust = -0.5, size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.8), inherit.aes = T) +
          labs(x = "", y = ifelse(is.element(measure, c("Odds ratio", "Ratio of means")), paste(measure, "(in logarithmic scale)"), measure), colour = "Analysis") +
          scale_x_discrete(breaks = as.factor(1:length(obs.comp)), labels = prepare$comparison[1:length(obs.comp)]) +
          #scale_y_continuous(trans = ifelse(measure != "Odds ratio" & measure != "Ratio of means", "identity", "log10")) +
          scale_y_continuous(trans = "identity") +
          scale_color_manual(breaks = c("Network meta-analysis", "Paiwise meta-analysis"), values = c("#009E73", "#D55E00")) +
          coord_flip() +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "none", legend.justification = c(0.13, 0),
                legend.text =  element_text(color = "black", size = 12), legend.title =  element_text(color = "black", face = "bold", size = 12))



  ## Forest plots of comparisons-specific between-trial standard deviation
  p2 <- if (model == "RE") {
    ggplot(data = prepare.tau, aes(x = as.factor(1:length(obs.comp)), y = as.numeric(median), ymin = as.numeric(lower), ymax = as.numeric(upper))) +
      geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
      geom_hline(yintercept = tau.full[3], lty = 2, size = 1.3, col = "#006CD1") +
      geom_hline(yintercept = tau.full[5], lty = 2, size = 1.3, col = "#006CD1") +
      geom_hline(yintercept = tau.full[7], lty = 2, size = 1.3, col = "#006CD1") +
      geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
      geom_text(aes(x = as.factor(1:length(obs.comp)), y = round(as.numeric(median), 2), label = round(as.numeric(median), 2)),
                color = "black", hjust = -0.1, vjust = -0.5, size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.8), inherit.aes = T) +
      scale_x_discrete(breaks = as.factor(1:length(obs.comp)), labels = prepare.tau$comparison[1:length(obs.comp)]) +
      labs(x = "", y = "Between-trial standard deviation") +
      coord_flip() +
      theme_classic() +
      theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
            axis.title.x = element_text(color = "black", face = "bold", size = 12))
  } else {
    NA
  }



  ## Bring together both forest-plots
  forest.plots <- if (model == "RE") {
    ggarrange(p1, p2, nrow = 1, ncol = 2, labels = c("A)", "B)"), common.legend = T, legend = "bottom")
  } else {
    p1
  }



  ## Write the table as .xlsx
  write_xlsx(EM.both.models, paste0(getwd(),"Table NMA vs MA.xlsx"))


  return(list(EM.both.models = EM.both.models, forest.plots = forest.plots))



}
