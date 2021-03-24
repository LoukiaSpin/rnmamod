#' Plot the results from the unrelated mean effects model
#'
#' @export
separate.meta.plot <- function(full, meta, drug.names) {



  ## The results on the following parameters will be used:
  # Posterior results on the effect estimates under NMA
  EM.full <- full$EM

  # Posterior results on the effect estimates under separate random-effect pairwise meta-analysis (RE-MAs)
  EM.meta <- meta$EM

  # Posterior results on between-trial standard deviation under NMA
  tau.full <- full$tau

  # Posterior results on between-trial standard deviation under RE-MAs
  tau.meta <- meta$tau

  # Comparisons with at leat two trials
  obs.comp <- paste0(EM.meta$t2, "vs", EM.meta$t1)


  ## Obtain all unique pairwise comparisons using the 'combn' functions
  nt <- (1 + sqrt(1 + 8*(length(EM.full[, 1]))))/2    # The quadratic formula for the roots of the general quadratic equation
  poss.pair.comp <- data.frame(t(combn(1:nt, 2))[, 2], t(combn(1:nt, 2))[, 1])
  colnames(poss.pair.comp) <- c("treat1", "treat2")


  poss.pair.comp$comp <- paste0(poss.pair.comp[, 1], "vs", poss.pair.comp[, 2])
  poss.pair.comp.clean <- poss.pair.comp[is.element(poss.pair.comp$comp, obs.comp), ]


  ## Replace intervention id with their original name
  # For treat1 (non-baseline arm)
  for(i in sort(unique(unlist(poss.pair.comp.clean[, 1])))) {
    poss.pair.comp.clean[poss.pair.comp.clean$treat1 == i, 1] <- drug.names[i]
  }

  # For treat2 (baseline arm)
  for(i in sort(unique(unlist(poss.pair.comp.clean[, 2])))) {
    poss.pair.comp.clean[poss.pair.comp.clean$treat2 == i, 2] <- drug.names[i]
  }



  ## Create a vector with the comparisons (non-baseline versus beseline) using the 'paste' function
  comparison <- paste(poss.pair.comp.clean[, 1], "vs", poss.pair.comp.clean[, 2])



  ## Keep only the effect estimates according to the 'poss.pair.comp.clean' - Consistency model
  EM.full.clean <- format(round(EM.full[is.element(poss.pair.comp$comp, obs.comp), c(1:3, 7)], 2), nsmall = 2)



  ## Effect estimate of separate RE-MAs
  EM.meta.clean <- format(round(EM.meta[, c(3:5, 9)], 2), nsmall = 2)



  ## Between-trial standard deviation of separate RE-MAs
  tau.meta.clean <- format(round(tau.meta[, c(3:5, 9)], 2), nsmall = 2)



  ## Keep only the 95% credible intervals (CrI) according to the 'poss.pair.comp.clean' - Consistency model
  CrI.full.clean <- paste0("(", EM.full.clean[, 3], ",", " ", EM.full.clean[, 4], ")", ifelse(as.numeric(EM.full.clean[, 3]) > 0 | as.numeric(EM.full.clean[, 4]) < 0, "*", " "))



  ## The 95% CrIs of the effect estimate of separate RE-MAs
  CrI.meta.clean <- paste0("(", EM.meta.clean[, 3], ",", " ", EM.meta.clean[, 4], ")", ifelse(as.numeric(EM.meta.clean[, 3]) > 0 | as.numeric(EM.meta.clean[, 4]) < 0, "*", " "))



  ## The 95% CrIs of the between-trial standard deviation of separate RE-MAs
  CrI.tau.meta <- paste0("(", tau.meta.clean[, 3], ",", " ", tau.meta.clean[, 4], ")")



  ## Create a data-frame with effect estimates on both models
  EM.both.models <- data.frame(comparison, EM.full.clean[, 1:2], CrI.full.clean, EM.meta.clean[, 1:2], CrI.meta.clean, tau.meta.clean[, 1:2], CrI.tau.meta)
  colnames(EM.both.models) <- c("Comparison", "Posterior mean NMA", "Posterior SD NMA", "95% CrI NMA", "Posterior mean MA",
                                "Posterior SD MA", "95% CrI MA", "Posterior median tau", "Posterior SD tau", "95% CrI tau")
  rownames(EM.both.models) <- NULL



  ## Prepare dataset for ggplot2
  # Effect estimate
  prepare <- data.frame(rep(1:length(obs.comp), 2), rep(comparison, 2), rbind(EM.full.clean[, -2], EM.meta.clean[, -2]), rep(c("Network meta-analysis", "Paiwise meta-analysis"), each = length(obs.comp)))
  colnames(prepare) <- c("order", "comparison", "mean", "lower", "upper", "analysis")
  rownames(prepare) <- NULL

  # Between-trial standard deviation
  prepare.tau <- data.frame(comparison, tau.meta.clean[, -2])
  colnames(prepare.tau) <- c("comparison", "median", "lower", "upper")



  ## Forest plots of comparisons on effect estimate
  p1 <- ggplot(data = prepare, aes(x = as.factor(order), y = as.numeric(mean), ymin = as.numeric(lower), ymax = as.numeric(upper), colour = analysis, group = analysis)) +
          geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
          geom_hline(yintercept = 0, lty = 2, size = 1.3, col = "grey53") +
          geom_point(size = 1.5,  colour = "black", stroke = 0.3, position = position_dodge(width = 0.5)) +
          geom_text(aes(x = as.factor(order), y = round(as.numeric(mean), 2), label = round(as.numeric(mean), 2)),
                    color = "black", hjust = -0.1, vjust = -0.5, size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.8), inherit.aes = T) +
          labs(x = "", y = "Effect estimate", colour = "Analysis") +
          scale_x_discrete(breaks = as.factor(1:length(obs.comp)), labels = prepare$comparison[1:length(obs.comp)]) +
          scale_color_manual(breaks = c("Network meta-analysis", "Paiwise meta-analysis"), values = c("#009E73", "#D55E00")) +
          coord_flip() +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                axis.title.x = element_text(color = "black", face = "bold", size = 12), legend.position = "none", legend.justification = c(0.13, 0),
                legend.text =  element_text(color = "black", size = 12), legend.title =  element_text(color = "black", face = "bold", size = 12))



  ## Forest plots of comparisons-specific between-trial standard deviation
  p2 <-  ggplot(data = prepare.tau, aes(x = as.factor(comparison), y = as.numeric(median), ymin = as.numeric(lower), ymax = as.numeric(upper))) +
           geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
           geom_hline(yintercept = tau.full[3], lty = 2, size = 1.3, col = "#006CD1") +
           geom_hline(yintercept = tau.full[5], lty = 2, size = 1.3, col = "#006CD1") +
           geom_hline(yintercept = tau.full[7], lty = 2, size = 1.3, col = "#006CD1") +
           geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
           geom_text(aes(x = as.factor(comparison), y = round(as.numeric(median), 2), label = round(as.numeric(median), 2)),
                     color = "black", hjust = -0.1, vjust = -0.5, size = 4.0, check_overlap = F, parse = F, position = position_dodge(width = 0.8), inherit.aes = T) +
           labs(x = "", y = "Between-trial standard deviation") +
           coord_flip() +
           theme_classic() +
           theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                 axis.title.x = element_text(color = "black", face = "bold", size = 12))



  ## Bring together both forest-plots
  forest.plots <- ggarrange(p1, p2, nrow = 1, ncol = 2, labels = c("A)", "B)"), common.legend = T, legend = "bottom")



  ## Write the table as .xlsx
  write_xlsx(EM.both.models, paste0(getwd(),"Table NMA vs MA.xlsx"))


  return(list(EM.both.models = EM.both.models, forest.plots = forest.plots))



}
