#' Plot the results from the node-splitting approach
#'
#' @export
nodesplit.plot <- function(net, drug.names, tau.values) {


  ## Keep results on 'direct evidence', 'indirect evidence', 'inconsistency factor', 'between-trial standard deviation',
  ## and model assessment measures (i.e., DIC, posterior mean of refisual deviance, and pD)
  direct0 <-net$direct; indirect0 <-net$EM; IF0 <-net$diff; tau0 <- net$tau; model.assess <- net$model.assessment


  ## Sort 'direct evidence', 'indirect evidence', 'inconsistency factor', and 'between-trial standard deviation' by DIC in ascending order
  direct <- direct0[order(model.assess$DIC), ]
  indirect <- indirect0[order(model.assess$DIC), ]
  IF <- IF0[order(model.assess$DIC), ]
  tau <- tau0[order(model.assess$DIC), ]


  ## Interventions' name: Replace code with original names
  # For treat1 (non-baseline arm)
  for(i in sort(unique(unlist(direct[, 1])))) {
    direct[direct$treat1 == i, 1] <- drug.names[i]
  }

  # For treat2 (baseline arm)
  for(i in sort(unique(unlist(direct[, 2])))) {
    direct[direct$treat2 == i, 2] <- drug.names[i]
  }


  ## Prepare the dataset to create the panel of forest-plot on the 'direct evidence', 'indirect evidence', and 'inconsistency factor' for each split node
  comp <- paste(direct[, 1], "vs", direct[, 2])
  prepare <- data.frame(rep(comp, 3), rbind(direct[, c(3, 5:6)], indirect[, c(3, 5:6)], IF[, c(3, 5:6)]), rep(c("direct", "indirect", "IF"), each = length(direct[, 1])))
  colnames(prepare) <- c("node", "mean", "lower", "upper", "evidence")
  prepare$stat.signif <- ifelse(prepare$lower > 0 | prepare$upper < 0  , "statistically significant", "statistically non-significant")
  prepare$stat.signif <- ifelse(prepare$evidence != "IF", NA, prepare$stat.signif)
  prepare$DIC <- sort(model.assess$DIC)


  ## Create the panel of forest-plots
  if(length(unique(comp)) <= 16) {

    p1 <- ggplot(data = prepare, aes(x = factor(evidence, levels = c("IF", "indirect", "direct")), y = mean, ymin = lower, ymax = upper, colour = stat.signif) ) +
            geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
            geom_hline(yintercept = 0, lty = 2, size = 1.5, col = "grey") +
            geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
            geom_text(aes(x = as.factor(evidence), y = round(mean, 2), label = round(mean, 2)), color = "black", hjust = -0.2, vjust = -0.3, size = 4.0,
                      check_overlap = F, parse = F, position = position_dodge(width = 0.8),  inherit.aes = T) +
            geom_label(aes(x = 3.5, y = -Inf, hjust = 0, vjust = 1, label = round(DIC, 0)), fill = "beige", colour = "black", fontface = "plain", size = 3.1) +
            facet_wrap(vars(factor(node, levels = unique(prepare$node))), scales = "free_x") +
            coord_flip() +
            labs(x = "", y = "", colour = "") +
            scale_color_manual(breaks = c("statistically significant", "statistically non-significant"), values = c("green4", "red"), na.value = "black") +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), legend.position = "none",
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
            geom_text(aes(x = as.factor(evidence), y = round(mean, 2), label = round(mean, 2)), color = "black", hjust = -0.2, vjust = -0.3, size = 4.0,
                      check_overlap = F, parse = F, position = position_dodge(width = 0.8),  inherit.aes = T) +
            geom_label(aes(x = 3.5, y = -Inf, hjust = 0, vjust = 1, label = round(DIC, 0)), fill = "beige", colour = "black", fontface = "plain", size = 3.1) +
            facet_wrap(vars(factor(node, levels = unique(prepare$node))), scales = "free_x") +
            coord_flip() +
            labs(x = "", y = "", colour = "") +
            scale_color_manual(breaks = c("statistically significant", "statistically non-significant"), values = c("green4", "red"), na.value = "black") +
            theme_classic() +
            theme(axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), legend.position = "none",
                  strip.text = element_text(size = 11))

  }



  ## Prepare the dataset to create the panel of forest-plot on the 'between-trial standard deviation' for each split node
  prepare.tau <- data.frame(comp, tau[, c(3, 5:6)], sort(model.assess$DIC))
  colnames(prepare.tau) <- c("node", "median", "lower", "upper", "DIC")


  ## Create the forest-plot for 'between-trial standard deviation'
  p2 <- ggplot(data = prepare.tau, aes(x = factor(node, levels = unique(prepare$node)), y = median, ymin = lower, ymax = upper) ) +
          geom_linerange(size = 2, position = position_dodge(width = 0.5)) +
          geom_hline(yintercept = tau.values[1], lty = 2, size = 1.2, col = "red") +
          geom_hline(yintercept = tau.values[2], lty = 2, size = 1.2, col = "red") +
          geom_hline(yintercept = tau.values[3], lty = 2, size = 1.2, col = "red") +
          geom_point(size = 1.5,  colour = "white", stroke = 0.3, position = position_dodge(width = 0.5)) +
          geom_text(aes(x = factor(node, levels = unique(prepare$node)), y = round(median, 2), label = round(median, 2)), color = "black", hjust = -0.2, vjust = -0.3, size = 4.0,
                    check_overlap = F, parse = F, position = position_dodge(width = 0.8),  inherit.aes = T) +
          geom_label(aes(x = factor(node, levels = unique(prepare$node)), y = upper, label = round(DIC, 0)), fill = "beige", colour = "black", fontface = "plain",  size = 3.1) +
          labs(x = "Split nodes (sorted by DIC in ascending order)", y = "Between-trial standard deviation") +
          scale_color_manual(breaks = c("statistically significant", "statistically non-significant"), values = c("green4", "red"), na.value = "black") +
          theme_classic() +
          theme(axis.text.x = element_text(color = "black", size = 12, angle = 45, hjust = 1), axis.text.y = element_text(color = "black", size = 12), legend.position = "none")


  return(list(node.split.forestplot = p1, tau.forestplot = p2))

}
