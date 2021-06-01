#' The heatmap of Kullback-Leibler divergence values of the observed comparisons under the network meta-analysis (NMA) and refined unrelated mean effects (UME) models
#'
#' @param robust The component \code{"robust"} of class \code{robustness.index}.
#' @param drug.names A vector of characteristics with name of the interventions as appear in the functions \code{run.model}, \code{run.UME}, or \code{run.UME.Dias}.
#' @param obs.comp A vector of the observed comparisons of interventions in the network obtained vias the \code{run.UME} function.
#' @param frail.comp A vector of the omitted comparisons of interventions among the observed comparisons in the network obtained vias the \code{run.UME} function.
#' @param threshold A number indicating the threshold of similarity.
#'
#' @return A heatmap on the Kullback-Leibler divergence (KLD) measure from the refined UME model to the NMA model for each observed comparison in the network.
#' The observed pairwise comparisons are read from left to right and are highlighted either with a green or red colour. Green and red colours imply high
#' and poor similarity of the posterior distribution of the treatment effects in the compared models. KLD values in white correspond to omitted comparisons.
similarity.index <- function(full, ume, threshold){


  ## Function for the Kullback-Leibler Divergence (comparing two univariate normal distributions)
  KLD.measure.univ <- function(mean.y, sd.y, mean.x, sd.x){

    # x is the 'truth' (e.g. the MAR assumption)
    KLD.xy <- 0.5*(((sd.x/sd.y)^2) + ((mean.y - mean.x)^2)/(sd.y^2) - 1 + 2*log(sd.y/sd.x))

    return(KLD.xy = KLD.xy)
  }


  kldxy <- rep(NA, length(full[, 1]))

  for(i in 1:length(full[, 1])){ ## We are interested in all observed comparisons of the network

    ## Returns the KLD of UME when compared with NMA (NMA as 'true') for comparison i
    kldxy[i] <- KLD.measure.univ(ume[i, 1], ume[i, 2], full[i, 1], full[i, 2])

  }


  robust <- ifelse(kldxy < threshold, "robust", "frail")

  return(list(robust = robust, KLD = kldxy))
}


heatmap.similarity.UME <- function(full, ume, drug.names, threshold){


  options(warn = -1)

  ## Possible and observed comparisons (with names)
  possible.comp <- possible.observed.frail.comparisons(drug.names, ume$obs.comp)


  ## Assign the robustness index to the observed comparisons
  EM.full <- full$EM
  ume.post <- ume$EM
  robust <- similarity.index(EM.full[is.element(possible.comp$poss.comp[, 4], possible.comp$obs.comp[, 3]), 1:2], ume.post[, 1:2], threshold)$KLD

  poss.comp <- cbind(possible.comp$poss.comp, rep(NA, length(possible.comp$poss.comp[, 1])))
  poss.comp[is.element(possible.comp$poss.comp[, 5], possible.comp$obs.comp[, 4]) == T, 6] <- robust
  colnames(poss.comp) <- c("ID", "treat1", "treat2", "comp", "comp.name", "RI")


  ## Create a dummy vector on whether a comparison is frail
  frail <- ifelse(is.element(possible.comp$poss.comp[, 4], ume$frail.comp), 1, 0)


  ## Lower triangular heatmap matrix - Comparisons are read from the left to the right
  ## CAREFUL: The interventions in the drug.names should follow the order you considered to run NMA pattern-mixture model!
  mat <- matrix(NA, nrow = length(drug.names) - 1, ncol = length(drug.names) - 1)
  mat[lower.tri(mat, diag = T)] <- round(poss.comp$RI, 2)
  colnames(mat) <- drug.names[1:(length(drug.names) - 1)]; rownames(mat) <- drug.names[2:length(drug.names)]
  mat.new <- melt(mat, na.rm = F)



  ## Create another 'mat' that indicate the frail comparisons with '1'
  mat2 <- matrix(NA, nrow = length(drug.names) - 1, ncol = length(drug.names) - 1)
  mat2[lower.tri(mat2, diag = T)] <- frail
  colnames(mat2) <- drug.names[1:(length(drug.names) - 1)]; rownames(mat2) <- drug.names[2:length(drug.names)]
  mat.new2 <- melt(mat2, na.rm = F)



  ## Create the heatmap for one network of interventions
  p <- ggplot(mat.new, aes(Var2, factor(Var1, level = drug.names[length(drug.names):2]), fill = ifelse(value < threshold, "high", "poor"))) +
         geom_tile(colour = "white") +
         geom_text(aes(Var2, Var1, label = value, fontface = "bold"), colour = ifelse(mat.new2$value < 1, "black", "white"), size = rel(4.5)) +
         scale_fill_manual(breaks = c("high", "poor"), values = c("#009E73", "#D55E00"), na.value = "white") +
         scale_x_discrete(position = "top") +
         labs(x = "", y = "", fill = "Similarity") +
         theme_bw() +
         theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
               legend.position = "bottom", legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12))

  return(p)
}

