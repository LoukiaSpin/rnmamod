#' Heatmap of Kullback-Leibler divergence values: consistency model versus unrelated mean effects model
#'
#' @description This functions facilitates the detection of observed pairwise comparisons that are sensible
#'   to applying the consistency model or the unrelated mean effects model. The heatmap is based on the Kullback-Leibler divergence (KLD) measure
#'   (Kullback and Leibler, 1951) in the summary effect size from the unrelated mean effects model to the consistency model (see \code{similarity.index.UME}).
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param ume An object of S3 class \code{\link{run.UME}}. See 'Value' in \code{\link{run.UME}}.
#' @param threshold A number indicating the threshold of similarity. See 'Details' below.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @details \code{heatmap.similarity.UME} is integrated in the \code{UME.plot} function. The heatmap illustrates the KLD values for the observed pairwise comparisons in the network.
#'   The observed pairwise comparisons are read from left to right and are highlighted either with a green or red colour.
#'   Comparisons highlighted with green or red colour imply high or poor similarity of the posterior distribution of the corresponding summary effect size in the compared models.
#'   Cells with KLD values in white correspond to frail comparisons as detected by the \code{run.UME} function (see 'Details' in \code{run.UME}).
#'
#'   The user may consider the values 0.28 and 0.17 as \code{threshold} for binary and continuous outcome data, respectively.
#'   These thresholds have been originally developed by Spineli et al. (2021) and considered also by Spineli (2021) in the proposed
#'   framework of global evaluation of the consistency assumption.
#'
#' @return A lower triangular heatmap matrix on the KLD measure in the summary effect size from the unrelated mean effects model
#'   to the consistency model (See, 'Description' in \code{similarity.index.UME}).
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{run.UME}}, \code{\link{UME.plot}}, \code{\link{similarity.index.UME}}
#'
#' @references
#' Spineli LM. A novel framework to evaluate the consistency assumption globally in a network of interventions. \emph{submitted} 2021.
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of primary analysis results: A case study on missing outcome data in pairwise and network meta-analysis. \emph{Res Synth Methods} 2021;\bold{12}(4):475--490. [\doi{10.1002/jrsm.1478}]
#'
#' Kullback S, Leibler RA. On information and sufficiency. \emph{Ann Math Stat} 1951;\bold{22}(1):79--86. [\doi{10.1214/aoms/1177729694}]
#'
#' @export
heatmap.similarity.UME <- function(full, ume, drug.names, threshold){


  options(warn = -1)

  ## Possible and observed comparisons (with names)
  possible.comp <- possible.observed.frail.comparisons(drug.names, ume$obs.comp)


  ## Assign the robustness index to the observed comparisons
  EM.full <- full$EM
  ume.post <- ume$EM
  robust <- similarity.index.UME(EM.full[is.element(possible.comp$poss.comp[, 4], possible.comp$obs.comp[, 3]), 1:2], ume.post[, 1:2], threshold)$KLD

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
  p <- ggplot(mat.new, aes(factor(Var2, level = drug.names[1:(length(drug.names) - 1)]), factor(Var1, level = drug.names[length(drug.names):2]), fill = ifelse(value < threshold, "high", "poor"))) +
         geom_tile(colour = "white") +
         geom_text(aes(factor(Var2, level = drug.names[1:(length(drug.names) - 1)]), factor(Var1, level = drug.names[length(drug.names):2]), label = value, fontface = "bold"), colour = ifelse(mat.new2$value < 1, "black", "white"), size = rel(4.5)) +
         scale_fill_manual(breaks = c("high", "poor"), values = c("#009E73", "#D55E00"), na.value = "white") +
         scale_x_discrete(position = "top") +
         labs(x = "", y = "", fill = "Similarity") +
         theme_bw() +
         theme(axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
               legend.position = "bottom", legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12))

  return(p)
}

