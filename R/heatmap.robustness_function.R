#' Heatmap of robsutness index: Investigate the impact of missing participant outcome data
#'
#' @description This functions facilitates the detection of comparisons that are sensible to different assumptions about the missingness mechanism in the compared interventions.
#'   The heatmap is based on the robustness index of each possible pairwise comparison of the investigated network (see \code{robustness.index})
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
#' @seealso \code{\link{robustness.index}}
#'
#' @references
#' Spineli LM. A novel framework to evaluate the consistency assumption globally in a network of interventions. \emph{submitted} 2021.
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of primary analysis results: A case study on missing outcome data in pairwise and network meta-analysis. \emph{Res Synth Methods} 2021;\bold{12}(4):475--490. [\doi{10.1002/jrsm.1478}]
#'
#' Kullback S, Leibler RA. On information and sufficiency. \emph{Ann Math Stat} 1951;\bold{22}(1):79--86. [\doi{10.1214/aoms/1177729694}]
#'
#' @export
heatmap.robustness <- function(robust, drug.names){


  if (is.na(robust)) {
    stop("Missing participant outcome data have *not* been collected. This function cannot be used.", call. = F)
  }


  RI <- robust$RI; threshold <- robust$threshold

  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The  argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- (1 + sqrt(1 + 8*length(robust$robust)))/2
    as.character(1:nt)
  } else {
    drug.names
  }


  if(length(drug.names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }


  if (missing(threshold) & is.element(robust$measure, "OR")) {
    threshold <- 0.28
    #message("The value 0.28 was assigned on 'threshold' by default")
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The value 0.28 was assigned on 'threshold' by default", "\033[0m", "\n")))
  } else if (missing(threshold) & is.element(robust$measure, c("MD", "SMD", "ROM"))) {
    threshold <- 0.17
    #message("The value 0.17 was assigned on 'threshold' by default")
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The value 0.17 was assigned on 'threshold' by default", "\033[0m", "\n")))
  } else {
    threshold <- threshold
    #message(paste("The value", threshold, "was assigned on 'threshold' for", effect.measure.name(full$measure)))
    message(cat(paste0("\033[0;", col = 32, "m", txt = paste("The value", threshold, "was assigned on 'threshold' for", effect.measure.name(robust$measure)), "\033[0m", "\n")))
  }


  ## Lower triangular heatmap matrix - Comparisons are read from the left to the right
  ## CAREFUL: The interventions in the drug.names should follow the order you considered to run NMA pattern-mixture model!
  mat <- matrix(NA, nrow = length(drug.names) - 1, ncol = length(drug.names) - 1)
  mat[lower.tri(mat, diag = T)] <- round(RI, 2)
  colnames(mat) <- drug.names[1:(length(drug.names) - 1)]; rownames(mat) <- drug.names[2:length(drug.names)]
  mat.new <- melt(mat, na.rm = T)


  ## Create the heatmap for one network of interventions
  p <- ggplot(mat.new, aes(factor(Var2, level = drug.names[1:(length(drug.names) - 1)]), factor(Var1, level = drug.names[length(drug.names):2]), fill = ifelse(value < threshold, "high", "poor"))) +
         geom_tile(colour = "white") +
         geom_text(aes(factor(Var2, level = drug.names[1:(length(drug.names) - 1)]), factor(Var1, level = drug.names[length(drug.names):2]), label = value, fontface = "bold"), size = rel(4.5)) +
         scale_fill_manual(breaks = c("high", "poor"), values = c("#009E73", "#D55E00")) +
         scale_x_discrete(position = "top") +
         labs(x = "", y = "") +
         theme(legend.position = "none", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12))

  return(p)
}

