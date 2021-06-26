#' Heatmap of robutness for NMA
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

