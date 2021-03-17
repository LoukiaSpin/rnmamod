#' Integrated Rankograms and SUCRA Curves
#'
#' @param net An object of S3 class \code{nma.continuous.full.model}.
#' @param drug.names A vector of characters with the names of the interventions in the order they appear in the function \code{nma.continuous.full.model}.
#'
#' @return Integrated rankograms and SUCRA curves for each intervention sorted in descending order by the SUCRA value.
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: an overview and tutorial. J Clin Epidemiol. 2011;64(2):163-71.
#'
#' \dontshow{load("netmodr/data/NMA_results.RData")}
#' @examples
#' drug.names <- sapply(1:14, function(x) letters[x])
#' rankogram.sucra.plots(net = res, drug.names = drug.names)
#'
#' @export
rankosucra.plot <- function(net, drug.names){

  SUCRA <- net$SUCRA; effectiveness <- net$effectiveness

  # Order techniques according to their SUCRA value (from best to worst)
  drug.names.order <- drug.names[order(-SUCRA[, 1])]

  nt <- length(drug.names)


  # Note: row is the drug, column is the order
  prob.rank0 <- matrix(effectiveness[, 1], nrow = nt, ncol = nt, byrow = T)
  prob.rank <- prob.rank0[, order(-SUCRA[, 1])]
  colnames(prob.rank) <- drug.names.order


  # Obtain cumulative rank probabilities for each intervention and order
  cum.prob.rank <- matrix(NA, nrow = nt, ncol = nt)
  for(i in 1:nt){

    cum.prob.rank[, i] <- cumsum(prob.rank[, i])
  }
  colnames(cum.prob.rank) <- drug.names.order


  # Merge ranking and cumulative raking probabilities (to proceed with ggplot2)
  rank.data <- cbind(melt(prob.rank*100), melt(cum.prob.rank*100)[, 3])
  colnames(rank.data) <- c("Order", "Intervention", "value.rank", "value.cum")

  dat_text <- data.frame(label = paste0(round(sort(SUCRA[, 1]*100, decreasing = T), 1), "%"), Intervention = unique(rank.data$Intervention), x = rep(1.2, nt), y = rep(98, nt))

  # Bars for the ranking probabilities and line for the SUCRA
  p <- ggplot(rank.data, aes(x = as.factor(Order), y = value.rank)) +
         geom_bar(stat = "identity", color = "brown2", fill = "brown2") +
         geom_line(aes(x = Order, y = value.cum), size = 1, color = "blue") +
         facet_wrap(vars(Intervention)) +
         geom_text(data = dat_text, aes(x = x, y = y, label = label), fontface = "bold") +
         labs(x = "Rank", y = "Probability (%)") +
         theme_classic() +
         theme(axis.title.y = element_text(color = "black", size = 12, face = "bold"), axis.title.x = element_text(color = "black", size = 12, face = "bold"),
               axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
               strip.text = element_text(size = 12))
  return(p)

}
