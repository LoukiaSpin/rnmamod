#' Integrated rankograms and SUCRA curves
#'
#' @description \code{rankosucra.plot} returns a panel of rankograms with integrated SUCRA curves for each intervention of the network.
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @return A panel of rankograms with integrated SUCRA curves for each intervention of the network (Salanti et al., 2011). The x-axis of each panel refers to the ranking, and the y-axis refers to the ranking probability in percent.
#'
#' @details Interventions are sorted in the descending order of their SUCRA value.
#'   The SUCRA value in percent appears on the top left corner of each panel.
#'
#'   \code{rankosucra.plot} can be used only for a network of interventions. In the case of two interventions, the execution of the function will be stopped and an error message will be printed in the R console.
#'
#' @author {Loukia M. Spineli}, {Chrysostomos Kalyvas}, {Katerina Papadimitropoulou}
#'
#' @seealso \code{\link{run.model}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: an overview and tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71. [\doi{10.1016/j.jclinepi.2010.03.016}]
#'
#' @examples
#' data("nma.baker2009")
#'
#' # Perform a random-effects network meta-analysis
#' res1 <- run.model(data = nma.baker2009,
#'                   measure = "OR",
#'                   model = "RE",
#'                   assumption = "IDE-ARM",
#'                   heter.prior = list("halfnormal", 0, 1),
#'                   mean.misspar = 0,
#'                   var.misspar = 1,
#'                   D = 1,
#'                   n.chains = 3,
#'                   n.iter = 10000,
#'                   n.burnin = 1000,
#'                   n.thin = 1)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("budesodine", "budesodine plus formoterol", "fluticasone", "fluticasone plus
#'                   salmeterol", "formoterol", "salmeterol", "tiotropium", "placebo")
#'
#' # Create the league heatmap
#' rankosucra.plot(full = res1, drug.names = interv.names)
#'
#' @export
rankosucra.plot <- function(full, drug.names){


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

  SUCRA <- full$SUCRA; effectiveness <- full$effectiveness; nt <- length(drug.names)

  # Order techniques according to their SUCRA value (from best to worst)
  drug.names.order <- drug.names[order(-SUCRA[, 1])]

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
         geom_bar(stat = "identity", color = "#D55E00", fill = "#D55E00") +
         geom_line(aes(x = Order, y = value.cum), size = 1, color = "blue") +
         #facet_wrap(vars(Intervention)) +
         facet_wrap(vars(factor(Intervention, levels = drug.names.order[1:length(drug.names.order)]))) +
         geom_text(data = dat_text, aes(x = x, y = y, label = label), fontface = "bold") +
         labs(x = "Rank", y = "Probability (%)") +
         theme_classic() +
         theme(axis.title.y = element_text(color = "black", size = 12, face = "bold"), axis.title.x = element_text(color = "black", size = 12, face = "bold"),
               axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
               strip.text = element_text(size = 12))
  return(p)

}
