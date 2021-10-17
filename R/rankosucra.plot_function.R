#' Rankograms and SUCRA curves
#'
#' @description Returns a panel of rankograms with integrated SUCRA curves for
#'   each intervention in the network.
#'
#' @param full An object of S3 class \code{\link{run_model}}. See 'Value' in
#'   \code{\link{run_model}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If the argument \code{drug_names} is not defined,
#'   the interventions are ordered as they appear in \code{data}.
#'
#' @return A panel of rankograms with integrated SUCRA curves for each
#'   intervention of the network (Salanti et al., 2011). The x-axis of each
#'   panel refers to the ranking, and the y-axis refers to the ranking
#'   probability expressed in percentage.
#'
#' @details Interventions are sorted in the descending order of their SUCRA
#'   value. The SUCRA value expressed in percentage appears on the top left
#'   corner of each panel.
#'
#'   \code{rankosucra_plot} can be used only for a network of interventions.
#'
#' @author {Loukia M. Spineli}, {Chrysostomos Kalyvas},
#'   {Katerina Papadimitropoulou}
#'
#' @seealso \code{\link{run_model}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries
#' for presenting results from multiple-treatment meta-analysis: an overview and
#' tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71.
#' [\doi{10.1016/j.jclinepi.2010.03.016}]
#'
#' @examples
#' data("nma.liu2013")
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' res <- run_model(data = nma.liu2013,
#'                  measure = "OR",
#'                  model = "RE",
#'                  assumption = "IDE-ARM",
#'                  heter_prior = list("halfnormal", 0, 1),
#'                  mean_misspar = c(0, 0),
#'                  var_misspar = 1,
#'                  D = 1,
#'                  n_chains = 3,
#'                  n_iter = 10000,
#'                  n_burnin = 1000,
#'                  n_thin = 1)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "pramipexole", "serotonin-norepinephrine
#'                   reuptake inhibitor", "serotonin reuptake inhibitor",
#'                   "tricyclic antidepressant", "pergolide")
#'
#' # Create the integrated rankograms and SUCRA curves
#' rankosucra_plot(full = res,
#'                 drug_names = interv_names)
#' }
#' @export
rankosucra_plot <- function(full, drug_names) {

  drug_names <- if (missing(drug_names)) {
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt = "The argument 'drug_names' has not been defined.
                       The intervention ID, as specified in 'data' is used as
                       intervention names", "\033[0m", "\n")))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug_names
  }

  if (length(drug_names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis",
         call. = F)
  }

  sucra <- full$SUCRA
  effectiveness <- full$effectiveness
  nt <- length(drug_names)

  # Order techniques according to their SUCRA value (from best to worst)
  drug_names_order <- drug_names[order(-sucra[, 1])]

  # Note: row is the drug, column is the order
  prob_rank0 <- matrix(effectiveness[, 1], nrow = nt, ncol = nt, byrow = T)
  prob_rank <- prob_rank0[, order(-sucra[, 1])]
  colnames(prob_rank) <- drug_names_order

  # Obtain cumulative rank probabilities for each intervention and order
  cum_prob_rank <- matrix(NA, nrow = nt, ncol = nt)
  for (i in 1:nt) {

    cum_prob_rank[, i] <- cumsum(prob_rank[, i])
  }
  colnames(cum_prob_rank) <- drug_names_order

  # Merge ranking and cumulative raking probabilities (to proceed with ggplot2)
  rank_data <- cbind(melt(prob_rank * 100), melt(cum_prob_rank * 100)[, 3])
  colnames(rank_data) <- c("order", "intervention", "value_rank", "value_cum")
  dat_text <- data.frame(label = paste0(round(sort(sucra[, 1] * 100,
                                                   decreasing = T), 1), "%"),
                         intervention = unique(rank_data$intervention),
                         x = rep(1.2, nt),
                         y = rep(98, nt))

  # Bars for the ranking probabilities and line for the SUCRA
  p <- ggplot(rank_data, aes(x = as.factor(order), y = value_rank)) +
         geom_bar(stat = "identity", color = "#D55E00", fill = "#D55E00") +
         geom_line(aes(x = order, y = value_cum), size = 1, color = "blue") +
         facet_wrap(vars(factor(intervention,
                                levels = drug_names_order[
                                  seq_len(length(drug_names_order))]))) +
         geom_text(data = dat_text,
                   aes(x = x, y = y, label = label),
                   fontface = "bold") +
         labs(x = "Rank", y = "Probability (%)") +
         theme_classic() +
         theme(axis.title.y = element_text(color = "black", size = 12,
                                           face = "bold"),
               axis.title.x = element_text(color = "black", size = 12,
                                           face = "bold"),
               axis.text.x = element_text(color = "black", size = 12),
               axis.text.y = element_text(color = "black", size = 12),
               strip.text = element_text(size = 12))
  return(p)
}
