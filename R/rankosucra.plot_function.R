#' Rankograms and SUCRA curves
#'
#' @description It returns a panel of rankograms with integrated SUCRA curves
#'   for each intervention in the network. The function can illustrate the
#'   results of a single or two outcomes simultaneously.
#'
#' @param full1 An object of S3 class \code{\link{run_model}} for network
#'   meta-analysis. See 'Value' in \code{\link{run_model}}.
#' @param full2 An object of S3 class \code{\link{run_model}} for network
#'   meta-analysis of a second outcome. See 'Value' in \code{\link{run_model}}.
#' @param drug_names1 A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}} for \code{full1}.
#' @param drug_names2 A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}} for \code{full2}. The elements must be a subset of
#'   \code{drug_names1}.
#' @param name1 The text for the title of the results that refer to
#'   the outcome under \code{full1}.
#' @param name2 The text for the title of the results that refer to
#'   the outcome under \code{full2}.
#'
#' @return A panel of rankograms (red bars) with integrated blue SUCRA curves
#'   for each intervention in the network (Salanti et al., 2011). The x-axis of
#'   each panel refers to the ranking, and the y-axis refers to the ranking
#'   probability expressed in percentage.
#'
#' @details Interventions are sorted in the descending order of their SUCRA
#'   value. The SUCRA value expressed in percentage appears on the top left
#'   corner of each panel. In the case of two outcomes, the SUCRA values of the
#'   outcome under the argument \code{full1} are considered to sort the
#'   interventions from the best to the worst.
#'
#'   When a second outcome is also considered, different colours are used to
#'   draw the corresponding SUCRA curves and the rankograms: green for the
#'   outcome under \code{full1}, and red for the outcome under \code{full2}.
#'
#'   \code{rankosucra_plot} can be used only for a network of interventions.
#'   Otherwise, the execution of the function will be stopped and an error
#'   message will be printed on the R console.
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
#' doi: 10.1016/j.jclinepi.2010.03.016
#'
#' @examples
#' data("nma.liu2013")
#'
#' # Read results from 'run_model' (using the default arguments)
#' res <- readRDS(system.file('extdata/res_liu.rds', package = 'rnmamod'))
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "pramipexole",
#'                   "serotonin-norepinephrine reuptake inhibitor",
#'                   "serotonin reuptake inhibitor",
#'                   "tricyclic antidepressant", "pergolide")
#'
#' # Create the integrated rankograms and SUCRA curves
#' rankosucra_plot(full1 = res,
#'                 drug_names1 = interv_names)
#'
#' @export
rankosucra_plot <- function(full1,
                            full2 = NULL,
                            drug_names1,
                            drug_names2 = NULL,
                            name1 = NULL,
                            name2 = NULL) {

  if (full1$type != "nma" || is.null(full1$type)) {
    stop("'full1' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  if (!is.null(full2) & (full2$type != "nma" || is.null(full2$type))) {
    stop("'full2' must be an object of S3 class 'run_model'.",
         call. = FALSE)
  }

  # Forcing to define 'drug_names1' & 'drug_names2'
  drug_names1 <- if (missing(drug_names1)) {
    stop("The argument 'drug_names1' has not been defined.", call. = FALSE)
  } else {
    drug_names1
  }

  drug_names2 <- if (!is.null(full2) & is.null(drug_names2)) {
    stop("The argument 'drug_names2' has not been defined.", call. = FALSE)
  } else if (!is.null(full2) & !is.null(drug_names2)) {
    drug_names2
  }

  if (length(unique(is.element(drug_names2, drug_names1))) > 1) {
    stop("The argument 'drug_names2' must be a subset of 'drug_names1'.",
         call. = FALSE)
  }

  drug_names <- if (is.null(full2) ||
                     (!is.null(full2) &
                      length(drug_names1) >= length(drug_names2))) {
    drug_names1
  } else if (!is.null(full2) & length(drug_names1) < length(drug_names2)) {
    stop("'drug_names1' must have greater length than 'drug_names2'.",
         call. = FALSE)
  }

  if (length(drug_names1) < 3 || (!is.null(full2) &length(drug_names2) < 3)) {
    stop("This function is *not* relevant for a pairwise meta-analysis.",
         call. = FALSE)
  }

  name1 <- if (!is.null(full2) & is.null(name1)) {
    stop("The argument 'name1' has not been defined.", call. = FALSE)
  } else if (!is.null(full2) & !is.null(name1)) {
    name1
  }

  name2 <- if (!is.null(full2) & is.null(name2)) {
    stop("The argument 'name2' has not been defined.", call. = FALSE)
  } else if (!is.null(full2) & !is.null(name2)) {
    name2
  }

  # Prepare first outcome
  sucra <- full1$SUCRA
  effectiveness <- full1$effectiveness
  nt <- length(drug_names)

  # Order techniques according to their SUCRA value (from best to worst)
  drug_names_order <- drug_names[order(-sucra[, 1])]

  # Note: row is the drug, column is the order
  prob_rank0 <- matrix(effectiveness[, 1], nrow = nt, ncol = nt, byrow = TRUE)
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
                                                   decreasing = TRUE), 1), "%"),
                         intervention = unique(rank_data$intervention),
                         x = rep(1.2, nt),
                         y = rep(98, nt))

  # Prepare second outcome
  if (!is.null(full2)) {
    sucra2 <- full2$SUCRA
    effectiveness2 <- full2$effectiveness

    prob_rank2 <- matrix(effectiveness2[, 1],
                         nrow = length(drug_names2),
                         ncol =  length(drug_names2),
                         byrow = TRUE)
    #prob_rank2 <- prob_rank20[, order(-sucra[, 1])]
    colnames(prob_rank2) <- drug_names2

    # Obtain cumulative rank probabilities for each intervention and order
    cum_prob_rank2 <- matrix(NA,
                             nrow = length(drug_names2),
                             ncol = length(drug_names2))
    for (i in 1:length(drug_names2)) {
      cum_prob_rank2[, i] <- cumsum(prob_rank2[, i])
    }
    colnames(cum_prob_rank2) <- drug_names2

    # Merge ranking and cumulative raking probabilities
    rank_data2 <- cbind(melt(prob_rank2 * 100),
                        melt(cum_prob_rank2 * 100)[, 3])
    colnames(rank_data2) <- c("order",
                              "intervention",
                              "value_rank",
                              "value_cum")

    rank_data_fin <- cbind(rbind(rank_data, rank_data2),
                           rep(c(name1, name2),
                               c(dim(rank_data)[1], dim(rank_data2)[1])))
    colnames(rank_data_fin)[5] <- "outcome"
  }

  # Bars for the ranking probabilities and line for the SUCRA
  p <- if (is.null(full2)) {
    ggplot(rank_data, aes(x = as.factor(order), y = value_rank)) +
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
  } else {
    ggplot(rank_data_fin, aes(x = as.factor(order),
                              y = value_rank,
                              group = outcome,
                              colour = outcome,
                              fill = outcome)) +
      geom_bar(stat = "identity", position = "identity", alpha=.3) +
      #geom_line(aes(x = as.factor(order), y = value_rank), size = 1, lty = 2) +
      geom_line(aes(x = as.factor(order), y = value_cum), size = 1) +
      facet_wrap(vars(factor(intervention,
                             levels = drug_names_order[
                               seq_len(length(drug_names_order))]))) +
      scale_color_manual(name = "Outcome",
                         breaks = c(name1, name2),
                         values = c("#D55E00", "#009E73")) +
      scale_fill_manual(name = "Outcome",
                         breaks = c(name1, name2),
                         values = c("#D55E00", "#009E73")) +
      labs(x = "Rank", y = "Probability (%)") +
      theme_classic() +
      theme(axis.title.y = element_text(color = "black", size = 12,
                                        face = "bold"),
            axis.title.x = element_text(color = "black", size = 12,
                                        face = "bold"),
            axis.text.x = element_text(color = "black", size = 12),
            axis.text.y = element_text(color = "black", size = 12),
            strip.text = element_text(size = 12),
            legend.position = "bottom",
            legend.text =  element_text(color = "black", size = 12),
            #legend.title = element_text(color = "black", face = "bold",
            #                             size = 12),
            legend.title = element_blank())
  }

  return(p)
}
