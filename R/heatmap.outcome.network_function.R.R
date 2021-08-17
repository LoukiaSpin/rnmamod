#' Heatmap of the distribution of the outcome in the network
#'
#' @description This function illustrates the distribution of the investigated outcome for each intervention and observed comparison of the network.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models.
#'   See 'Format' in \code{\link[rnmamod]{run.model}} function for the specification of the columns.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in the argument \code{data} is used, instead.
#'
#' @return A heatmap on the summarised outcome in each intervention and observed comparison of the network.
#'   Each cell presents the median, minimum and maximum (in parenthesis) outcome across the corresponding trials.
#'   When the outcome is binary, first, the proportion of observed events in each trial-arm of the network is calculated to obtain the aforementioned summary statistics for each intervention.
#'   Then, the total proportion of observed events in each pairwise comparison of each trial is calculated to obtain the aforementioned summary statistics for each observed comparison of the network.
#'
#'   In the case of a continuous outcome, first, the t-statistic in each trial-arm of the network is calculated to obtain the aforementioned summary statistics for each intervention.
#'   The trial-arm t-statistic is the ratio of the corresponding extracted mean outcome, \code{y}, to the standard error, \code{se} (see the \code{data.preparation} function).
#'   Then, the t-statistic in each pairwise comparison of each trial is calculated to obtain the aforementioned summary statistics for each observed comparison of the network.
#'   For a pairwise comparison of a trial, the t-statistic refers to the standardised mean difference defined as the ratio of mean difference to the standard error of the mean difference.
#'
#'   The summary statistics in each intervention and observed comparison are depicted using white and black colour in the main diagonal and lower off-diagonal, respectively, of the heatmap.
#'   The pairwise comparisons are read from left to right. The function is redundant for a pairwise meta-analysis.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link[rnmamod]{data.preparation}}.
#'
#' @examples
#' data("nma.stowe2011")
#'
#' # Return the first six trials of the dataset
#' head(nma.stowe2011)
#' #              study t1 t2    y1    y2  sd1  sd2 m1 m2  n1  n2
#' #   DA (B): Interntl  1  2 -0.30 -1.20 4.36 4.32  7  3  83  84
#' #      DA (C): Spain  1  2 -2.47 -3.33 3.91 3.48  8  9  20  23
#' #         DA (C): UK  1  2 -0.70 -2.00 2.24 2.33  2  1  18  19
#' #      DA (C): USA 1  1  2 -0.77 -2.08 3.32 3.21 19 34  65 123
#' # DA (Pe): N America  1  2 -0.20 -1.80 4.79 4.81  0  0 187 189
#' # DA (Pr): CLEOPATRA  1  2 -0.90 -2.80 5.00 2.83  1  1 101 201
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("PBO+LD", "DA+LD", "COMTI+LD", "MAOBI+LD")
#'
#' # Create the heatmap
#' heatmap.outcome.network(data = nma.stowe2011, drug.names = interv.names)
#'
#' @export
heatmap.outcome.network <- function(data, drug.names){

  options(warn = -1)

  if(dim(data %>% dplyr::select(starts_with("r")))[2] > 0) {
    measure <- "OR"
  } else {
    measure <- "MD"
  }

  # Use the 'data.preparation' function
  dat <- data.preparation(data, measure)

  # Condition when 'drug.names' is not defined
  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in argument 'data' is used as intervention names", "\033[0m", "\n")))
    as.character(1:dat$nt)
  } else {
    drug.names
  }

  if (measure == "OR") {
    # Proportion of observed events per trial-arm
    arm.risk <- dat$r/(dat$N - dat$m)

    # Minimum proportion of observed events per intervention in % (across the corresponding trials)
    min.out.interv <- round(aggregate(unlist(arm.risk), by = list(unlist(dat$t)), min)[, 2], 2)*100

    # Median proportion of observed events per intervention in % (across the corresponding trials)
    median.out.interv <- round(aggregate(unlist(arm.risk), by = list(unlist(dat$t)), median)[, 2], 2)*100

    # Max proportion of observed events per intervention in % (across the corresponding trials)
    max.out.interv <- round(aggregate(unlist(arm.risk), by = list(unlist(dat$t)), max)[, 2], 2)*100

    # Turn into long format using the 'pairwise' function (netmeta): binary outcome
    pair.out <- pairwise(as.list(dat$t),
                         event = as.list(dat$r),
                         n = as.list(dat$N),
                         data = cbind(dat$t, dat$r, dat$N),
                         studlab = 1:dat$ns)[, c(3:6, 8, 7, 9)]
    colnames(pair.out) <- c("study", "t1", "t2", "r1", "r2", "n1", "n2")

    # The comparison between the second and first arms of each trial
    comp <- paste(pair.out[, "t2"], "vs", pair.out[, "t1"])

    # Turn into long format using the 'pairwise' function (netmeta): missing participant outcome data
    pair.mod <- pairwise(as.list(dat$t),
                         event = as.list(dat$m),
                         n = as.list(dat$N),
                         data = cbind(dat$t, dat$m, dat$N),
                         studlab = 1:dat$ns)[, c(3:6, 8, 7, 9)]
    colnames(pair.mod) <- c("study", "t1", "t2", "m1", "m2", "n1", "n2")

    # Proportion of observed events per trial-comparison
    trial.risk <- apply(pair.out[, c("r1", "r2")], 1, sum)/(apply(pair.out[, c("n1", "n2")], 1, sum) - apply(pair.mod[, c("m1", "m2")], 1, sum))

    # Minimum proportion of observed events per observed comparison in % (across the corresponding trials)
    min.out.comp <- round(aggregate(trial.risk, by = list(comp), min)[, 2], 2)*100

    # Median proportion of observed events per observed comparison in % (across the corresponding trials)
    median.out.comp <- round(aggregate(trial.risk, by = list(comp), median)[, 2], 2)*100

    # Maximum proportion of observed events per observed comparison in % (across the corresponding trials)
    max.out.comp <- round(aggregate(trial.risk, by = list(comp), max)[, 2], 2)*100
  } else {
    # Minimum t-statistic per intervention (across the corresponding trials)
    min.out.interv <- round(aggregate(unlist(dat$y0)/unlist(dat$se0), by = list(unlist(dat$t)), min)[, 2], 2)

    # Median t-statistic per intervention (across the corresponding trials)
    median.out.interv <- round(aggregate(unlist(dat$y0)/unlist(dat$se0), by = list(unlist(dat$t)), median)[, 2], 2)

    # Maximum t-statistic per intervention (across the corresponding trials)
    max.out.interv <- round(aggregate(unlist(dat$y0)/unlist(dat$se0), by = list(unlist(dat$t)), max)[, 2], 2)

    # Turn into long format using the 'pairwise' function (netmeta): continuous outcome
    pair.out <- pairwise(as.list(dat$t),
                         n = as.list(dat$N),
                         mean = as.list(dat$y0),
                         sd = as.list(dat$se0),
                         data = cbind(dat$t, dat$N, dat$y0, dat$se0),
                         studlab = 1:dat$ns)[, c(3:5, 7, 10, 8, 11, 6, 9)]
    colnames(pair.out) <- c("study", "t1", "t2", "y1", "y2", "se1", "se2","n1", "n2")

    # The comparison between the second and first arms of each trial
    comp <- paste(pair.out[, "t2"], "vs", pair.out[, "t1"])

    # t-statistic (t2 versus t1) per trial-comparison
    t.stat.trial <- (pair.out$y2 - pair.out$y1)/sqrt((pair.out$se2)^2 + (pair.out$se1)^2)

    # Minimum t-statistic per observed comparison (across the corresponding trials)
    min.out.comp <- round(aggregate(t.stat.trial, by = list(comp), min)[, 2], 2)

    # Median t-statistic per observed comparison (across the corresponding trials)
    median.out.comp <- round(aggregate(t.stat.trial, by = list(comp), median)[, 2], 2)

    # Maximum t-statistic per observed comparison (across the corresponding trials)
    max.out.comp <- round(aggregate(t.stat.trial, by = list(comp), max)[, 2], 2)
  }

  # Lower triangular heatmap matrix - Comparisons are read from the left to the right
  median.mat <- min.mat <- max.mat <- matrix(NA, nrow = length(drug.names), ncol = length(drug.names))
  median.mat[cbind(pair.out$t2, pair.out$t1)] <- median.out.comp
  min.mat[cbind(pair.out$t2, pair.out$t1)] <- min.out.comp
  max.mat[cbind(pair.out$t2, pair.out$t1)] <- max.out.comp

  # 'Paste' the median with the minimum and maximum (the latter two in a parenthesis)
  final <- matrix(paste0(median.mat,  "\n", "(", min.mat, ",", " ", max.mat, ")"), nrow = length(drug.names), ncol = length(drug.names))
  diag(final) <- paste0(median.out.interv,  "\n", "(", min.out.interv, ",", " ", max.out.interv, ")")
  colnames(final) <- rownames(final) <- drug.names

  # Preparing the dataset for the ggplot2
  mat <- median.mat
  diag(mat) <- median.out.interv
  mat.new <- cbind(melt(final, na.rm = T), melt(mat, na.rm = F)[, 3])
  colnames(mat.new) <- c("Var1", "Var2", "value", "value2")

  # Keep only the rows without 'NA' and use for ggplot2
  mat.final <- mat.new[complete.cases(mat.new), ]
  mat.final[, 5] <- ifelse(mat.final$Var1 == mat.final$Var2, 1, 0)
  colnames(mat.final) <- c("Var1", "Var2", "value", "value2", "value.interv")

  # Create the heatmap for one network of interventions
  ggplot(mat.final, aes(as.factor(Var2), factor(Var1, levels = drug.names[length(drug.names):1]), fill = value2)) +
    geom_tile(colour = "white") +
    geom_text(aes(Var2, Var1, label = value, fontface = "bold"), colour = ifelse(mat.final$value.interv < 1, "black", "white"), size = rel(4.5)) +
    scale_fill_gradient(low = "#009E73", high = "#D55E00", limits = c(ifelse(measure == "OR", 0, min(mat.final$value2)),ifelse(measure == "OR", 100, max(mat.final$value2)) )) +
    scale_x_discrete(position = "top") +
    labs(x = "", y = "", caption = paste("Summary of", ifelse(measure == "OR", "observed events (%)", "t-statistic"), ": median (minimum, maximum)"),
         fill = "Median value") +
    theme_classic() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12))

}
