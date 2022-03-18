#' Heatmap of proportion of missing participants in the network
#'
#' @description Illustrates the distribution of missing participants and the
#'   associated risk of bias for each intervention and observed comparison in
#'   the network.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level
#'   data of each trial. See 'Format' in \code{\link{run_model}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data}. If the argument
#'   \code{drug_names} is not defined, the order of the interventions as they
#'   appear in \code{data} is used, instead.
#'
#' @return A heatmap with the proportion of missing participants in each
#'   intervention and observed comparison in the network. Each cell annotates
#'   the median, minimum and maximum (the latter two in parenthesis) proportion
#'   of missing participants across the corresponding trials.
#'   The proportion of missing participants in each intervention and observed
#'   comparison are depicted in the main diagonal and lower off-diagonal with
#'   white and black colour, respectively. The pairwise comparisons are read
#'   from left to right.
#'
#'   The 'five-and-twenty' rule of Sackett and colleagues (1997) is used to
#'   characterise the \bold{median} proportion of missing participants as being
#'   associated with low (up to 5\%), moderate (more than 5\% and up to 20\%),
#'   and high risk of bias (more than 20\%). Low, moderate, and high risk of
#'   bias associated with missing participants is indicated using green, orange,
#'   and red colour, respectively. If missing participants have not been
#'   reported for an intervention or comparison, the corresponding cell is
#'   indicated in grey.
#'
#'   The summary statistics (median, minimum and maximum) for each intervention
#'   (main diagonal; white font) result from calculating the proportion of
#'   missing participants in each arm of every trial and then summarising across
#'   the corresponding trial-arms. Similarly, the summary statistics for each
#'   observed comparison (lower off-diagonal; black font) result from
#'   calculating the proportion of total missing participants in each trial and
#'   then summarising across the corresponding trials.
#'
#'   \code{heatmap_missing_network} can be used only for a network of
#'   interventions. Otherwise, the execution of the function will be stopped and
#'   an error message will be printed on the R console. Likewise, when the
#'   number of missing participants has not been extracted for any arm of the
#'   trials.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run_model}}
#'
#' @references
#' Sackett DL, Richardson WS, Rosenberg WM, Haynes RB. Evidence-based medicine:
#' how to practice and teach EBM. New York: Churchill Livingstone 1997.
#' ISBN: 0-443-05686-2.
#'
#' @examples
#' data("nma.stowe2011")
#'
#' # Return the first six trials of the dataset
#' head(nma.stowe2011)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("PBO+LD", "DA+LD", "COMTI+LD", "MAOBI+LD")
#'
#' # Create the heatmap
#' heatmap_missing_network(data = nma.stowe2011,
#'                         drug_names = interv_names)
#'
#' @export
heatmap_missing_network <- function(data, drug_names) {


  if (dim(data %>% select(starts_with("m")))[2] == 0) {
    aa <- "Missing participant outcome data have *not* been collected."
    stop(paste(aa, "This function cannot be used."), call. = FALSE)
  }

  if (dim(data %>% dplyr::select(starts_with("r")))[2] > 0) {
    measure <- "OR"
  } else {
    measure <- "MD"
  }

  # Use the 'data.preparation' function
  dat <- data_preparation(data, measure)

  # Condition when 'drug_names' is not defined
  drug_names <- if (missing(drug_names)) {
    aa <- "The argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as specified in argument 'data' is"
    cc <- "used, instead."
    message(paste(aa, bb, cc))
    as.character(1:dat$nt)
  } else {
    drug_names
  }
  len_drugs <- length(drug_names)

  # Unique comparisons with the baseline intervention
  # A function to extract numbers from a character.
  # Source: http://stla.github.io/stlapblog/posts/numextract.html
  numextract <- function(string) {
    unlist(regmatches(string, gregexpr("[[:digit:]]+\\.*[[:digit:]]*", string)))
  }

  # Proportion of MOD per trial-arm
  arm_mod <- dat$m_pseudo / dat$N
  arm_mod[arm_mod < 0] <- NA

  # Minimum proportion of MOD per intervention in % (across the trials)
  min_mod_interv <- round(aggregate(unlist(arm_mod),
                                    by = list(unlist(dat$t)), min)[, 2],
                          2) * 100
  min_mod_interv[is.na(min_mod_interv)] <- -1

  # Median proportion of MOD per intervention in % (across the trials)
  median_mod_interv <- round(aggregate(unlist(arm_mod),
                                       by = list(unlist(dat$t)), median)[, 2],
                             2) * 100
  median_mod_interv[is.na(median_mod_interv)] <- -1

  # Max proportion of MOD per intervention in % (across the trials)
  max_mod_interv <- round(aggregate(unlist(arm_mod),
                                    by = list(unlist(dat$t)), max)[, 2],
                          2) * 100
  max_mod_interv[is.na(max_mod_interv)] <- -1

  # Turn into long format using the 'pairwise' function (netmeta)
  pair_mod <- pairwise(as.list(dat$t),
                        event = as.list(dat$m_pseudo),
                        n = as.list(dat$N),
                        data = cbind(dat$t, dat$m_pseudo, dat$N),
                        studlab = 1:dat$ns)[, c(3:6, 8, 7, 9)]
  colnames(pair_mod) <- c("study", "t1", "t2", "m1", "m2", "n1", "n2")
  pair_mod[pair_mod < 0] <- NA

  # The comparison between the second and first arms of each trial
  comp <- paste(pair_mod[, "t2"], "vs", pair_mod[, "t1"])

  # Proportion of MOD per trial-comparison
  trial_mod <- apply(pair_mod[, c("m2", "m1")], 1, sum, na.rm = TRUE) /
    apply(pair_mod[, c("n2", "n1")], 1, sum)
  # Keep as 'NA' trial-arms without reported missing participants
  cond <- apply(pair_mod[, c("m2", "m1")], 1, sum, na.rm = FALSE)
  trial_mod[is.na(cond) & trial_mod == 0] <- NA

  # Minimum proportion of MOD per observed comparison in % (across the trials)
  min_mod_comp <- round(aggregate(trial_mod,
                                  by = list(comp), min)[, 2], 2) * 100
  min_mod_comp[is.na(min_mod_comp)] <- -1

  # Median proportion of MOD per observed comparison in % (across the trials)
  median_mod_comp <- round(aggregate(trial_mod,
                                     by = list(comp), median)[, 2], 2) * 100
  median_mod_comp[is.na(median_mod_comp)] <- -1

  # Maximum proportion of MOD per observed comparison in % (across the trials)
  max_mod_comp <- round(aggregate(trial_mod,
                                  by = list(comp), max)[, 2], 2) * 100
  max_mod_comp[is.na(max_mod_comp)] <- -1

  # Lower triangular matrix: comparisons read from the left to the right
  unique_comp0 <- aggregate(trial_mod, by = list(comp), max)[, 1]
  unique_comp <- matrix(as.numeric(numextract(unique_comp0)),
                        nrow = length(unique_comp0),
                        ncol = 2,
                        byrow = TRUE)
  median_mat <- matrix(NA, nrow = len_drugs, ncol = len_drugs)
  min_mat <- max_mat <- median_mat
  median_mat[unique_comp] <- median_mod_comp
  min_mat[unique_comp] <- min_mod_comp
  max_mat[unique_comp] <- max_mod_comp

  # 'Paste' the median with the minimum and maximum
  final <- matrix(paste0(median_mat, "\n",
                         "(", min_mat, ",", " ", max_mat, ")"),
                  nrow = len_drugs,
                  ncol = len_drugs)
  diag(final) <- paste0(median_mod_interv,  "\n",
                        "(", min_mod_interv, ",", " ", max_mod_interv, ")")
  colnames(final) <- rownames(final) <- drug_names

  # Preparing the dataset for the ggplot2
  mat <- median_mat
  diag(mat) <- median_mod_interv
  mat_new <- cbind(melt(final, na.rm = TRUE), melt(mat, na.rm = FALSE)[, 3]) #HERE
  colnames(mat_new) <- c("Var1", "Var2", "value", "value2")

  # Keep only the rows without 'NA' and use for ggplot2
  mat_final <- mat_new[complete.cases(mat_new), ]
  value2_cat <- value_interv <- NULL
  mat_final$value2_cat <-
    ifelse(mat_final$value2 < 0, "not reported",
           ifelse(mat_final$value2 >= 0 & mat_final$value2 <= 5, "low",
                  ifelse(mat_final$value2 > 20, "high", "moderate")))
  value_interv <- ifelse(mat_final$Var1 == mat_final$Var2, 1, 0)
  mat_final$value[mat_final$value2 < 0] <- " "

  # Create the heatmap for one network of interventions
  ggplot(mat_final,
         aes(as.factor(Var2),
             factor(Var1, levels = drug_names[rev(seq_len(len_drugs))]),
             fill = value2_cat)) +
    geom_tile(colour = "white") +
    geom_text(aes(as.factor(Var2),
                  factor(Var1, levels = drug_names[rev(seq_len(len_drugs))]),
                  label = value,
                  fontface = "bold"),
              colour = ifelse(value_interv < 1, "black", "white"),
              size = rel(4.5)) +
    scale_fill_manual(breaks = c("low", "moderate", "high", "not reported"),
                      values = c("#009E73", "orange", "#D55E00", "grey")) +
    scale_x_discrete(position = "top") +
    labs(x = "",
         y = "",
         caption = "Proportion of missing participants:
         median (minimum, maximum)",
         fill = "Risk of bias due to missing participants") +
    theme_classic() +
    theme(legend.position = "bottom",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 12))
}
