#' Heatmap of the propoertion of missing participant outcome data in the network
#'
#' @description Illustrates the distribution of and the risk of bias associated
#'   with missing participant outcome data (MOD) for each intervention and
#'   observed comparison in the network.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level
#'   data of each trial. See 'Format' in \code{\link[rnmamod]{run_model}}
#'   function for the specification of the columns.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data}. If the argument
#'   \code{drug_names} is not defined, the interventions are ordered
#'   as they appear in the argument \code{data}.
#'
#' @return A heatmap on the proportion of MOD in each intervention and observed
#'   comparison in the network. Each cell annotates the median, minimum and
#'   maximum (in parenthesis) proportion of MOD across the corresponding trials.
#'   The proportion of MOD in each intervention and observed comparison are
#'   depicted in white and black colour in the main diagonal and lower
#'   off-diagonal, of the heatmap, respectively, The pairwise comparisons are
#'   read from left to right.  The 'five-and-twenty' rule of Sackett and
#'   colleagues is used to characterise the \bold{median} proportion of MOD as
#'   being associated with low (up to 5\%), moderate (more than 5\% and up to
#'   20\%), and high risk of attrition bias (more than 20\%). Low, moderate, and
#'   high risk of bias associated with MOD is indicated using green, orange, and
#'   red colour, respectively. The function is redundant for a pairwise
#'   meta-analysis.
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
#' #              study t1 t2    y1    y2  sd1  sd2 m1 m2  n1  n2
#' #   DA (B): Interntl  1  2 -0.30 -1.20 4.36 4.32  7  3  83  84
#' #      DA (C): Spain  1  2 -2.47 -3.33 3.91 3.48  8  9  20  23
#' #         DA (C): UK  1  2 -0.70 -2.00 2.24 2.33  2  1  18  19
#' #      DA (C): USA 1  1  2 -0.77 -2.08 3.32 3.21 19 34  65 123
#' # DA (Pe): N America  1  2 -0.20 -1.80 4.79 4.81  0  0 187 189
#' # DA (Pr): CLEOPATRA  1  2 -0.90 -2.80 5.00 2.83  1  1 101 201
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

  options(warn = -1)

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
    cc <- "used as intervention names"
    message(cat(paste0("\033[0;", col = 32, "m", aa, " ", bb, " ", cc,
                       "\033[0m", "\n")))
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
  arm_mod <- dat$m / dat$N

  # Minimum proportion of MOD per intervention in % (across the trials)
  min_mod_interv <- round(aggregate(unlist(arm_mod),
                                    by = list(unlist(dat$t)), min)[, 2],
                          2) * 100

  # Median proportion of MOD per intervention in % (across the trials)
  median_mod_interv <- round(aggregate(unlist(arm_mod),
                                       by = list(unlist(dat$t)), median)[, 2],
                             2) * 100

  # Max proportion of MOD per intervention in % (across the trials)
  max_mod_interv <- round(aggregate(unlist(arm_mod),
                                    by = list(unlist(dat$t)), max)[, 2],
                          2) * 100

  # Turn into long format using the 'pairwise' function (netmeta)
  pair_mod0 <- pairwise(as.list(dat$t),
                        event = as.list(dat$m),
                        n = as.list(dat$N),
                        data = cbind(dat$t, dat$m, dat$N),
                        studlab = 1:dat$ns)[, c(3:6, 8, 7, 9)]
  colnames(pair_mod0) <- c("study", "t1", "t2", "m1", "m2", "n1", "n2")

  # Ensure that t1 < t2 and correspondingly for the other elements
  treat <- treat0 <- pair_mod0[, 2:3]
  miss <- miss0 <- pair_mod0[, 4:5]
  rand <- rand0 <- pair_mod0[, 6:7]
  for (i in seq_len(length(pair_mod0[, 1]))) {
    treat[i, ] <- treat0[i, order(treat0[i, ], na.last = T)]
    miss[i, ] <- miss0[i, order(treat0[i, ], na.last = T)]
    rand[i, ] <- rand0[i, order(treat0[i, ], na.last = T)]
  }
  pair_mod <- data.frame(study = pair_mod0$study, treat, miss, rand)

  # The comparison between the second and first arms of each trial
  comp <- paste(pair_mod[, "t2"], "vs", pair_mod[, "t1"])

  # Proportion of MOD per trial-comparison
  trial_mod <- apply(pair_mod[, c("m2", "m1")], 1, sum) /
    apply(pair_mod[, c("n2", "n1")], 1, sum)

  # Minimum proportion of MOD per observed comparison in % (across the trials)
  min_mod_comp <- round(aggregate(trial_mod,
                                  by = list(comp), min)[, 2], 2) * 100

  # Median proportion of MOD per observed comparison in % (across the trials)
  median_mod_comp <- round(aggregate(trial_mod,
                                     by = list(comp), median)[, 2], 2) * 100

  # Maximum proportion of MOD per observed comparison in % (across the trials)
  max_mod_comp <- round(aggregate(trial_mod,
                                  by = list(comp), max)[, 2], 2) * 100

  # Lower triangular matrix: comparisons read from the left to the right
  unique_comp0 <- aggregate(trial_mod, by = list(comp), max)[, 1]
  unique_comp <- matrix(as.numeric(numextract(unique_comp0)),
                        nrow = length(unique_comp0),
                        ncol = 2,
                        byrow = T)
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
  mat_new <- cbind(melt(final, na.rm = T), melt(mat, na.rm = F)[, 3])
  colnames(mat_new) <- c("Var1", "Var2", "value", "value2")

  # Keep only the rows without 'NA' and use for ggplot2
  mat_final <- mat_new[complete.cases(mat_new), ]
  value2_cat <- value_interv <- NULL
  mat_final$value2_cat <-
    ifelse(mat_final$value2 <= 5, "low",
           ifelse(mat_final$value2 > 20, "high", "moderate"))
  value_interv <- ifelse(mat_final$Var1 == mat_final$Var2, 1, 0)

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
    scale_fill_manual(breaks = c("low", "moderate", "high"),
                      values = c("#009E73", "orange", "#D55E00")) +
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
