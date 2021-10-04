#' Heatmap of the distribution of missing participant outcome data in the network
#'
#' @description Illustrates the distribution of and the risk of bias associated with missing participant outcome data (MOD) for each intervention and observed comparison in the network.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial.
#'   See 'Format' in \code{\link[rnmamod]{run.model}} function for the specification of the columns.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data}. If the argument \code{drug.names} is not defined, the interventions are ordered
#'   as they appear in the argument \code{data}.
#'
#' @return A heatmap on the proportion of MOD in each intervention and observed comparison in the network.
#'   Each cell annotates the median, minimum and maximum (in parenthesis) proportion of MOD across the corresponding trials.
#'   The proportion of MOD in each intervention and observed comparison are depicted in white and black colour in the main diagonal and lower off-diagonal, of the heatmap, respectively, 
#'   The pairwise comparisons are read from left to right.  The 'five-and-twenty' rule of Sackett and colleagues is used to characterise the \bold{median} proportion of MOD as being associated with low (up to 5\%),
#'   moderate (more than 5\% and up to 20\%), and high risk of attrition bias (more than 20\%). Low, moderate, and high risk of bias associated with MOD is indicated using green, orange, and red colour, respectively.
#'   The function is redundant for a pairwise meta-analysis.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}
#'
#' @references
#' Sackett DL, Richardson WS, Rosenberg WM, Haynes RB. Evidence-based medicine: how to practice and teach EBM.
#' New York: Churchill Livingstone 1997. ISBN: 0-443-05686-2.
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
#' heatmap.missing.network(data = nma.stowe2011, drug.names = interv.names)
#'
#' @export
heatmap.missing.network <- function(data, drug.names){

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

  ## Unique comparisons with the baseline intervention
  # A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
  Numextract <- function(string){
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
  }

  # Proportion of MOD per trial-arm
  arm.mod <- dat$m/dat$N

  # Minimum proportion of MOD per intervention in % (across the corresponding trials)
  min.mod.interv <- round(aggregate(unlist(arm.mod), by = list(unlist(dat$t)), min)[, 2], 2)*100

  # Median proportion of MOD per intervention in % (across the corresponding trials)
  median.mod.interv <- round(aggregate(unlist(arm.mod), by = list(unlist(dat$t)), median)[, 2], 2)*100

  # Max proportion of MOD per intervention in % (across the corresponding trials)
  max.mod.interv <- round(aggregate(unlist(arm.mod), by = list(unlist(dat$t)), max)[, 2], 2)*100

  # Turn into long format using the 'pairwise' function (netmeta): MOD
  pair.mod <- pairwise(as.list(dat$t),
                       event = as.list(dat$m),
                       n = as.list(dat$N),
                       data = cbind(dat$t, dat$m, dat$N),
                       studlab = 1:dat$ns)[, c(3:6, 8, 7, 9)]
  colnames(pair.mod) <- c("study", "t1", "t2", "m1", "m2", "n1", "n2")

  # The comparison between the second and first arms of each trial
  comp <- paste(pair.mod[, "t2"], "vs", pair.mod[, "t1"])

  # Proportion of MOD per trial-comparison
  trial.mod <- apply(pair.mod[, c("m1", "m2")], 1, sum)/apply(pair.mod[, c("n1", "n2")], 1, sum)

  # Minimum proportion of MOD per observed comparison in % (across the corresponding trials)
  min.mod.comp <- round(aggregate(trial.mod, by = list(comp), min)[, 2], 2)*100

  # Median proportion of MOD per observed comparison in % (across the corresponding trials)
  median.mod.comp <- round(aggregate(trial.mod, by = list(comp), median)[, 2], 2)*100

  # Maximum proportion of MOD per observed comparison in % (across the corresponding trials)
  max.mod.comp <- round(aggregate(trial.mod, by = list(comp), max)[, 2], 2)*100

  # Lower triangular heatmap matrix - Comparisons are read from the left to the right
  unique.comp0 <- aggregate(trial.mod, by = list(comp), max)[, 1]
  unique.comp <- matrix(as.numeric(Numextract(unique.comp0)), nrow = length(unique.comp0), ncol = 2, byrow = T)
  median.mat <- min.mat <- max.mat <- matrix(NA, nrow = length(drug.names), ncol = length(drug.names))
  median.mat[unique.comp] <- median.mod.comp
  min.mat[unique.comp] <- min.mod.comp
  max.mat[unique.comp] <- max.mod.comp


  # 'Paste' the median with the minimum and maximum (the latter two in a parenthesis)
  final <- matrix(paste0(median.mat,  "\n", "(", min.mat, ",", " ", max.mat, ")"), nrow = length(drug.names), ncol = length(drug.names))
  diag(final) <- paste0(median.mod.interv,  "\n", "(", min.mod.interv, ",", " ", max.mod.interv, ")")
  colnames(final) <- rownames(final) <- drug.names

  # Preparing the dataset for the ggplot2
  mat <- median.mat
  diag(mat) <- median.mod.interv
  mat.new <- cbind(melt(final, na.rm = T), melt(mat, na.rm = F)[, 3])
  colnames(mat.new) <- c("Var1", "Var2", "value", "value2")

  # Keep only the rows without 'NA' and use for ggplot2
  mat.final <- mat.new[complete.cases(mat.new), ]
  value2.cat <- value.interv <- NULL
  mat.final$value2.cat <- ifelse(mat.final$value2 <= 5, "low", ifelse(mat.final$value2 > 20, "high", "moderate"))
  mat.final$value.interv <- ifelse(mat.final$Var1 == mat.final$Var2, 1, 0)

  # Create the heatmap for one network of interventions
  ggplot(mat.final, aes(as.factor(Var2), factor(Var1, levels = drug.names[length(drug.names):1]), fill = value2.cat)) +
    geom_tile(colour = "white") +
    geom_text(aes(Var2, Var1, label = value, fontface = "bold"), colour = ifelse(mat.final$value.interv < 1, "black", "white"), size = rel(4.5)) +
    scale_fill_manual(breaks = c("low", "moderate", "high"), values = c("#009E73", "orange", "#D55E00")) +
    scale_x_discrete(position = "top") +
    labs(x = "", y = "", caption = "Proportion of missing participants: median (minimum, maximum)",
         fill = "Risk of bias due to missing participants") +
    theme_classic() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12))

}
