#' Heatmap of the risk of bias associated with missing participant outcome data in the network
#'
#' @description This function illustrates the risk of bias associated with missing participant outcome data (MOD) for each intervention and observed comparison of the network.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models.
#'   See 'Format' in \code{\link[rnmamod]{run.model}} function for the specification of the columns.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in the argument \code{data} is used, instead.
#'
#' @return The percentage of MOD in each intervention and observed comparison are depicted using white and black colour in the main diagonal and lower off-diagonal, respectively, of the heatmap.
#'   The percentage of MOD in an intervention is defined as the ratio of the sum of MOD for that intervention across the corresponding trials to the sum of the
#'   randomised participants in that intervention. Similarly, the percentage of MOD in an observed comparison is the ratio of the sum of MOD for that comparison
#'   across the corresponding trials to the sum of the randomised participants in that comparison.
#'   The pairwise comparisons are read from left to right. We used the 'five-and-twenty' rule of Sackett and colleagues to characterise the percentage of MOD as being associated with low (up to 5\%),
#'   moderate (more than 5\% and up to 20\%), and high risk of bias (more than 20\%). Low, moderate, and high risk of bias associated with MOD is indicated using green, orange, and red colour, respectively.
#'   The function is also applicable for a pairwise meta-analysis.
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
#' data("nma.liu2013")
#'
#' # Return the first six trials of the dataset
#' head(nma.liu2013)
#' #            study t1 t2 t3 r1 r2 r3 m1 m2 m3  n1  n2 n3
#' #    Richard, 2012  1  3  4 15 16 23  6  8  4  39  42 34
#' #     Barone, 2010  1  2 NA 27 38 NA 19 20 NA 152 144 NA
#' # Weinbtraub, 2010  1  3 NA  2  5 NA  6  6 NA  27  28 NA
#' #      Menza, 2009  1  4  5  4  2  9  6  7  5  17  18 17
#' #      Devos, 2008  1  4  5  4  8 11  0  2  1  16  15 17
#' #   Antonini, 2006  4  5 NA 10  8 NA  4  4 NA  16  15 NA
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("placebo", "pramipexole", "SNRI", "SSRI", "TCA", "pergolide")
#'
#' # Create the heatmap
#' heatmap.mod.network(data = nma.liu2013, drug.names = interv.names)
#'
#' @export
heatmap.mod.network <- function(data, drug.names){


  ## Obtain dataset
  m <- if (dim(data %>% select(starts_with("m")))[2] == 0) {
    stop("Missing participant outcome data have *not* been collected. This function cannot be used.", call. = F)
  } else {
    data %>% select(starts_with("m"))                                    # Number of missing participants in each arm of every trial
  }
  n <- data %>% select(starts_with("n"))
  t <- data %>% select(starts_with("t"))
  nt <- length(table(as.matrix(t)))
  ns <- length(m[, 1])


  ## Turn arm-level to contrast-level dataset
  (pair <- pairwise(as.list(t), event = as.list(m), n = as.list(n), data = cbind(t, m, n), studlab = 1:ns)[, c(3:6, 8, 7, 9)])
  colnames(pair) <- c("study", "t1", "t2", "m1", "m2", "n1", "n2")


  ## Calculate summary of %MOD in each intervention
  risk.drug <- round(aggregate(unlist(m), by = list(unlist(t)), sum)[, 2]/aggregate(unlist(n), by = list(unlist(t)), sum)[, 2], 2)


  ## Calculate %total MOD per observed comparison
  trial.mod <- apply(pair[, 4:5], 1, sum)
  trial.size <- apply(pair[, 6:7], 1, sum)
  comp <- paste0(pair[, 3], "vs", pair[, 2])
  risk.comp <- round(aggregate(trial.mod, by = list(comp), sum)[, 2]/aggregate(trial.size, by = list(comp), sum)[, 2], 2)


  ## Keep the unique pairwise comparisons and order by the baseline arm
  ## A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
  Numextract <- function(string) {
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*", string)))
  }

  comp.mat0 <- aggregate(trial.mod, by = list(comp), sum)[, 1]
  observed.comp <- matrix(Numextract(comp.mat0), nrow = length(comp.mat0), ncol = 2, byrow = T)
  t1.obs.com <- as.numeric(as.character(observed.comp[, 1]))
  t2.obs.com <- as.numeric(as.character(observed.comp[, 2]))
  comp.mat <- data.frame(t1.obs.com, t2.obs.com)



  ## Lower triangular heatmap matrix - Comparisons are read from the left to the right
  ## CAREFUL: The interventions in the drug.names should follow the order you considered to run NMA!
  mat <- matrix(NA, nrow = length(drug.names), ncol = length(drug.names))
  for (i in 1:length(comp.mat[, 1])) {
    mat[as.matrix(comp.mat[i, ])] <- round(risk.comp[i], 2)
  }
  diag(mat) <- round(risk.drug, 2)
  colnames(mat) <- drug.names[1:length(drug.names)]; rownames(mat) <- drug.names[1:length(drug.names)]
  mat.new <- melt(mat, na.rm = T)
  mat.new$value <- round(mat.new$value*100, 0)
  mat.new$risk <- ifelse(mat.new$value <= 5, "low", ifelse(mat.new$value > 20, "high", "moderate"))
  mat.new$value2 <- ifelse(mat.new$Var1 == mat.new$Var2, 1, 0)


  ## Create the heatmap for one network of interventions
  ggplot(mat.new, aes(Var2, factor(Var1, levels = drug.names[length(drug.names):1]), fill = risk)) +
    geom_tile(colour = "white") +
    geom_text(aes(Var2, Var1, label = paste0(value, "%"), fontface = "bold"), colour = ifelse(mat.new$value2 < 1, "black", "white"), size = rel(4.5)) +
    scale_fill_manual(breaks = c("low", "moderate", "high"), values = c("#009E73", "orange", "#D55E00")) +
    scale_x_discrete(position = "top") +
    labs(x = "", y = "", fill = "Risk due to missing outcome data") +
    theme_classic() +
    theme(legend.position = "bottom", axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
          legend.title = element_text(size = 12, face = "bold"), legend.text = element_text(size = 12))

}

