#' Heatmap of the risk of bias associated with missing participant outcome data in each trial-arm
#'
#' @description This function illustrates the risk of bias associated with missing participant outcome data (MOD) in each arm of every trial in the dataset.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models.
#'   See 'Format' in \code{\link[rnmamod]{run.model}} function for the specification of the columns.
#' @param trial.names A vector of labels with the name of the trials in the order they appear in the argument \code{data}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data}. If \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @return A heatmap presenting the percentage of MOD in each trial-arm of the dataset. The columns and the rows of the heatmap refer to the investigated interventions and trials, respectively.
#'   We used the 'five-and-twenty' rule of Sackett and colleagues to characterise the percentage of MOD as being associated with low (up to 5\%), moderate (more than 5\% and up to 20\%),
#'   and high risk of bias (more than 20\%). Low, moderate, high risk of bias due to MOD are indicated using green, orange, and red colour, respectively.
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
#' heatmap.mod.trial(data = nma.liu2013, trial.names = nma.liu2013$study, drug.names = interv.names)
#'
#' @export
heatmap.mod.trial <- function(data, trial.names, drug.names) {


  m <- if (dim(data %>% dplyr::select(starts_with("m")))[2] == 0) {
    stop("Missing participant outcome data have *not* been collected. This function cannot be used.", call. = F)
  } else {
    data %>% dplyr::select(starts_with("m"))                                    # Number of missing participants in each arm of every trial
  }
  n <- data %>% dplyr::select(starts_with("n"))   # Number randomised
  t <- data %>% dplyr::select(starts_with("t"))
  nt <- length(table(as.matrix(t)))
  ns <- length(m[, 1])
  na..  <- rep(0, length(m[, 1]))
  for(i in 1:length(m[, 1])){
    na..[i] <- table(!is.na(t[i, ]))["TRUE"]
  }


  trial.names <- if (missing(trial.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'trial.names' has not been defined. The trial ID, as specified in the argument 'data' is used as trial names", "\033[0m", "\n")))
    as.character(1:ns)
  } else {
    trial.names
  }


  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in argument 'data' is used as intervention names", "\033[0m", "\n")))
    as.character(1:nt)
  } else {
    drug.names
  }


  ## Rename properly to use gemtc
  names(m) <- paste0("m..",1:length(m[1, ]),".")
  names(n) <- paste0("n..",1:length(n[1, ]),".")
  names(t) <- paste0("t..",1:length(t[1, ]),".")


  ## Turn one row per trial to one row per trial-arm
  (transform <- mtc.data.studyrow(cbind(t, m, n, na..), armVars = c('treatment'= 't', 'response'='m', 'sampleSize'='n'), nArmsVar='na'))


  ## Turn all columns into numeric
  for(i in 1:ncol(transform)){
    transform[, i] <- as.numeric(transform[, i])
  }


  ## Rename interventions
  oldvals <- sort(unique(transform$treatment))
  for(i in 1:length(oldvals)) {
    transform[transform$treatment == oldvals[i], 2] <- drug.names[i]
  }


  ## Rename trials
  for(i in 1:length(unique(transform$study))){

    transform[transform$study == i, 1] <- trial.names[i]
  }


  ## Calculate percentage of MOD in trial-arm
  transform$m.prop <- round((transform$response/transform$sampleSize)*100, 0)


  ## For more than 80 trialsm do not show text on tiles
  if(ns < 80){

    #ggplot(transform, aes(treatment, factor(study, level = 1:ns), fill = ifelse(m.prop <= 5, "low", ifelse(m.prop > 20, "high", "moderate") ))) +
    ggplot(transform, aes(treatment, study, fill = ifelse(m.prop <= 5, "low", ifelse(m.prop > 20, "high", "moderate") ))) +
      geom_tile(colour = "white") +
      geom_text(aes(treatment, study, label = paste0(m.prop, "%"), fontface = "plain"), size = rel(3.8)) +
      scale_fill_manual(breaks = c("low", "moderate", "high"), values = c("#009E73", "orange", "#D55E00")) +
      scale_x_discrete(position = "top") +
      labs(x = "", y = "", fill = "Risk of bias due to missingness") +
      theme_bw() +
      theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
            legend.position = "bottom", legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 11))

  } else {

    ggplot(transform, aes(treatment, study, fill = ifelse(m.prop <= 5, "low", ifelse(m.prop > 20, "high", "moderate") ))) +
      geom_tile(colour = "white") +
      scale_fill_manual(breaks = c("low", "moderate", "high"), values = c("green3", "orange", "firebrick1")) +
      scale_x_discrete(position = "top") +
      labs(x = "", y = "", fill = "Risk of bias due to missingness") +
      theme_bw() +
      theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
            legend.position = "bottom", legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 11))
  }


}
