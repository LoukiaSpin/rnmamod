#' Heatmap of proportion of missing participants in the dataset
#'
#' @description Illustrates the proportion of missing participants and the
#'   associated risk of bias in each arm of every trial in the dataset.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level
#'   data of each trial. See 'Format' in \code{\link{run_model}}.
#' @param trial_names A vector of labels with the name of the trials in the
#'   order they appear in the argument \code{data}. If \code{trial_names} is not
#'   defined, the order of the trials as they appear in \code{data} is used,
#'   instead.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data}. If \code{drug_names} is
#'   not defined, the order of the interventions as they appear in \code{data}
#'   is used, instead.
#'
#' @return A heatmap presenting the proportion of missing participants in each
#'   trial-arm of the dataset. The columns and the rows of the heatmap
#'   correspond to the interventions and trials, respectively.
#'   The 'five-and-twenty' rule of Sackett and colleagues (1997) is used to
#'   characterise the proportion of missing participants as being associated
#'   with low (up to 5\%), moderate (more than 5\% and up to 20\%), and high
#'   risk of bias (more than 20\%). Low, moderate, and high risk of bias due to
#'   missing participants are indicated using green, orange, and red colour,
#'   respectively. The function is also applicable for a pairwise meta-analysis.
#'   If missing participants have not been reported for a trial-arm, the
#'   corresponding cell is indicated in grey.
#'
#' @details When the number of missing participants has not been extracted for
#'   any arm of the trials, the execution of the function will be stopped and
#'   an error message will be printed on the R console.
#'
#'   When there are more than 80 trials, the value on the proportion of missing
#'   participants will not appear on the heatmap. This is because the number on
#'   the cells will not be distinguishable.
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
#' data("nma.schwingshackl2014")
#'
#' # Return the first six trials of the dataset
#' head(nma.schwingshackl2014)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("aerobic", "resistance", "combined training")
#'
#' # Create the heatmap
#' heatmap_missing_dataset(data = nma.schwingshackl2014,
#'                         trial_names = nma.schwingshackl2014$study,
#'                         drug_names = interv_names)
#'
#' @export
heatmap_missing_dataset <- function(data, trial_names, drug_names) {


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
  m <- dat$m_pseudo
  n <- dat$N
  t <- dat$t
  nt <- dat$nt
  ns <- dat$ns
  na..  <- dat$na
  for (i in seq_len(length(m[, 1]))) {
    na..[i] <- table(!is.na(t[i, ]))["TRUE"]
  }

  trial_names <- if (missing(trial_names)) {
    aa <- "The argument 'trial_names' has not been defined."
    bb <- "The trial ID, as specified in the argument 'data' is used"
    cc <- "as trial names"
    message(cat(paste0("\033[0;", col = 32, "m", aa, " ", bb, " ", cc,
                       "\033[0m", "\n")))
    as.character(1:ns)
  } else {
    trial_names
  }

  drug_names <- if (missing(drug_names)) {
    aa <- "The argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as specified in argument 'data' is"
    cc <- "used as intervention names"
    message(cat(paste0("\033[0;", col = 32, "m", aa, " ", bb, " ", cc,
                       "\033[0m", "\n")))
    as.character(1:nt)
  } else {
    drug_names
  }

  # Rename properly to use gemtc
  names(m) <- paste0("m..", seq_len(max(na..)), ".")
  names(n) <- paste0("n..", seq_len(max(na..)), ".")
  names(t) <- paste0("t..", seq_len(max(na..)), ".")

  # Turn one row per trial to one row per trial-arm
  transform <- mtc.data.studyrow(cbind(t, m, n, na..),
                                 armVars = c("treatment" = "t",
                                             "response" = "m",
                                             "sampleSize" = "n"),
                                 nArmsVar = "na")

  # Turn all columns into numeric
  for (i in 1:dim(transform)[2]) {
    transform[, i] <- as.numeric(transform[, i])
  }

  # Rename interventions
  oldvals <- sort(unique(transform$treatment))
  for (i in seq_len(length(oldvals))) {
    transform[transform$treatment == oldvals[i], 2] <- drug_names[i]
  }

  # Rename trials
  for (i in seq_len(length(unique(transform$study)))) {
    transform[transform$study == i, 1] <- trial_names[i]
  }

  # Calculate proportion of MOD in trial-arm
  transform$m_prop <-
    round((transform$response / transform$sampleSize) * 100, 0)

  # Some necessary transformations
  labeling <-
    ifelse(transform$m_prop >= 0, paste0(transform$m_prop, "%"), " ")

  # For more than 80 trials do not show text on tiles
  if (ns < 80) {
    ggplot(transform,
           aes(factor(treatment, levels = drug_names),
               factor(study, levels = trial_names),
               fill = ifelse(m_prop < 0, "not reported",
                             ifelse(m_prop >= 0 & m_prop <= 5, "low",
                                    ifelse(m_prop > 20, "high", "moderate"))))
           ) +
      geom_tile(colour = "white") +
      geom_text(aes(factor(treatment, levels = drug_names),
                    factor(study, levels = trial_names),
                    label = labeling,
                    fontface = "plain"),
                size = rel(3.8)) +
      scale_fill_manual(breaks = c("low", "moderate", "high", "not reported"),
                        values = c("#009E73", "orange", "#D55E00", "grey")) +
      scale_x_discrete(position = "top") +
      labs(x = "", y = "", fill = "Risk of bias due to missing participants") +
      theme_classic() +
      theme(axis.text.x = element_text(size = 11),
            axis.text.y = element_text(size = 11),
            legend.position = "bottom",
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 11))
  } else {
    ggplot(transform,
           aes(factor(treatment, levels = drug_names),
               factor(study, levels = trial_names),
               fill = ifelse(m_prop < 0, "not reported",
                             ifelse(m_prop >= 0 & m_prop <= 5, "low",
                                    ifelse(m_prop > 20, "high", "moderate"))))
           ) +
      geom_tile(colour = "white") +
      scale_fill_manual(breaks = c("low", "moderate", "high", "not reported"),
                        values = c("#009E73", "orange", "#D55E00", "greyk1")) +
      scale_x_discrete(position = "top") +
      labs(x = "", y = "", fill = "Risk of bias due to missing participants") +
      theme_classic() +
      theme(axis.text.x = element_text(size = 11),
            axis.text.y = element_text(size = 11),
            legend.position = "bottom",
            legend.title = element_text(size = 11, face = "bold"),
            legend.text = element_text(size = 11))
  }
}
