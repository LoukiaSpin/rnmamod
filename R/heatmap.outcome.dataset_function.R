#' Heatmap of the distribution of the outcome in the dataset
#'
#' @description This function illustrates the the distribution of the investigated outcome in each arm of every trial in the dataset.
#'
#' @param data A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models.
#'   See 'Format' in \code{\link[rnmamod]{run.model}} function for the specification of the columns.
#' @param trial.names A vector of labels with the name of the trials in the order they appear in the argument \code{data}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data}. If \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' ##CORRECT
#' @return A heatmap presenting the outcome in each arm of every trial. IN case of a binary outcome, the heatmap illustrates the proportion of
#'   observed events in each arm of every trial. In the case of a continuous outcome, the heatmap illustrates the t-statistic in each arm of
#'   every trial. The trial-arm t-statistic is the ratio of the corresponding extracted mean outcome, \code{y}, to the standard error, \code{se} (see the \code{data.preparation} function).
#' The columns and the rows of the heatmap refer to the investigated interventions and trials, respectively.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{run.model}}, \code{\link{data.preparation}}.
#'
#' @examples
#' data("nma.schwingshackl2014")
#'
#' # Return the first six trials of the dataset
#' head(nma.schwingshackl2014)
#' #          study t1 t2 t3    y1    y2 y3  sd1  sd2 sd3 m1 m2 m3 n1 n2 n3
#' #   Bacchi, 2012  1  2 NA -0.40 -0.35 NA 0.44 0.48  NA  1  1 NA 20 20 NA
#' #       Ku, 2010  1  2 NA -0.60 -0.30 NA 1.20 0.90  NA  0  0 NA 15 13 NA
#' #      Moe, 2011  1  2 NA -0.53 -0.35 NA 0.45 0.40  NA  1  2 NA 13 13 NA
#' #       Ng, 2010  1  2 NA -0.30 -0.40 NA 0.88 0.60  NA  0  0 NA 30 30 NA
#' #   Sukala, 2012  1  2 NA -0.10 -0.10 NA 0.51 0.93  NA  4  4 NA 13 13 NA
#' # Balducci, 2009  1  3 NA  6.34  6.65 NA 0.94 1.08  NA  0  0 NA 20 22 NA
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("aerobic", "resistance", "combined training")
#'
#' # Create the heatmap
#' heatmap.outcome.dataset(data = nma.schwingshackl2014,
#'                         trial.names = nma.schwingshackl2014$study,
#'                         drug.names = interv.names)
#'
#' @export
heatmap.outcome.dataset <- function(data, trial.names, drug.names) {

  options(warn = -1)

  if(dim(data %>% dplyr::select(starts_with("r")))[2] > 0) {
    measure <- "OR"
  } else {
    measure <- "MD"
  }

  # Use the 'data.preparation' function
  dat <- data.preparation(data, measure)

  trial.names <- if (missing(trial.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'trial.names' has not been defined. The trial ID, as specified in the argument 'data' is used as trial names", "\033[0m", "\n")))
    as.character(1:dat$ns)
  } else {
    trial.names
  }

  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in argument 'data' is used as intervention names", "\033[0m", "\n")))
    as.character(1:dat$nt)
  } else {
    drug.names
  }

  ## Rename properly to use gemtc
  names(dat$t) <- paste0("t..",1:length(dat$t[1, ]),".")
  na.. <- dat$na

  outcome <- NULL
  if (measure == "OR") {
    names(dat$r) <- paste0("r..",1:length(dat$r[1, ]),".")

    n <- dat$N - dat$m  # Completers
    names(n) <- paste0("n..",1:length(n[1, ]),".")

    ## Turn one row per trial to one row per trial-arm
    transform <- mtc.data.studyrow(cbind(dat$t, dat$r, n, na..), armVars = c('treatment'= 't', 'response'='r', 'sampleSize'='n'), nArmsVar='na')

    ## Calculate the risk of observed event in each trial-arm
    transform$outcome <- round((transform$response/transform$sampleSize), 2)*100
  } else {
    names(dat$y0) <- paste0("y..",1:length(dat$y0[1, ]),".")

    names(dat$se0) <- paste0("se..",1:length(dat$se0[1, ]),".")

    transform <- mtc.data.studyrow(cbind(dat$t, dat$y0, dat$se0, na..), armVars = c('treatment'='t', 'mean'='y', 'std.err'='se'), nArmsVar='na')

    ## Calculate the t-statistic in each trial-arm
    transform$outcome <- round(transform$mean/transform$std.err, 2)
  }

  ## Turn all columns into numeric
  #for(i in 1:ncol(transform)){
  #  transform[, i] <- as.numeric(transform[, i])
  #}

  ## Rename interventions
  oldvals <- sort(unique(transform$treatment))
  for(i in 1:length(oldvals)) {
    transform[transform$treatment == oldvals[i], 2] <- drug.names[i]
  }

  ## Rename trials
  for(i in 1:length(unique(transform$study))){
    transform[transform$study == i, 1] <- trial.names[i]
  }

  ## For more than 80 trials do not show text on tiles
  if(dat$ns < 80){

    ggplot(transform, aes(factor(treatment, levels = drug.names), factor(study, levels = trial.names), fill = outcome)) +
      geom_tile(colour = "white") +
      geom_text(aes(factor(treatment, levels = drug.names), factor(study, levels = trial.names), label = outcome, fontface = "plain"), size = rel(3.8)) +
      scale_fill_gradient(low = "#009E73", high = "#D55E00", limits = c(ifelse(measure == "OR", 0, min(transform$outcome)),ifelse(measure == "OR", 100, max(transform$outcome)) )) +
      scale_x_discrete(position = "top") +
      labs(x = "", y = "", fill = ifelse(measure == "OR", "Proportion of observed events (%)", "t-statistic")) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
            legend.position = "bottom", legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 11))

  } else {

    ggplot(transform, aes(factor(treatment, levels = drug.names), factor(study, levels = trial.names), fill = outcome)) +
      geom_tile(colour = "white") +
      scale_fill_gradient(low = "#009E73", high = "#D55E00", limits = c(ifelse(measure == "OR", 0, min(transform$outcome)),ifelse(measure == "OR", 100, max(transform$outcome)) )) +
      scale_x_discrete(position = "top") +
      labs(x = "", y = "", fill = ifelse(measure == "OR", "Proportion of observed events (%)", "t-statistic")) +
      theme_classic() +
      theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
            legend.position = "bottom", legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 11))
  }


}
