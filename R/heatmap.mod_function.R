#' Heatmap to illustrate the proportion of MOD in each intervention of every trial
#'
#' @export
heatmap.mod <- function(data, trial.names, drug.names) {


  m <- data %>% dplyr::select(starts_with("m"))
  c <- data %>% dplyr::select(starts_with("c"))   # Number of completers in each arm of every trial
  n0 <- data %>% dplyr::select(starts_with("n"))  # Number randomised
  if (dim(c)[2] == 0) {
    n <- n0
  } else {
    n <- m + c
  }
  t <- data %>% dplyr::select(starts_with("t"))
  nt <- length(table(as.matrix(t)))
  ns <- length(m[, 1])
  na..  <- rep(0, length(m[, 1]))
  for(i in 1:length(m[, 1])){
    na..[i] <- table(!is.na(t[i, ]))["TRUE"]
  }


  ## Rename properly to use gemtc
  names(m) <- paste0("m..",1:length(m[1, ]),".")
  names(n) <- paste0("n..",1:length(n[1, ]),".")
  names(t) <- paste0("t..",1:length(t[1, ]),".")


  ## Turn one row per trial to one row per trial-arm
  (transform <- mtc.data.studyrow(cbind(t, m, n, na..), armVars = c('treatment'= 't', 'response'='m', 'sampleSize'='n'), nArmsVar='na'))
  for(i in 1:length(unique(transform$study))){
    transform[transform$treatment == i, 2] <- drug.names[i]
    transform[transform$study == i, 1] <- trial.names[i]
  }
  transform$m.prop <- round((transform$response/transform$sampleSize)*100, 0)


  ## For more than 80 trialsm do not show text on tiles
  if(ns < 80){

    ggplot(transform, aes(treatment, factor(study, level = 1:ns), fill = ifelse(m.prop <= 5, "low", ifelse(m.prop > 20, "high", "moderate") ))) +
      geom_tile(colour = "white") +
      geom_text(aes(treatment, study, label = paste0(m.prop, "%"), fontface = "plain"), size = rel(3.8)) +
      scale_fill_manual(breaks = c("low", "moderate", "high"), values = c("green3", "orange", "firebrick1")) +
      scale_x_discrete(position = "top") +
      labs(x = "", y = "", fill = "Risk of bias due to missingness") +
      theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
            legend.position = "bottom", legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 11))

  } else {

    ggplot(transform, aes(treatment, factor(study, level = 1:ns), fill = ifelse(m.prop <= 5, "low", ifelse(m.prop > 20, "high", "moderate") ))) +
      geom_tile(colour = "white") +
      scale_fill_manual(breaks = c("low", "moderate", "high"), values = c("green3", "orange", "firebrick1")) +
      scale_x_discrete(position = "top") +
      labs(x = "", y = "", fill = "Risk of bias due to missingness") +
      theme(axis.text.x = element_text(size = 11), axis.text.y = element_text(size = 11),
            legend.position = "bottom", legend.title = element_text(size = 11, face = "bold"), legend.text = element_text(size = 11))
  }


}
