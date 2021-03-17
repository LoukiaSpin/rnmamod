#' Scatterplots on the posterior mean of deviance between consistency and unrelated mean effects (UME) models
#'
#' @export
scatterplots.dev <- function(full, ume, drug.names) {

  #full <- res1; ume <- ume1; drug.names = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M", "N")

  ## The results on the following parameters will be used:
  # Posterior mean of deviance contribution (missing outcomes) - Consistency model
  dev.m.full <- full$dev.m

  # Posterior mean of deviance contribution (observed outcomes) - Consistency model
  dev.o.full <- full$dev.o

  # Posterior mean of deviance contribution (missing outcomes) - UME model
  dev.m.ume <- ume$dev.m

  # Posterior mean of deviance contribution (observed outcomes) - UME model
  dev.o.ume <- ume$dev.o



  ## A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
  Numextract <- function(string){
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
  }



  ## A matriX with the trial-arms of the analysed network: first column is the trial id, second column is the arm id (i.e. 2, 3, and so on)
  (trial.arm0 <- matrix(as.numeric(Numextract(names(dev.o.full[, 1]))), nrow = length(dev.m.full[, 1]), ncol = 2, byrow = T))



  ## use the matrix above to create a vector of trial-arms using the 'paste0' function: 'trial id', 'arm id'
  (trial.arm <- paste0(trial.arm0[, 1], ",", trial.arm0[, 2]))



  ## Prepare the data-frame for the scatterplot (using ggplot2)
  (prepare.dev <- data.frame(trial.arm, dev.o.full[, 1], dev.o.ume[, 1], dev.m.full[, 1], dev.m.ume[, 1]))
  colnames(prepare.dev) <- c("trial_arm", "dev.full.o", "dev.ume.o", "dev.full.m", "dev.ume.m")



  ## Scatterplot on the observed outcomes
  observed <- ggplot(data = prepare.dev, aes(x = dev.full.o, y = dev.ume.o)) +
                geom_point(size = 2, colour = "red") +
                geom_abline(slope = 1,lty = 2, size = 1) +
                geom_text(aes(x = dev.full.o, y = dev.ume.o, label = trial_arm), color = "black", fontface = "bold",
                         hjust = -0.2, vjust = -0.3, size = 4.0, check_overlap = T, inherit.aes = T) +
                labs(x = "Consistency model", y = "Unrelated mean effects model") +
                xlim(min(prepare.dev$dev.full.o, prepare.dev$dev.ume.o), max(prepare.dev$dev.full.o, prepare.dev$dev.ume.o)) +
                ylim(min(prepare.dev$dev.full.o, prepare.dev$dev.ume.o), max(prepare.dev$dev.full.o, prepare.dev$dev.ume.o)) +
                ggtitle("Observed outcome given the completers") +
                theme_classic() +
                theme(axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12),
                      axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))



  ## Scatterplot on the missing outcomes
  missing <- ggplot(data = prepare.dev, aes(x = dev.full.m, y = dev.ume.m)) +
               geom_point(size = 2, colour = "green4") +
               geom_abline(slope = 1,lty = 2, size = 1) +
               geom_text(aes(x = dev.full.m, y = dev.ume.m, label = trial_arm), color = "black", fontface = "bold",
                         hjust = -0.2, vjust = -0.3, size = 4.0, check_overlap = T, inherit.aes = T) +
               labs(x = "Consistency model", y = "Unrelated mean effects model") +
               xlim(min(prepare.dev$dev.full.m, prepare.dev$dev.ume.m), max(prepare.dev$dev.full.m, prepare.dev$dev.ume.m)) +
               ylim(min(prepare.dev$dev.full.m, prepare.dev$dev.ume.m), max(prepare.dev$dev.full.m, prepare.dev$dev.ume.m)) +
               ggtitle("Number of missing outcome data") +
               theme_classic() +
               theme(axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12),
                     axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))


  scatterplots <- ggarrange(observed, missing, nrow = 1, ncol = 2, labels = c("A)", "B)"))


  return(scatterplots = scatterplots)

}
