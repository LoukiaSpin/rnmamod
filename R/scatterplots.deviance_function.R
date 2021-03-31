#' Scatterplots on the posterior mean of deviance between consistency and unrelated mean effects (UME) models
#'
#' @export
scatterplots.dev <- function(full, ume, colour) {



  ## A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
  Numextract <- function(string){
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
  }



  ## A matriX with the trial-arms of the analysed network: first column is the trial id, second column is the arm id (i.e. 2, 3, and so on)
  (trial.arm0 <- matrix(as.numeric(Numextract(names(full))), nrow = length(full), ncol = 2, byrow = T))



  ## Use the matrix above to create a vector of trial-arms using the 'paste0' function: 'trial id', 'arm id'
  (trial.arm <- paste0(trial.arm0[, 1], ",", trial.arm0[, 2]))



  ## Prepare the data-frame for the scatterplot (using ggplot2)
  (prepare.dev <- data.frame(trial.arm, full, ume))



  ## Round to the second decimal
  scaleFUN <- function(x) sprintf("%.2f", x)



  ## Scatterplot on the observed outcomes
  ggplot(data = prepare.dev, aes(x = full, y = ume)) +
    geom_point(size = 2, colour = colour) +
    geom_abline(slope = 1,lty = 2, size = 1) +
    geom_text(aes(x = full, y = ume, label = trial.arm), color = "black", fontface = "bold",
              hjust = -0.2, vjust = -0.3, size = 4.0, check_overlap = T, inherit.aes = T) +
    labs(x = "Consistency model", y = "Unrelated mean effects model") +
    #xlim(min(prepare.dev$full, prepare.dev$ume), max(prepare.dev$full, prepare.dev$ume)) +
    #ylim(min(prepare.dev$full, prepare.dev$ume), max(prepare.dev$full, prepare.dev$ume)) +
    scale_y_continuous(labels = scaleFUN) +
    scale_x_continuous(labels = scaleFUN) +
    theme_classic() +
    theme(axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12),
          axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12))


}
