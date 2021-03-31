#' Bland-Altman plot on the posterior mean of deviance between consistency and unrelated mean effects (UME) models
#'
#' @export
BlandAltman.plot <- function(group1, group2, colour){



  ## A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
  Numextract <- function(string){
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
  }



  ## A matriX with the trial-arms of the analysed network: first column is the trial id, second column is the arm id (i.e. 2, 3, and so on)
  (trial.arm0 <- matrix(as.numeric(Numextract(names(group1))), nrow = length(group1), ncol = 2, byrow = T))



  ## Use the matrix above to create a vector of trial-arms using the 'paste0' function: 'trial id', 'arm id'
  (trial.arm <- paste0(trial.arm0[, 1], ",", trial.arm0[, 2]))


  ## Round to the second decimal
  scaleFUN <- function(x) sprintf("%.2f", x)


  difference <- group1 - group2
  average <- apply(cbind(group1, group2), 1, mean)

  ## Dataset to be plotted in ggplot2
  data <- data.frame(round(average, 4), round(difference, 4), trial.arm)

  ## Bias and Limits of Agreement (LoA)
  bias <- mean(difference)
  lowerLoA <- bias - 1.96*sd(difference)
  upperLoA <- bias + 1.96*sd(difference)

  ggplot(data, aes(x = data[, 1], y = data[, 2])) +
    geom_point(size = 2, colour = colour) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_hline(yintercept = bias, linetype = 1) +
    geom_hline(yintercept = lowerLoA, linetype = 1) +
    geom_hline(yintercept = upperLoA, linetype = 1) +
    geom_text(aes(x = data[, 1], y = data[, 2], label = data[, 3]), color = "black", fontface = "bold",
              hjust = -0.2, vjust = -0.3, size = 4.0, check_overlap = T, inherit.aes = T) +
    scale_y_continuous(labels = scaleFUN) +
    scale_x_continuous(labels = scaleFUN) +
    theme_classic() +
    labs(x = paste("Average posterior mean deviance contribution")) +
    labs(y = paste("Difference in posterior mean deviance contribution")) +
    theme(axis.text = element_text(size = 12), axis.text.x = element_text(size = 12),
            axis.text.y = element_text(size = 12), axis.title = element_text(size = 12))


}
