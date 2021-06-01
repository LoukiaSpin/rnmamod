#' Create the leverage plot
#'
#' @export
leverage.plot <- function(net, drug.names, title.o, title.m) {



  ## The results on the following parameters will be used:
  # Posterior mean of deviance contribution (missing outcomes)
  dev.m <- net$dev.m

  # Posterior mean of deviance contribution (observed outcomes)
  dev.o <- net$dev.o

  # Leverage (missing outcomes)
  lev.m <- net$leverage.m

  # Leverage (observed outcomes)
  lev.o <- net$leverage.o

  # The sign of the difference between observed and fitted number of missing outcome data
  sign.m <- net$sign.dev.m

  # The sign of the difference between observed and fitted number of observed outcome
  sign.o <- net$sign.dev.o;



  ## A function to extract numbers from a character. Source: http://stla.github.io/stlapblog/posts/Numextract.html
  Numextract <- function(string){
    unlist(regmatches(string,gregexpr("[[:digit:]]+\\.*[[:digit:]]*",string)))
  }



  ## A matriX with the trial-arms of the analysed network: first column is the trial id, second column is the arm id (i.e. 2, 3, and so on)
  (trial.arm0 <- matrix(as.numeric(Numextract(names(dev.o[, 1]))), nrow = length(dev.o[, 1]), ncol = 2, byrow = T))



  ## use the matrix above to create a vector of trial-arms using the 'paste0' function: 'trial id', 'arm id'
  (trial.arm <- paste0(trial.arm0[, 1], ",", trial.arm0[, 2]))



  ## Prepare the data-frame for the leverage plot (using ggplot2)
  prepare.lev <- round(data.frame(sign.o*sqrt(as.vector(dev.o[, 1])), lev.o, sign.m*sqrt(as.vector(dev.m[, 1])), lev.m), 2)
  colnames(prepare.lev) <- c("signed.dev.o", "lev.o", "signed.dev.m", "lev.m")



  ## Keep only trial-arms that exceed the parabola y + (x)^2 = c at c = 2.
  # Observed outcomes
  poor.o <- ifelse(prepare.lev$lev.o > 2 - (prepare.lev$signed.dev.o^2) | prepare.lev$lev.o < -(2 - (prepare.lev$signed.dev.o^2)), trial.arm, NA)
  poor.fit.o <- data.frame(prepare.lev[!is.na(poor.o), 1:2], poor.o[!is.na(poor.o)])
  colnames(poor.fit.o) <- c("signed.dev", "leverage", "poor")

  # Missing outcomes
  poor.m <- ifelse(prepare.lev$lev.m > 2 - (prepare.lev$signed.dev.m^2) | prepare.lev$lev.m < -(2 - (prepare.lev$signed.dev.m^2)), trial.arm, NA)
  poor.fit.m <- data.frame(prepare.lev[!is.na(poor.m), 3:4], poor.m[!is.na(poor.m)])
  colnames(poor.fit.m) <- c("signed.dev", "leverage", "poor")



  ## Leverage plot for observed outcomes
  observed <- ggplot(data = prepare.lev, aes(x = signed.dev.o, y = lev.o)) +
               geom_point(size = 2, colour = "black") +
               geom_smooth(aes(x = signed.dev.o, y = 1 - (signed.dev.o^2)), method = 'loess', formula = 'y ~ x',  colour = "#009E73", linetype = 2) +
               geom_smooth(aes(x = signed.dev.o, y = 2 - (signed.dev.o^2)), method = 'loess', formula = 'y ~ x',  colour = "orange", linetype = 2) +
               geom_smooth(aes(x = signed.dev.o, y = 3 - (signed.dev.o^2)), method = 'loess', formula = 'y ~ x',  colour = "#D55E00", linetype = 2) +
               geom_text_repel(data = poor.fit.o, aes(x = signed.dev, y = leverage, label = poor),
                         color = "blue", fontface = "bold", hjust = "right", size = 3.8, max.overlaps = Inf, nudge_x = -0.1, direction = "y") +
               labs(x = expression(""%+-% sqrt("Posterior mean of the residual deviance")), y = "Leverage (each data point)") +
               coord_cartesian(xlim = c(min(prepare.lev$signed.dev.o), max(prepare.lev$signed.dev.o)),
                               ylim = c(0, max(3 - (prepare.lev$signed.dev.o^2))), expand = T) +
               ggtitle(title.o) +
               theme_classic() +
               theme(axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12),
                     axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                     plot.title = element_text(color = "black", size = 11, face = "bold"))



  ## Leverage plot for missing outcomes
  missing <- ggplot(data = prepare.lev, aes(x = signed.dev.m, y = lev.m)) +
               geom_point(size = 2, colour = "black") +
               geom_smooth(aes(x = signed.dev.m, y = 1 - (signed.dev.m^2)), method = 'loess', formula = 'y ~ x', colour = "#009E73", linetype = 2) +
               geom_smooth(aes(x = signed.dev.m, y = 2 - (signed.dev.m^2)), method = 'loess', formula = 'y ~ x', colour = "orange", linetype = 2) +
               geom_smooth(aes(x = signed.dev.m, y = 3 - (signed.dev.m^2)), method = 'loess', formula = 'y ~ x', colour = "#D55E00", linetype = 2) +
               geom_text_repel(data = poor.fit.m, aes(x = signed.dev, y = leverage, label = poor),
                               color = "blue", fontface = "bold", hjust = "right", size = 3.8, max.overlaps = Inf, nudge_x = -0.1, direction = "y") +
               labs(x = expression(""%+-% sqrt("Posterior mean of the residual deviance")), y = "Leverage (each data point)") +
               coord_cartesian(xlim = c(min(prepare.lev$signed.dev.m), max(prepare.lev$signed.dev.m)),
                               ylim = c(0, max(3 - (prepare.lev$signed.dev.m^2))), expand = T) +
               ggtitle(title.m) +
               theme_classic() +
               theme(axis.title.x = element_text(color = "black", size = 12), axis.title.y = element_text(color = "black", size = 12),
                     axis.text.x = element_text(color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12),
                     plot.title = element_text(color = "black", size = 11, face = "bold"))


  return(list(leverage.plot.observed = observed, leverage.plot.missing = missing))

}
