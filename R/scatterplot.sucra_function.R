#' Scatterplot of surface under the cumulative ranking curve (SUCRA) values for two models
#'
#' @description XX
#'
#' @param full An object of S3 class \code{\link{run.model}}. See 'Value' in \code{\link{run.model}}.
#' @param reg An object of S3 class \code{\link{run.metareg}}. See 'Value' in \code{\link{run.metareg}}.
#' @param cov.value A vector of two elements in the following order: a number that corresponds to a value of the covariate considered in \code{\link{run.metareg}},
#'   and a character object to indicate the name of the covariate.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @export
scatteplot.sucra <- function(full, reg, cov.value, drug.names) {


  options(warn = -1)

  if (length(unique(reg$covariate)) < 3 & !is.element(cov.value[1], reg$covariate)) {
    stop("The first element of the argument 'cov.value' is out of the value range of the analysed covariate", call. = F)
  } else if (length(unique(reg$covariate)) > 2 & (cov.value[1] < min(reg$covariate) | cov.value[1] > max(reg$covariate))) {
    stop("The first element of the argument 'cov.value' is out of the value range of the analysed covariate", call. = F)
  }

  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    nt <- length(full$SUCRA[, 1])
    as.character(1:nt)
  } else {
    drug.names
  }

  cov.value <- if (!is.null(reg$beta.all) & missing(cov.value)) {
    stop("The argument 'cov.value' has not been defined", call. = F)
  } else if (!is.null(reg$beta.all) & length(cov.value) < 2) {
    stop("The argument 'cov.value' must be a vector with elements a number and a character", call. = F)
  } else if (!is.null(reg$beta.all) & length(cov.value) == 2) {
    cov.value
  }

  covar <- if (length(unique(reg$covariate)) < 3) {
    as.numeric(cov.value[1])
  } else {
    as.numeric(cov.value[1]) - mean(reg$covariate)
  }


  ## Source: https://rdrr.io/github/nfultz/stackoverflow/man/reflect_triangle.html
  reflect_triangle <- function(m, from=c("lower", "upper")) {
    ix <- switch(match.arg(from), lower=upper.tri, upper=lower.tri)(m, diag=FALSE)
    m[ix] <- t(-m)[ix]
    m
  }

  sucra.nma <- round(full$SUCRA[, 1]*100, 0)
  par.mean <- reg$EM[, 1] + reg$beta.all[, 1]*covar
  par.sd <- sqrt(((reg$EM[, 2])^2) + ((reg$beta.all[, 2]*covar)^2))
  par.lower <- par.mean - 1.96*par.sd
  par.upper <- par.mean + 1.96*par.sd
  par <- data.frame(par.mean, par.lower, par.upper)
  z.test <- par.mean/par.sd
  z.test.mat <- matrix(NA, nrow = length(drug.names), ncol = length(drug.names))
  z.test.mat[lower.tri(z.test.mat, diag = F)] <- z.test*(-1)
  z.test.mat <- reflect_triangle(z.test.mat, from = "lower")
  prob.diff <- if (full$D == 0) {
    pnorm(z.test.mat)
  } else {
    1 - pnorm(z.test.mat)
  }
  sucra.nmr <- round(apply(prob.diff, 1, sum, na.rm = T)/(length(drug.names) - 1)*100, 0) # The p-scores per intervention

  dataset <- data.frame(sucra.nma, sucra.nmr, drug.names)

  ggplot(dataset, aes(x = sucra.nmr, y = sucra.nma)) +
    geom_rect(aes(xmin = 80, xmax = 100, ymin = 80, ymax = 100), fill = "#009E73", alpha = 0.1) +
    geom_rect(aes(xmin = 50, xmax = 80, ymin = 50, ymax = 80), fill = "orange", alpha = 0.1) +
    geom_rect(aes(xmin = 0, xmax = 50, ymin = 0, ymax = 50), fill = "#D55E00", alpha = 0.1) +
    geom_point(size = 2, color = "blue") +
    geom_abline(lty = 1, size = 1, col = "grey60") +
    geom_text_repel(data = dataset, aes(x = sucra.nmr, y = sucra.nma, label = drug.names),
                    color = "black", fontface = "bold", hjust = "right", size = 3.8, max.overlaps = Inf, nudge_x = -0.1, direction = "both") +
    scale_x_continuous(expand = c(0, 0), limits = c(0, 100)) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100)) +
    labs(x = "Network meta-regression SUCRA values (%)", y = "Network meta-analysis SUCRA values (%)",
         caption = paste("Results for", cov.value[2], cov.value[1])) +
    theme_classic() +
    theme(axis.title = element_text(color = "black", face = "bold", size = 12),
          axis.text = element_text(color = "black", size = 12),
          plot.caption = element_text(color = "black", size = 11))
}


