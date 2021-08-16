.onLoad <- function(libname, pkgname) {
  invisible(suppressPackageStartupMessages(
    sapply(c("mcmcplots", "R2jags"),
           requireNamespace, quietly = TRUE)
  ))
}


