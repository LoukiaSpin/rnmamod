.onLoad <- function(libname, pkgname) {
  invisible(suppressPackageStartupMessages(
    sapply(c("rjags", "coda", "R2jags"),
           requireNamespace, quietly = TRUE)
  ))
}


