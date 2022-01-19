effect_measure_name <- function(x) {

  abbrev <- c("OR", "MD", "SMD", "ROM")
  name <- c("Odds ratio",
            "Mean difference",
            "Standardised mean difference",
            "Ratio of means")

  if (is.element(x, abbrev)) {
    name[which(x == abbrev)]
  } else {
    message("Insert 'OR', 'MD', 'SMD', or 'ROM'")
  }

}
