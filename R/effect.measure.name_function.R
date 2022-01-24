effect_measure_name <- function(x) {

  abbrev <- c("OR", "RR", "RD", "MD", "SMD", "ROM")
  name <- c("Odds ratio",
            "Relative risk",
            "Risk difference",
            "Mean difference",
            "Standardised mean difference",
            "Ratio of means")

  if (is.element(x, abbrev)) {
    name[which(x == abbrev)]
  } else {
    message("Insert 'OR', 'RR', 'RD', 'MD', 'SMD', or 'ROM'")
  }

}
