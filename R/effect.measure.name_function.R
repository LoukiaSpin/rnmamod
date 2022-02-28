effect_measure_name <- function(x, lower) {

  abbrev <- c("OR", "RR", "RD", "MD", "SMD", "ROM")
  name1 <- c("Odds ratio",
            "Relative risk",
            "Risk difference",
            "Mean difference",
            "Standardised mean difference",
            "Ratio of means")

  name2 <- c("odds ratio",
             "relative risk",
             "risk difference",
             "mean difference",
             "standardised mean difference",
             "ratio of means")

  if (is.element(x, abbrev) & lower == FALSE) {
    name1[which(x == abbrev)]
  } else if (is.element(x, abbrev) & lower == TRUE) {
    name2[which(x == abbrev)]
  } else if (!is.element(x, abbrev))  {
    message("Insert 'OR', 'RR', 'RD', 'MD', 'SMD', or 'ROM'")
  }

}
