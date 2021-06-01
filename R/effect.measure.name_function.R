effect.measure.name <- function(name) {

  if (name == "OR") {
    spell.out <- "Odds ratio"
  } else if (name == "MD") {
    spell.out <- "Mean difference"
  } else if (name == "SMD") {
    spell.out <- "Standardised mean difference"
  } else if (name == "ROM") {
    spell.out <- "Ratio of means"
  }

}
