effect_measure_name <- function(name) {

  spell_out <- if (name == "OR") {
    "Odds ratio"
  } else if (name == "MD") {
    "Mean difference"
  } else if (name == "SMD") {
    "Standardised mean difference"
  } else if (name == "ROM") {
    "Ratio of means"
  }

  return(spell_out)
}
