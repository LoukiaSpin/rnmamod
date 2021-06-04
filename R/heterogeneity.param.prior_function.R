#' Determine the prior distribution for the heterogeneity parameter
#'
#' @param measure The effect measure
#' @param model The analysis model
#' @param heter.prior The prior
#'
#' @export
heterogeneity.param.prior <- function(measure, model, heter.prior) {


  ## Specification of the prior distribution for the between-trial parameter
  if (model == "RE" & missing(heter.prior)) {
    stop("The 'heter.prior' needs to be defined")
  } else if (model == "FE" & missing(heter.prior)) {
    list(NA, NA, NA)
  } else if (model == "FE") {
    message("The argument 'heter.prior' has been ignored")
    list(NA, NA, NA)
  } else if (model == "RE" & measure == "OR" & heter.prior[[1]] != "halfnormal" & heter.prior[[1]] != "uniform" & heter.prior[[1]] != "lognormal") {
    stop("Insert 'halfnormal', 'uniform', or 'lognormal'")
  } else if (model == "RE" & measure == "SMD" & heter.prior[[1]] != "halfnormal" & heter.prior[[1]] != "uniform" & heter.prior[[1]] != "logt") {
    stop("Insert 'halfnormal', 'uniform', or 'logt'")
  } else if (model == "RE" & (measure == "MD" || measure == "ROM") & heter.prior[[1]] != "halfnormal" & heter.prior[[1]] != "uniform") {
    stop("Insert 'halfnormal', or 'uniform'")
  } else if (model == "RE" & heter.prior[[1]] == "halfnormal") {
    as.numeric(c(0, heter.prior[[3]], 1))
  } else if (model == "RE" & heter.prior[[1]] == "uniform") {
    as.numeric(c(0, heter.prior[[3]], 2))
  } else if (model == "RE" & measure == "OR" & heter.prior[[1]] == "lognormal")  {
    as.numeric(c(heter.prior[[2]], heter.prior[[3]], 3))
  } else if (model == "RE" & measure != "OR" & heter.prior[[1]] == "lognormal") {
    stop("This is not the proper prior distribution when the effect measure refers to continuous outcome data")
  } else if (model == "RE" & measure == "SMD" & heter.prior[[1]] == "logt") {
    as.numeric(c(heter.prior[[2]], heter.prior[[3]], 3))
  } else if (model == "RE" & (measure == "MD" || measure == "ROM") & heter.prior[[1]] == "logt") {
    stop("There are currently no empirically-based prior distributions for MD and ROM. Choose a half-normal or a uniform prior distribution, instead")
  } else if (model == "RE" & measure == "OR" & heter.prior[[1]] == "logt") {
    stop("This is not the proper prior distribution when the effect measure refers to binary outcome data")
  }

}
