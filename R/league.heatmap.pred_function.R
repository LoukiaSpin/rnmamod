#' League heatmap of all possible comparisons: prediction
#'
#' @description
#' A function to create a heatmap with the predicted effects of all possible comparisons of interventions in the network.
#' \code{league.heatmap.pred} can be used only for a random-effects network meta-analysis and network meta-regression.
#' \code{league.heatmap.pred} is applied for one outcome only.
#'
#' @param full An object of S3 class \code{\link{run.model}} for network meta-analysis or \code{\link{run.metareg}} for network meta-regression. See 'Value' in \code{\link{run.model}} and \code{\link{run.metareg}}.
#' @param cov.value A vector of two elements in the following order: a number for the covariate value of interest and a character for the name of the covariate.
#' @param drug.names A vector of labels with the name of the interventions in the order they appear in the argument \code{data} of \code{\link{run.model}}. If the argument \code{drug.names} is not defined, the order of the interventions
#'   as they appear in \code{data} is used, instead.
#'
#' @return A league heatmap of the posterior mean, 95\% predictive interval of the effects for all possible comparisons in the off-diagonals, and the posterior mean of the SUCRA values in the diagonal.
#'
#' @details The rows and columns of the heatmap display the names of interventions which are sorted by decreasing order from the best to the worst
#'   based on their SUCRA value (Salanti et al., 2011). The main diagonal contains the SUCRA values of the corresponding interventions.
#'
#'   Results in the lower triangle refer to comparisons in the opposite direction after converting negative values into positive values (in absolute or logarithmic scale), and vice versa.
#'   Darker shades of red and green correspond to larger treatment effects in the upper and lower triangle, respectively, for a beneficial outcome, and the other way around for a harmful outcome.
#'   Odds ratio and ratio of means are reported in the original scale after exponentiation of the logarithmic scale.
#'
#'   Comparisons between interventions should be read from left to right.
#'   Therefore, each cell refers to the corresponding row-defining intervention against the column-defining intervention. Results that indicate strong
#'   evidence in favor of the row-defining intervention (i.e. the respective 95\% predictive interval does not include the null value) are indicated in bold.
#'
#'   \code{league.heatmap.pred} can be used only for a network of interventions. In the case of two interventions, the execution of the function will be stopped and an error message will be printed in the R console.
#'   Similarly, when the function is executed for a fixed-effect network meta-analysis.
#'
#' @author {Loukia M. Spineli}, {Chrysostomos Kalyvas}, {Katerina Papadimitropoulou}
#'
#' @seealso \code{\link{run.model}}, \code{\link{run.metareg}}
#'
#' @references
#' Salanti G, Ades AE, Ioannidis JP. Graphical methods and numerical summaries for presenting results from multiple-treatment meta-analysis: an overview and tutorial. \emph{J Clin Epidemiol} 2011;\bold{64}(2):163--71. [\doi{10.1016/j.jclinepi.2010.03.016}]
#'
#' @examples
#' data("nma.liu2013")
#'
#' \dontrun{
#' # Perform a random-effects network meta-analysis
#' res <- run.model(data = nma.liu2013,
#'                  measure = "OR",
#'                  model = "RE",
#'                  assumption = "IDE-ARM",
#'                  heter.prior = list("halfnormal", 0, 1),
#'                  mean.misspar = c(0, 0),
#'                  var.misspar = 1,
#'                  D = 1,
#'                  n.chains = 3,
#'                  n.iter = 10000,
#'                  n.burnin = 1000,
#'                  n.thin = 1)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv.names <- c("placebo", "pramipexole", "serotonin–norepinephrine reuptake inhibitor",
#'                   "serotonin reuptake inhibitor", "tricyclic antidepressant", "pergolide")
#'
#' # Create the league heatmap
#' league.heatmap.pred(full = res, drug.names = interv.names)
#' }
#'
#' @export
league.heatmap.pred <- function(full, cov.value = NULL, drug.names){


  drug.names <- if (missing(drug.names)) {
    message(cat(paste0("\033[0;", col = 32, "m", txt = "The argument 'drug.names' has not been defined. The intervention ID, as specified in 'data' is used as intervention names", "\033[0m", "\n")))
    as.character(1:length(full$SUCRA[, 1]))
  } else {
    drug.names
  }

  if(length(drug.names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }


  cov.value <- if (!is.null(full$beta.all) & missing(cov.value)) {
    stop("The argument 'cov.value' has not been defined", call. = F)
  } else if (!is.null(full$beta.all) & length(cov.value) < 2) {
    stop("The argument 'cov.value' must be a vector with elements a number and a character", call. = F)
  } else if (!is.null(full$beta.all) & length(cov.value) == 2) {
    cov.value
  }

  if (is.null(full$beta.all)) {
    par <- full$EM
    sucra <- full$SUCRA
    measure <- full$measure
  } else {
    par <- full$EM + full$beta.all*ifelse(length(unique(full$covariate)) < 3, as.numeric(cov.value[1]), as.numeric(cov.value[1]) - mean(full$covariate))
    sucra <- full$SUCRA
    measure <- full$measure
  }


  ## Source: https://rdrr.io/github/nfultz/stackoverflow/man/reflect_triangle.html
  reflect_triangle <- function(m, from=c("lower", "upper")) {
    ix <- switch(match.arg(from), lower=upper.tri, upper=lower.tri)(m, diag=FALSE)
    m[ix] <- t(-m)[ix]
    m
  }


  ## Matrix of effect measure for all possible comparisons
  # Lower triangle
  point0 <- lower0 <- upper0 <- matrix(NA, nrow = length(drug.names), ncol = length(drug.names))
  point0[lower.tri(point0, diag = F)] <- round(par[, 1], 2)
  # Incorporate upper triangle
  point1 <- reflect_triangle(point0, from = "lower")


  ## Matrix of lower and upper bound of effect measure for all possible comparisons
  # Lower triangle
  lower0[lower.tri(lower0, diag = F)] <- round(par[, 3], 2)
  upper0[lower.tri(upper0, diag = F)] <- round(par[, 7], 2)
  # Incorporate upper triangle
  lower1 <- reflect_triangle(upper0, from = "lower")
  lower1[lower.tri(lower1, diag = F)] <- round(par[, 3], 2)
  upper1 <- reflect_triangle(lower0, from = "lower")
  upper1[lower.tri(upper1, diag = F)] <- round(par[, 7], 2)


  ## Interventions order according to their SUCRA value (from best to worst)
  drug.order.col <- drug.order.row <- order(-sucra[, 1])


  ## Order interventions according to their SUCRA value (from the best to the worst)
  order.drug <- drug.names[order(-sucra[, 1])]


  ## Symmetric matrix for effect measure and its bounds after ordering rows and columns from the best to the worst intervention
  if(measure != "OR" & measure != "ROM") {
    point <- point1[drug.order.col, drug.order.row]
    lower <- lower1[drug.order.col, drug.order.row]
    upper <- upper1[drug.order.col, drug.order.row]


    ## Spot the statistically significant comparisons (i.e. the 95% CrI does not include the value of no difference)
    (signif.status <- melt(ifelse(upper < 0 | lower > 0, "significant", "non-significant"), na.rm = F)[3])
    signif.status[is.na(signif.status)] <- 0
  } else {
    point <- round(exp(point1[drug.order.col, drug.order.row]), 2)
    lower <- round(exp(lower1[drug.order.col, drug.order.row]), 2)
    upper <- round(exp(upper1[drug.order.col, drug.order.row]), 2)


    ## Spot the statistically significant comparisons (i.e. the 95% CrI does not include the value of no difference)
    signif.status <- melt(ifelse(upper < 1 | lower > 1, "significant", "non-significant"), na.rm = F)[3]
    signif.status[is.na(signif.status)] <- 1
  }



  ## Merge point estimate with 95% credible interval in a new symmetric matric
  #(final <- matrix(paste0(point, signif.status, "\n", "(", lower, ",", " ", upper, ")"), nrow = length(drug.names), ncol = length(drug.names)))
  final <- matrix(paste0(sprintf("%.2f", point),  "\n", "(", sprintf("%.2f", lower), ",", " ", sprintf("%.2f", upper), ")"), nrow = length(drug.names), ncol = length(drug.names))
  colnames(final) <- order.drug; rownames(final) <- order.drug


  ## Include SUCRA values in the diagonal of the new matrix
  diag(final) <- paste0(round(sort(sucra[, 1]*100, decreasing = T), 1), "%")


  ## Preparing the dataset for the ggplot2
  mat.new1 <- melt(final, na.rm = F)



  ## Merge both datasets to be used for ggplot2
  mat <- point
  mat.new <- cbind(mat.new1, melt(mat, na.rm = F)[, 3])
  colnames(mat.new) <- c("Var1", "Var2", "value", "value2")


  ## The final dataset for ggplot2
  diag(mat) <- NA
  final_col <- melt(mat)
  mat.new$value.SUCRA <- final_col$value


  ## Hooray, the precious league table as a heatmap!
  p <- ggplot(mat.new, aes(factor(Var2, levels = order.drug[1:length(order.drug)]), factor(Var1, levels = order.drug[length(order.drug):1]), fill = value2)) +
    geom_tile(aes(fill = value.SUCRA)) +
    geom_fit_text(aes(factor(Var2, levels = order.drug[1:length(order.drug)]), factor(Var1, levels = order.drug[length(order.drug):1]), label = value), reflow = T) +
    geom_fit_text(aes(factor(Var2, levels = order.drug[1:length(order.drug)]), factor(Var1, levels = order.drug[length(order.drug):1]), label = value, fontface = ifelse(signif.status == "significant", "bold", "plain")), reflow = T) +
    scale_fill_gradientn(colours = c("#009E73", "white", "#D55E00"),
                         values = rescale(c(min(mat.new$value2, na.rm = T), ifelse(measure != "OR" & measure != "ROM", 0, 1), max(mat.new$value2, na.rm = T))),
                         limits = c(min(mat.new$value2, na.rm = T), max(mat.new$value2, na.rm = T))) +
    scale_x_discrete(position = "top") +
    labs(x = "", y = "", caption = ifelse(!is.null(full$beta.all), paste("Posterior mean (95% predible interval) for", cov.value[2], "equal to", cov.value[1]),
                                          "Posterior mean (95% predible interval)") ) +
    theme_bw() +
    theme(legend.position = "none", axis.text.x = element_text(size = 12, angle = 50, hjust = 0.0), axis.text.y = element_text(size = 12),
          plot.caption = element_text(hjust = 0.01))


  return(p)

}


