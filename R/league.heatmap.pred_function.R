#' A heatmap league table with ordered interventions for the predictions
#'
#' @description
#' A function to create a heatmap with the treatment effects of all possible comparisons of interventions.
#' The rows and columns of the heatmap refer to the interventions which are sorted by decreasing order from the best to the worst
#' by their SUCRA value. The off-diagonals contain the posterior mean and the 95\% credible interval
#' of the corresponding comparisons. The main diagonal contains the SUCRA values of the corresponding interventions. Results in the lower
#' triangle refer to comparisons in the opposite direction after converting negative values into positive values, and vice versa.
#' Darker shades of red correspond to larger treatment effects. Comparisons between interventions should be read from left to right and
#' the estimate in the cell refers to the row-defining intervention against the column-defining intervention. Results that indicate strong
#' evidence in favor of the row-defining intervention (i.e. the respective 95\% credible interval does not include the zero value of no difference)
#' are indicated with a double asterisk.
#'
#' @param net An object of S3 class \code{nma.continuous.full.model}.
#' @param drug.names A vector of characteristics with name of the interventions as appear in the function \code{nma.continuous.full.model}.
#'
#' @return A heatmap of the treatment effects of all possible comparisons in the off-diagonals, and the SUCRA values in the diagonals
#'
#' \dontshow{load("netmodr/data/NMA_results.RData")}
#' @examples
#' drug.names <- sapply(1:14, function(x) letters[x])
#' league.heatmap(net = res, drug.names = drug.names)
#'
#' @export
league.heatmap.pred <- function(net, drug.names){


  if(length(drug.names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis", call. = F)
  }


  model <- if (net$model == "FE") {
    stop("Prediction is *not* relevant in the fixed-effect model")
  }

  par.pred <- net$EM.pred; par <- net$EM; sucra <- net$SUCRA; measure <- net$measure


  ## Source: https://rdrr.io/github/nfultz/stackoverflow/man/reflect_triangle.html
  reflect_triangle <- function(m, from=c("lower", "upper")) {
    ix <- switch(match.arg(from), lower=upper.tri, upper=lower.tri)(m, diag=FALSE)
    m[ix] <- t(-m)[ix]
    m
  }


  ## Matrix of effect measure for all possible comparisons
  # Lower triangle
  point0 <- lower0 <- upper0 <- lower.pred0 <- upper.pred0 <- matrix(NA, nrow = length(drug.names), ncol = length(drug.names))
  point0[lower.tri(point0, diag = F)] <- round(par[, 1], 2)
  # Incorporate upper triangle
  point1 <- reflect_triangle(point0, from = "lower")


  ## Matrix of lower and upper bound of effect measure for all possible comparisons
  # Lower triangle
  lower0[lower.tri(lower0, diag = F)] <- round(par[, 3], 2)
  upper0[lower.tri(upper0, diag = F)] <- round(par[, 4], 2)
  # Incorporate upper triangle
  lower1 <- reflect_triangle(upper0, from = "lower")
  lower1[lower.tri(lower1, diag = F)] <- round(par[, 3], 2)
  upper1 <- reflect_triangle(lower0, from = "lower")
  upper1[lower.tri(upper1, diag = F)] <- round(par[, 4], 2)


  ## Matrix of lower and upper bound of predictions for all possible comparisons
  # Lower triangle
  lower.pred0[lower.tri(lower.pred0, diag = F)] <- round(par.pred[, 3], 2)
  upper.pred0[lower.tri(upper.pred0, diag = F)] <- round(par.pred[, 4], 2)
  # Incorporate upper triangle
  lower.pred1 <- reflect_triangle(upper.pred0, from = "lower")
  lower.pred1[lower.tri(lower.pred1, diag = F)] <- round(par.pred[, 3], 2)
  upper.pred1 <- reflect_triangle(lower.pred0, from = "lower")
  upper.pred1[lower.tri(upper.pred1, diag = F)] <- round(par.pred[, 4], 2)


  ## Interventions order according to their SUCRA value (from best to worst)
  drug.order.col <- drug.order.row <- order(drug.names[order(-sucra[, 1])])


  ## Order interventions according to their SUCRA value (from the best to the worst)
  order.drug <- drug.names[order(-sucra[, 1])]


  ## Symmetric matrix for effect measure and its bounds after ordering rows and columns from the best to the worst intervention
  if (measure != "OR" & measure != "ROM") {
    # Effect estimates
    point <- point1[order(drug.order.col), order(drug.order.row)]
    lower <- lower1[order(drug.order.col), order(drug.order.row)]
    upper <- upper1[order(drug.order.col), order(drug.order.row)]
    # Predictions
    lower.pred <- lower.pred1[order(drug.order.col), order(drug.order.row)]
    upper.pred  <- upper.pred1[order(drug.order.col), order(drug.order.row)]
  } else {
    # Effect estimates
    point <- round(exp(point1[order(drug.order.col), order(drug.order.row)]), 2)
    lower <- round(exp(lower1[order(drug.order.col), order(drug.order.row)]), 2)
    upper <- round(exp(upper1[order(drug.order.col), order(drug.order.row)]), 2)
    # Predictions
    lower.pred <- round(exp(lower.pred1[order(drug.order.col), order(drug.order.row)]), 2)
    upper.pred  <- round(exp(upper.pred1[order(drug.order.col), order(drug.order.row)]), 2)
  }



  ## Merge point estimate with 95% credible interval in a new symmetric matric
  # Effect estimates
  (final <- matrix(paste0(point,  "\n", "(", lower, ",", " ", upper, ")"), nrow = length(drug.names), ncol = length(drug.names)))
  colnames(final) <- order.drug; rownames(final) <- order.drug
  # Predictions
  (final.pred <- matrix(paste0(point,  "\n", "(", lower.pred , ",", " ", upper.pred , ")"), nrow = length(drug.names), ncol = length(drug.names)))
  colnames(final.pred) <- order.drug; rownames(final.pred) <- order.drug


  ## Include SUCRA values in the diagonal of the new matrix
  diag(final) <- diag(final.pred) <- paste0(round(sort(sucra[, 1]*100, decreasing = T), 1), "%")


  ## Preparing the dataset for the ggplot2
  mat.new1 <- melt(final, na.rm = F); mat.new.pred1 <- melt(final.pred, na.rm = F)


  ## Upper triangle for effect estimates, lower triangle fro predictions
  dummy <- melt(lower.tri(matrix(1:(length(drug.names)*length(drug.names)), nrow = length(drug.names), ncol = length(drug.names)), diag = T))
  mat.final <- mat.new1
  mat.final[, 3] <- ifelse(dummy[, 3] == T, mat.new.pred1[, 3], mat.new1[, 3])


  ## Spot the statistically significant comparisons (i.e. the 95% CrI does not include the value of no difference)
  if (measure != "OR" & measure != "ROM") {
    # Effect estimate
    (signif.status <- melt(ifelse(upper < 0 | lower > 0, "significant", "non-significant"), na.rm = F)[3])
    signif.status[is.na(signif.status)] <- 0
    # Prediction
    (signif.status.pred <- melt(ifelse(upper.pred < 0 | lower.pred > 0, "significant", "non-significant"), na.rm = F)[3])
    signif.status.pred[is.na(signif.status.pred)] <- 0
  } else {
    # Effect estimate
    (signif.status <- melt(ifelse(upper < 1 | lower > 1, "significant", "non-significant"), na.rm = F)[3])
    signif.status[is.na(signif.status)] <- 1
    # Prediction
    (signif.status.pred <- melt(ifelse(upper.pred < 1 | lower.pred > 1, "significant", "non-significant"), na.rm = F)[3])
    signif.status.pred[is.na(signif.status.pred)] <- 1
  }



  ## Separate statistical significance for effect estiamte and prediction results
  signif.status.final <- dummy[, 3]
  signif.status.final<- ifelse(dummy[, 3] == T, signif.status.pred[, 1], signif.status[, 1])


  ## Necessary sub-dataset to color the cells according to the value of the effect measure
  # Reflect the upper triangle of the effect measure to the lower triangle
  mat <- point


  ## Merge both datasets to be used for ggplot2
  mat.new <- cbind(mat.final, melt(mat, na.rm = F)[, 3])
  colnames(mat.new) <- c("Var1", "Var2", "value", "value2")


  diag(mat) <- NA
  final_col <- melt(mat)
  mat.new$value.SUCRA <- final_col$value


  ## Hooray, the precious league table as a heatmap!
  p <- ggplot(mat.new, aes(Var2, factor(Var1, level = order.drug[length(order.drug):1]), fill = value2)) +
         geom_tile(aes(fill = value.SUCRA)) +
         geom_fit_text(aes(Var2, Var1, label = value), reflow = T) +
         geom_fit_text(aes(Var2, Var1, label = value, fontface = ifelse(signif.status.final == "significant", "bold", "plain")), reflow = T) +
         scale_fill_gradientn(colours = c("blue", "white", "red"),
                              values = rescale(c(min(mat.new$value2, na.rm = T), ifelse(measure != "OR" & measure != "ROM", 0, 1), max(mat.new$value2, na.rm = T))),
                              limits = c(min(mat.new$value2, na.rm = T), max(mat.new$value2, na.rm = T))) +
         scale_x_discrete(position = "top") +
         labs(x = "", y = "") +
         theme_bw() +
         theme(legend.position = "none", axis.text.x = element_text(size = 12, angle = 50, hjust = 0.0), axis.text.y = element_text(size = 12))

  return(p)

}


