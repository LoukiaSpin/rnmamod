#' A heatmap league table with ordered interventions
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
league.heatmap <- function(net, drug.names, expon){


  par <- net$EM; sucra <- net$SUCRA


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
  upper0[lower.tri(upper0, diag = F)] <- round(par[, 4], 2)
  # Incorporate upper triangle
  lower1 <- reflect_triangle(upper0, from = "lower")
  lower1[lower.tri(lower1, diag = F)] <- round(par[, 3], 2)
  upper1 <- reflect_triangle(lower0, from = "lower")
  upper1[lower.tri(upper1, diag = F)] <- round(par[, 4], 2)


  ## Interventions order according to their SUCRA value (from best to worst)
  drug.order.col <- drug.order.row <- order(drug.names[order(-sucra[, 1])])


  ## Order interventions according to their SUCRA value (from the best to the worst)
  order.drug <- drug.names[order(-sucra[, 1])]


  ## Symmetric matrix for effect measure and its bounds after ordering rows and columns from the best to the worst intervention
  if(missing(expon) || expon == F) {

    point <- point1[order(drug.order.col), order(drug.order.row)]
    lower <- lower1[order(drug.order.col), order(drug.order.row)]
    upper <- upper1[order(drug.order.col), order(drug.order.row)]


    ## Spot the statistically significant comparisons (i.e. the 95% CrI does not include the value of no difference)
    (signif.status <- melt(ifelse(upper < 0 | lower > 0, "significant", "non-significant"), na.rm = F)[3])
    signif.status[is.na(signif.status)] <- 0

  } else {

    point <- round(exp(point1[order(drug.order.col), order(drug.order.row)]), 2)
    lower <- round(exp(lower1[order(drug.order.col), order(drug.order.row)]), 2)
    upper <- round(exp(upper1[order(drug.order.col), order(drug.order.row)]), 2)


    ## Spot the statistically significant comparisons (i.e. the 95% CrI does not include the value of no difference)
    (signif.status <- melt(ifelse(upper < 1 | lower > 1, "significant", "non-significant"), na.rm = F)[3])
    signif.status[is.na(signif.status)] <- 1

  }



  ## Merge point estimate with 95% credible interval in a new symmetric matric
  #(final <- matrix(paste0(point, signif.status, "\n", "(", lower, ",", " ", upper, ")"), nrow = length(drug.names), ncol = length(drug.names)))
  (final <- matrix(paste0(point,  "\n", "(", lower, ",", " ", upper, ")"), nrow = length(drug.names), ncol = length(drug.names)))
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
  p <- ggplot(mat.new, aes(Var2, factor(Var1, level = order.drug[length(order.drug):1]), fill = value2)) +
         geom_tile(aes(fill = value.SUCRA)) +
         geom_fit_text(aes(Var2, Var1, label = value), reflow = T) +
         geom_fit_text(aes(Var2, Var1, label = value, fontface = ifelse(signif.status == "significant", "bold", "plain")), reflow = T) +
         scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = ifelse(missing(expon) || expon == F, 0, 1), na.value = "grey70") +
         scale_x_discrete(position = "top") +
         labs(x = "", y = "") +
         theme_bw() +
         theme(legend.position = "none", axis.text.x = element_text(size = 12, angle = 50, hjust = 0.0), axis.text.y = element_text(size = 12))


  return(p)

}


