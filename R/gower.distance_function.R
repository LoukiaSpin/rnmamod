#' Gower's dissimilarity measure
#' (Trials' comparability for transitivity evaluation)
#'
#' @description
#'   \code{gower_distance} calculate the Gower's dissimilarity coefficient for 
#'   all pairs of trials included in a network of interventions, considering 
#'   several characteristics measured at trial level. 
#'   It takes values from 0 to 1, with 0 implying complete similarity and 1
#'   complete dissimilarity.
#'
#' @param input A data-frame in the long arm-based format. Two-arm trials occupy
#'   one row in the data-frame. Multi-arm trials occupy as many rows as the
#'   number of possible comparisons among the interventions. The first two
#'   columns refer to the trial name, and the pairwise comparison,
#'   respectively. The remaining columns refer to summary characteristics. See
#'   'Details' for the specification of the columns.
#'
#' @return
#'   \code{gower_distance} returns the following list of elements:
#'   \item{Dissimilarity_table}{A lower off-diagonal matrix of 'dist' class
#'   with the dissimilarities of all pairs of trials.}
#'   \item{Types_used}{A data-frame with type mode (i.e., double or integer) of
#'   each characteristic.}
#'   \item{Total_missing}{The percentage of missing cases in the comparison,
#'   calculated as the ratio of total missing cases to the product of the number
#'   of studies with the number of characteristics.}
#'
#' @details
#'   The correct type mode of columns in \code{input} must be ensured to use
#'   the function \code{gower_distance}. The first two columns referring to
#'   the trial name, and pairwise comparison, respectively, must
#'   be \strong{character}. The remaining columns referring to the
#'   characteristics must be \strong{double} or \strong{integer} depending on
#'   whether the corresponding characteristic refers to a quantitative or
#'   qualitative variable. The type mode of each column is assessed by
#'   \code{gower_distance} using the base function \code{typeof}. Note that
#'   \code{gower_distance} invites unordered and ordered variables; for the
#'   latter, add the argument \code{ordered = TRUE} in the base function
#'   \bold{factor()}.
#'
#'  \code{gower_distance} is integrated in the function
#'   \code{\link{comp_clustering}}.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso
#'   \code{\link{comp_clustering}}
#'
#' @references
#' Gower J. General Coefficient of Similarity and Some of Its Properties.
#' \emph{Biometrics} 1971;\bold{27}(4):857--71.
#' doi: 10.2307/2528823
#'
#' @export
gower_distance <- function (input) {


  ## Check if the dataset is correct
  input <- if (any(sapply(input, typeof)[1:2] != "character")) {
    stop("The first two columns (trial and comparison) must be 'characters').",
         call. = FALSE)
  } else if (any(sapply(input, typeof)[-c(1, 2)] == "character")) {
    stop("The characteristics must be 'double' or 'integer'.", call. = FALSE)
  } else {
    input
  }


  ## Remove the first two columns
  data <- input[, -c(1:2)]


  ## Table with the variable type
  char_type <- data.frame(characteristic = colnames(data),
                          type = sapply(data, typeof))
  rownames(char_type) <- 1:dim(data)[2]


  ## Set the necessary matrices
  data_dist0 <- data_dist <- delta_dummy <-
    matrix(NA, nrow = dim(combn(dim(data)[1], 2))[2], ncol = dim(data)[2])


  ## Obtain the distances and deltas per variable
  for (i in 1:dim(data)[2]) {

    ## Then, get the dissimilarities per variable
    data_dist0[, i] <-
      # Numeric
      if (typeof(data[, i]) == "double" &
          length(na.omit(unique(data[, i]))) > 1) {
        abs(apply(combn(data[, i], 2), 2, diff)) / diff(range(data[, i],
                                                              na.rm = TRUE))

        # Numeric (Same value for all trials)
      } else if (typeof(data[, i]) == "double" &
                 length(na.omit(unique(data[, i]))) == 1) {
        0

        # Nominal unordered (if binary, it is considered symmetric)
      } else if (typeof(data[, i]) != "double" &
                 length(class(data[, i])) == 1) {
        ifelse(apply(combn(data[, i], 2), 2, diff) == 0, 0, 1)

        # Nominal ordered
      } else if (typeof(data[, i]) != "double" &
                 length(class(data[, i])) == 2) {
        abs(apply(combn(rank(data[, i]), 2), 2, diff)) /
          diff(range(rank(data[, i]), na.rm = TRUE))
      }

    ## Finally, get the deltas per variable
    delta_dummy[, i] <-
      ifelse(!is.na(apply(combn(data[, i], 2), 2, diff)), 1, 0)
  }


  ## Gower's formula (weighted average of distances)
  data_dist <- if (max(data_dist0, na.rm = TRUE) == 0) {
    #stop("Dissimilarity matrix is zero for all comparisons.", call. = FALSE)
    0
  } else {
    (apply(data_dist0 * delta_dummy, 1, sum, na.rm = TRUE) /
       apply(delta_dummy, 1, sum))
  }


  ## Turn Gower dissimilarity into a lower triangle
  dist_mat <- matrix(NA, nrow = dim(data)[1], ncol = dim(data)[1])
  dist_mat[lower.tri(dist_mat, diag = FALSE)] <- data_dist
  diag(dist_mat) <- 0


  ## Remove redundant row and column
  #dist_mat <- as.data.frame(dist_mat0[-1, -dim(data)[1]])
  #rownames(dist_mat) <- input[-1, 1]
  #colnames(dist_mat) <- input[-dim(input)[1], 1]


  ## Percentage total missing data
  total_mod <-
    round((sum(is.na(input[, -c(1, 2)]) == TRUE) /
             (dim(input[, -c(1, 2)])[1] * dim(input[, -c(1, 2)])[2])) * 100, 2)


  ## Collect the results
  results <- list(Dissimilarity_table = dist_mat,
                  Types_used = char_type,
                  Total_missing = paste0(total_mod, "%"))

  return(results)
}
