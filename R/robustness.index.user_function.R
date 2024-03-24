#' Robustness index when 'metafor' or 'netmeta' are used
#'
#' @description
#'   Calculates the robustness index for a sensitivity analysis
#'   (Spineli et al., 2021) performed using the results of the analysis
#'   performed via the R-package
#'   \href{https://CRAN.R-project.org/package=netmeta}{netmeta} or
#'   \href{https://CRAN.R-project.org/package=metafor}{metafor}.
#'   The user defines the input and the function returns the robustness index.
#'
#' @param sens A list of R objects of class
#'   \code{\link[netmeta:netmeta]{netmeta}},
#'   \code{\link[netmeta:netmetabin]{netmetabin}} (see
#'   \href{https://CRAN.R-project.org/package=netmeta}{netmeta}) or
#'   \code{\link[metafor:rma]{rma}}, \code{\link[metafor:rma.glmm]{rma.glmm}},
#'   \code{\link[metafor:rma.mh]{rma.mh}}, \code{\link[metafor:rma.mv]{rma.mv}},
#'   \code{\link[metafor:rma.peto]{rma.peto}}, and
#'   \code{\link[metafor:rma.uni]{rma.uni}}
#'   (see \href{https://CRAN.R-project.org/package=metafor}{metafor}).
#'   The number of elements equals the number of analyses using the same dataset
#'   and the same R-package. The first element should refer to the primary
#'   analysis. Hence, the list should include at least two elements
#'   (see 'Details').
#' @param pkg Character string indicating the R-package with values
#'   \code{"netmeta"}, or \code{"metafor"}.
#' @param attribute This is relevant only for
#'   \href{https://CRAN.R-project.org/package=netmeta}{netmeta}. A vector of
#'   at least two characters with values \code{"TE.common"} or
#'   \code{"TE.random"}. See 'Values' in \code{\link[netmeta:netmeta]{netmeta}}
#'   or \code{\link[netmeta:netmetabin]{netmetabin}}.
#' @param threshold A number indicating the threshold of robustness, that is,
#'   the minimally allowed deviation between the primary analysis (the first
#'   element in \code{sens}) and re-analysis results. See 'Details' below.
#'
#' @return \code{robustness_index_user} prints on the R console a message in
#'   red text on the threshold of robustness determined by the user.
#'   Then, the function returns the following list of elements:
#'   \item{robust_index}{A numeric scalar or vector on the robustness
#'   index values. In the case of a pairwise meta-analysis,
#'   \code{robust_index} is scalar as only one summary effect size is obtained.
#'   In the case of network meta-analysis, \code{robust_index} is a vector with
#'   length equal to the number of possible pairwise comparisons;
#'   one robustness index per pairwise comparison.}
#'   \item{robust}{A character or character vector (of same length with
#'   \code{robust_index}) on whether the primary analysis results are
#'   \emph{robust} or \emph{frail} to the different re-analyses.}
#'   \item{kld}{A vector or matrix on the Kullback-Leibler divergence
#'   (KLD) measure in the summary effect size from a subsequent re-analysis to
#'   the primary analysis. In the case of a pairwise meta-analysis, \code{kld}
#'   is a vector with length equal to the number of total analyses (one KLD
#'   value is obtained per analysis). The number of total analyses equals the
#'   length of \code{sens}.
#'   In the case of network meta-analysis, \code{robust_index} is a matrix with
#'   number of rows equal to the number of total analyses and number of columns
#'   equal to the number of  possible pairwise comparisons; one KLD
#'   value per analysis and possible comparison.}
#'   \item{attribute}{The attributes considered.}
#'   \item{threshold}{The threshold used to be inherited by the
#'   \code{\link{heatmap_robustness}} function. See 'Details'.}
#'
#' @details Thresholds of robustness have been proposed only for the odds ratio
#'   and standardised mean difference (Spineli et al., 2021).
#'   The user may consider the values 0.28 and 0.17 in the argument
#'   \code{threshold} for the odds ratio and standardised mean difference effect
#'   measures (the default values), respectively, or consider other plausible
#'   values. When the argument \code{threshold} has not been defined,
#'   \code{robustness_index} considers the default values 0.28 and 0.17 as
#'   threshold for robustness for binary and continuous outcome, respectively,
#'   regardless of the effect measure (the default thresholds may not be proper
#'   choices for other effect measures; hence, use these threshold with great
#'   caution in this case). Spineli et al. (2021) offers a discussion on
#'   specifying the \code{threshold} of robustness.
#'
#'   When other effect measure is used (other than odds ratio or standardised
#'   mean difference) or the elements in \code{sens} refer to different effect
#'   measures, the execution of the function will be stopped and an error
#'   message will be printed in the R console.
#'
#'   In \code{robust}, the value \code{"robust"} appears when the calculated
#'   \code{robust_index} is less than \code{threshold}; otherwise, the value
#'   \code{"frail"} appears.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[metafor:rma]{rma}},
#'   \code{\link[metafor:rma.glmm]{rma.glmm}},
#'   \code{\link[metafor:rma.mh]{rma.mh}}, \code{\link[metafor:rma.mv]{rma.mv}},
#'   \code{\link[metafor:rma.peto]{rma.peto}},
#'   \code{\link[metafor:rma.uni]{rma.uni}},
#'   \code{\link[netmeta:netmeta]{netmeta}},
#'   \code{\link[netmeta:netmetabin]{netmetabin}},
#'   \code{\link{heatmap_robustness}}
#'
#' @references
#' Kullback S, Leibler RA. On information and sufficiency.
#' \emph{Ann Math Stat} 1951;\bold{22}(1):79--86. doi: 10.1214/aoms/1177729694
#'
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of
#' primary analysis results: A case study on missing outcome data in pairwise
#' and network meta-analysis.
#' \emph{Res Synth Methods} 2021;\bold{12}(4):475--90. doi: 10.1002/jrsm.1478
#'
#' @examples
#'
#' \dontrun{
#' library(netmeta)
#'
#' data(Baker2009)
#'
#' # Transform from arm-based to contrast-based format
#' p1 <- pairwise(treatment, exac, total, studlab = paste(study, year),
#' data = Baker2009, sm = "OR")
#'
#' # Conduct standard network meta-analysis
#' net1 <- netmeta(p1, ref = "Placebo")
#'
#' # Calculate the robustness index (random-effects versus fixed-effect)
#' robustness_index_user(sens = list(net1, net1),
#'                       pkg = "netmeta",
#'                       attribute = c("TE.random", "TE.common"),
#'                       threshold = 0.28)
#' }
#'
#' @export
robustness_index_user <- function(sens, pkg, attribute, threshold) {

  pkg <- if (missing(pkg) || !is.element(pkg, c("netmeta", "metafor"))) {
    stop("Insert 'netmeta' or 'metafor'", call. = FALSE)
  } else {
    pkg
  }

  # Check the class of the elements
  metafor_classes <-
    c("rma", "rma.glmm", "rma.mh", "rma.mv", "rma.peto", "rma.uni")
  netmeta_classes <- c("netmeta", "netmetabin")
  get_class <- unique(unlist(lapply(sens, class)))
  class_check <- if (pkg == "netmeta" &
                     any(!is.element(get_class, netmeta_classes))) {
    a <- "'netmeta', or 'netmetabin'."
    stop(paste("The elements must be objects of class", a), call. = FALSE)
  } else if (pkg == "metafor" & any(!is.element(get_class, metafor_classes))) {
    a <- "'rma', 'rma.glmm', 'rma.mh', 'rma.mv', 'rma.peto', or 'rma.uni'."
    stop(paste("The elements must be objects of class", a), call. = FALSE)
  }

  # Check whether 'sens' is a list with elements from the correct class
  sens <- if(!is.list(sens) || length(sens) < 2) {
    a <- "'metafor', or 'netmeta'."
    stop(paste("A list of at least two objects from", a), call. = FALSE)
  } else if (is.list(sens) & length(sens) > 1) {
    sens
  }

  # Get to know the class of the elements in 'attribute'
  attribute <- if (missing(attribute) & pkg == "netmeta") {
    stop("The argument 'attribute' needs to be defined.", call. = FALSE)
  } else if (pkg == "netmeta" & length(attribute) != length(sens)) {
    a <- "the number of elements in 'sens'."
    stop(paste("The length of 'attribute' does not equal", a), call. = FALSE)
  } else if (pkg == "netmeta" &
             any(!is.element(attribute, c("TE.common", "TE.random")))) {
    stop("Possible values are 'TE.common' and 'TE.random'", call. = FALSE)
  } else {
    attribute
  }

  # Number of analyses (i.e., primary + (length(sens) - 1))
  n_scenar <- length(sens)

  # Check the summary effect in element of 'sens'
  measure_check <- rep(NA, length(sens))
  for (i in 1:length(sens)) {
    if (pkg == "netmeta") {
      measure_check[i] <- sens[[i]]$sm
    } else if (pkg == "metafor") {
      measure_check[i] <- sens[[i]]$measure
    }
  }

  # Check whether the same summary effect is used in 'sens'
  measure <- if (length(unique(measure_check)) > 1 ||
                 any(!is.element(unique(measure_check), c("OR", "SMD")))) {
    stop("Use odds ratio or standardised mean difference as effect measure.",
         call. = FALSE)
  } else {
    unique(measure_check)
  }

  # Obtain the dataset (summary effect and std errors) for all scenarios
  treat_effect0 <- treat_effect <- std_error0 <- std_error <- list()
  num_elem_te <- num_elem_se <- rep(NA, length(sens))
  vect_effect <- vect_stderr <- rep(NA, length(sens))
  if(pkg == "netmeta") {
    for (i in 1:length(sens)) {
      # Find the position of attribute for the corresponding element in 'sens'
      num_elem_te[i] <- match(attribute[i], attributes(sens[[i]])$names)
      num_elem_se[i] <- match(paste0("se", attribute[i]),
                              attributes(sens[[i]])$names)

      # A list of summary effects per scenario
      treat_effect0[[i]] <- as.data.frame(sens[[i]][num_elem_te[i]])
      treat_effect[[i]] <- treat_effect0[[i]][lower.tri(treat_effect0[[i]])]

      # A list of summary standard errors per scenario
      std_error0[[i]] <- as.data.frame(sens[[i]][num_elem_se[i]])
      std_error[[i]] <- std_error0[[i]][lower.tri(std_error0[[i]])]
    }
    # A k-by-r matrix of results for all scenarios (k) and comparisons (r)
    mean_mat <- do.call(rbind, treat_effect) # Summary effects
    se_mat <- do.call(rbind, std_error)      # Summary standard errors
  } else if (pkg == "metafor") {
    for (i in 1:length(sens)) {
      # Find the position of attribute for the corresponding element in 'sens'
      num_elem_te[i] <- match(attribute[i], attributes(sens[[i]])$names)
      num_elem_se[i] <- match("se", attributes(sens[[i]])$names)

      # A vector of results for all scenarios
      vect_effect[i] <- sens[[i]][num_elem_te[i]] # Summary effects
      vect_stderr[i] <- sens[[i]][num_elem_se[i]] # Summary standard errors
    }
    # Turn into a k-by-1 matrix
    mean_mat <- as.matrix(vect_effect, nrow = length(sens), ncol = 1)
    se_mat <- as.matrix(vect_stderr, nrow = length(sens), ncol = 1)
  }

  # Check the number of comparisons per element. Relevant only for 'netmeta'
  poss_comp_check <- rep(NA, length(sens))
  if (pkg == "netmeta") {
    for (i in 1:length(sens)) {
      poss_comp_check[i] <- length(treat_effect[[i]])
    }
  }

  # Check whether the same network is used (based on the number of comparisons).
  # Relevant only for 'netmeta'
  poss_comp0 <- if (pkg == "netmeta" & length(unique(poss_comp_check)) > 1) {
    stop("Use the same network (same set of interventions) for all analyses.",
         call. = FALSE)
  } else if (pkg == "netmeta" & length(unique(poss_comp_check)) == 1) {
    unique(poss_comp_check)
  }

  # Number of possible pairwise comparisons (one in pairwise meta-analysis)
  poss_comp <- if(pkg == "netmeta") {
    poss_comp0
  } else if (pkg == "metafor") {
    1
  }

  # The first scenario is *always* the primary scenario
  primary_scenar <- 1

  # Check the selected 'threshold'
  if (missing(threshold) & measure == "OR") {
    threshold <- 0.28
    message("The value 0.28 was assigned as 'threshold' by default.")
  } else if (missing(threshold) & measure == "SMD") {
    threshold <- 0.17
    message("The value 0.17 was assigned as 'threshold' by default.")
  } else {
    threshold <- threshold
    aa <- "was assigned as 'threshold' for"
    effect_measure <- effect_measure_name(measure, lower = TRUE)
    message(paste("The value", threshold, aa, paste0(effect_measure, ".")))
  }

  # Function for the Kullback-Leibler Divergence (two normal distributions)
  kld_measure_univ <- function(mean_y, sd_y, mean_x, sd_x) {
    # x is the 'truth' (e.g. the MAR assumption)
    kld_xy <- 0.5 * (((sd_x / sd_y)^2) + ((mean_y - mean_x)^2)
                     / (sd_y^2) - 1 + 2 * log(sd_y / sd_x))

    return(kld_xy)
  }

  # Calculate the Kullback-Leibler Divergence and robustness index
  kldxy <- list()
  robust_index <- rep(NA, poss_comp)
  for (i in 1:poss_comp) {
    kldxy[[i]] <- rep(NA, n_scenar)
    for (j in 1:n_scenar) {
      # Returns the KLD of scenario j compared with primary analysis
      # for comparison i
      kldxy[[i]][j] <- kld_measure_univ(mean_mat[j, i],
                                        se_mat[j, i],
                                        mean_mat[primary_scenar, i],
                                        se_mat[primary_scenar, i])
    }
    # This refers to the primary analysis
    kldxy[[i]][primary_scenar] <- 0

    # Returns the Robustness Index of comparison i across all scenarios
    robust_index[i] <- round(sqrt(t(kldxy[[i]]) %*% kldxy[[i]]), 3)
  }
  kld <- matrix(unlist(kldxy), ncol = n_scenar, byrow = TRUE)
  robust <- ifelse(robust_index < threshold, "robust", "frail")

  # Collect results in a list
  results <- list(robust_index = robust_index,
                  robust = robust,
                  kld = kld,
                  measure = measure,
                  threshold = threshold,
                  attribute = attribute)

  class(results) <- "robustness_index_user"

  return(results)
}
