#' Rainbow of Gower's similarity values for each study (Transitivity evaluation)
#'
#' @description
#' Illustrating the range of Gower's similarity values for each study in the
#' network.
#'
#' @param results An object of S3 class \code{\link{comp_clustering}}.
#'   See 'Value' in \code{\link{comp_clustering}}.
#' @param axis_title_size A positive integer for the font size of axis title
#'   (both axes). \code{axis_title_size} determines the axis.title argument
#'   found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param axis_text_size A positive integer for the font size of axis text (both
#'   axes). \code{axis_text_size} determines the axis.text argument found in the
#'   theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param label_size A positive integer for the font size of labels appearing on
#'   each study-specific segment. \code{label_size} determines the size argument
#'   found in the geom's aesthetic properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'
#' @return A horizontal bar plot illustrating the range of Gower's similarity
#' values for each study with those found in other comparisons using shades of
#' red and green to indicate low and substantial between-comparison similarity,
#' respectively: the darker the red, the lower the between-comparison similarity
#' (corresponding to values close to 0), whilst the darker the green, the higher
#' the between-comparison similarity (corresponding to values close to 1). The
#' study names appear on the y-axis in the order they appear in \code{results}
#' and the similarity values appear on the x-axis. Red and blue points refer to
#' the (average) within-comparison and between-comparison similarity,
#' respectively, for each study.
#'
#' @details
#' The range of Gower's similarity values for each study result from calculating
#' the Gower's dissimilarity of a study versus the remaining studies in the
#' network for a set of clinical and methodological characteristics that may act
#' as effect modifiers. Then, the Gower's dissimilarities are transformed into
#' similarities by subtracting each value from 1: Gower's dissimilarities take
#' values from 0 to 1, with 0 and 1 implying perfect similarity and perfect
#' dissimilarity, respectively.
#'
#' The unique similarity values appear as dotted, vertical, black lines on each
#' bar.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{comp_clustering}}
#'
#' @references
#' Gower J. General Coefficient of Similarity and Some of Its Properties.
#' \emph{Biometrics} 1971;\bold{27}(4):857--71.
#' doi: 10.2307/2528823
#'
#' @export
rainbow_similarities <- function(results,
                                 axis_title_size = 12,
                                 axis_text_size = 12,
                                 label_size = 3.5) {


  ## Check default
  if (inherits(results, "comp_clustering") == FALSE) {
    stop("The argument must be an object of S3 class 'comp_clustering'",
         call. = FALSE)
  }


  ## Consider the matrix of study dissimilarities
  diss <- results$Trials_diss_table


  ## Copy-paste the lower off-diagonal elements to the corresponding upper
  diss[upper.tri(diss)] <- t(diss)[upper.tri(diss)]


  ## Turn into weights (1 - diss)
  weight <- 1 - diss


  ## Turn diagonal into NA
  diag(weight) <- NA


  ## Distinguish between two-arm and multi-arm studies
  # Get the unique study ID (remove the compared treatments)
  index <- gsub("^\\s+|\\s+$", "",
                sub("\\(.*", "", gsub('.{3}$', '', rownames(weight))))

  # Split dataset by 'index'
  split_multi_arms <- split(weight, factor(index, levels = unique(index)))

  # Get the comparison for each study
  comp_index <- gsub("^\\s+|\\s+$", "",
                     substring(rownames(weight), nchar(rownames(weight)) - 3))

  # Split 'weight' further by 'rownames(weight)'
  split_comp <- split(weight, factor(rownames(weight), levels = unique(rownames(weight))))

  # Split 'split_comp' further by 'comp_index'
  split_study_comp <- lapply(split_comp, function(x) split(x, comp_index))

  # Set of within-comparison similarities per study
  within_set <-
    lapply(1:length(split_study_comp),
           function(x)
             unlist(split_study_comp[[x]][is.element(names(split_study_comp[[x]]),
                                                     comp_index[x])]))

  # Within-comparison similarity per study
  within_comp <-
    sqrt(unlist(lapply(within_set, function(x) {sum(na.omit(x)^2) /
        (length(x) - sum(is.na(x)))})))

  # Set of between-comparison similarities per study
  between_set <-
    lapply(1:length(split_study_comp),
           function(x)
             unlist(split_study_comp[[x]][!is.element(names(split_study_comp[[x]]),
                                                      comp_index[x])]))

  # Between-comparison similarity per study
  between_comp <-
    sqrt(unlist(lapply(between_set, function(x) {sum(x^2) / length(x)})))


  ## Transform the weights
  # Get the lower bound for the uniform distribution
  lower <- aggregate(unlist(lapply(between_set, function(x) min(x))),
                     by = list(factor(index, levels = unique(index))), min)[, 2]

  # Get the upper bound for the uniform distribution
  upper <- aggregate(unlist(lapply(between_set, function(x) max(x))),
                     by = list(factor(index, levels = unique(index))), min)[, 2]

  # Bring together
  tranform_result <- data.frame(lower, upper); rownames(tranform_result) <- unique(index)


  ## Create a sequence of values to define the colour shades
  # Prepare dataset
  vals <- lapply(1:dim(tranform_result)[1],
                 function(x) seq(tranform_result$lower[x],
                                 tranform_result$upper[x],
                                 0.001))

  # Repeat the study id
  mid <- rep(1:length(split_multi_arms), lengths(vals))

  # Bring all together
  data_set1 <-
    data.frame(y = mid - 0.4,
               yend = mid + 0.4,
               x = unlist(vals),
               xend = unlist(vals))


  ## Raw data to be plotted in 'geom_crossbar'
  # Prepare dataset
  data_raw <- melt(split_multi_arms); colnames(data_raw)[2] <- "study"

  # Define the following variables
  study_id <- xend <- yend <- NULL

  # Add the study id
  data_raw$study_id <-
    unlist(lapply(1:length(split_multi_arms),
                  function(x) rep(x, length(split_multi_arms[[x]]))))

  # Add the bounds
  data_raw$min <- data_raw$value; data_raw$max <- data_raw$value


  ## Prepare extra dataset for the 'within_comp' and 'between_comp'
  # Include 'within-comparison' similarities
  data_rms_within <- data.frame(var = unique(names(split_multi_arms)),
                                melt(within_comp),
                                study_id = 1:length(within_comp))

  # Include 'between-comparison' similarities
  data_rms_between <- data.frame(var = unique(names(split_multi_arms)),
                                 melt(between_comp),
                                 study_id = 1:length(between_comp))


  ## Get the plot! :-)
  p1 <-
    ggplot(data = data_set1,
           aes(x = x,
               xend = xend,
               y = y,
               yend = yend,
               color = x)) +
      geom_segment(linewidth = 2) +
      scale_color_gradient2(low = "#F98866", mid = "white", high = "#A1BE95",
                            midpoint = max(data_set1$x)/2) +
      geom_crossbar(data = na.omit(data_raw),
                    aes(x = value,
                        y = study_id,
                        xmin = min,
                        xmax = max),
                    colour = "black",
                    linetype = "dotted",
                    width = 0.79,
                    inherit.aes = FALSE) +
    scale_y_reverse(breaks = 1:dim(tranform_result)[1],
                    labels = rownames(tranform_result),
                    expand = c(0, 0)) +
    geom_point(data = data_rms_within,
               aes(x = value,
                   y = study_id,
                   fill = "Within-comparison"),
               color = "red",
               size = 3.0,
               shape = "diamond",
               inherit.aes = FALSE) +
    geom_point(data = data_rms_between,
               aes(x = value,
                   y = study_id,
                   fill = "Between-comparison"),
               color = "blue",
               size = 3.0,
               shape = "diamond",
               inherit.aes = FALSE) +
    geom_text(data = data_rms_within,
              aes(x = value,
                  y = study_id,
                  label = sprintf("%.2f", value)),
              size = label_size,
              vjust = 0.5, # -0.85
              hjust = 1.25,
              inherit.aes = FALSE) +
    geom_text(data = data_rms_between,
              aes(x = value,
                  y = study_id,
                  label = sprintf("%.2f", value)),
              size = label_size,
              vjust = 0.5, # -0.85
              hjust = 1.25,
              inherit.aes = FALSE) +
    labs(x = "Range of Gower's similarity values",
         y = "",
         fill = "Similarity type") +
    scale_x_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.1),
                       labels = sprintf("%.2f", seq(0, 1, 0.1)),
                       #labels = scales::percent,
                       expand = c(0.01, 0)) +
    scale_fill_manual(values = c("Within-comparison" = "red",
                                 "Between-comparison" = "blue")) +
    guides(colour = "none") +
    theme_classic() +
    theme(axis.title = element_text(size = axis_title_size,
                                    face = "bold",
                                    colour = "black"),
          axis.text = element_text(size = axis_text_size),
          panel.border = element_blank(),
          legend.position = "bottom",
          legend.title = element_text(size = axis_title_size,
                                      face = "bold",
                                      colour = "black"),
          legend.text = element_text(size = axis_text_size))

  return(p1)
}
