#' Plot Gower's disimilarity values for each study (Transitivity evaluation)
#'
#' @description
#' Illustrating the range of Gower's dissimilarity values for each study in the
#' network, as well as their between- and within-comparison dissimilarities
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
#' @param strip_text_size A positive integer for the font size of facet labels.
#'   \code{strip_text_size} determines the strip.text argument found in the
#'   theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param label_size A positive integer for the font size of labels appearing on
#'   each study-specific segment. \code{label_size} determines the size argument
#'   found in the geom's aesthetic properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'
#' @return A horizontal bar plot illustrating the range of Gower's dissimilarity
#' values for each study with those found in other comparisons. The
#' study names appear on the y-axis in the order they appear in \code{results}
#' and the dissimilarity values appear on the x-axis. Red and blue points refer
#' to the (average) within-comparison and between-comparison dissimilarity,
#' respectively, for each study.
#'
#' @details
#' The range of Gower's dissimilarity values for each study versus the remaining
#' studies in the network for a set of clinical and methodological
#' characteristics that may act as effect modifiers. Gower's dissimilarities take
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
#' @examples
#' \donttest{
#' # Fictional dataset
#' data_set <- data.frame(Trial_name = as.character(1:7),
#'                       arm1 = c("1", "1", "1", "1", "1", "2", "2"),
#'                       arm2 = c("2", "2", "2", "3", "3", "3", "3"),
#'                       sample = c(140, 145, 150, 40, 45, 75, 80),
#'                       age = c(18, 18, 18, 48, 48, 35, 35),
#'                       blinding = factor(c("yes", "yes", "yes", "no", "no", "no", "no")))
#'
#' # Obtain comparison dissimilarities (informative = TRUE)
#' res <- comp_clustering(input = data_set,
#'                        drug_names = c("A", "B", "C"),
#'                        threshold = 0.13,  # General research setting
#'                        informative = TRUE,
#'                        get_plots = TRUE)
#'
#' plot_study_dissimilarities(results = res,
#'                            axis_title_size = 12,
#'                            axis_text_size = 12,
#'                            strip_text_size = 11,
#'                            label_size = 3.5)
#' }
#'
#' @export
plot_study_dissimilarities <- function(results,
                                       axis_title_size = 12,
                                       axis_text_size = 12,
                                       strip_text_size = 11,
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


  ## Turn diagonal into NA
  diag(diss) <- NA


  ## Distinguish between two-arm and multi-arm studies
  # Get the unique study ID (remove the compared treatments)
  #index <- gsub("^\\s+|\\s+$", "",
  #              sub("\\(.*", "", gsub('.{0}$', '', rownames(diss))))
  index <- sub("\\s+[^ ]+$", "", rownames(diss)) #gsub( " .*$", "", rownames(diss))

  # Split dataset by 'index'
  split_multi_arms <- split(diss, factor(index, levels = unique(index)))

  # Get the comparison for each study
  #comp_index <- gsub("^\\s+|\\s+$", "",
  #                   substring(rownames(diss), nchar(rownames(diss)) - 3))
  comp_index <- sub(".*\\s", "", rownames(diss)) # sub(".* ", "", rownames(diss))

  # Split 'diss' further by 'rownames(diss)'
  split_comp <- split(diss, factor(rownames(diss), levels = unique(rownames(diss))))

  # Split 'split_comp' further by 'comp_index'
  split_study_comp <- lapply(split_comp, function(x) split(x, comp_index))

  # Set of within-comparison similarities per study
  within_set <-
    lapply(1:length(split_study_comp),
           function(x)
             unlist(split_study_comp[[x]][is.element(names(split_study_comp[[x]]),
                                                     comp_index[x])]))

  # Within-comparison similarity per study
  within_comp0 <-
    sqrt(unlist(lapply(within_set, function(x) {sum(na.omit(x)^2) /
        (length(x) - sum(is.na(x)))})))

  # Use the minimum 'within_comp' among the comparisons of each multi-arm trial
  within_comp <-
    unlist(lapply(split(within_comp0, factor(index, levels = unique(index))),
                  function(x) min(na.omit(x))))

  # Set of between-comparison similarities per study
  between_set <-
    lapply(1:length(split_study_comp),
           function(x)
             unlist(split_study_comp[[x]][!is.element(names(split_study_comp[[x]]),
                                                      comp_index[x])]))

  # Between-comparison similarity per study
  between_comp0 <-
    sqrt(unlist(lapply(between_set, function(x) {sum(x^2) / length(x)})))

  # Use the minimum 'between_comp' among the comparisons of each multi-arm trial
  between_comp <-
    unlist(lapply(split(between_comp0, factor(index, levels = unique(index))),
                  function(x) min(x)))


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

  # Add comparison for each study
  data_raw$comp <- rep(comp_index, each =  dim(diss)[1])


  ## Prepare extra dataset for the 'within_comp' and 'between_comp'
  # Include 'within-comparison' similarities
  data_rms_within <- data.frame(var = unique(names(split_multi_arms)),
                                melt(within_comp),
                                study_id = 1:length(within_comp),
                                comp = comp_index)

  # Include 'between-comparison' similarities
  data_rms_between <- data.frame(var = unique(names(split_multi_arms)),
                                 melt(between_comp),
                                 study_id = 1:length(between_comp),
                                 comp = comp_index)


  ## Get the plot! :-)
  p1 <-
    ggplot(na.omit(data_raw),
           aes(x = value,
               y = study_id,
               xmin = min,
               xmax = max)) +
    facet_grid(comp ~ .,
               scales = "free",
               space = "free") +
    geom_crossbar(colour = "black",
                  linetype = "dotted",
                  width = 0.79) +
    scale_y_reverse(breaks = 1:length(unique(index)),
                    labels = unique(index),
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
    labs(x = "Range of Gower's dissimilarity values",
         y = "",
         fill = "Dissimilarity measure") +
    scale_x_continuous(limits = c(0, 1),
                       breaks = seq(0, 1, 0.1),
                       labels = sprintf("%.2f", seq(0, 1, 0.1)),
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
          strip.text = element_text(size = strip_text_size),
          legend.position = "bottom",
          legend.title = element_text(size = axis_title_size,
                                      face = "bold",
                                      colour = "black"),
          legend.text = element_text(size = axis_text_size))

  return(p1)
}
