#' Visualising the distribution of characteristics
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{distr_characteristics} uses violin and bar plots to visualise the
#'   distribution of each characteristic in the dataset either per comparison
#'   or cluster of comparisons.
#'
#' @param input A data-frame in the long arm-based format. Two-arm trials occupy
#'   one row in the data-frame. Multi-arm trials occupy as many rows as the
#'   number of possible comparisons among the interventions. The first three
#'   columns refer to the trial name, first and second arm of the comparison
#'   (their identifier number), respectively. The remaining columns refer to
#'   summary characteristics. See 'Details' for specifying the columns.
#' @param drug_names A vector of labels with the name of the interventions
#'   in the order they appear in the argument \code{input}.
#' @param rename_char A list of two elements: (i) a numeric vector with the
#'   position of the characteristics in \code{input}, and (ii) a character
#'   vector with the names of the characteristics, as they are wished to
#'   appear in the title of the plots. This argument is optional, in case the
#'   user wants to control the appearance of the titles.
#' @param cluster An object of S3 class \code{\link{comp_clustering}} that has
#'   information on the cluster of each comparison. See 'Value' in
#'   \code{\link{comp_clustering}}. If \code{cluster} is not provided, the
#'   function presents the distribution of characteristics per comparison;
#'   otherwise per cluster. In the latter, the function prints a table with the
#'   comparisons and the corresponding cluster.
#' @param label_size A positive integer for the font size of labels in the
#'   plots. \code{label_size} determines the size argument found in the geom's
#'   aesthetic properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param title_size A positive integer for the font size of legend title in
#'   the stacked barplot on the percentage studies of each comparison found in
#'   the clusters. \code{title_size} determines the title argument
#'   found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param axis_title_size A positive integer for the font size of axis title in
#'   the plots. \code{axis_title_size} determines the axis.title argument found
#'   in the theme's properties in the
#'   R-package \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param axis_text_size A positive integer for the font size of axis text in
#'   the plots. \code{axis_text_size} determines the axis.text argument found in
#'   the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param axis_x_text_angle A positive integer for the angle of axis text in
#'   the plots. \code{axis_text_angle} determines the axis.text.x argument found
#'   in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param legend_text_size A positive integer for the font size of legend text
#'   in the plots. \code{legend_text_size} determines the legend.text argument
#'   found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'
#' @return
#'   \code{distr_characteristics} returns a list of plots using the proper plot
#'   (violin or bar plot) for each characteristic. The size of the dots in the
#'   violin plot (with amalgamated box plots and dots) are proportional to the
#'   total sample size of the study: the large the sample size of the study, the
#'   larger the size of the corresponding point.
#'
#' @details
#'   The correct type mode of columns in \code{input} must be ensured to use
#'   the function \code{distr_characteristics}. The first three columns
#'   referring to the trial name, first and second arm of the comparison,
#'   respectively, must be \strong{character}. The remaining columns referring
#'   to the characteristics must be \strong{double} or \strong{integer}
#'   depending on whether the corresponding characteristic refers to a
#'   quantitative or qualitative variable. The type mode of each column is
#'   assessed by \code{distr_characteristics} using the base function
#'   \code{typeof}.
#'
#'   The interventions should be sorted in an ascending order of their
#'   identifier number within the trials so that the first treatment column
#'   (second column in \code{input}) is the control arm for every pairwise
#'   comparison. This is important to ensure consistency in the order of
#'   interventions within the comparisons obtained from the other related
#'   functions.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{comp_clustering}}
#'
#' @export
distr_characteristics <- function (input,
                                   drug_names,
                                   rename_char = NULL,
                                   cluster = NULL,
                                   label_size = 4,
                                   title_size = 14,
                                   axis_title_size = 14,
                                   axis_text_size = 14,
                                   axis_x_text_angle = 0,
                                   legend_text_size = 13) {


  ## Check defaults
  # Dataset
  input0 <- if (any(sapply(input, typeof)[1:3] != "character")) {
    stop("The first three columns (trial and arms) must be 'characters'.",
         call. = FALSE)
  } else if (any(sapply(input, typeof)[-c(1:3)] == "character")) {
    stop("The characteristics must be 'double' or 'integer'.", call. = FALSE)
  } else {
    input
  }
  colnames(input0)[1:3] <- c("Trial_name", "Arm1", "Arm2")

  # Intervention names
  drug_names <- if (missing(drug_names)) {
    as.character(1:length(unique(unlist(input0[, 2:3]))))
  } else {
    drug_names
  }

  # Clustered comparisons
  if (!is.null(cluster) & !inherits(cluster, "comp_clustering")) {
    stop("'cluster' must be an object of S3 class 'comp_clustering'.",
         call. = FALSE)
  }


  ## Create new input and name the treatments
  input <- input0
  input[, 2:3] <- matrix(drug_names[as.numeric(unlist(input[, 2:3]))],
                         nrow = dim(input)[1],
                         ncol = 2)


  ## Insert 'Comparison' in the dataset (control appears second in the compar.)
  input$Comparison <- as.character(paste0(input$Arm2, "-", input$Arm1))


  ## Reduce dataset to trial, comparison & characteristics
  input_new0 <- input[, c(1, dim(input)[2], 4:(dim(input)[2] - 1))]


  ## Split 'dataset' by 'Comparison'
  split_dataset0 <- split(input_new0, f = input$Comparison)


  ## Find the completely missing columns in all non-single-study comparisons
  col_all_miss <-
    unique(unlist(
      lapply(split_dataset0, function(x) if (dim(x)[1] > 1)
        as.vector(which(colSums(is.na(x)) == nrow(x) |
                          colSums(is.na(x)) == nrow(x) - 1))))) #as.vector(which(colSums(is.na(x)) == nrow(x))))))


  ## Keep the names of the completely missing columns in all comparisons
  col_all_miss_names <-
    unique(unlist(lapply(split_dataset0,
                         function(x) colnames(x)[col_all_miss])))


  ## Message on the dropped characteristics
  dropped_chars <- ifelse(length(col_all_miss_names) == 0,
                          "none",
                          col_all_miss_names)
  message(paste("- Dropped characteristics:", paste(dropped_chars,
                                                    collapse = ", ")))


  ## Now remove these columns for *all* comparisons!
  split_dataset <- min_split <- max_split <-
    lapply(split_dataset0, function(x) x[!names(x) %in% col_all_miss_names])


  ## Remove these columns also from the dataset for the moment
  input_new <-
    if (length(col_all_miss) > 0) {
      subset(input_new0, select = -col_all_miss)
    } else {
      input_new0
    }


  ## Bind all lists into a data-frame
  dataset_new <- do.call(rbind, split_dataset)


  ## Rename columns if indicated
  if (!is.null(rename_char)) {
    colnames(dataset_new)[rename_char[[1]] - 1] <- rename_char[[2]]
  }


  ## Variable on sample size
  colnames(dataset_new)[with(dataset_new,
                             startsWith(names(dataset_new),
                                        c("sample", "Sample")))] <- "Sample size"


  ## Function for first letter capital (Source: https://stackoverflow.com/questions/18509527/first-letter-to-upper-case)
  firstup <- function(x) {
    substr(x, 1, 1) <- toupper(substr(x, 1, 1))
    x
  }


  ## Visualise characteristics per cluster (!is.null(cluster)) or comparison
  if (!is.null(cluster)) {

    ## Comparisons with their cluster
    clustered_comp <- cluster$Cluster_comp


    ## Include a column with the cluster of the comparisons
    # Copy-paste the trials to a new column
    dataset_new$`Studies cluster` <- dataset_new$Trial_name

    # Match the study with the cluster
    for (i in 1:dim(dataset_new)[1]) {
      dataset_new$`Studies cluster`[
        dataset_new$`Studies cluster` == clustered_comp[i, 1]] <-
        clustered_comp[i, 2]
    }


    ## Visualise distribution according to characteristic type
    # Double type
    double_type <- function (yvar) {
      ggplot(subset(dataset_new, !is.na(dataset_new[, yvar])),
             aes_(x = ~`Studies cluster`,
                  y = as.name(yvar))) +
        geom_violin(fill = "#99CCFF",
                    trim = FALSE,
                    alpha = 0.3) +
        geom_boxplot(outlier.alpha = 0.3,
                     fill = "white",
                     colour = "black",
                     varwidth = TRUE) +
        geom_point(aes_(size = ~`Sample size`)) +
        stat_boxplot(geom = 'errorbar',
                     width = 0.2,
                     linetype = "dashed") +
        labs(x = "Clusters",
             y = " ") +
        guides(size = FALSE,
               colour = guide_legend(override.aes = list(size = 3.5))) +
        theme_classic() +
        ggtitle(as.name(yvar)) +
        theme(plot.title = element_text(size = title_size, face = "bold"),
              axis.title = element_text(size = axis_title_size, face = "bold"),
              axis.text = element_text(size = axis_text_size),
              axis.text.x = element_text(angle = axis_x_text_angle,
                                         hjust =
                                           ifelse(axis_x_text_angle == 0, 0.5, 1)
              ),
              legend.position = "bottom",
              legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = title_size, face = "bold"),
              plot.caption = element_text(size = 10, hjust = 0.0))
    }

    # Integer type
    factor_type <- function (yvar) {

      # Turn fist letter into capital
      levels(dataset_new[, yvar]) <-
        firstup(as.character(sort(unique(dataset_new[, yvar]))))

      # Get the bar plot
      ggplot(subset(dataset_new, !is.na(dataset_new[, yvar])),
             aes_(x = ~`Studies cluster`)) +
        geom_bar(stat = "count",
                 aes_(fill = as.name(yvar))) +
        geom_text(data = data.frame(prop.table(table(dataset_new[, "Studies cluster"],
                                                     dataset_new[, yvar]),
                                               margin = 1),
                                    count = table(dataset_new[, "Studies cluster"],
                                                  dataset_new[, yvar])),
                  inherit.aes = FALSE,
                  aes_(x = ~Var1,
                       y = ~count.Freq,
                       group = ~Var2,
                       label = ~ifelse(Freq != 0,
                                       paste0(round(Freq * 100, 1), "% (",
                                              count.Freq,")"), " ")),
                  hjust = 0.5,
                  vjust = 1.0,
                  size = label_size,
                  position = "stack",
                  colour = "white") +
        labs(x = "Clusters",
             y = "Count",
             fill = "Categories") +
        theme_classic() +
        ggtitle(as.name(yvar)) +
        theme(plot.title = element_text(size = title_size, face = "bold"),
              axis.title = element_text(size = axis_title_size, face = "bold"),
              axis.text = element_text(size = axis_text_size),
              axis.text.x = element_text(angle = axis_x_text_angle,
                                         hjust =
                                           ifelse(axis_x_text_angle == 0, 0.5, 1)
              ),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(legend_text_size),
              legend.title = element_text(size = title_size, face = "bold"))
    }


    ## Remind the users which comparisons belong to each cluster
    message(paste0(capture.output(
      knitr::kable(clustered_comp[order(clustered_comp$cluster),],
                   align = "ll",
                   caption = "Comparisons with their clusters")),
      collapse = "\n"))

    ## Report the size of the clusters
    message(" ")
    message(do.call(cbind,
            lapply(sort(unique(clustered_comp$cluster)),
                   function(x)
                     paste0("Cluster", " ", x, ": ",
                            round(
                              prop.table(
                                table(dataset_new$`Studies cluster`))[x] *
                                100, 1), "%", " "))))

  } else {

    ## Visualise distribution according to characteristic type
    # Double type
    double_type <- function (yvar) {
      ggplot(subset(dataset_new, !is.na(dataset_new[, yvar])),
             aes_(x = ~Comparison,
                  y = as.name(yvar))) +
        geom_violin(fill = "#99CCFF",
                    trim = FALSE,
                    alpha = 0.3) +
        geom_boxplot(outlier.alpha = 0.3,
                     fill = "white",
                     colour = "black",
                     varwidth = TRUE) +
        geom_point(aes_(size = ~`Sample size`)) +
        stat_boxplot(geom = 'errorbar',
                     width = 0.2,
                     linetype = "dashed") +
        labs(x = " ",
             y = " ") +
        guides(size = FALSE,
               colour = guide_legend(override.aes = list(size = 3.5))) +
        theme_classic() +
        ggtitle(as.name(yvar)) +
        theme(plot.title = element_text(size = title_size, face = "bold"),
              axis.title = element_text(size = axis_title_size, face = "bold"),
              axis.text = element_text(size = axis_text_size),
              axis.text.x = element_text(angle = axis_x_text_angle,
                                         hjust =
                                           ifelse(axis_x_text_angle == 0, 0.5, 1)
              ),
              legend.position = "bottom",
              legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = title_size, face = "bold"),
              plot.caption = element_text(size = 10, hjust = 0.0))
    }

    # Integer type
    factor_type <- function (yvar) {

      # Turn fist letter into capital
      levels(dataset_new[, yvar]) <-
        firstup(as.character(sort(unique(dataset_new[, yvar]))))

      # Get the bar plot
      ggplot(subset(dataset_new, !is.na(dataset_new[, yvar])),
             aes_(x = ~Comparison)) +
        geom_bar(stat = "count",
                 aes_(fill = as.name(yvar))) +
        geom_text(data = data.frame(prop.table(table(dataset_new[, 2],
                                                     dataset_new[, yvar]),
                                               margin = 1),
                                    count = table(dataset_new[, 2],
                                                  dataset_new[, yvar])),
                  inherit.aes = FALSE,
                  aes_(x = ~Var1,
                       y = ~count.Freq,
                       group = ~Var2,
                       label = ~ifelse(Freq != 0,
                                       paste0(round(Freq * 100, 1), "% (",
                                              count.Freq,")"), " ")),
                  hjust = 0.5,
                  vjust = 1.0,
                  size = label_size,
                  position = "stack",
                  colour = "white") +
        labs(x = " ",
             y = "Count",
             fill = "Categories") +
        theme_classic() +
        ggtitle(as.name(yvar)) +
        theme(plot.title = element_text(face = "bold"),
              axis.title = element_text(size = axis_title_size, face = "bold"),
              axis.text = element_text(size = axis_text_size),
              axis.text.x = element_text(angle = axis_x_text_angle,
                                         hjust =
                                           ifelse(axis_x_text_angle == 0, 0.5, 1)
              ),
              legend.position = "bottom",
              legend.direction = "horizontal",
              legend.text = element_text(size = legend_text_size),
              legend.title = element_text(size = title_size, face = "bold"))
    }

  }

  ## List of graphs and suppressing warning
  suppressWarnings({

    suppressWarnings(
      plots <-
        lapply(names(dataset_new[-c(1, 2)]), function(x)
          if(typeof(dataset_new[, x]) == "double") double_type(x) else
            factor_type(x)))
    names(plots) <- colnames(dataset_new)[-c(1, 2)]

    return(plots)
 })

}
