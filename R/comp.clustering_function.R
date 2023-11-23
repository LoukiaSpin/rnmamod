#' End-user-ready results for informative and heuristic clustering of comparisons
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{comp_clustering} hosts a toolkit of functions that facilitates
#'   conducting, visualising and evaluating informative and heuristic (hierarchical
#'   agglomerative) clustering of observed comparisons of interventions for a
#'   specific network and set of characteristics.
#'
#' @param input A data-frame in the long arm-based format. Two-arm trials occupy
#'   one row in the data-frame. Multi-arm trials occupy as many rows as the
#'   number of possible comparisons among the interventions. The first three
#'   columns refer to the trial name, first and second arm of the comparison,
#'   respectively. The remaining columns refer to summary characteristics. See
#'   'Details' for the specification of the columns.
#' @param drug_names A vector of labels with the name of the interventions
#'   in the order they have been defined in the argument \code{input}.
#' @param rename_char A list of two elements: (i) a numeric vector with the
#'   position of the characteristics in \code{input}, and (ii) a character
#'   vector with the names of the characteristics, as they are wished to
#'   appear in the title of the plots. This argument is optional, in case the
#'   user wants to control the appearance of the titles.
#' @param height Logical with \code{TRUE} for performing informative clustering
#'   and \code{FALSE} for performing heuristic clustering, thus, allowing the
#'   user to define the number of clusters via the argument \code{optimal_clusters}.
#'   The default argument is \code{FALSE}.
#' @param num_neighb A positive integer for the number of neighbouring
#'   comparisons that is found in the \code{\link{connectivity_index}} and
#'   \code{\link{internal_measures_plot}} functions. It takes values from two to
#'   the number of comparisons minus one. The default argument equals half the
#'   number of observed comparisons.
#' @param optimal_clusters A positive integer for the optimal number of clusters
#'   based on three internal validation measures. The default argument is 2,
#'   but the user \bold{must} adjust the value, if needed, after inspecting the
#'   results of the internal validation measures. It takes values from two to
#'   the number of comparisons minus one.
#' @param get_plots Logical with values \code{TRUE} for returning all plots and
#'   \code{FALSE} for concealing the plots. The default argument is
#'   \code{FALSE}.
#' @param label_size A positive integer for the font size of labels in the
#'   bubble plot for the characteristics contribution, the violin plot for the
#'   study dissimilarities per comparison, the bar plot for total dissimilarity
#'   and the silhouette widths plot. \code{label_size} determines the size
#'   argument found in the geom's aesthetic properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param title_size A positive integer for the font size of legend title in
#'   the bubble plot for the characteristics contribution, the violin plot for
#'   the study dissimilarities per comparison, the bar plot for total
#'   dissimilarity and the silhouette widths plot. \code{title_size} determines
#'   the title argument found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param axis_title_size A positive integer for the font size of axis title in
#'   the bubble plot for the characteristics contribution, the violin plot for
#'   the study dissimilarities per comparison, the bar plot for total
#'   dissimilarity and the silhouette widths plot. \code{axis_title_size}
#'   determines the axis.title argument found in the theme's properties in the
#'   R-package \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param axis_text_size A positive integer for the font size of axis text in
#'   the bubble plot for the characteristics contribution, the violin plot for
#'   the study dissimilarities per comparison, the bar plot for total
#'   dissimilarity and the silhouette widths plot. \code{axis_text_size}
#'   determines the axis.text argument found in the theme's properties in the
#'   R-package \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param axis_x_text_angle A positive integer for the angle of axis text in
#'   the bubble plot for the characteristics contribution, the violin plot for
#'   the study dissimilarities per comparison, and the bar plot for total
#'   dissimilarity. \code{axis_x_text_angle} determines the axis.text.x argument
#'   found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param legend_text_size A positive integer for the font size of legend text
#'   in the bubble plot for the characteristics contribution, the violin plot
#'   for the study dissimilarities per comparison, the bar plot for total
#'   dissimilarity and the silhouette widths plot. \code{legend_text_size}
#'   determines the legend.text argument found in the theme's properties in the
#'   R-package \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#'
#' @return
#'   Initially, \code{comp_clustering} prints on the console the following
#'   messages: the number of observed comparisons (and number of single-study
#'   comparisons, if any); the number of dropped characteristics due to many
#'   missing data; the maximum value of the cophenetic correlation coefficient;
#'   the optimal dissimilarity measure and optimal linkage method selected
#'   based on the cophenetic correlation coefficient. Then, the function
#'   returns the following list of elements:
#'   \item{Total_dissimilarity}{A data-frame on the total dissimilarity for each
#'   observed comparison. The data-frame has been sorted in decreasing order of
#'   the total dissimilarity.}
#'   \item{Types_used}{A data-frame with type mode (i.e., double or integer) of
#'   each characteristic.}
#'   \item{Total_missing}{The percentage of missing cases in the dataset,
#'   calculated as the ratio of total missing cases to the product of the number
#'   of studies with the number of characteristics.}
#'   \item{Cluster_color}{A data-frame on the number and color of the cluster
#'   assigned to each comparison. This can be used in the
#'   \code{\link{network_comparisons}} function to colour the edges according to
#'   the cluster assigned. The data-frame has been sorted in decreasing order of
#'   the total dissimilarity.}
#'   \item{Dissimilarity_table}{A lower off-diagonal matrix of 'dist' class
#'   with the dissimilarities of all pairs of comparisons. The row and column
#'   names has been sorted in decreasing order of the total dissimilarity.}
#'   \item{Table_internal_measures}{A data-frame with the connectivity index,
#'   silhouette width, and Dunn index for a range of 2 to P-1 clusters, with P
#'   being the number of comparisons.}
#'   \item{Table_cophenetic_coefficient}{A data-frame on the cophenetic
#'   correlation coefficient for all pairwise combinations of two dissimilarity
#'   measures (Euclidean, and Canberra) with eight linkage methods (Ward's two
#'   versions, single, complete, average, Mcquitty, median and centroid). The
#'   data-frame has been sorted in decreasing order of the cophenetic correlation
#'   coefficient.}
#'   \item{Optimal_dist}{The optimal dissimilarity measure (Euclidean or Canberra)
#'   based on the cophenetic correlation coefficient.}
#'   \item{Optimal_link}{The optimal linkage method (ward.D, ward.D2, single,
#'   complete, average, mcquitty, median, or centroid) based on the cophenetic
#'   correlation coefficient.}
#'
#'   If \code{get_plots = FALSE} only the list of elements mentioned above is
#'   returned. If \code{get_plots = TRUE}, \code{comp_clustering} returns a
#'   series of plots in addition to the list of elements mentioned above:
#'   \item{Dissimilarity_comparison}{A violin plot with integrated box plots and
#'   dots on the study dissimilarities per comparison (x-axis). Violins are
#'   sorted in descending order of the total dissimilarity (red point).}
#'   \item{Characteristics_contribution}{A bubble plot on the percentage of
#'   average contribution of each characteristic (x-axis) to dissimilarities of
#'   pairs of studies in the corresponding comparison (y-axis).}
#'   \item{Total_dissimilarity_plot}{A barplot on the total dissimilarity of
#'   comparisons. Bars are sorted in descending order of the total dissimilarity
#'   and have been coloured based on the cluster they belong, with the clusters
#'   referring to the optimal partitioning determined by the argument
#'   \code{optimal_clusters}, when heuristic clustering is performed, or based
#'   on the informative clustering.}
#'   \item{Internal_measures_panel}{A panel of profile plots on the connectivity
#'   index, silhouette width, and Dunn index for a range of 2 to P-1 clusters,
#'   with P being the number of comparisons. The candidate optimal number of
#'   clusters is indicated with a red point directly on the line.}
#'   \item{Silhouette_comparisons}{A silhouette plot illustrating the silhouette
#'   width for each comparison. The comparisons are sorted in descending order
#'   of the silhouette width and have been coloured based on the cluster they
#'   belong, with the clusters referring to the optimal partitioning. For
#'   comparisons with zero silhouette width, 0.01 is added (an arbitrary very
#'   small number) to make visible the colour of the corresponding bars.}
#'
#' @details
#'   The correct type mode of columns in \code{input} must be ensured to use
#'   the function \code{comp_clustering}. The first three columns referring to
#'   the trial name, first and second arm of the comparison, respectively, must
#'   be \strong{character}. The remaining columns referring to the
#'   characteristics must be \strong{double} or \strong{integer} depending on
#'   whether the corresponding characteristic refers to a quantitative or
#'   qualitative variable. The type mode of each column is assessed by
#'   \code{comp_clustering} using the base function \code{typeof}. Note that
#'   \code{comp_clustering} invites unordered and ordered variables; for the
#'   latter, add the argument \code{ordered = TRUE} in the base function
#'   \bold{factor()}.
#'
#'   The interventions should be sorted in an ascending order of their
#'   identifier number within the trials so that the first intervention column
#'   (second column in \code{input}) is the control arm for every pairwise
#'   comparison. This is important to ensure consistency in the intervention
#'   order within the comparisons obtained from the other related functions.
#'
#'   \code{comp_clustering} excludes from the dataset the following type of
#'   characteristics: (i) completely missing characteristics and
#'   (ii) characteristics with missing values in all but one studies for at
#'   least one non-single-stufy comparison. Then it proceeds with the clustering
#'   process.
#'
#'   The cophenetic correlation coefficient is calculated using the
#'   \code{\link[stats:cophenetic]{cophenetic}} function alongside the
#'   \code{\link[stats:dist]{dist}} function for selected dissimilarity measures
#'   found in the latter and the \code{\link[stats:hclust]{hclust}} function for
#'   selected linkage methods.
#'
#'  \code{comp_clustering} can be used only for a network with at least three
#'   comparisons. Otherwise, the execution of the function will be stopped and
#'   an error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso
#'  \code{\link{connectivity_index}},
#'  \code{\link[stats:cophenetic]{cophenetic}}, \code{\link[stats:dist]{dist}},
#'  \code{\link[stats:hclust]{hclust}}, \code{\link{internal_measures_plot}},
#'  \code{\link{network_comparisons}}
#'
#' @references
#' Gower J. General Coefficient of Similarity and Some of Its Properties.
#' \emph{Biometrics} 1971;\bold{27}(4):857--71.
#' doi: 10.2307/2528823
#'
#' Sokal R, Rohlf F. The Comparison of Dendrograms by Objective Methods.
#' \emph{Int Assoc Plant Taxon} 1962;\bold{11}(2):33--40.
#' doi: 10.2307/1217208
#'
#' Handl J, Knowles J, Kell DB. Computational cluster validation in post-genomic
#' data analysis. \emph{Biometrics} 2005;\bold{21}(15):3201--120.
#' doi: 10.1093/bioinformatics/bti517
#'
#' Rousseeuw PJ. Silhouettes: A graphical aid to the interpretation and
#' validation of cluster analysis.
#' \emph{J Comput Appl Math} 1987;\bold{20}:53--65.
#'
#' Dunn J. Well-separated clusters and optimal fuzzy partitions.
#' \emph{J Cybern} 1974;\bold{4}(1):95--104.
#'
#' @export
comp_clustering <- function (input,
                             drug_names,
                             rename_char = NULL,
                             height = FALSE,
                             num_neighb,
                             optimal_clusters,
                             get_plots = "none",
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


  ## Create new input and name the treatments
  input <- input0
  input[, 2:3] <- matrix(drug_names[as.numeric(unlist(input[, 2:3]))],
                         nrow = dim(input)[1],
                         ncol = 2)


  ## Defaults with 'get_plots'
  get_plots <- if (missing(get_plots)) {
    FALSE
  } else if (!is.element(get_plots, c(TRUE, FALSE))) {
    stop("'get_plots is logical.", call. = FALSE)
  } else {
    get_plots
  }


  ## Table with the variable type
  char_type <- data.frame(characteristic = colnames(input[, -c(1:3)]),
                          type = sapply(input[, -c(1:3)], typeof))
  rownames(char_type) <- 1:(dim(input)[2] - 3)


  ## Insert 'Comparison' in the dataset (control appears second in the compar.)
  input$Comparison <- as.character(paste(input$Arm2, "vs", input$Arm1))


  ## Number of unique observed comparisons
  unique_comp <- length(unique(input$Comparison))


  ## Number of single-study comparisons
  single_study_comp <- length(which(table(input$Comparison) == 1))


  ## Message on the number of comparisons and single-study comparisons
  message(paste0("-", " ", unique_comp, " ", "observed comparisons (",
                 single_study_comp, " ", "single-study comparisons)"))


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


  ## Address single-trial comparisons
  for (i in 1:length(split_dataset)) {
    if (dim(split_dataset[[i]])[1] < 2) {
      min_split[[i]] <-
        sapply(subset(input_new,
                      Comparison != unique(split_dataset[[i]][, "Comparison"]))
               [-c(1, 2)], function(x) if (typeof(x) == "double") {
                 min(x, na.rm = TRUE)
                 } else if (typeof(x) == "integer") {
                   names(which.min(table(x))) # as.integer(names(which.min(table(x))))
                   })
      max_split[[i]] <-
        sapply(subset(input_new,
                      Comparison != unique(split_dataset[[i]][, "Comparison"]))
               [-c(1, 2)], function(x) if (typeof(x) == "double") {
                 max(x, na.rm = TRUE)
                 } else if (typeof(x) == "integer") {
                   names(which.max(table(x))) # as.integer(names(which.max(table(x))))
                   })
      split_dataset[[i]] <- rbind(split_dataset[[i]],
                                  data.frame(Trial_name =
                                               c(paste0("new_", i, "_min"),
                                                 paste0("new_", i, "_max")),
                                             Comparison =
                                               rep(unique(split_dataset[[i]]
                                                          [, "Comparison"]), 2),
                                             rbind(min_split[[i]],
                                                   max_split[[i]])))
    }
  }


  ## Ensure that the characteristics have the correct type
  for (i in 1:length(split_dataset)) {
    for (j in 1:dim(input_new)[2]) { # input_new0
      if (is.character(input_new[, j])) { # input_new0
        split_dataset[[i]][, j] <- as.character(split_dataset[[i]][, j])
      } else if (is.double(input_new[, j])) { # input_new0
        split_dataset[[i]][, j] <- as.double(split_dataset[[i]][, j])
      } else if (is.integer(input_new[, j])) { # input_new0
        split_dataset[[i]][, j] <- as.integer(split_dataset[[i]][, j])
      }
    }
  }


  ## Stop for networks with two comparisons only (clustering is redundant)
  if (length(split_dataset) < 3) {
    stop(paste0("Clustering is redundant for two comparisons only!"),
         call. = FALSE)
  }


  ## Calculate the Gower dissimilarity among trials by comparison
  comparison_gower <- lapply(split_dataset,
                             function(x) {gower_distance(input = x)})


  ## Rename columns if indicated
  if (!is.null(rename_char)) {
    colnames(input_new)[rename_char[[1]] - 1] <- rename_char[[2]]
  }


  ## Variable on sample size
  colnames(input_new)[with(input_new,
                           startsWith(names(input_new),
                                      c("sample", "Sample")))] <- "Sample size"


  ## Prepare dataset for characteristic average contribution
  # Get the character by comparison data-frame
  # zero contribution means that all studies had the same value for the corresponding characteristic
  contr_dataset0 <- as.data.frame(sapply(comparison_gower, function(x) x[[4]]))
  rownames(contr_dataset0) <- names(input_new[, -c(1:2)])


  # Turn into long format for ggplot2
  suppressMessages({
  contr_dataset <-
    data.frame(melt(contr_dataset0),
               char = rep(names(input_new[, -c(1:2)]), dim(contr_dataset0)[2]))
  contr_dataset$size <-
    ifelse(contr_dataset$value < 0.25, "low",
           ifelse(contr_dataset$value >= 0.25 &
                    contr_dataset$value < 0.50, "moderate",
                  ifelse(contr_dataset$value >= 0.50 &
                           contr_dataset$value < 0.74, "high", "very high")))
  colnames(contr_dataset)[c(1, 3)] <- c("comparison", "characteristic")
  })

  # Get the bubble-plot
  char_contr <-
    ggplot(contr_dataset,
           aes(x = characteristic,
               y = comparison)) +
    geom_point(aes(size = ((value - min(value)) / (max(value) - min(value)))*15,
                   fill = size),
               shape = 21,
               color = "black",
               stroke = 1.2,
               alpha = 0.5) +
    geom_text(aes(label = ifelse(value > 0, round(value * 100, 1), " ")),
              size = label_size,
              color="black",
              fontface = "bold") +
    labs(x = " ",
         y = " ") +
    scale_size_identity() +
    scale_x_discrete(position = "top",
                     labels = function(x) str_wrap(x, width = 2)) +
    scale_fill_manual(name = "Contribution",
                      breaks = c("low", "moderate", "high", "very high"),
                      limits = c("low", "moderate", "high", "very high"),
                      values = c("#A6D854", "#E6AB02", "#D95F02", "#E31A1C"),
                      labels = c("Low", "Moderate", "High", "Very high")) +
    theme_bw() +
    guides(fill = guide_legend(override.aes = list(size = 6))) +
    theme(panel.grid.major = element_line(linetype = 2, color = "grey"),
          title = element_text(size = title_size, face = "bold"),
          axis.title = element_text(size = axis_title_size, face = "bold"),
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = axis_x_text_angle,
                                     hjust =
                                       ifelse(axis_x_text_angle == 0, 0.5, 1)),
          legend.position = "bottom",
          legend.text = element_text(size = legend_text_size))


  ## Data-frame of total dissimilarity
  total_diss0 <- data.frame(as.character(names(comparison_gower)),
                            as.character(1:length(comparison_gower)),
                            round(
                              sapply(comparison_gower,
                                     function(x)
                                       sqrt(mean(na.omit(unlist(x[[1]]))^2))),
                              3),
                            stringsAsFactors = FALSE)
  colnames(total_diss0) <- c("comparison", "id", "total_dissimilarity")
  rownames(total_diss0) <- as.character(names(comparison_gower))


  ## Data-frame on number of trials per comparison
  comp_size <-
    data.frame(total_diss0[, c(1, 3)],
               unlist(lapply(split_dataset0, function(x) dim(x)[1])))
  colnames(comp_size)[3] <- "num_trials"


  ## Prepare dataset on comparison dissimilarities and total dissimilarities
  diss_dataset <-
    data.frame(diss = unlist(sapply(comparison_gower,
                                    function(x) na.omit(unlist(x[[1]])))),
               comp = rep(names(sapply(comparison_gower,
                                       function(x) na.omit(unlist(x[[1]])))),
                          sapply(comparison_gower,
                                 function(x) length(na.omit(unlist(x[[1]]))))),
               total = rep(total_diss0$total_dissimilarity,
                           sapply(comparison_gower,
                                  function(x) length(na.omit(unlist(x[[1]]))))),
               comp_size =
                 rep(comp_size$num_trials,
                     sapply(comparison_gower,
                            function(x) length(na.omit(unlist(x[[1]])))))
               )
  diss_dataset$col <- ifelse(diss_dataset$comp_size > 1, "No", "Yes")


  ## Violin plot on dissimilarity distribution per comparison
  comp_diss_plot <-
    ggplot(diss_dataset,
           aes(x = reorder(comp, total, decreasing = TRUE),
               y = diss)) +
    geom_violin(aes(fill = col),
                trim = TRUE, #FALSE
                alpha = 0.3) +
    geom_boxplot(outlier.alpha = 0.3,
                 fill = "white",
                 colour = "black",
                 varwidth = TRUE) +
    geom_point() +
    geom_point(data = total_diss0,
               aes(x = comparison,
                   y = total_dissimilarity),
               color = "red",
               size = 2.5,
               shape = 21,
               stroke = 1.5) +
    geom_text(data = total_diss0,
              aes(x = comparison,
                  y = total_dissimilarity,
                  label = sprintf("%0.2f", round(total_dissimilarity, 2))),
              hjust = 1.3, #1.2
              vjust = 0.2,
              size = label_size,
              fontface = "bold",
              colour = "blue") +
    geom_text(data = comp_size,
              aes(x = comparison,
                  y = 0,
                  label = paste("n =", num_trials)),
              hjust = 0.5,
              vjust = 2.8,
              size = label_size,
              fontface = "plain",
              colour = "black") +
    stat_boxplot(geom = 'errorbar',
                 width = 0.2,
                 linetype = "dashed") +
    scale_fill_manual(name = "Includes pseudostudies",
                      breaks = c("Yes", "No"),
                      limits = c("Yes", "No"),
                      values = c("red", "#99CCFF"),
                      labels = c("Yes", "No")) +
    labs(x = "Comparisons",
         y = "Gower's dissimilarity") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic() +
    theme(title = element_text(size = title_size, face = "bold"),
          axis.title = element_text(size = axis_title_size , face = "bold"),
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = axis_x_text_angle,
                                     hjust =
                                       ifelse(axis_x_text_angle == 0, 0.5, 1)),
          legend.position = "bottom",
          legend.text = element_text(size = legend_text_size))

  ## Sort 'total_diss' in decreasing order of total dissimilarity
  total_diss <- total_diss0[order(total_diss0$total_dissimilarity, decreasing = TRUE),]


  ## Linkage methods of the 'hclust' function
  linkage_methods <- c("ward.D", "ward.D2", "single", "complete", "average",
                       "mcquitty", "median", "centroid")


  ## Different route depending on whether we choose informative or heuristic clustering
  if (height == TRUE) { # Informative clustering

    ## Get the clusters per comparison
    clusters0 <-
      do.call(rbind,
              lapply(total_diss[, 3],
                     function(x)
                       ifelse(x <= 0.25, 1,
                              ifelse(x > 0.25 & x <= 0.50, 2,
                                     ifelse(x > 0.50 & x <= 0.75, 3, 4)))))


    ## Data-frame of comparisons and corresponding cluster
    clusters <- data.frame(comparison = total_diss[, 1],
                           cluster = match(clusters0, unique(clusters0)))


    ## Clusters obtained
    optimal_clusters <- length(unique(clusters0))


    ## Barplot on total dissimilarity (Prepare the dataset with the clusters)
    total_diss_new <- data.frame(total_diss[, c(1, 3)],
                                 cluster = clusters[, 2])

  } else { # Heuristic clustering

    ## Checking further defaults
    # Number of 'optimal' clusters (based on the internal measures)
    optimal_clusters <- if (missing(optimal_clusters)) {
      2
      #stop("The argument 'optimal_clusters' must be defined", call. = FALSE)
    } else if ((optimal_clusters > length(unique(input$Comparison)) - 1 ||
                optimal_clusters < 2) & length(unique(input$Comparison)) > 3) {
      stop(paste0("'optimal_clusters' must range from 2 to", " ",
                  length(unique(input$Comparison)) - 1, "."), call. = FALSE)
    } else if ((optimal_clusters > length(unique(input$Comparison)) - 1 ||
                optimal_clusters < 2) & length(unique(input$Comparison)) == 3) {
      stop(paste0("'optimal_clusters' must equal exactly 2."), call. = FALSE)
    } else {
      optimal_clusters
    }


    ## Default (to be used in 'connectivity_index')
    num_neighb <- if (missing(num_neighb)) {
      round(length(split_dataset) / 2, 0)
    } else if (num_neighb > length(split_dataset) || num_neighb < 2) {
      stop(paste0("'num_neighb' must range from 2 to", " ",
                  length(split_dataset) - 1, "."), call. = FALSE)
    } else {
      num_neighb
    }


    ## Distance methods of the 'dist' function
    distance_methods <- c("euclidean", "canberra")


    ## Possible combinations of 'distance_methods' with 'linkage_methods'
    # Data-frame of possible combinations
    poss_comb <- expand.grid(distance = distance_methods,
                             linkage = linkage_methods,
                             stringsAsFactors = FALSE)

    # Obtain results on cophenetic correlation coefficient
    table_coph <-
      data.frame(poss_comb,
                 results =
                   mapply(function(x, y)
                     round(cor(dist(total_diss[, 3], method = x, p = 3),
                               cophenetic(hclust(dist(total_diss[, 3],
                                                      method = x, p = 3),
                                                 method = y))), 3),
                     x = as.character(poss_comb$distance),
                     y = as.character(poss_comb$linkage)))


    # Bring both tables together
    table_cophenetic <-
      table_coph[order(table_coph$results, decreasing = TRUE), ]


    ## Select the distance and linkage methods for the max cophenetic coefficient
    optimal_dist_link <- subset(table_cophenetic, results == max(results))


    ## When more distances or linkages are proper for the same cophenetic coeff.
    if (length(unique(optimal_dist_link[, 1])) > 1) {
      optimal_dist <- optimal_dist_link[1, 1]
      optimal_link <-
        optimal_dist_link[optimal_dist_link$distance == optimal_dist, 2][1]
    } else if (length(unique(optimal_dist_link[, 1])) == 1) {
      optimal_dist <- unique(optimal_dist_link[, 1])
      optimal_link <- optimal_dist_link[1, 2]
    } else if (dim(optimal_dist_link[, 1])[1] == 1) {
      optimal_dist <- optimal_dist_link[1]
      optimal_link <- optimal_dist_link[2]
    }


    ## Report the optimal dissimilarity measure and linkage method
    message(paste("- Cophenetic coefficient:", max(table_cophenetic[, 3])))
    message(paste("- Optimal dissimilarity measure:", optimal_dist))
    message(paste("- Optimal linkage method:", optimal_link))


    ## Dissimilarity matrix of comparisons for the optimal dissimilarity measure
    data_cluster <-
      dist(matrix(total_diss[, 3], dimnames = list(rownames(total_diss))),
           method = optimal_dist)


    ## Table on internal measures results for all combinations
    table_internal_measures <-
      internal_measures_plot(input = data_cluster,
                             num_neighb = num_neighb,
                             optimal_link = optimal_link)$Table_internal_measures


    ## Panel of internal measures
    internal_measures_panel <- if (dim(table_internal_measures)[1] > 1) {
      internal_measures_plot(input = data_cluster,
                             num_neighb = num_neighb,
                             optimal_link = optimal_link)$Internal_measures_panel
    } else {
      a <- "At least four comparisons are needed to create a panel with plots"
      b <- "on internal measures for a range of clusters!"
      message(paste(a, b))
    }


    ## Silhouette width per comparison for selected cluster and linkage method
    silhouette_comp <-
      silhouette_index(input = data_cluster,
                       method = optimal_link,
                       num_clusters = optimal_clusters)$silhoutte_comp

    ## Add 0.02 to comparison with 0 silhouette so that the bar is coloured
    silhouette_comp$silhouette_new <-
      ifelse(silhouette_comp$silhouette == 0, 0.01, silhouette_comp$silhouette)


    ## Average silhouette width
    average_silhouette <-
      silhouette_index(input = data_cluster,
                       method = optimal_link,
                       num_clusters = optimal_clusters)$silhoutte_width


    ## Barplot on total dissimilarity
    # Prepare the dataset with the clusters
    total_diss_new <- data.frame(total_diss[, c(1, 3)],
                                 cluster = silhouette_comp[, 2])


    ## Plot silhouette by comparison and cluster
    # SOS: Because the bars are order by silhouette value, the colour order is
    # not the same with that in 'Total diss bar plot', but both plots share the
    # same colours for the same interventions (which is the desired)!
    plot_comp_silhouette <-
      ggplot(silhouette_comp,
             aes(x = silhouette_new,
                 y = reorder(comp, silhouette_new),
                 group = reorder(factor(cluster), silhouette_new),
                 fill = factor(cluster))) +
      geom_bar(stat = "identity") +
      geom_vline(xintercept = average_silhouette,
                 colour = "black",
                 linewidth = 0.6,
                 linetype = 3) +
      geom_text(aes(label = sprintf("%0.2f",round(silhouette, 2))),
                hjust = 1.1,
                vjust = 0.2,
                size = label_size,
                colour = "black") +
      geom_text(aes(x = average_silhouette,
                    y = 0.44,
                    label = sprintf("%0.2f",round(average_silhouette, 2))),
                hjust = 0.5,
                vjust = 0.0,
                colour = "blue",
                size = label_size) +
      scale_x_continuous(limits = c(-1, 1)) +
      labs(x = "Silhouette width",
           y = " ",
           fill = "Cluster") +
      theme_classic() +
      guides(fill = guide_legend(nrow = 1)) +
      scale_fill_discrete(limits = levels(factor(total_diss_new$cluster)),
                          labels = factor(1:max(total_diss_new$cluster))) +
      theme(title = element_text(size = title_size, face = "bold"),
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_text_size),
            legend.position = "bottom",
            legend.text = element_text(size = legend_text_size),
            plot.caption = element_text(size = 10, hjust = 0.0))
  }


  ## Barplot on total dissimilarity
  # Create the plot
  total_diss_plot <-
    ggplot(total_diss_new,
           aes(x = reorder(factor(comparison),
                           total_dissimilarity,
                           decreasing = TRUE),
               y = round(total_dissimilarity, 2),
               fill = factor(cluster))) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = sprintf("%0.2f", round(total_dissimilarity, 2))),
              hjust = 0.5,
              vjust = -0.5,
              size = label_size,
              colour = "black") +
    labs(x = "Comparisons",
         y = "Total dissimilarity",
         fill = "Cluster") +
    theme_classic() +
    coord_cartesian(ylim = c(0, 1)) +
    guides(fill = guide_legend(nrow = 1)) +
    scale_fill_discrete(limits = levels(factor(total_diss_new$cluster)),
                        labels = factor(1:max(total_diss_new$cluster))) +
    theme(title = element_text(size = title_size, face = "bold"),
          axis.title = element_text(size = axis_title_size),
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = axis_x_text_angle,
                                     hjust =
                                       ifelse(axis_x_text_angle == 0, 0.5, 1)),
          legend.position = "bottom",
          legend.text = element_text(size = legend_text_size))


  ## Percentage total missing data
  total_mod <-
    round((sum(is.na(input_new0[, -c(1, 2)]) == TRUE) /
             (dim(input_new0[, -c(1, 2)])[1] *
                dim(input_new0[, -c(1, 2)])[2])) * 100, 2)


  ## Data-frame with the colour and cluster per comparison
  cluster_color <- data.frame(comparison = total_diss_new[, 1],
                              cluster = total_diss_new[, 3],
                              colour = scales::hue_pal()(optimal_clusters)[
                                total_diss_new[, 3]])


  ## Collect the results
  # First without the table with the internal measure results
  collect0 <- list(Total_dissimilarity = total_diss,
                   Types_used = char_type,
                   Total_missing = paste0(total_mod, "%"),
                   Cluster_color = cluster_color)

  # Define the results based on the argument 'height'
  collect <- if (height == TRUE) {
    collect0
  } else {
    append(collect0, list(Dissimilarity_table = round(data_cluster, 3),
                          Table_internal_measures = table_internal_measures,
                          Table_cophenetic_coefficient = table_cophenetic,
                          Optimal_dist = optimal_dist,
                          Optimal_link = optimal_link))
  }


  ## Report results based on 'get_plots'
  results <- if (get_plots == FALSE) {
    collect
  } else if (height == FALSE & get_plots == TRUE) { #& !is.null(internal_measures_panel)) {
    append(collect, list(Dissimilarity_comparison = comp_diss_plot,
                         Characteristics_contribution = char_contr,
                         Total_dissimilarity_plot = total_diss_plot,
                         Internal_measures_panel = internal_measures_panel,
                         Silhouette_comparisons = plot_comp_silhouette))
  } else if (height == TRUE & get_plots == TRUE) {
    append(collect, list(Dissimilarity_comparison = comp_diss_plot,
                         Characteristics_contribution = char_contr,
                         Total_dissimilarity_plot = total_diss_plot))
  }

  class(results) <- "comp_clustering"

  return(suppressWarnings({print(results)}))
}
