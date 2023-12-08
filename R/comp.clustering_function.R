#' End-user-ready results for informative and heuristic clustering of comparisons
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{comp_clustering} hosts a toolkit of functions that facilitates
#'   conducting, visualising and evaluating informative decision and
#'   hierarchical agglomerative of observed comparisons of interventions for a
#'   specific network and set of characteristics that act as effect modifiers.
#'
#' @param input A data-frame in the long arm-based format. Two-arm trials occupy
#'   one row in the data-frame. Multi-arm trials occupy as many rows as the
#'   number of possible comparisons among the interventions. The first three
#'   columns refer to the trial name, first and second arm of the comparison,
#'   respectively. The remaining columns refer to summary characteristics. See
#'   'Details' for the specification of the columns.
#' @param drug_names A vector of labels with the name of the interventions
#'   in the order they have been defined in the argument \code{input}.
#' @param threshold A scalar in the interval [0, 1) to indicate the cut-off of
#'   low dissimilarity of two comparisons. The value much be low. The default
#'   value is 0.13.
#' @param informative Logical with \code{TRUE} for performing informative
#'   decision and \code{FALSE} for performing hierarchical agglomerative
#'   clustering, thus, allowing the user to define the number of clusters via
#'   the argument \code{optimal_clusters}. The default argument is \code{TRUE}.
#' @param optimal_clusters A positive integer for the optimal number of clusters
#'   based on the silhouette width, an internal validation measure. The user
#'   \bold{must} define the value, after inspecting the profile plot on the
#'   average silhouette width. It takes values from two to the number of
#'   comparisons minus one.
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
#'   and the optimal linkage method selected based on the cophenetic correlation
#'   coefficient. Then, the function returns the following list of elements:
#'   \item{Trials_diss_table}{A lower off-diagonal matrix of 'dist' class
#'   with the Gower dissimilarities of all pairs of studies in the network.}
#'   \item{Comparisons_diss_table}{A lower off-diagonal matrix of 'dist' class
#'   with the total across-comparison dissimilarities of all pairs of observed
#'   intervention comparisons in the network.}
#'   \item{Informative_heatmap}{A heatmap on total within-comparison and total
#'   across-comparison dissimilarities when the 'informative decision' is
#'   applied (\code{informative = TRUE}). Diagonal elements refer to total
#'   within-comparison dissimilarity, and off-diagonal elements refer to
#'   across-comparisons dissimilarity. Using a threshold of high similarity
#'   (default at 0.13), cells exceeding this threshold are highlighted in red;
#'   otherwise, in green. This heatmap aids in finding 'hot spots' of
#'   comparisons that may violate the plausibility of transitivity in the
#'   network.}
#'   \item{Total_dissimilarity}{A data-frame on the observed comparisons and
#'   comparisons between comparisons, alongside the corresponding total
#'   within-comparison and total across-comparisons dissimilarity.
#'   The data-frame has been sorted in decreasing order of the total (within- or
#'   across-comparison) dissimilarity.}
#'   \item{Types_used}{A data-frame with type mode (i.e., double or integer) of
#'   each characteristic.}
#'   \item{Total_missing}{The percentage of missing cases in the dataset,
#'   calculated as the ratio of total missing cases to the product of the number
#'   of studies with the number of characteristics.}
#'   \item{Cluster_color}{A data-frame on the number and color of the cluster
#'   assigned to each comparison.}
#'   \item{Table_average_silhouette_width}{A data-frame with the average
#'   silhouette width for a range of 2 to P-1 clusters, with P being the number
#'   of comparisons.}
#'   \item{Table_cophenetic_coefficient}{A data-frame on the cophenetic
#'   correlation coefficient for eight linkage methods (Ward's two
#'   versions, single, complete, average, Mcquitty, median and centroid). The
#'   data-frame has been sorted in decreasing order of the cophenetic correlation
#'   coefficient.}
#'   \item{Optimal_link}{The optimal linkage method (ward.D, ward.D2, single,
#'   complete, average, mcquitty, median, or centroid) based on the cophenetic
#'   correlation coefficient.}
#'
#'   If \code{get_plots = FALSE} only the list of elements mentioned above is
#'   returned. If \code{get_plots = TRUE}, \code{comp_clustering} returns a
#'   series of plots in addition to the list of elements mentioned above:
#'   \item{Dissimilarity_comparison}{A violin plot with integrated box plots and
#'   dots on the study dissimilarities per comparison and comparison between
#'   comparisons (x-axis). Violins are sorted in descending order of the total
#'   dissimilarity (red point) and are split into those for within-comparison
#'   dissimilarities and across-comparison dissimilarities.}
#'   \item{Profile_plot}{A profile plot on the average silhouette width for a
#'   range of 2 to P-1 clusters, with P being the number of comparisons.
#'   The candidate optimal number of  clusters is indicated with a red point
#'   directly on the line.}
#'   \item{Silhouette_width_plot}{A silhouette plot illustrating the silhouette
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
#'   \code{\link[stats:hclust]{hclust}} function for selected linkage methods.
#'
#'  \code{comp_clustering} can be used only for a network with at least three
#'   comparisons. Otherwise, the execution of the function will be stopped and
#'   an error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso
#'  \code{\link[stats:cophenetic]{cophenetic}},
#'  \code{\link[stats:hclust]{hclust}}, \code{\link{internal_measures_plot}}
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
#' @export
comp_clustering <- function (input,
                             drug_names,
                             threshold = 0.13,
                             informative = TRUE,
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

  # To print or not plots
  get_plots <- if (missing(get_plots)) {
    FALSE
  } else if (!is.element(get_plots, c(TRUE, FALSE))) {
    stop("'get_plots is logical.", call. = FALSE)
  } else {
    get_plots
  }


  ## Create new input and name the treatments
  input <- input0
  input[, 2:3] <- matrix(drug_names[as.numeric(unlist(input[, 2:3]))],
                         nrow = dim(input)[1],
                         ncol = 2)


  ## Table with the variable type
  char_type <- data.frame(characteristic = colnames(input[, -c(1:3)]),
                          type = sapply(input[, -c(1:3)], typeof))
  rownames(char_type) <- 1:(dim(input)[2] - 3)


  ## Insert 'Comparison' in the dataset (control appears second in the compar.)
  input$Comparison <- as.character(paste0(input$Arm2, "-", input$Arm1))


  ## Unique comparisons
  unique_comp <- unique(input$Comparison)


  ## Number of unique observed comparisons
  num_unique_comp <- length(unique_comp)


  ## Single-study comparisons
  single_study_comp <- names(which(table(input$Comparison) == 1))


  ## Number of single-study comparisons
  num_single_study_comp <- length(single_study_comp)


  ## Message on the number of comparisons and single-study comparisons
  message(paste0("-", " ", num_unique_comp, " ", "observed comparisons (",
                 num_single_study_comp, " ", "single-study comparisons)"))


  ## Reduce dataset to trial, comparison & characteristics
  input_new0 <- input[, c(1, dim(input)[2], 4:(dim(input)[2] - 1))]


  ## Split 'dataset' by 'Comparison'
  split_dataset0 <- split(input_new0, f = input$Comparison)


  ## Stop for networks with two comparisons only (clustering is redundant)
  if (length(split_dataset0) < 3) {
    stop(paste0("Clustering is redundant for two comparisons only!"),
         call. = FALSE)
  }


  ## Find the completely missing columns in all non-single-study comparisons
  col_all_miss <-
    unique(unlist(
      lapply(split_dataset0, function(x) if (dim(x)[1] > 1)
        as.vector(which(colSums(is.na(x)) == nrow(x) |
                          colSums(is.na(x)) == nrow(x) - 1)))))


  ## Keep the names of the completely missing columns in all comparisons
  col_all_miss_names <-
    unique(unlist(lapply(split_dataset0,
                         function(x) colnames(x)[col_all_miss])))


  ## Message on the dropped characteristics
  dropped_chars <- if (length(col_all_miss_names) == 0) {
    "none"
  } else {
    col_all_miss_names
  }
  message(paste("- Dropped characteristics:", paste(dropped_chars,
                                                    collapse = ", ")))


  ## Remove these columns also from the dataset for the moment
  input_new <-
    if (length(col_all_miss) > 0) {
      subset(input_new0, select = -col_all_miss)
    } else {
      input_new0
    }


  ## 'Re-name' the multi-arm trials as their name is repeated!
  input_new$Trial_name <-
    ave(input_new$Trial_name, input_new$Trial_name,
        FUN = function(x) if (length(x) > 1) paste0(x[1], "(", seq_along(x), ")") else x[1])


  ## Calculate the Gower dissimilarity of all study pairs in the network
  gower_res <- gower_distance(input = input_new)


  ## Re-name the columns/rows with the corresponding comparisons
  gower_diss_mat <- as.matrix(gower_res$Dissimilarity_table)
  colnames(gower_diss_mat) <- input_new[1:(dim(input_new)[1] - 1), 2]
  rownames(gower_diss_mat) <- input_new[2:(dim(input_new)[1]), 2]


  ## Data-frame on compared comparisons, and corresponding Gower value
  # First turn 'gower_diss_mat' into data.frame with 'melt'
  dataset_diss <- na.omit(melt(gower_diss_mat))

  # Re-order the comparisons within based on the order in unique_comp
  dataset_diss[, 1:2] <-
    t(apply(dataset_diss[, 1:2], 1,
            function(x) x[order(match(x, unique_comp))]))

  # Create the comparison of comparisons using 'paste'
  dataset_diss$comp <- apply(dataset_diss[, 1:2], 1, paste, collapse = " vs ")


  ## Split 'dataset' by 'Comparison of comparisons'
  split_dataset <- split(dataset_diss, f = dataset_diss$comp)


  ## Calculate within & across comparison total dissimilarity (Dp)
  d_p <- round(sapply(split_dataset,
                      function(x) sqrt(mean(na.omit(unlist(x[[3]]))^2))), 2)


  ## Append the single-study comparisons (NA value)
  if (num_single_study_comp > 0) {

    # Comparisons of single-study comparison with itself
    name_single_study_comp <-
      apply(t(single_study_comp), 1, function(x) paste(x, " vs ", x))

    # Append the single-study comparison
    d_p <- append(d_p, rep(0, length(name_single_study_comp)))
    names(d_p)[length(split_dataset) + 1] <- name_single_study_comp
  }


  ## Lower triangular matrix of within & between comparisons total dissimilarity
  dist_mat <- matrix(NA, nrow = num_unique_comp, ncol = num_unique_comp)
  dist_mat[lower.tri(dist_mat, diag = TRUE)] <- d_p
  colnames(dist_mat) <- unique_comp
  rownames(dist_mat) <- unique_comp


  ## Prepare dataset on comparison dissimilarities and total dissimilarities
  # Set index for 'comparison' and 'comparison of comparisons'
  index_type <-
    apply(dataset_diss[!duplicated(dataset_diss[, 1:2]), ], 1,
          function(x) ifelse(setequal(x[1], x[2]),
                             "Within-comparison", "Across-comparison"))

  # Select name based on the 'index_type'
  name_type <- ifelse(index_type == "Within-comparison",
                      dataset_diss[!duplicated(dataset_diss[, 1:2]), 1],
                      dataset_diss[!duplicated(dataset_diss[, 1:2]), 4])

  # Create the data.frame
  diss_dataset <-
    data.frame(diss = unlist(lapply(split_dataset,
                                    function(x) na.omit(unlist(x[[3]])))),
               comp = rep(name_type,
                          lapply(split_dataset, function(x) dim(x)[1])),
               index = rep(index_type,
                           lapply(split_dataset, function(x) dim(x)[1])),
               total = rep(na.omit(d_p),
                           lapply(split_dataset,
                                  function(x) length(na.omit(unlist(x[[1]]))))))
  rownames(diss_dataset) <- NULL


  ## Violin plot on dissimilarity distribution per dissimilarity type
  # within-comparison & across-comparison dissimilarity
  suppressWarnings({
  comp_diss_plot <-
    ggplot(diss_dataset,
           aes(x = reorder(comp, total, decreasing = TRUE),
               y = diss)) +
    geom_violin(fill = "#99CCFF",
                trim = TRUE, #FALSE
                alpha = 0.2) +
    geom_boxplot(outlier.alpha = 0.3,
                 fill = "white",
                 colour = "black",
                 varwidth = TRUE) +
    stat_boxplot(geom = 'errorbar',
                 width = 0.2,
                 linetype = "dashed") +
    geom_point() +
    geom_point(aes(x = reorder(comp, total, decreasing = TRUE),
                   y = total),
               color = "red",
               size = 2.5,
               shape = 21,
               stroke = 1.5) +
    geom_text(aes(x = reorder(comp, total, decreasing = TRUE),
                  y = total,
                  label = sprintf("%0.2f", total)),
              hjust = 1.3, #1.2
              vjust = 0.2,
              size = label_size,
              fontface = "bold",
              colour = "blue",
              inherit.aes = FALSE) +
    facet_wrap(.~ index, scales = "free_x") +
    labs(x = "Comparisons",
         y = "Gower's dissimilarity") +
    coord_cartesian(ylim = c(0, 1)) +
    theme_classic() +
    theme(title = element_text(size = title_size, face = "bold"),
          axis.title = element_text(size = axis_title_size, face = "bold"),
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = axis_x_text_angle,
                                     hjust =
                                       ifelse(axis_x_text_angle == 0, 0.5, 1)),
          strip.text = element_text(size = axis_text_size, face = "bold"))
  })


  ## Data-frame of total dissimilarity
  total_diss0 <- data.frame(name_type,
                            na.omit(d_p),
                            index_type,
                            stringsAsFactors = FALSE)
  colnames(total_diss0)[1:2] <- c("comparison", "total_dissimilarity")
  rownames(total_diss0) <- NULL


  ## Sort to bring all 'within-comparison' at the beginning
  total_diss <- total_diss0[order(total_diss0$index_type, decreasing = TRUE), ]


  ## Different route depending on whether we choose informative decision or hierarchical clustering
  if (informative == TRUE) { # Informative decision

    ## Prepare dataset for dissimilarity heatmap
    # Keep single-study comparisons (NA values)
    mat_new <- melt(dist_mat, na.rm = FALSE)


    ## To create the orders of the lower diagonal
    xmin1 <- rep(seq(0.5, num_unique_comp - 0.5, 1), each = num_unique_comp)
    xmax1 <- rep(seq(1.5, num_unique_comp + 0.5, 1), each = num_unique_comp)
    ymin1 <- rep(seq(num_unique_comp - 0.5, 0.5, -1), each = num_unique_comp)
    ymax1 <- ymin1


    ## Create the heatmap for one network of interventions
    informative_heatmap <-
      ggplot(mat_new,
             aes(factor(Var2, levels = unique_comp[1:num_unique_comp]),
                 factor(Var1, levels = unique_comp[num_unique_comp:1]),
                 fill = ifelse(value < threshold, "high", "poor"))) +
      geom_tile(colour = "white",
                alpha = 0.5) +
      geom_text(aes(factor(Var2,
                           levels = unique_comp[1:num_unique_comp]),
                    factor(Var1, levels = unique_comp[num_unique_comp:1]),
                    label = ifelse(is.na(value), "", value),
                    fontface = "bold"),
                size = rel(label_size)) +
      geom_rect(aes(xmin = xmin1, xmax = xmax1, ymin = ymin1, ymax = ymax1),
                color = "black", linewidth = 1) +
      geom_rect(aes(xmin = ymin1, xmax = ymax1, ymin = xmin1, ymax = xmax1),
                color = "black", linewidth = 1) +
      scale_fill_manual(breaks = c("high", "poor"),
                        values = c("#009E73", "#D55E00"),
                        na.value = "white") +
      scale_x_discrete(position = "top") +
      labs(x = "", y = "") +
      theme_bw() +
      theme(legend.position = "none",
            axis.text = element_text(size = axis_text_size))

  } else { # Hierarchical agglomerative clustering


    ## Stop if triangle inequality is violated
    if (rdist::triangle_inequality(as.dist(dist_mat)) == FALSE) {
      stop(paste0("Triangle inequality is violated!"), call. = FALSE)
    }


    ## Number of 'optimal' clusters (based on the internal measures)
    optimal_clusters <- if (missing(optimal_clusters)) {
      stop("The argument 'optimal_clusters' must be defined", call. = FALSE)
    } else if ((optimal_clusters > num_unique_comp - 1 ||
                optimal_clusters < 2) & num_unique_comp > 3) {
      stop(paste0("'optimal_clusters' must range from 2 to", " ",
                  num_unique_comp - 1, "."), call. = FALSE)
    } else if ((optimal_clusters > num_unique_comp - 1 ||
                optimal_clusters < 2) & num_unique_comp == 3) {
      stop(paste0("Maximum two clusters are possible."), call. = FALSE)
    } else {
      optimal_clusters
    }


    ## Linkage methods of the 'hclust' function
    linkage_methods <- c("ward.D", "ward.D2", "single", "complete", "average",
                         "mcquitty", "median", "centroid")


    ## Obtain results on cophenetic correlation coefficient
    table_coph <-
      data.frame(linkage_methods,
                 results =
                   sapply(linkage_methods, function(x)
                     round(cor(as.dist(dist_mat),
                               cophenetic(hclust(as.dist(dist_mat),
                                                 method = x))), 3)))


    ## Sort in decreasing order
    table_cophenetic <-
      table_coph[order(table_coph$results, decreasing = TRUE), ]


    ## Select the linkage method for the maximum cophenetic coefficient
    optimal_link <- if (length(table_cophenetic[, 1]) > 1) {
      subset(table_cophenetic, results == max(results))[1, 1]
    } else {
      subset(table_cophenetic, results == max(results))
    }


    ## Report the optimal dissimilarity measure and linkage method
    message(paste("- Cophenetic coefficient:", max(table_cophenetic[, 2])))
    message(paste("- Optimal linkage method:", optimal_link))


    ## Table on average silhouette width results for all combinations
    table_internal_measures <-
      internal_measures_plot(input = as.dist(dist_mat),
                             optimal_link = optimal_link)$Table_internal_measures


    ## Panel of internal measures
    internal_measures_panel <- if (dim(table_internal_measures)[1] > 1) {
    internal_measures_plot(input = as.dist(dist_mat),
                           optimal_link = optimal_link)$Internal_measures_panel
    } else {
      a <- "At least four comparisons are needed to create the profile plot"
      b <- "for a range of clusters!"
      message(paste(a, b))
    }


    ## Silhouette width per comparison for selected cluster and linkage method
    silhouette_comp <-
      silhouette_index(input = as.dist(dist_mat),
                       method = optimal_link,
                       num_clusters = optimal_clusters)$silhoutte_comp

    ## Add 0.02 to comparison with 0 silhouette so that the bar is coloured
    silhouette_comp$silhouette_new <-
      ifelse(silhouette_comp$silhouette == 0, 0.01, silhouette_comp$silhouette)


    ## Average silhouette width
    average_silhouette <-
      silhouette_index(input = as.dist(dist_mat),
                       method = optimal_link,
                       num_clusters = optimal_clusters)$silhoutte_width


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
      geom_text(aes(label = sprintf("%0.2f", round(silhouette, 2))),
                hjust = 1.1,
                vjust = 0.2,
                size = label_size,
                colour = "black") +
      geom_text(aes(x = average_silhouette,
                    y = 0.44,
                    label = sprintf("%0.2f", round(average_silhouette, 2))),
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
      scale_fill_discrete(limits = levels(factor(silhouette_comp$cluster)),
                          labels = factor(1:max(silhouette_comp$cluster))) +
      theme(title = element_text(size = title_size, face = "bold"),
            axis.title = element_text(size = axis_title_size),
            axis.text = element_text(size = axis_text_size),
            legend.position = "bottom",
            legend.text = element_text(size = legend_text_size),
            plot.caption = element_text(size = 10, hjust = 0.0))


    ## Data-frame with the cluster per comparison
    cluster_comp0 <- data.frame(silhouette_comp[, 1:3])
    colnames(cluster_comp0)[1] <- c("comparison")


    ## Sort 'cluster' in decreasing order of cluster
    cluster_comp <- cluster_comp0[order(cluster_comp0$cluster,
                                        decreasing = FALSE), ]
  }


  ## Percentage total missing data
  total_mod <-
    round((sum(is.na(input_new0[, -c(1, 2)]) == TRUE) /
             (dim(input_new0[, -c(1, 2)])[1] *
                dim(input_new0[, -c(1, 2)])[2])) * 100, 2)


  ## Collect the results
  # First without the table with the internal measure results
  collect0 <-
    list(Trials_diss_table = round(gower_res$Dissimilarity_table, 3),
         Comparisons_diss_table = as.dist(dist_mat),
         Total_dissimilarity = total_diss,
         Types_used = char_type,
         Total_missing = paste0(total_mod, "%"))

  # Define the results based on the argument 'informative'
  collect <- if (informative == TRUE) {
    collect0
  } else {
    append(collect0,
           list(Table_average_silhouette_width = table_internal_measures,
                Table_cophenetic_coefficient = table_cophenetic,
                Optimal_link = optimal_link,
                Cluster_comp = cluster_comp))
  }


  ## Report results based on 'get_plots'
  results <- if (get_plots == FALSE) {
    collect
  } else if (informative == FALSE & get_plots == TRUE) {
    append(collect, list(Dissimilarity_comparison = comp_diss_plot,
                         Profile_plot = internal_measures_panel,
                         Silhouette_width_plot = plot_comp_silhouette))
  } else if (informative == TRUE & get_plots == TRUE) {
    append(collect, list(Dissimilarity_comparison = comp_diss_plot,
                         Informative_heatmap = informative_heatmap))
  }

  class(results) <- "comp_clustering"

  return(suppressWarnings({print(results)}))
}
