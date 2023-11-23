#' Visualising missing data in characteristics
#' (Comparisons' comparability for transitivity evaluation)
#'
#' @description
#'   \code{miss_characteristics} hosts a set of visualisation tools to assess
#'   the size and pattern of missing characteristics values in the dataset.
#'
#' @param input A data-frame in the long arm-based format. Two-arm trials occupy
#'   one row in the data-frame. Multi-arm trials occupy as many rows as the
#'   number of possible comparisons among the interventions. The first two
#'   columns refer to the trial name, and the pairwise comparison,
#'   respectively. The remaining columns refer to summary characteristics. See
#'   'Details' for the specification of the columns.
#' @param drug_names A vector of labels with the name of the interventions
#'   in the order they appear in the argument \code{input}.
#' @param rename_char A list of two elements: (i) a numeric vector with the
#'   position of the characteristics in \code{input}, and (ii) a character
#'   vector with the names of the characteristics, as they are wished to
#'   appear in the title of the plots. This argument is optional, in case the
#'   user wants to control the appearance of the titles.
#' @param label_size A positive integer for the font size of labels in the
#'   plots. \code{label_size} determines the size argument found in the geom's
#'   aesthetic properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param axis_title_size A positive integer for the font size of axis titles in
#'   the plots. \code{axis_title_size} determines the axis.title argument found
#'   in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param axis_text_size A positive integer for the font size of axis text in
#'   the plots. \code{axis_text_size} determines the axis.text argument found in
#'   the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param axis_x_text_angle A positive integer for the angle of axis text in
#'   plots related to missing data. \code{axis_text_angle} determines the
#'   axis.text.x argument found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param legend_text_size A positive integer for the font size of legend text
#'   in the plots. \code{legend_text_size} determines the legend.text argument
#'   found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param legend_title_size A positive integer for the font size of legend title
#'   in the plots. \code{legend_title_size} determines the legend.title argument
#'   found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param strip_text_size A positive integer for the font size of strip text
#'   in the plots. \code{strip_text_size} determines the strip.text argument
#'   found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#' @param strip_text_angle A positive integer for the angle of strip text
#'   in the plots. \code{strip_text_angle} determines the strip.text argument
#'   found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}).
#'
#' @return
#'   \code{miss_characteristics} returns the following list of elements:
#'   \item{Barplot_missing_combined}{A panel of barplots on the percentage of
#'   missing and observed cases for each comparison and characteristic.}
#'   \item{Barplot_missing_characteristics}{A barplot on the percentage of
#'   missing and observed cases for each comparison.}
#'   \item{Tileplot_missing}{A plot that illustrates the position of missing
#'   cases for each trial, comparison and characteristic.}
#'
#' @details
#'   The correct type mode of columns in \code{input} must be ensured to use
#'   the function \code{miss_characteristics}. The first two columns referring
#'   to the trial name, and pairwise comparison, respectively, must be
#'   \strong{character}. The remaining columns referring to the characteristics
#'   must be \strong{double} or \strong{integer} depending on whether the
#'   corresponding characteristic refers to a quantitative or qualitative
#'   variable. The type mode of each column is assessed by
#'   \code{miss_characteristics} using the base function \code{typeof}.
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
#' @export
miss_characteristics <- function (input,
                                  drug_names,
                                  rename_char = NULL,
                                  label_size = 4,
                                  axis_title_size = 14,
                                  axis_text_size = 14,
                                  axis_x_text_angle = 0,
                                  legend_text_size = 14,
                                  legend_title_size = 14,
                                  strip_text_size = 14,
                                  strip_text_angle = 0) {


  ## Check if the dataset is correct
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


  ## Assign the intervention names (if applicable)
  input <- input0
  input[, 2:3] <- matrix(drug_names[as.numeric(unlist(input[, 2:3]))],
                         nrow = dim(input)[1],
                         ncol = 2)


  ## Insert 'Comparison' in the dataset (control appears second in the compar.)
  input$Comparison <- as.character(paste(input$Arm2, "vs", input$Arm1))


  ## Reduce dataset to trial, comparison & characteristics
  input_new <- input[, c(1, dim(input)[2], 4:(dim(input)[2] - 1))]


  ## Rename columns if indicated
  if (!is.null(rename_char)) {
    colnames(input_new)[rename_char[[1]] - 1] <- rename_char[[2]]
  }


  ## Variable on sample size
  colnames(input_new)[with(input_new,
                           startsWith(names(input_new),
                                      c("sample", "Sample")))] <- "Sample size"


  ## Dataframe with missing and observed data per characteristic
  yes <- colSums(is.na(input_new[, -c(1, 2)]))
  no <- colSums(!is.na(input_new[, -c(1, 2)]))
  data_mod0 <- melt(rbind(yes, no))
  data_mod0$perc <- round((data_mod0$value/ dim(input_new)[1]) * 100, 2)
  data_mod0$value2 <- ifelse(data_mod0$value == 0, NA, data_mod0$value)
  data_mod0$perc2 <- ifelse(data_mod0$value == 0, NA, data_mod0$perc)
  colnames(data_mod0)[1:4] <- c("missing", "char", "value", "perc")


  ## If dataset is large, present only baplots where missing data exist
  data_mod <- if (length(unique(data_mod0$char)) > 15) {
    subset(data_mod0, perc > 0 & perc < 100)
  } else {
    data_mod0
  }


  ## Stacked barplot with observed and missing data
  barplot_mod_char <-
    ggplot(data_mod,
           aes(x = char,
               y = value,
               fill = missing)) +
    geom_bar(position = "stack",
             stat = "identity") +
    geom_text(aes(label = ifelse(!is.na(perc2),
                                 paste0(perc2, "%", " ", "(",
                                        value2, ")"), " ")),
              hjust = 0.5,
              vjust = 1.0, # -0.2
              size = label_size,
              position = "stack",
              colour = "white") +
    scale_fill_manual(breaks = c("yes", "no"),
                      labels = c("Yes", "No"),
                      values = c("black", "grey70")) +
    scale_y_continuous(breaks = seq(1, dim(input_new)[1],
                                    ceiling(sqrt(dim(input_new)[1]))),
                       labels = seq(1, dim(input_new)[1],
                                    ceiling(sqrt(dim(input_new)[1]))),
                       limits = c(0, dim(input_new)[1])) +
    labs(x = "Characteristics",
         y = "Counts",
         fill = "Missing") +
    coord_cartesian(expand = FALSE) +
    theme_classic() +
    theme(axis.title = element_text(size = axis_title_size, face = "bold"),
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = axis_x_text_angle,
                                     hjust =
                                       ifelse(axis_x_text_angle == 0, 0.5, 1)),
          legend.position = "bottom",
          legend.text = element_text(size = legend_text_size),
          legend.title = element_text(size = legend_title_size, face = "bold"))


  ## Split 'dataset' by 'Comparison'
  split_dataset <- split(input_new, f = input_new$Comparison)


  ## Dataframe with missing and observed data per comparison & characteristic
  # Absolute frequencies
  yes_char_count_comp <-
    do.call(rbind,
            lapply(split_dataset,
                   function(x) {colSums(is.na(x))}))[, -c(1, 2)]
  no_char_count_comp <-
    do.call(rbind,
            lapply(split_dataset,
                   function(x) {colSums(!is.na(x))}))[, -c(1, 2)]
  data_mod_char_count_comp <-
    melt(rbind(yes_char_count_comp, no_char_count_comp))

  # Relative frequencies
  yes_char_comp <-
    do.call(rbind,
            lapply(split_dataset,
                   function(x) {colSums(is.na(x)) / dim(x)[1]}))[, -c(1, 2)]
  no_char_comp <-
    do.call(rbind,
            lapply(split_dataset,
                   function(x) {colSums(!is.na(x)) / dim(x)[1]}))[, -c(1, 2)]
  data_mod_char_comp <- melt(rbind(yes_char_comp, no_char_comp))
  data_mod_char_comp$value2 <-
    ifelse(data_mod_char_comp$value == 0, NA, data_mod_char_count_comp$value)
  data_mod_char_comp$missing <-
    rep(c("yes", "no"), each = dim(yes_char_comp)[1])
  data_mod_char_comp$value3 <-
    ifelse(data_mod_char_comp$value == 0, NA,
           round(data_mod_char_comp$value * 100, 2))
  colnames(data_mod_char_comp)[1:6] <-
    c("compar", "char", "perc", "value2", "missing", "perc2")


  ## Stacked barplot with observed & missing data (Characteristic by comparison)
  barplot_mod_comp_char <-
    ggplot(data_mod_char_comp,
           aes(x = "",
               y = perc * 100,
               fill = factor(missing, levels = c("yes", "no")))) +
    geom_bar(position = "stack",
             stat = "identity") +
    geom_text(aes(label =
                    ifelse(!is.na(perc2) & missing == "no",
                           paste0(perc2, "%", " ", "(",
                                  value2, ")"), " ")),
              hjust = 0.5,
              vjust = 1.0, # -0.2
              size = label_size,
              position = "stack",
              colour = "white") +
    facet_grid(vars(compar), vars(char)) +
    scale_fill_manual(breaks = c("yes", "no"),
                      labels = c("Yes", "No"),
                      values = c("black", "grey70")) +
    labs(x = " ",
         y = "Percentage (%)",
         fill = "Missing") +
    coord_cartesian(expand = FALSE) +
    theme_classic() +
    theme(axis.title = element_text(size = axis_title_size, face = "bold"),
          axis.text = element_text(size = axis_text_size),
          legend.position = "bottom",
          legend.text = element_text(size = legend_text_size),
          legend.title = element_text(size = legend_title_size, face = "bold"),
          strip.text = element_text(size = strip_text_size, face = "bold"),
          strip.text.y.right = element_text(angle = strip_text_angle),
          axis.ticks.x = element_blank())


  ## Indicate missing cases for each trial and characteristic
  data_mod_dummy0 <- ifelse(is.na(input_new[, -c(1, 2)]), "yes", "no")
  rownames(data_mod_dummy0) <- input_new[, 1]
  data_mod_dummy <- melt(data_mod_dummy0)
  data_mod_dummy$comparison <- rep(input_new[, 2], dim(input_new)[2] - 2)
  colnames(data_mod_dummy)[1:3] <- c("trial", "char", "missing")


  ## Percentage total missing data
  total_mod <-
    round((sum(is.na(input_new[, -c(1, 2)]) == TRUE) /
             (dim(input_new[, -c(1, 2)])[1] *
                dim(input_new[, -c(1, 2)])[2])) * 100, 2)


  ## Tiles with missing and observed data for each trial and characteristic
  tileplot_mod <-
    ggplot(data_mod_dummy,
           aes(x = char,
               y = trial)) +
    geom_tile(aes(fill = missing)) +
    scale_fill_manual(breaks = c("yes", "no"),
                      labels =
                        c(paste0("Yes", " ", "(", total_mod, "%", ")"),
                          paste0("No", " ", "(", 100 - total_mod, "%", ")")),
                      values = c("black", "grey70"),
                      limits = c("yes", "no"),
                      drop = FALSE) +
    #scale_x_discrete(labels = function(x) str_wrap(x, width = 2)) +
    scale_y_discrete(guide = guide_axis(check.overlap = TRUE))+
    facet_grid(vars(comparison),
               scales = "free",
               space = "free_y") +
    labs(x = "Characteristics",
         y = "Study",
         fill = "Missing") +
    theme_classic() +
    theme(axis.title = element_text(size = axis_title_size, face = "bold"),
          axis.text = element_text(size = axis_text_size),
          axis.text.x = element_text(angle = axis_x_text_angle,
                                     hjust =
                                       ifelse(axis_x_text_angle == 0, 0.5, 1)),
          legend.position = "bottom",
          legend.text = element_text(size = legend_text_size),
          legend.title = element_text(size = legend_title_size, face = "bold"),
          strip.text = element_text(size = strip_text_size, face = "bold"),
          strip.text.y.right = element_text(angle = strip_text_angle))


  ## Collect results
  results <- list(Barplot_characteristics = barplot_mod_char,
                  Barplot_combined = barplot_mod_comp_char,
                  Tileplot = tileplot_mod)

  return(results)
}
