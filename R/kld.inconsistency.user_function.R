#'  Density plots of local inconsistency results and Kullback-Leibler divergence
#'  (When dataset is created by the user)
#'
#' @description
#' When the user has extracted results obtained from a method of local
#' inconsistency evaluation (e.g., loop-specific, back-calculation or
#' node-splitting approaches) as reported in publication, this function
#' provides the same output with the function \code{\link{kld_inconsistency}}.
#' A panel of density plots on the direct and indirect estimates of the
#' selected comparisons based on approach for local inconsistency evaluation,
#' such as back-calculation and node-splitting approaches (Dias et al., 2010;
#' van Valkenhoef et al., 2016) and loop-specific approach (Bucher et al., 1997)
#' accompanied by the average Kullback-Leibler divergence. Additionally, stacked
#' bar plots on the percentage contribution of either Kullback-Leibler
#' divergence (from direct to indirect, and vice-versa) to the total information
#' loss for each selected comparison are presented.
#'
#' @param dataset A data-frame of seven columns and as many rows as the split
#'   nodes. The first column contains the names of the split nodes, and the
#'   remaining columns have the point estimate and standard error of the direct,
#'   indirect and inconsistency parameter in that order.
#' @param threshold A positive number indicating the threshold of not concerning
#'   inconsistency, that is, the minimally allowed deviation between the direct
#'   and indirect estimates for a split node that does raise concerns for
#'   material inconsistency. The argument is optional.
#' @param level A number indicating the significance level. Suggested values
#'   are 0.05 and 0.10. The default value is 0.05.
#' @param outcome Optional argument to describe the effect measure used (the
#'   x-axis of the plots).
#' @param scales A character on whether both axes should be fixed
#'   (\code{"fixed"}) or free (\code{"free"}) or only one of them be free
#'   (\code{"free_x"} or \code{"free_y"}). \code{scales} determines the scales
#'   argument found in function (\code{\link[ggplot2:facet_wrap]{facet_wrap}})
#'   in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}. The default is
#'   (\code{"free"}).
#' @param show_incons Logical to indicate whether to present the point estimate
#'   and 95% interval of the inconsistency parameter. The default is \code{TRUE}
#'   (report).
#' @param y_axis_name Logical to indicate whether to present the title of y-axis
#'   ('Density'). The default is \code{TRUE} (report).
#' @param title_name Text for the title of the plot. \code{title_name}
#'   determines the labs argument of the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param axis_title_size A positive integer for the font size of axis title.
#'   \code{axis_title_size} determines the axis.title argument found in the
#'   theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'   The default option is 13.
#' @param axis_text_size A positive integer for the font size of axis text.
#'   \code{axis_text_size} determines the axis.text argument found in the
#'   theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'   The default option is 13.
#' @param text_size A positive integer for the font size of labels.
#'   \code{text_size} determines the size argument found in the geom_text
#'   function in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'   The default option is 3.5.
#' @param strip_text_size A positive integer for the font size of facet labels.
#'   \code{legend_text_size} determines the legend.text argument found in
#'   the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'   The default option is 13.
#' @param legend_title_size A positive integer for the font size of legend
#'   title. \code{legend_text_size} determines the legend.text argument found in
#'   the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'   The default option is 13.
#' @param legend_text_size A positive integer for the font size of legend text.
#'   \code{legend_text_size} determines the legend.text argument found in the
#'   theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'   The default option is 13.
#' @param str_wrap_width A positive integer for wrapping the axis labels in the
#'   percent stacked bar-plot. \code{str_wrap_width} determines the
#'   \code{\link[stringr:str_wrap]{str_wrap}} function of the R-package
#'   \href{https://CRAN.R-project.org/package=stringr}{stringr}.
#'
#' @return The first plot is a panel of density plots for each split node sorted
#' in ascending order of the Kullback-Leibler divergence value. Blue and black
#' lines refer to the direct and indirect estimates, respectively. The grey
#' segment refers to the (1 - \code{level})\% 'pseudo' confidence interval of
#' the inconsistency parameter based on the corresponding normal z-scores,
#' with a darker grey line  referring to the point estimate. The names of the
#' selected comparisons appear at the top of each plot. The mean estimate on
#' the scale of the selected effect measure appears at the top of each density
#' curve.
#'
#' The Kullback-Leibler divergence value appears at the top left of each plot
#' in three colours: black, if no threshold has been defined (the default),
#' green, if the Kullback-Leibler divergence is below the specified
#' \code{threshold} (not concerning inconsistency) and red, if the
#' Kullback-Leibler divergence is at least the specified \code{threshold}
#' (substantial inconsistency).
#'
#' The second plot is a percent stacked bar plot on the percentage contribution
#' of approximating direct with indirect estimate (and vice-versa) to the total
#' information loss for each target comparison. Total information loss is
#' defined as the sum of the Kullback-Leibler divergence value
#' when approximating the direct with indirect estimate (blue bars), and the
#' Kullback-Leibler divergence value when approximating the indirect
#' with direct estimate (black bars). Values parentheses refer to the
#' corresponding Kullback-Leibler divergence value. Bars are sorted in ascending
#' order of the average Kullback-Leibler divergence value.
#'
#' The function also returns the data-frame \code{average_KLD} that includes the
#' split comparisons and the corresponding average Kullback-Leibler divergence
#' value.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[ggplot2:facet_wrap]{facet_wrap}},
#'   \code{\link{kld_inconsistency}},
#'   \code{\link{kld_measure}}
#'
#' @references
#' Bucher HC, Guyatt GH, Griffith LE, Walter SD. The results of direct and
#' indirect treatment comparisons in meta-analysis of randomized controlled
#' trials. \emph{J Clin Epidemiol} 1997;\bold{50}(6):683--91.
#'
#' Dias S, Welton NJ, Caldwell DM, Ades AE. Checking consistency in mixed
#' treatment comparison meta-analysis.
#' \emph{Stat Med} 2010;\bold{29}(7-8):932--44.
#' doi: 10.1002/sim.3767
#'
#' Kullback S, Leibler RA. On information and sufficiency.
#' \emph{Ann Math Stat} 1951;\bold{22}(1):79--86. doi: 10.1214/aoms/1177729694
#'
#' van Valkenhoef G, Dias S, Ades AE, Welton NJ. Automated generation of
#' node-splitting models for assessment of inconsistency in network
#' meta-analysis. \emph{Res Synth Methods} 2016;\bold{7}(1):80--93.
#' doi: 10.1002/jrsm.1167
#'
#' @examples
#'
#' \dontrun{
#' ## Data are taken from Table II in Dias et al. (2010)
#' # Treatments compared
#' treat <-
#' c("SK", "t-PA", "Acc t-PA", "SK+t-PA", "r-PA", "TNK", "PTCA", "UK", "ASPAC")
#'
#' # Baseline arm (from each selected comparison)
#' base <- rep(1:3, c(6, 3, 5))
#'
#' # Non-baseline arm (from each selected comparison)
#' nonbase <- c(2, 3, 5, 7, 8, 9, 7, 8, 9, 4, 5, 7, 8, 9)
#'
#' # Compared treatments with their names
#' treat_comp <-
#' mapply(function(x, y) paste(treat[x], "vs", treat[y]), base, nonbase)
#'
#' # Direct results
#' direct_mean <- c(0.000, -0.158, -0.060, -0.666, -0.369, 0.009, -0.545,
#'                  -0.295, 0.006, 0.126, 0.019, -0.216, 0.143, 1.409)
#' direct_sd <- c(0.030, 0.048, 0.089, 0.185, 0.518, 0.037, 0.417, 0.347, 0.037,
#'                0.054, 0.066, 0.118, 0.356, 0.415)
#'
#' # Indirect results
#' indirect_mean <- c(0.189, -0.247, -0.175, -0.393, -0.168, 0.424, -0.475,
#'                    -0.144, 0.471, 0.630, 0.135, -0.477, -0.136, 0.165)
#' indirect_sd <- c(0.235, 0.092, 0.081, 0.120, 0.244, 0.252, 0.108, 0.290,
#'                  0.241, 0.697, 0.101, 0.174, 0.288, 0.057)
#'
#' # Inconsistency
#' incons_mean <- c(-0.190, 0.088, 0.115, -0.272, -0.207, -0.413, -0.073,
#'                  -0.155, -0.468, -0.506, -0.116, 0.260, 0.277, 1.239)
#' incons_sd <- c(0.236, 0.104, 0.121, 0.222, 0.575, 0.253, 0.432, 0.452, 0.241,
#'                0.696, 0.120, 0.211, 0.461, 0.420)
#'
#' # Collect results in a data-frame (exactly as required from the function)
#' dias_results <- data.frame(treat_comp, direct_mean, direct_sd, indirect_mean,
#'                            indirect_sd, incons_mean, incons_sd)
#'
#' # Apply the function
#' kld_inconsistency_user(dataset = dias_results,
#'                        threshold = 0.13,
#'                        outcome = "Odds ratio (logarithmic scale)")
#' }
#'
#' @export
kld_inconsistency_user <- function(dataset,
                                   threshold = 0.00001,
                                   level = 0.05,
                                   outcome = NULL,
                                   scales = "free",
                                   show_incons = TRUE,
                                   y_axis_name = TRUE,
                                   title_name = NULL,
                                   axis_title_size = 13,
                                   axis_text_size = 13,
                                   text_size = 3.5,
                                   strip_text_size = 13,
                                   legend_title_size = 13,
                                   legend_text_size = 13,
                                   str_wrap_width = 10) {


  # General message
  a1 <- "Note: Make sure that you have created the dataset"
  b1 <- "according to the description of the argument 'dataset'."
  message(paste(a1, b1))

  # Default arguments
  dataset <- if (dim(dataset)[2] != 7) {
    aa <- "The argument 'dataset' must have 7 columns; the first column must"
    bb <- "be a character vector. See parameter description and example."
    stop(paste(aa, bb), call. = FALSE)
  } else {
    dataset
  }
  threshold <- if (threshold <= 0) {
    stop("The argument 'threshold' must be a positive number.", call. = FALSE)
  } else if (0 < threshold & threshold <= 0.00001) {
    threshold
  } else if (threshold > 0.00001) {
    message(paste0("Threshold specified at ", threshold, "."))
    threshold
  }
  if (missing(level)) {
    message("Significance level specified at 0.05 (the default).")
  } else {
    level <- level
    message("Significance level specified at ", level)
  }
  scales <- if (missing(scales)) {
    "free"
  } else if (!is.element(scales, c("fixed", "free", "free_x", "free_y"))) {
    stop("Insert one of the following: 'fixed', 'free', 'free_x', or 'free_y'.",
         call. = FALSE)
  } else if (is.element(scales, c("fixed", "free", "free_x", "free_y"))) {
    scales
  }
  y_axis_name <- if (y_axis_name == TRUE) {
    "Density"
  } else {
    ""
  }

  # Function for the Kullback-Leibler Divergence (two normal distributions)
  #kld_measure <- function(mean_y, sd_y, mean_x, sd_x) {
  #  # x is the 'truth' (e.g. the direct estimate)
  #  kld_xy <- 0.5 * (((sd_x / sd_y)^2) + ((mean_y - mean_x)^2)
  #                   / (sd_y^2) - 1 + 2 * log(sd_y / sd_x))
  #
  #  # y is the 'truth' (e.g. the indirect estimate)
  #  kld_yx <- 0.5 * (((sd_y / sd_x)^2) + ((mean_x - mean_y)^2)
  #                   / (sd_x^2) - 1 + 2 * log(sd_x / sd_y))
  #
  #  # Symmetric KLD, also known as J-divergence
  #  sym_kld <- (kld_xy + kld_yx) / 2
  #
  #  return(list(kld_sym = sym_kld,
  #              kld_dir = kld_xy,
  #              kld_ind = kld_yx))
  #}

  # Obtain the Kullback-Leibler Divergence values for each selected comparison
  kld_value <-
    unlist(lapply(1:dim(dataset)[1],
                  function(x) kld_measure(mean_y = dataset[x, 4],
                                          sd_y = dataset[x, 5],
                                          mean_x = dataset[x, 2],
                                          sd_x = dataset[x, 3])$kld_sym))

  # Kullback-Leibler Divergence by approximating direct with indirect
  kld_dir <-
    unlist(lapply(1:dim(dataset)[1],
                  function(x) kld_measure(mean_y = dataset[x, 4],
                                          sd_y = dataset[x, 5],
                                          mean_x = dataset[x, 2],
                                          sd_x = dataset[x, 3])$kld_x_true))

  # Kullback-Leibler Divergence by approximating indirect with direct
  kld_ind <-
    unlist(lapply(1:dim(dataset)[1],
                  function(x) kld_measure(mean_y = dataset[x, 4],
                                          sd_y = dataset[x, 5],
                                          mean_x = dataset[x, 2],
                                          sd_x = dataset[x, 3])$kld_y_true))

  # Bring both divergences together per target comparison
  kld_dataset <-
    data.frame(value = c(kld_dir, kld_ind),
               perc = c(kld_dir, kld_ind) / (kld_dir + kld_ind),
               approx = rep(c("Direct estimate", "Indirect estimate"),
                            each = length(kld_dir)),
               compar = rep(dataset[, 1], 2))

  # Obtain the 0.1th and 99.9th percentile of direct estimates per comparison
  range_dir <-
    lapply(1:dim(dataset)[1],
           function(x) c(qnorm(1 - 0.999, dataset[x, 2], dataset[x, 3]),
                         qnorm(0.999, dataset[x, 2], dataset[x, 3])))

  # Obtain the 0.1th and 99.9th percentile of indirect estimates per comparison
  range_ind <-
    lapply(1:dim(dataset)[1],
           function(x) c(qnorm(1 - 0.999, dataset[x, 4], dataset[x, 5]),
                         qnorm(0.999, dataset[x, 4], dataset[x, 5])))

  # The range based on direct and indirect estimates in a vector per comparison
  range_x0 <- lapply(1:dim(dataset)[1],
                     function(x) range(c(range_dir[[x]], range_ind[[x]])))

  # Using the range create a sequence of values per selected comparison
  range_x <-
    lapply(1:dim(dataset)[1],
           function(x) seq(from = min(range_x0[[x]]),
                           to = max(range_x0[[x]]),
                           by = 0.01))

  # Probability density of direct estimates per selected comparison
  prob_dir <-
    lapply(1:dim(dataset)[1],
           function(x) dnorm(range_x[[x]], dataset[x, 2], dataset[x, 3]))

  # Probability density of indirect estimates per selected comparison
  prob_ind <-
    lapply(1:dim(dataset)[1],
           function(x) dnorm(range_x[[x]], dataset[x, 4], dataset[x, 5]))

  # Kullback-Leibler Divergence between the distributions for each percentile
  diff <-
    lapply(1:dim(dataset)[1],
           function(x) prob_dir[[x]] * log(prob_dir[[x]] / prob_ind[[x]]))

  # The KLD value per selected comparison
  # (repeat as many times as the length 'range_x' of that comparison)
  KLD <- lapply(1:dim(dataset)[1],
                function(x) rep(kld_value[x], length(range_x[[x]])))

  # The names of selected comparison
  # (repeat as many times as the length 'range_x' of that comparison)
  compar <- lapply(1:dim(dataset)[1],
                   function(x) rep(dataset[x, 1], length(range_x[[x]])))

  # Based on argument 'threshold' return a decision regarding (in)consistency
  decision <-
    lapply(1:dim(dataset)[1],
           function(x)
             ifelse(threshold == 0.00001, "No threshold defined",
                    ifelse(KLD[[x]] < threshold & threshold > 0.00001,
                           "Consistency","Inconsistency")))

  # Bring all together per selected comparison
  output0 <-
    lapply(1:dim(dataset)[1],
           function(x) data.frame(time = unlist(range_x[[x]]),
                                  prob_dir = unlist(prob_dir[[x]]),
                                  prob_ind = unlist(prob_ind[[x]]),
                                  diff = unlist(diff[[x]]),
                                  KLD = unlist(KLD[[x]]),
                                  compar = unlist(compar[[x]]),
                                  decision = unlist(decision[[x]])))

  # Bind by row all comparison-specific results
  output <- do.call(rbind, output0)

  # Fake dummy (to add the estimate type, direct and indirect, in the legend)
  estimate <- NULL
  output$estimate <-
    rep(c("Direct estimate", "Indirect estimate"), c(20, dim(output)[1] - 20))

  # Create a new dataset that generates 'pseudo' CrI based on the posterior mean
  # and SD of inconsistency factor
  dataset_new <-
    data.frame(compar = dataset[, 1],
               incons_cri =
                 do.call(rbind,
                         lapply(1:dim(dataset)[1],
                                function(x)
                                  c(dataset[x, 6] +
                                      qnorm(level/2, 0, 1) * dataset[x, 7],
                                    dataset[x, 6] +
                                      qnorm(1 - (level/2), 0, 1) * dataset[x, 7]))),
               mean = dataset[, 6],
               KLD = kld_value)
  colnames(dataset_new)[2:3] <- c("lower", "upper")

  # Maximum density probability of direct and indirect estimates
  max_prob_dir <- unlist(lapply(1:dim(dataset)[1], function(x) max(prob_dir[[x]])))
  max_prob_ind <- unlist(lapply(1:dim(dataset)[1], function(x) max(prob_ind[[x]])))

  # Data-frame with the mean direct and indirect estimates per comparison
  prob <- NULL
  data_mean <- data.frame(mean = c(dataset[, 2], dataset[, 4]),
                          prob = c(max_prob_dir, max_prob_ind),
                          compar = dataset[, 1],
                          source = rep(c("direct", "indirect"),
                                       each = length(dataset[, 1])))

  # Get the panel of density plots
  plot <-
    ggplot() +
    {if (show_incons)
    geom_rect(data = dataset_new,
              aes(xmin = lower,
                  xmax = upper,
                  ymin = 0,
                  ymax = Inf),
              fill = "grey90")} +
    {if (show_incons)
    geom_rect(data = dataset_new,
              aes(xmin = mean,
                  xmax = mean,
                  ymin = 0,
                  ymax = Inf),
              col = "grey75")} +
    geom_area(data = output,
              aes(x = time,
                  y = prob_ind),
              linewidth = 1.3,
              colour = "black",
              fill = "white",
              alpha = 0) +
    geom_area(data = output,
              aes(x = time,
                  y = prob_dir),
              linewidth = 1.3,
              colour = "#0072B2",
              fill = "white",
              alpha = 0) +
    geom_vline(xintercept = 0,
               linetype = 2) +
    geom_hline(yintercept = 0) +
    geom_text(data = output,
              x = -Inf,
              y = Inf,
              aes(label = paste0("D=", sprintf("%.2f", KLD)),
                  col = decision),
              fontface = "bold",
              size = text_size,
              hjust = -0.2,
              vjust = 1.4,
              show.legend = FALSE) +
    geom_text(data = subset(data_mean, source == "direct"),
              aes(x = mean,
                  y = prob,
                  label = sprintf("%.2f", mean),
                  vjust = -0.5),
              fontface = "bold",
              size = text_size,
              show.legend = FALSE,
              inherit.aes = FALSE) +
    geom_text(data = subset(data_mean, source == "indirect"),
              aes(x = mean,
                  y = prob,
                  label = sprintf("%.2f", mean),
                  vjust = -0.5),
              fontface = "bold",
              size = text_size,
              show.legend = FALSE,
              inherit.aes = FALSE) +
    geom_point(data = output,
               aes(x = time,
                   y = prob_dir,
                   fill = estimate),
               alpha = 0) +
    {if (dim(dataset)[1] > 1)
    facet_wrap(~factor(compar,
                       levels = dataset[, 1][order(kld_value,
                                                   decreasing = FALSE)]),
               scales = scales)} +
    scale_colour_manual(values = c("No threshold defined" = "black",
                                   "Consistency" = "#009E73",
                                   "Inconsistency" = "#D55E00")) +
    scale_fill_manual(values = c("Direct estimate" = "#0072B2",
                                 "Indirect estimate" = "black")) +
    labs(x = outcome,
         y = y_axis_name,
         title_name = title_name,
         fill = " ") +
    guides(colour = "none",
           fill = guide_legend(override.aes = list(size = 3,
                                                   alpha = 1,
                                                   colour = c("#0072B2",
                                                              "black")))) +
    scale_y_continuous(expand = c(0.20, 0)) +
    theme_classic() +
    theme(plot.title = element_text(size = axis_title_size, face = "bold"),
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size, face = "bold"),
          strip.text = element_text(size = strip_text_size, face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = legend_text_size),
          legend.title = element_text(size = legend_title_size, face = "bold"))

  # Percent stacked barplot of comparison divergence when approximating direct
  # with indirect estimate and vice-versa
  barplot <-
    ggplot(kld_dataset,
           aes(x = compar,
               y = perc,
               fill = approx)) +
    geom_bar(position = "fill",
             stat = "identity") +
    geom_hline(yintercept = 0.5,
               linetype = "dashed",
               linewidth = 0.8,
               colour = "white") +
    geom_text(aes(label = paste0(sprintf("%.0f", perc * 100),"%", " ",
                                 "(", round(value, 2), ")")),
              hjust = 0.5,
              vjust = 1.0,
              size = text_size,
              position = "stack",
              colour = "white") +
    labs(x = "Target comparisons",
         y = "% contribution to total information loss",
         fill = "Approximating") +
    scale_y_continuous(labels = scales::label_percent(suffix = " ")) +
    scale_x_discrete(labels = function(x) str_wrap(x,
                                                   width = str_wrap_width),
                     limits = dataset[, 1][order(kld_value,
                                                 decreasing = FALSE)]) +
    scale_fill_manual(values = c("Direct estimate" = "#0072B2",
                                 "Indirect estimate" = "black")) +
    theme_classic() +
    theme(plot.title = element_text(size = axis_title_size, face = "bold"),
          axis.text = element_text(size = axis_text_size),
          axis.title = element_text(size = axis_title_size, face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = legend_text_size),
          legend.title = element_text(size = legend_title_size, face = "bold"))

  return(list(Density_plot = plot,
              Barplot = barplot,
              average_KLD = data.frame(comparison, kld_value)))
}
