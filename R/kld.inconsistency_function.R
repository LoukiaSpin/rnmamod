#'  Density plots of local inconsistency results and Kullback-Leibler divergence
#'  when 'rnmamod', 'netmeta' or 'gemtc' R packages are used
#'
#' @description
#' A panel of density plots on the direct and indirect estimates of the
#' selected comparisons based on approach for local inconsistency evaluation,
#' such as back-calculation and node-splitting approaches (Dias et al., 2010;
#' van Valkenhoef et al., 2016) and loop-specific approach (Bucher et al., 1997)
#' accompanied by the average Kullback-Leibler divergence. Additionally, stacked
#' bar plots on the percentage contribution of either Kullback-Leibler
#' divergence (from direct to indirect, and vice-versa) to the total information
#' loss for each selected comparison are presented.
#' The function handles results also from the R-packages
#'   \href{https://CRAN.R-project.org/package=gemtc}{gemtc} and
#'   \href{https://CRAN.R-project.org/package=netmeta}{netmeta}.
#'
#' @param node An object of S3 class \code{\link{run_nodesplit}} or class
#'   \code{\link[gemtc:mtc.nodesplit]{mtc.nodesplit}}
#'   (see \href{https://CRAN.R-project.org/package=gemtc}{gemtc}) or
#'   class \code{\link[netmeta:netsplit]{netsplit}} (see
#'   \href{https://CRAN.R-project.org/package=netmeta}{netmeta}).
#' @param threshold A positive number indicating the threshold of not concerning
#'   inconsistency, that is, the minimally allowed deviation between the direct
#'   and indirect estimates for a split node that does raise concerns for
#'   material inconsistency. The argument is optional.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data}. It is not relevant for
#'   \href{https://CRAN.R-project.org/package=gemtc}{gemtc} and
#'   \href{https://CRAN.R-project.org/package=netmeta}{netmeta}.
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
#' segment refers to the 95\% credible (confidence) interval of the
#' inconsistency parameter, when \code{\link{run_nodesplit}}
#' (\code{\link[netmeta:netsplit]{netsplit}}) has been applied, with a darker
#' grey line referring to the point estimate.
#' When \code{\link[gemtc:mtc.nodesplit]{mtc.nodesplit}} has been employed, the
#' 95\% confidence interval has been approximated using the Bucher's approach
#' based on the corresponding direct and indirect results. This was necessary
#' because \code{\link[gemtc:mtc.nodesplit]{mtc.nodesplit}} (version 1.0-2)
#' returns only the inconsistency p-values rather than the posterior results on
#' the inconsistency parameters. The mean estimate on
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
#' defined as the sum of the KLD value when approximating the direct with
#' indirect estimate (blue bars), and the KLD when approximating the indirect
#' with direct estimate (black bars). Values parentheses refer to the
#' corresponding KLD value.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link[ggplot2:facet_wrap]{facet_wrap}},
#'   \code{\link[gemtc:mtc.nodesplit]{mtc.nodesplit}},
#'   \code{\link{kld_measure}},
#'   \code{\link[netmeta:netsplit]{netsplit}},
#'   \code{\link{run_nodesplit}},
#'   \code{\link[stringr:str_wrap]{str_wrap}}
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
#' data("nma.baker2009")
#'
#' # Read results from 'run_nodesplit' (using the default arguments)
#' node <- readRDS(system.file('extdata/node_baker.rds', package = 'rnmamod'))
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "budesonide", "budesonide plus formoterol",
#'                   "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium")
#'
#' # Apply the function
#' kld_inconsistency(node = node,
#'                   threshold = 0.23,
#'                   drug_names = interv_names,
#'                   outcome = "Odds ratio (logarithmic scale)",
#'                   str_wrap_width = 15)
#' }
#'
#' @export
kld_inconsistency <- function(node,
                              threshold = 0.00001,
                              drug_names = NULL,
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

  if (all(c(inherits(node, "run_nodesplit"), inherits(node, "mtc.nodesplit"),
            inherits(node, "netsplit")) == FALSE)) {
    aa <- "'node' must be an object of S3 class 'run_nodesplit' (rnmamod R package)"
    bb <- "'mtc.nodesplit' (gemtc R package) or 'netsplit' (netmeta R package)"
    stop(paste(aa, bb), call. = FALSE)
  }

  # Default arguments
  threshold <- if (threshold <= 0) {
    stop("The argument 'threshold' must be a positive number.", call. = FALSE)
  } else if (0 < threshold & threshold <= 0.00001) {
    threshold
  } else if (threshold > 0.00001) {
    message(paste0("Threshold specified at ", threshold, "."))
    threshold
  }
  drug_names <- if (is.null(drug_names) & inherits(node, "run_nodesplit")) {
    aa <- "The argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as specified in 'data' is used, instead."
    message(paste(aa, bb))
    as.character(1:max(unlist(node$direct[, 1:2])))
  } else if (inherits(node, "run_nodesplit")) {
    drug_names
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

  # Extract results based on the class (and hence, the R package)
  if (inherits(node, "mtc.nodesplit")) { # R package: gemtc

    # Get direct estimates
    direct0 <- summary(node)[[1]]

    # Restrict to direct mean
    direct_mean <- direct0[, 3]

    # Get indirect estimates
    indirect0 <- summary(node)[[2]]

    # Restrict to indirect mean
    indirect_mean <- indirect0[, 3]

    # Number of split nodes
    num_nodes <- length(names(node))

    # Remove 'consistency' node
    node_new <- node[-(num_nodes)]

    # Get mcmc summaries per split node
    summary_res <- lapply(node_new, function(x) summary(coda::as.mcmc.list(x)))

    # Extract the posterior SD for direct estimate
    direct_se <-
      unlist(lapply(summary_res,
                    function(x)
                      as.data.frame(x$statistics)[dim(x$statistics)[1] - 3, 2]))

    # Extract the posterior SD for indirect estimate
    indirect_se <-
      unlist(lapply(summary_res,
                    function(x)
                      as.data.frame(x$statistics)[dim(x$statistics)[1] - 2, 2]))

    # Bring direct mean and standard error in a data-frame
    direct <- data.frame(direct_mean, direct_se)

    # Bring indirect mean and standard error in a data-frame
    indirect <- data.frame(indirect_mean, indirect_se)

    # Calculate (approximately) the inconsistency using Bucher's approach
    # 'gemtc' returns only the inconsistency p-value
    incons_mean <- direct_mean - indirect_mean

    # Calculate (approximately) the inconsistency standard error
    incons_se <- (direct_se^2) + (indirect_se^2)

    # Calculate (approximately) the inconsistency lower bound (95%)
    incons_lower <- incons_mean - 1.96 * incons_se

    # Calculate (approximately) the inconsistency upper bound (95%)
    incons_upper <- incons_mean + 1.96 * incons_se

    # Bring inconsistency mean and CI bounds in a data-frame
    inconsistency <- data.frame(incons_mean, incons_lower, incons_upper)

    # Vector of comparison names
    comparison <- paste(indirect0$t2, "vs", indirect0$t1)

    message("Results on inconsistency parameters are calculated approximately.")

  } else if (inherits(node, "netsplit")) { # R package: netmeta

    # Direct results
    direct_res <- node$direct.random

    # Bring direct mean and standard error in a data-frame
    direct0 <- data.frame(direct_res$TE, direct_res$seTE)

    # Indirect results
    indirect_res <- node$indirect.random

    # Bring indirect mean and standard error in a data-frame
    indirect0 <- data.frame(indirect_res$TE, indirect_res$seTE)

    # Inconsistency results
    incons_res <- node$compare.random

    # Bring inconsistency mean and CI bounds in a data-frame
    incons0 <- data.frame(incons_res$TE, incons_res$lower, incons_res$upper)

    # Find the rows that correspond to non-split treatments
    row_na <- which(is.na(incons0[, 1]))

    # Data-frame with direct results after removing non-split treatments
    direct <- direct0[-row_na, ]

    # Data-frame with indirect results after removing non-split treatments
    indirect <- indirect0[-row_na, ]

    # Data-frame with inconsistency results after removing non-split treatments
    inconsistency <- incons0[-row_na, ]

    # Vector of comparison names (after removing non-split treatments)
    comparison0 <- node$comparison[-row_na]
    comparison <- gsub(":"," vs ", comparison0)

  } else if (inherits(node, "run_nodesplit")) { # R package: rnmamod

    # Direct effects (node split)
    direct <- node$direct[, 3:4]

    # Inirect effects (node split)
    indirect <- node$indirect[, 3:4]

    # Inconsistency factor
    inconsistency <- node$diff[, c(3, 5:6)]

    # Comparisons
    split_nodes0 <- node$direct[, 1:2]

    # Interventions' name: Replace code with original names
    # (only when the argument 'drug_names' has been defined)
    first_arm <- unlist(lapply(1:dim(split_nodes0)[1],
                               function(x) drug_names[split_nodes0[x, 1]]))
    second_arm <- unlist(lapply(1:dim(split_nodes0)[1],
                                function(x) drug_names[split_nodes0[x, 2]]))
    split_nodes <- cbind(first_arm, second_arm)

    # Vector of comparison names
    comparison <- paste(split_nodes[, 1], "vs", split_nodes[, 2])
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
    unlist(lapply(1:dim(direct)[1],
                  function(x) kld_measure(mean_y = indirect[x, 1],
                                          sd_y = indirect[x, 2],
                                          mean_x = direct[x, 1],
                                          sd_x = direct[x, 2])$kld_sym))

  # Kullback-Leibler Divergence by approximating direct with indirect
  kld_dir <-
    unlist(lapply(1:dim(direct)[1],
                  function(x) kld_measure(mean_y = indirect[x, 1],
                                          sd_y = indirect[x, 2],
                                          mean_x = direct[x, 1],
                                          sd_x = direct[x, 2])$kld_x_true))

  # Kullback-Leibler Divergence by approximating indirect with direct
  kld_ind <-
    unlist(lapply(1:dim(direct)[1],
                  function(x) kld_measure(mean_y = indirect[x, 1],
                                          sd_y = indirect[x, 2],
                                          mean_x = direct[x, 1],
                                          sd_x = direct[x, 2])$kld_y_true))

  # Bring both divergences together per target comparison
  kld_dataset <-
    data.frame(value = c(kld_dir, kld_ind),
               perc = c(kld_dir, kld_ind) / (kld_dir + kld_ind),
               approx = rep(c("Direct estimate", "Indirect estimate"),
                            each = length(kld_dir)),
               compar = rep(comparison, 2))

  # Obtain the 0.1th and 99.9th percentile of direct estimates per comparison
  range_dir <-
    lapply(1:dim(direct)[1],
           function(x) c(qnorm(1 - 0.999, direct[x, 1], direct[x, 2]),
                         qnorm(0.999, direct[x, 1], direct[x, 2])))

  # Obtain the 0.1th and 99.9th percentile of indirect estimates per comparison
  range_ind <-
    lapply(1:dim(indirect)[1],
           function(x) c(qnorm(1 - 0.999, indirect[x, 1], indirect[x, 2]),
                         qnorm(0.999, indirect[x, 1], indirect[x, 2])))

  # The range based on direct and indirect estimates in a vector per comparison
  range_x0 <- lapply(1:dim(direct)[1],
                     function(x) range(c(range_dir[[x]], range_ind[[x]])))

  # Using the range create a sequence of values per selected comparison
  range_x <-
    lapply(1:dim(direct)[1],
           function(x) seq(from = min(range_x0[[x]]),
                           to = max(range_x0[[x]]),
                           by = 0.01))

  # Probability density of direct estimates per selected comparison
  prob_dir <-
    lapply(1:dim(direct)[1],
           function(x) dnorm(range_x[[x]], direct[x, 1], direct[x, 2]))

  # Probability density of indirect estimates per selected comparison
  prob_ind <- lapply(1:dim(indirect)[1],
                     function(x) dnorm(range_x[[x]], indirect[x, 1], indirect[x, 2]))

  # Kullback-Leibler Divergence between the distributions for each percentile
  diff <- lapply(1:dim(direct)[1], function(x) prob_ind[[x]] - prob_dir[[x]])

  # The KLD value per selected comparison
  # (repeat as many times as the length 'range_x' of that comparison)
  KLD <- lapply(1:dim(direct)[1],
                function(x) rep(kld_value[x], length(range_x[[x]])))

  # The names of split node
  # (repeat as many times as the length 'range_x' of that comparison)
  compar <- lapply(1:dim(direct)[1],
                   function(x) rep(comparison[x],
                                   length(range_x[[x]])))

  # Based on argument 'threshold' return a decision regarding (in)consistency
  decision <-
    lapply(1:dim(direct)[1],
           function(x)
             ifelse(threshold == 0.00001, "No threshold defined",
                    ifelse(KLD[[x]] < threshold & threshold > 0.00001,
                           "Consistency","Inconsistency")))

  # Bring all together per selected comparison
  output0 <-
    lapply(1:dim(direct)[1],
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

  # Create a new dataset that generates 'pseudo' CrI based on the posterior mean and SD of inconsistency factor
  dataset_new <- data.frame(compar = comparison,
                            mean = inconsistency[, 1],
                            lower = inconsistency[, 2],
                            upper = inconsistency[, 3])

  # Maximum density probability of direct and indirect estimates
  max_prob_dir <- unlist(lapply(1:dim(direct)[1], function(x) max(prob_dir[[x]])))
  max_prob_ind <- unlist(lapply(1:dim(direct)[1], function(x) max(prob_ind[[x]])))

  # Data-frame with the mean direct and indirect estimates per comparison
  prob <- NULL
  data_mean <- data.frame(mean = c(direct[, 1], indirect[, 1]),
                          prob = c(max_prob_dir, max_prob_ind),
                          compar = unique(unlist(compar)),
                          source = rep(c("direct", "indirect"),
                                       each = length(direct[, 1])))

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
    {if (length(kld_value) > 1)
    facet_wrap(~factor(compar,
                       levels =
                         comparison[order(kld_value, decreasing = FALSE)]),
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
                     limits = comparison[order(kld_value,
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
              Barplot = barplot))
}
