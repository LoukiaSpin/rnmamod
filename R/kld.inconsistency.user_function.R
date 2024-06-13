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
#' accompanied by the Kullback-Leibler divergence from the indirect to direct
#' estimate.
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
#'
#' @return A panel of density plots for each split node sorted in ascending
#' order of the Kullback-Leibler divergence value. Blue and black lines refer to
#' the direct and indirect estimates, respectively. The grey segment refers to
#' the (1 - \code{level})\% 'pseudo' confidence interval of the inconsistency
#' parameter based on the corresponding normal z-scores, with a darker grey line
#' referring to the point estimate. The names of the selected comparisons appear
#' at the top of each plot.
#'
#' The Kullback-Leibler divergence value appears at the top left of each plot
#' in three colours: black, if no threshold has been defined (the default),
#' green, if the Kullback-Leibler divergence is below the specified
#' \code{threshold} (not concerning inconsistency) and red, if the
#' Kullback-Leibler divergence is at least the specified \code{threshold}
#' (substantial inconsistency).
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{kld_inconsistency}}
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
#'
#' @export
kld_inconsistency_user <- function(dataset,
                                   threshold = 0.00001,
                                   level = 0.05,
                                   outcome = NULL) {


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
  #threshold <- if (missing(threshold)) {
  #  stop("The argument 'threshold' has not been defined.", call. = FALSE)
  #} else {
  #  message(paste0("Threshold specified at ", threshold, "."))
  #  threshold
  #}
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

  # Function for the Kullback-Leibler Divergence (two normal distributions)
  kld_measure_univ <- function(mean_y, sd_y, mean_x, sd_x) {
    # x is the 'truth' (e.g. the direct estimate)
    kld_xy <- 0.5 * (((sd_x / sd_y)^2) + ((mean_y - mean_x)^2)
                     / (sd_y^2) - 1 + 2 * log(sd_y / sd_x))

    return(kld_xy)
  }

  # Obtain the Kullback-Leibler Divergence values for each selected comparison
  kld_value <-
    unlist(lapply(1:dim(dataset)[1],
                  function(x) kld_measure_univ(mean_y =  dataset[x, 4],
                                               sd_y =  dataset[x, 5],
                                               mean_x = dataset[x, 2],
                                               sd_x = dataset[x, 3])))

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

  # Get the panel of density plots
  plot <-
    ggplot() +
    geom_rect(data = dataset_new,
              aes(xmin = lower,
                  xmax = upper,
                  ymin = 0,
                  ymax = Inf),
              fill = "grey90") +
    geom_rect(data = dataset_new,
              aes(xmin = mean,
                  xmax = mean,
                  ymin = 0,
                  ymax = Inf),
              col = "grey75") +
    #geom_area(data = output_fin,
    #          aes(x = time,
    #              y = diff),
    #         linewidth = 1.0,
    #          fill = "lightblue",
    #          alpha = 0.75) +
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
              aes(label = paste0("KLD=", sprintf("%.2f", KLD)),
                  col = decision),
              fontface = "bold",
              size = 3.5,
              hjust = -0.2,
              vjust = 1.4,
              show.legend = FALSE) +
    geom_point(data = output,
               aes(x = time,
                   y = prob_dir,
                   fill = estimate),
               alpha = 0) +
    facet_wrap(~factor(compar,
                       levels = dataset[, 1][order(kld_value,
                                                   decreasing = FALSE)]),
               scales = "free") +
    scale_colour_manual(values = c("No threshold defined" = "black",
                                   "Consistency" = "#009E73",
                                   "Inconsistency" = "#D55E00")) +
    scale_fill_manual(values = c("Direct estimate" = "#0072B2",
                                 "Indirect estimate" = "black")) +
    labs(x = outcome,
         y = " ",
         fill = " ") +
    guides(colour = "none",
           fill = guide_legend(override.aes = list(size = 3,
                                                   alpha = 1,
                                                   colour = c("#0072B2",
                                                              "black")))) +
    theme_classic() +
    theme(axis.text = element_text(size = 13),
          axis.title = element_text(size = 13, face = "bold"),
          strip.text = element_text(size = 13, face = "bold"),
          legend.position = "bottom",
          legend.text = element_text(size = 13),
          legend.title = element_text(size = 13, face = "bold"))

  return(plot)
}
