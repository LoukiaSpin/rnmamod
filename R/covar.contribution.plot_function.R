#' Visualising study percentage contributions against a covariate
#'
#' @description
#' A scatter plot of the study percentage contributions against the values of a
#' continuous study-level covariate for the treatment effects of comparisons
#' referring to the basic parameters, functional parameters or both.
#' Contributions on the estimated regression coefficients are also presented.
#' Study percentage contributions are based on the proposed methodology of
#' Donegan and colleagues (2018).
#'
#' @param contr_res An object of S3 class \code{\link{study_perc_contrib}}. This
#'   object contains the study percentage contributions to the treatment effects
#'   (or regression coefficients, if relevant) of all possible comparisons in
#'   the network. See 'Value' in \code{\link{study_perc_contrib}}.
#' @param comparisons Character string indicating the type of comparisons to
#'   plot, with possible values: \code{"basic"}, \code{"functional"}, or
#'   \code{"all"} to consider only the basic parameters, only the functional
#'   parameters, or both, respectively. The default argument is \code{"basic"}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{contr_res}.
#'   If \code{drug_names} is not defined, the order of the interventions as
#'   they appear in \code{contr_res} is used, instead.
#' @param upper_limit A positive number to define the upper bound of range of
#'   percentage values for the y-axis. The default argument is 100.
#' @param name_x_axis Text for the x axis title through the \code{labs} function
#'   found in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param axis_title_size A positive integer for the font size of x axis title.
#'   \code{axis_title_size} determines the axis.title (and legend.title)
#'   arguments found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param axis_text_size A positive integer for the font size of axis text (both
#'   axes). \code{axis_text_size} determines the axis.text (and legend.text)
#'   arguments found in the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param strip_text_size A positive integer for the font size of strip text in
#'   facets. \code{strip_text_size} determines the strip.text argument found in
#'   the theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param subtitle_size A positive integer for the font size of subtitle.
#'   \code{subtitle_size} determines the plot.subtitle argument found in the
#'   theme's properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param label_size A positive integer for the font size of labels appearing on
#'   each data point. \code{label_size} determines the size argument found in
#'   the geom's aesthetic properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#' @param seq_by A positive integer for the sequence of values in the x-axis.
#'   \code{seq_by} appears in the arguments breaks and labels found in the
#'   scale_x_continuous aesthetic properties in the R-package
#'   \href{https://CRAN.R-project.org/package=ggplot2}{ggplot2}.
#'
#' @return If interest lies only on the study percentage contributions to the
#' summary treatment effects of all possible pairwise comparisons, the function
#' returns one plot named 'plot_treat'. If interest lies also on the study
#' percentage contributions to the regression coefficient(s), the function
#' returns also the plot named 'plot_reg'.
#'
#' @details
#' A panel of scatter plots is returned on the study percentage contributions to
#' the treatment effects (and also regression coefficients, if relevant) against
#' a continuous covariate for each comparison defined by the argument
#' \code{comparisons}; namely, only those referring to the basic or functional
#' parameters or all possible pairwise comparisons. Blue and red points indicate
#' the studies investigating the corresponding comparisons directly and
#' indirectly, respectively. Each point displays the number of the corresponding
#' study in the dataset.
#'
#' If interest also lies on the study percentage contributions to the regression
#' coefficients, the regression coefficients can be determined to be common
#' across the comparisons, independent or exchangeable and this assumption is
#' specified in the \code{\link{study_perc_contrib}} function.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{study_perc_contrib}}
#'
#' @references
#' Donegan S, Dias S, Tudur-Smith C, Marinho V, Welton NJ. Graphs of study
#' contributions and covariate distributions for network meta-regression.
#' \emph{Res Synth Methods} 2018;\bold{9}(2):243--60. doi: 10.1002/jrsm.1292
#'
#' @examples
#'
#' \dontrun{
#' data("nma.fluoride.donegan2018")
#'
#' # Get study contributions to random-effects network meta-regression
#' # results under the assumption of independent treatment-by-covariate
#' # interaction
#' res <- study_perc_contrib(study_name = nma.fluoride.donegan2018$study,
#'                           base_t = nma.fluoride.donegan2018$t1,
#'                           exp_t = nma.fluoride.donegan2018$t2,
#'                           ref_t = 1,
#'                           obs_se = nma.fluoride.donegan2018$SE,
#'                           obs_cov = nma.fluoride.donegan2018$Cov,
#'                           covar = nma.fluoride.donegan2018$year,
#'                           covar_assum = "independent",
#'                           model = "RE",
#'                           tau = sqrt(0.03))
#'
#' # Covariate-contribution plot on the basic parameters only
#' covar_contribution_plot(contr_res = res,
#'                         comparisons = "basic",
#'                         drug_names = c("NT", "PL", "DE", "RI", "GE", "VA"),
#'                         upper_limit = 15,
#'                         name_x_axis = "Randomisation year",
#'                         seq_by = 10)
#' }
#'
#' @export
covar_contribution_plot <- function (contr_res,
                                     comparisons = "basic",
                                     drug_names,
                                     upper_limit = 100,
                                     name_x_axis = NULL,
                                     axis_title_size = 14,
                                     axis_text_size = 14,
                                     strip_text_size = 14,
                                     subtitle_size = 14,
                                     label_size = 4,
                                     seq_by = 0.1) {


  ## Default arguments
  if (!inherits(contr_res, "study_perc_contrib")) {
    stop("'contr_res' must be an object of S3 class 'study_perc_contrib'.",
         call. = FALSE)
  }
  comparisons <- if (!is.element(comparisons,
                         c("basic", "functional", "all"))) {
    stop("Insert one of the following: 'basic', 'functional', or 'all'.",
         call. = FALSE)
  } else {
    comparisons
  }
  drug_names <- if (missing(drug_names)) {
    aa <- "The argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as provided in 'contr_res' is used, instead."
    message(paste(aa, bb))
    as.character(seq_len(max(unique(unlist(contr_res$perc_contribute[, 2:3])))))
  } else {
    drug_names
  }
  upper_limit <- if (upper_limit > 100 || upper_limit < 0) {
    stop("'upper_limit' must be a number from 0 to 100", call. = FALSE)
  } else {
    upper_limit
  }

  ## Capture results on study percentage contribution
  study_contr <- contr_res$perc_contribute

  ## Distinguish between treatment effects & regression coefficients
  # Treatment effects
  suppressMessages({
    treat_effect <- melt(study_contr[startsWith(colnames(study_contr), "d")])
    })

  # Regression coefficients
  suppressMessages({
    reg_coeff <- if (contr_res$covar_assumption == "common") {
      study_contr[startsWith(colnames(study_contr), "beta")]
    } else if (is.element(contr_res$covar_assumption,
                          c("independent", "exchangenable"))) {
      melt(study_contr[startsWith(colnames(study_contr), "beta")])
    }
  })


  ## Prepare datasets for ggplot2 (*All* treatment effects)
  # Number of comparisons (basic & functional)
  comp_length <- dim(study_contr[startsWith(colnames(study_contr), "d")])[2]

  # Numnber of basic parameters
  basic_comp_length <- dim(study_contr[startsWith(colnames(study_contr), "d1")])[2]

  # Indicate the comparison type (basic/functional)
  comp_indic <- rep(c("basic", "functional"),
                    c(basic_comp_length * dim(study_contr)[1],
                      (comp_length - basic_comp_length) * dim(study_contr)[1]))

  # Study specific comparison
  colour_study_basic0 <- rep(paste0(study_contr$comparator_arm,
                                    study_contr$experimental_arm), comp_length)

  # Indicate the studies providing direct evidence to the corresponding comparison (yes/no)
  colour_study_basic <- ifelse(colour_study_basic0 == sub("d", "", treat_effect$variable), "Yes", "No")

  # Bring all together
  study_id <- direct_evid <- covar <- NULL
  dataset_treat <- data.frame(study_id = rep(study_contr$study_name, comp_length),
                              treat_effect,
                              comp_indic,
                              direct_evid = factor(colour_study_basic, levels = c("Yes", "No")),
                              covar = rep(study_contr$covariate, comp_length))

  # Get the comparisons as second versus first arm
  obs_comp <- paste0(substr(sub("d", "", treat_effect$variable), 2, 2), "vs",
                     substr(sub("d", "", treat_effect$variable), 1, 1))

  # Attach names to the pairwise comparisons
  named_compar <- possible_observed_comparisons(drug_names, obs_comp)$poss_comp

  # Repeat as many times as the unique elements in 'treat_effect$variable'
  dataset_treat$variable <-
    rep(named_compar$comp_name, table(treat_effect$variable))

  # *All* regression coefficients
  dataset_reg <- if (contr_res$covar_assumption == "common") {
    data.frame(study_id = study_contr$study_name,
               reg_coeff,
               covar = study_contr$covariate)
  } else if (is.element(contr_res$covar_assumption,
                        c("independent", "exchangenable"))) {
    data.frame(study_id = rep(study_contr$study_name, comp_length),
               reg_coeff,
               comp_indic,
               direct_evid = factor(colour_study_basic, levels = c("Yes", "No")),
               covar = rep(study_contr$covariate, comp_length))
  }


  ## Panel of contribution plots for treatment effects
  plot_treat <-
    if (comparisons == "basic") {
      ggplot(subset(dataset_treat, comp_indic == "basic"),
             aes(x = covar,
                 y = value)) +
        geom_text_repel(aes(label = study_id),
                        size = label_size) +
        geom_hline(yintercept = 50,
                   colour = "grey") +
        geom_point(aes(col = direct_evid)) +
        facet_wrap(~variable) +
        scale_colour_manual(values = c("Yes" = "blue", "No" = "red")) +
        scale_y_continuous(limits = c(0, upper_limit), expand = c(0.03, 0)) +
        scale_x_continuous(limits = c(min(dataset_treat$covar), max(dataset_treat$covar)),
                           expand = c(0.02, 0),
                           breaks = seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by),
                           labels = sprintf("%.2f", seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by))) +
        labs(x = name_x_axis,
             y = "Study contributions (%)",
             colour = "Provides direct evidence",
             subtitle = "Treatment effects (basic parameters)") +
        guides(colour = guide_legend(override.aes = list(shape = 19, size = 3.5))) +
        theme_bw() +
        theme(axis.title = element_text(size = axis_title_size, face = "bold"),
              axis.text = element_text(size = axis_text_size),
              strip.text = element_text(size = strip_text_size, face = "bold"),
              plot.subtitle = element_text(size = subtitle_size, face = "bold"),
              legend.position = "bottom",
              legend.title = element_text(size = axis_title_size, face = "bold"),
              legend.text = element_text(size = axis_text_size))
    } else if (comparisons == "functional") {
      ggplot(subset(dataset_treat, comp_indic == "functional"),
             aes(x = covar,
                 y = value)) +
        geom_text_repel(aes(label = study_id),
                        size = label_size) +
        geom_hline(yintercept = 50,
                   colour = "grey") +
        geom_point(colour = "red") +
        facet_wrap(~variable) +
        scale_y_continuous(limits = c(0, upper_limit), expand = c(0.03, 0)) +
        scale_x_continuous(limits = c(min(dataset_treat$covar), max(dataset_treat$covar)),
                           expand = c(0.02, 0),
                           breaks = seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by),
                           labels = sprintf("%.2f", seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by))) +
        labs(x = name_x_axis,
             y = "Study contributions (%)",
             subtitle = "Treatment effects (functional parameters)") +
        theme_bw() +
        theme(axis.title = element_text(size = axis_title_size, face = "bold"),
              axis.text = element_text(size = axis_text_size),
              panel.background = element_rect(fill = "white"),
              strip.text = element_text(size = strip_text_size, face = "bold"),
              plot.subtitle = element_text(size = subtitle_size, face = "bold"))
    } else if (comparisons == "all") {
      ggplot(dataset_treat,
             aes(x = covar,
                 y = value)) +
        geom_text_repel(aes(label = study_id),
                        size = label_size) +
        geom_hline(yintercept = 50,
                   colour = "grey") +
        geom_point(aes(col = direct_evid)) +
        facet_wrap(~variable) +
        scale_colour_manual(values = c("Yes" = "blue", "No" = "red")) +
        scale_y_continuous(limits = c(0, upper_limit), expand = c(0.03, 0)) +
        scale_x_continuous(limits = c(min(dataset_treat$covar), max(dataset_treat$covar)),
                           expand = c(0.02, 0),
                           breaks = seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by),
                           labels = sprintf("%.2f", seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by))) +
        labs(x = name_x_axis,
             y = "Study contributions (%)",
             colour = "Provides direct evidence",
             subtitle = "Treatment effects") +
        guides(colour = guide_legend(override.aes = list(shape = 19, size = 3.5))) +
        theme_bw() +
        theme(axis.title = element_text(size = axis_title_size, face = "bold"),
              axis.text = element_text(size = axis_text_size),
              panel.background = element_rect(fill = "white"),
              strip.text = element_text(size = strip_text_size, face = "bold"),
              plot.subtitle = element_text(size = subtitle_size, face = "bold"),
              legend.position = "bottom",
              legend.title = element_text(size = axis_title_size, face = "bold"),
              legend.text = element_text(size = axis_text_size))
    }


  ## Panel of contribution plots for regression coefficients
  plot_reg <-
    if(contr_res$covar_assumption != "common" & comparisons == "basic") {
      ggplot(subset(dataset_reg, comp_indic == "basic"),
             aes(x = covar,
                 y = value)) +
        geom_text_repel(aes(label = study_id),
                        size = label_size) +
        geom_hline(yintercept = 50,
                   colour = "grey") +
        geom_point(aes(col = direct_evid)) +
        facet_wrap(~variable) +
        scale_colour_manual(values = c("Yes" = "blue", "No" = "red")) +
        scale_y_continuous(limits = c(0, upper_limit), expand = c(0.03, 0)) +
        scale_x_continuous(limits = c(min(dataset_treat$covar), max(dataset_treat$covar)),
                           expand = c(0.02, 0),
                           breaks = seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by),
                           labels = sprintf("%.2f", seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by))) +
        labs(x = name_x_axis,
             y = "Study contributions (%)",
             colour = "Provides direct evidence",
             subtitle = "Regression coefficients (basic parameters)") +
        guides(colour = guide_legend(override.aes = list(shape = 19, size = 3.5))) +
        theme_bw() +
        theme(axis.title = element_text(size = axis_title_size, face = "bold"),
              axis.text = element_text(size = axis_text_size),
              panel.background = element_rect(fill = "white"),
              strip.text = element_text(size = strip_text_size, face = "bold"),
              plot.subtitle = element_text(size = subtitle_size, face = "bold"),
              legend.position = "bottom",
              legend.title = element_text(size = axis_title_size, face = "bold"),
              legend.text = element_text(size = axis_text_size))
    } else if (contr_res$covar_assumption != "common" & comparisons == "functional") {
      ggplot(subset(dataset_reg, comp_indic == "functional"),
             aes(x = covar,
                 y = value)) +
        geom_text_repel(aes(label = study_id),
                        size = label_size) +
        geom_hline(yintercept = 50,
                   colour = "grey") +
        geom_point(colour = "red") +
        facet_wrap(~variable) +
        scale_y_continuous(limits = c(0, upper_limit), expand = c(0.03, 0)) +
        scale_x_continuous(limits = c(min(dataset_treat$covar), max(dataset_treat$covar)),
                           expand = c(0.02, 0),
                           breaks = seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by),
                           labels = sprintf("%.2f", seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by))) +
        labs(x = name_x_axis,
             y = "Study contributions (%)",
             subtitle = "Regression coefficients (functional parameters)") +
        guides(colour = guide_legend(override.aes = list(shape = 19, size = 3.5))) +
        theme_bw() +
        theme(axis.title = element_text(size = axis_title_size, face = "bold"),
              axis.text = element_text(size = axis_text_size),
              panel.background = element_rect(fill = "white"),
              strip.text = element_text(size = strip_text_size, face = "bold"),
              plot.subtitle = element_text(size = subtitle_size, face = "bold"))
    } else if (contr_res$covar_assumption != "common"  & comparisons == "all") {
      ggplot(dataset_reg,
             aes(x = covar,
                 y = value)) +
        geom_text_repel(aes(label = study_id),
                        size = label_size) +
        geom_hline(yintercept = 50,
                   colour = "grey") +
        geom_point(aes(col = direct_evid)) +
        facet_wrap(~variable) +
        scale_colour_manual(values = c("Yes" = "blue", "No" = "red")) +
        scale_y_continuous(limits = c(0, upper_limit), expand = c(0.03, 0)) +
        scale_x_continuous(limits = c(min(dataset_treat$covar), max(dataset_treat$covar)),
                           expand = c(0.02, 0),
                           breaks = seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by),
                           labels = sprintf("%.2f", seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by))) +
        labs(x = name_x_axis,
             y = "Study contributions (%)",
             colour = "Provides direct evidence",
             subtitle = "Regression coefficients") +
        guides(colour = guide_legend(override.aes = list(shape = 19, size = 3.5))) +
        theme_bw() +
        theme(axis.title = element_text(size = axis_title_size, face = "bold"),
              axis.text = element_text(size = axis_text_size),
              panel.background = element_rect(fill = "white"),
              strip.text = element_text(size = strip_text_size, face = "bold"),
              plot.subtitle = element_text(size = subtitle_size, face = "bold"),
              legend.position = "bottom",
              legend.title = element_text(size = axis_title_size, face = "bold"),
              legend.text = element_text(size = axis_text_size))
    } else if (contr_res$covar_assumption == "common") {
      ggplot(dataset_reg,
             aes(x = covar,
                 y = value)) +
        geom_text_repel(aes(label = study_id),
                        size = label_size) +
        geom_hline(yintercept = 50,
                   colour = "grey") +
        geom_point(colour = "red") +
        scale_y_continuous(limits = c(0, upper_limit), expand = c(0.03, 0)) +
        scale_x_continuous(limits = c(min(dataset_treat$covar), max(dataset_treat$covar)),
                           expand = c(0.02, 0),
                           breaks = seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by),
                           labels = sprintf("%.2f", seq(min(dataset_treat$covar), max(dataset_treat$covar), seq_by))) +
        labs(x = name_x_axis,
             y = "Study contributions (%)",
             subtitle = "Regression coefficient (common interaction)") +
        theme_bw() +
        theme(axis.title = element_text(size = axis_title_size, face = "bold"),
              axis.text = element_text(size = axis_text_size),
              panel.background = element_rect(fill = "white"),
              strip.text = element_text(size = strip_text_size, face = "bold"),
              plot.subtitle = element_text(size = subtitle_size, face = "bold"))
    } else if (contr_res$covar_assumption == "no") {
      NULL
    }


  ## Bring together
  results <- if (contr_res$covar_assumption == "no") plot_treat else
    list(plot_treat = plot_treat,
         plot_reg = plot_reg)

  return(suppressWarnings(print(results)))
}
