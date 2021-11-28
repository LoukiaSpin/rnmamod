#' Heatmap of robustness
#'
#' @description Facilitates the detection of comparisons that are associated
#'   with lack of robustness in the context of a sensitivity analysis.
#'
#' @param robust An object of S3 class \code{\link{robustness_index}}.
#'   See 'Value' in \code{\link{robustness_index}}.
#' @param drug_names A vector of labels with the name of the interventions in
#'   the order they appear in the argument \code{data} of
#'   \code{\link{run_model}}. If \code{drug_names} is not defined,
#'   the order of the interventions as they appear in \code{data} is used,
#'   instead.
#'
#' @return \code{heatmap_robustness} first prints on the R console a message on
#'   the threshold of robustness determined by the user in
#'   \code{\link{robustness_index}}. Then, it returns a lower triangular
#'   heatmap matrix with the robustness index value of all possible pairwise
#'   comparisons.
#'
#' @details The heatmap illustrates the robustness index for each possible
#'   pairwise comparison in the network. The pairwise comparisons are read from
#'   left to right. Comparisons highlighted with green or red colour imply
#'   robust or frail conclusions for the primary analysis, respectively.
#'   This corresponds to robustness index below or at least the selected
#'   threshold of robustness (see 'Details' in \code{\link{robustness_index}}).
#'   The robustness index value of each pairwise comparison also appears in the
#'   corresponding cell.
#'   When there is at least one comparison with frail conclusions, the primary
#'   analysis results may be questionable for the whole network
#'   (Spineli et al., 2021).
#'
#'   \code{heatmap_robustness} is \emph{not} restricted to the sensitivity
#'   analysis concerning the impact of missing participant outcome data.
#'
#'   \code{heatmap_robustness} uses the threshold of robustness selected in the
#'   \code{\link{robustness_index}} function.
#'
#'   \code{heatmap_robustness} can be used only for a network of interventions.
#'   Otherwise, the execution of the function will be stopped and an
#'   error message will be printed on the R console.
#'
#' @author {Loukia M. Spineli}
#'
#' @seealso \code{\link{robustness_index}}, \code{\link{run_model}}
#'
#' @references
#' Spineli LM, Kalyvas C, Papadimitropoulou K. Quantifying the robustness of
#' primary analysis results: A case study on missing outcome data in pairwise
#' and network meta-analysis.
#' \emph{Res Synth Methods} 2021;\bold{12}(4):475--490.
#' \doi{10.1002/jrsm.1478}
#'
#' @examples
#' data("nma.baker2009")
#'
#' # Read results from 'run_sensitivity' (using the default arguments)
#' res_sens <- readRDS(system.file('extdata/res_sens_baker.rds',
#'                     package = 'rnmamod'))
#'
#' # Calculate the robustness index
#' robust <- robustness_index(sens = res_sens,
#'                            threshold = 0.28)
#'
#' # The names of the interventions in the order they appear in the dataset
#' interv_names <- c("placebo", "budesonide", "budesonide plus formoterol",
#'                   "fluticasone", "fluticasone plus salmeterol",
#'                   "formoterol", "salmeterol", "tiotropium")
#'
#' # Create the heatmap of robustness
#' heatmap_robustness(robust = robust,
#'                    drug_names = interv_names)
#'
#' @export
heatmap_robustness <- function(robust, drug_names) {

  if (any(is.na(robust))) {
    aa <- "Missing participant outcome data have *not* been collected."
    stop(paste(aa, "This function cannot be used."), call. = FALSE)
  }

  drug_names <- if (missing(drug_names)) {
    aa <- "The  argument 'drug_names' has not been defined."
    bb <- "The intervention ID, as specified in 'data' is used as"
    cc <- "intervention names"
    message(cat(paste0("\033[0;", col = 32, "m", aa, bb, cc, "\033[0m", "\n")))
    nt <- (1 + sqrt(1 + 8 * length(robust$robust))) / 2
    as.character(1:nt)
  } else {
    drug_names
  }
  len_drugs <- length(drug_names)

  if (len_drugs < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis",
         call. = FALSE)
  }

  robust_index <- robust$robust_index
  threshold <- robust$threshold

  if (missing(threshold) & is.element(robust$measure, "OR")) {
    threshold <- 0.28
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt =
                       "The value 0.28 was assigned on 'threshold' by default",
                       "\033[0m", "\n")))
  } else if (missing(threshold) & is.element(robust$measure,
                                             c("MD", "SMD", "ROM"))) {
    threshold <- 0.17
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt =
                       "The value 0.17 was assigned on 'threshold' by default",
                       "\033[0m", "\n")))
  } else {
    threshold <- threshold
    message(cat(paste0("\033[0;",
                       col = 32,
                       "m",
                       txt = paste("The value", threshold,
                                   "was assigned on 'threshold' for",
                                   effect_measure_name(robust$measure)),
                       "\033[0m", "\n")))
  }

  # Lower triangular matrix: comparisons are read from the left to the right
  mat <- matrix(NA,
                nrow = len_drugs - 1,
                ncol = len_drugs - 1)
  mat[lower.tri(mat, diag = TRUE)] <- sprintf("%.2f", robust_index)
  colnames(mat) <- drug_names[1:(len_drugs - 1)]
  rownames(mat) <- drug_names[2:len_drugs]
  mat_new <- melt(mat, na.rm = TRUE)


  ## Create the heatmap for one network of interventions
  p <- ggplot(mat_new,
              aes(factor(Var2, levels = drug_names[1:(len_drugs - 1)]),
                  factor(Var1, levels = drug_names[len_drugs:2]),
                  fill = ifelse(value < threshold, "high", "poor"))) +
         geom_tile(colour = "white") +
         geom_text(aes(factor(Var2,
                              levels = drug_names[1:(len_drugs - 1)]),
                       factor(Var1, levels = drug_names[len_drugs:2]),
                       label = value,
                       fontface = "bold"),
                   size = rel(4.5)) +
         scale_fill_manual(breaks = c("high", "poor"),
                           values = c("#009E73", "#D55E00")) +
         scale_x_discrete(position = "top") +
         labs(x = "", y = "") +
         theme_bw() +
         theme(legend.position = "none",
               axis.text.x = element_text(size = 12, angle = 50, hjust = 0.0),
               axis.text.y = element_text(size = 12))

  return(p)
}
