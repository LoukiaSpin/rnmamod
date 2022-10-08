#'  Rating the quality of transitivity evaluation
#'
#' @description Classifies a systematic review with multiple interventions as
#'   having low, unclear or high quality regarding the transitivity evaluation.
#'
#' @param plan_protocol Character string that indicates how the systematic
#'   review planned in the protocol to assess the transitivity assumption. The
#'   following values can be considered: \code{"Both"}, \code{"No"},
#'   \code{"No protocol"}, \code{"Only direct methods"}, and
#'   \code{"Only indirect methods"}.
#' @param plan_methods Character string that indicates whether the systematic
#'   review described in the methods section a strategy to assess the
#'   transitivity assumption. The following values can be considered:
#'   \code{"Yes"}, and \code{"No"}.
#' @param report_results Character string that indicates whether the systematic
#'   review reported in the results section the transitivity evaluation and
#'   which strategy was employed. The following values can be considered:
#'   \code{"Both"}, \code{"No"}, \code{"Only direct methods"}, and
#'   \code{"Only indirect methods"}.
#' @param discuss_trans Character string that indicates whether the systematic
#'   review discussed the transitivity assumption and which model parameters
#'   where considered. The following values can be considered: \code{"Both"},
#'   \code{"No"}, \code{"Only treatment effects"}, \code{"Other parameter"},
#'   and \code{"NMA not conducted"}.
#' @param proper_table Character string that indicates whether the systematic
#'   review reported a proper table of characteristics. The following values can
#'   be considered: \code{"No"}, \code{"No table"}, and \code{"Yes"}.
#'
#' @return A character with value \code{"Low"}, \code{"Unclear"}, or
#'   \code{"High"} to indicate low, unclear, or high-quality of transitivity
#'   evaluation.
#'
#' @details A systematic review with \code{"Low"} quality of transitivity
#'   evaluation does not provide a protocol, nor describe the evaluation
#'   strategy in the methods section, does not report the evaluation results,
#'   nor discusses the transitivity evaluation and does not provide a proper
#'   table of characteristics. On the contrary, a systematic review with
#'   \code{"High"} quality of transitivity evaluation provides an evaluation
#'   plan in the protocol (including at least one direct method), describes the
#'   evaluation strategy in the methods section (including at least one direct
#'   method), reports the evaluation results in the results section, discusses
#'   the transitivity evaluation while considering at least one model parameter,
#'   if NMA has been conducted (i.e. \code{nma_abstain = FALSE}), and provides
#'   a proper table of characteristics. Otherwise, the systematic review is
#'   judged to have an \code{"Unclear"} quality of transitivity evaluation.
#'
#' @author {Loukia M. Spineli}
#'
#' @references
#' Spineli LM, Kalyvas C, Seide SE, Papadimitropoulou K. Low awareness of the
#' transitivity assumption in complex networks of interventions: empirical
#' evidence from 356 reviews. 2022 \emph{submitted}
#'
#' @export
trans_quality <- function(plan_protocol,
                          plan_methods,
                          report_results,
                          discuss_trans,
                          proper_table) {

  # Missing and default arguments
  prot_val <- c("Both", "No", "No protocol",
                "Only direct methods", "Only indirect methods")
  a <- "Insert 'Both', 'No', 'No protocol', 'Only direct methods',"
  b <- "or 'Only indirect methods'"
  prot_val_text1 <- paste(a, b)
  prot_val_text2 <- "for plan_protocol."
  plan_protocol <- if (missing(plan_protocol) ||
                       !is.element(plan_protocol, prot_val)) {
    stop(paste(prot_val_text1, prot_val_text2), call. = FALSE)
  } else {
    plan_protocol
  }
  plan_methods <- if (missing(plan_methods) ||
                      !is.element(plan_methods, c("No", "Yes"))) {
    stop( "Insert 'No', or 'Yes' for plan_methods.", call. = FALSE)
  } else {
    plan_methods
  }
  res_val <- c("Both", "No", "Only direct methods", "Only indirect methods")
  res_val_text1 <-
    "Insert 'Both', 'No', 'Only direct', or 'Only indirect'"
  res_val_text2 <- "for report_results."
  report_results <- if (missing(report_results) ||
                        !is.element(report_results, res_val)) {
    stop(paste(res_val_text1, res_val_text2), call. = FALSE)
  } else {
    report_results
  }
  disc_val <- c("Both", "No", "Only treatment effects", "Other parameter",
                "NMA not conducted")
  disc_val_text1 <-
    "Insert 'Both', 'No', 'Only treatment effects', 'Other parameter'"
  disc_val_text2 <- "or 'NMA not conducted' for discuss_trans."
  discuss_trans <- if (missing(discuss_trans) ||
                       !is.element(discuss_trans, disc_val)) {
    stop(paste(disc_val_text1, disc_val_text2), call. = FALSE)
  } else {
    discuss_trans
  }
  proper_table <- if (missing(proper_table) ||
                      !is.element(proper_table, c("No", "No table", "Yes"))) {
    stop( "Insert 'No', 'No table', or 'Yes' for proper_table.", call. = FALSE)
  } else {
    proper_table
  }

  # The function
  quality <-
    if (is.element(plan_protocol, c("Both", "Only direct methods")) &
        plan_methods == "Yes" &
        is.element(report_results, c("Both", "Only direct methods")) &
        discuss_trans != "No" &
        proper_table == "Yes") {
      "High"
    } else if (plan_protocol == "No protocol" &
               plan_methods == "No" &
               report_results == "No" &
               discuss_trans == "No" &
               proper_table == "No table") {
      "Low"
    } else {
      "Unclear"
    }

  return(quality)
}
