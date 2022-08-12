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
#'   \code{"No"}, \code{"Only treatment effects"}, and \code{"Other parameter"}.
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
#'   method), reports the evaluation results in the results section,  discusses
#'   the transitivity evaluation while considering at least one model parameter,
#'   and provides a proper table of characteristics. Otherwise, the systematic
#'   review is judged to have an \code{"Unclear"} quality of transitivity
#'   evaluation.
#'
#' @author {Loukia M. Spineli}
#'
#' @export
trans_quality <- function(plan_protocol,
                          plan_methods,
                          report_results,
                          discuss_trans,
                          proper_table) {

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
