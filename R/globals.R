#' @import
#'   ggplot2
#'   stats
#'
#' @importFrom dplyr select starts_with ends_with %>%
#' @importFrom fdrtool rhalfnorm
#' @importFrom MASS fractions
#' @importFrom gemtc mtc.network mtc.data.studyrow mtc.nodesplit.comparisons
#' @importFrom ggfittext geom_fit_text
#' @importFrom ggpubr ggarrange
#' @importFrom ggrepel geom_text_repel
#' @importFrom grDevices rainbow
#' @importFrom heatmaply heatmaply
#' @importFrom igraph E layout_in_circle make_graph plot.igraph V
#' @importFrom knitr kable
#' @importFrom netmeta netconnection
#' @importFrom reshape2 melt
#' @importFrom R2jags jags autojags
#' @importFrom rdist triangle_inequality
#' @importFrom scales hue_pal rescale percent
#' @importFrom writexl write_xlsx
#' @importFrom utils combn packageDescription

.onAttach <- function (libname, pkgname) {
  welcome <- paste0("Loading 'rnmamod', version",
                    " ", packageDescription("rnmamod")$Version,
                    ". Welcome to Network Meta-analysis with",
                    " " , "Missing Participants! ^_^")
  packageStartupMessage(welcome)
}

utils::globalVariables(c("active",
                         "analysis",
                         "capture.output",
                         "char",
                         "char1",
                         "char2",
                         "char3",
                         "characteristic",
                         "cluster",
                         "clusters",
                         "comp",
                         "compar",
                         "comparison",
                         "Comparison",
                         "DIC",
                         "diss",
                         "EM",
                         "em_ref00",
                         "evidence",
                         "frail",
                         "index",
                         "input0",
                         "interval",
                         "intervention",
                         "label",
                         "leverage",
                         "lower",
                         "m_prop",
                         "num_trials",
                         "Order",
                         "outcome",
                         "perc",
                         "perc2",
                         "point",
                         "poor",
                         "results",
                         "risk",
                         "sd_value",
                         "signed_dev",
                         "signed_dev_o",
                         "silhouette",
                         "silhouette_new",
                         "size",
                         "stat_sign",
                         "str_wrap",
                         "study",
                         "studlab",
                         "total",
                         "total_dissimilarity",
                         "treatment",
                         "treatment1",
                         "treat1",
                         "treat2",
                         "trial",
                         "type",
                         "upper",
                         "value",
                         "value2",
                         "value_cum",
                         "value_rank",
                         "value_sucra",
                         "Var1",
                         "Var2",
                         "x",
                         "y"))

#Sys.setenv(`_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_` = "false")
Sys.setenv(JAGS_HOME = "C:/Program Files/JAGS/JAGS-4.2.0")
Sys.setenv(JAGS_ROOT = "C:/Program Files/JAGS/JAGS-4.2.0")

