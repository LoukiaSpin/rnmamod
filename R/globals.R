#' @import
#'   ggplot2
#'   stats
#'
#' @importFrom coda as.mcmc
#' @importFrom dplyr select starts_with ends_with %>%
#' @importFrom fdrtool rhalfnorm
#' @importFrom MASS fractions
#' @importFrom gemtc mtc.network mtc.data.studyrow mtc.nodesplit.comparisons
#' @importFrom ggfittext geom_fit_text
#' @importFrom ggpubr ggarrange
#' @importFrom ggrepel geom_text_repel
#' @importFrom knitr kable
#' @importFrom mcmcplots mcmcplot
#' @importFrom netmeta pairwise netconnection
#' @importFrom pcnetmeta nma.networkplot
#' @importFrom reshape2 melt
#' @importFrom R2jags jags
#' @importFrom scales rescale percent
#' @importFrom writexl write_xlsx
#' @importFrom utils combn


utils::globalVariables(c("active",
                         "analysis",
                         "char1",
                         "char2",
                         "char3",
                         "comp",
                         "comparison",
                         "DIC",
                         "EM",
                         "em_ref00",
                         "evidence",
                         "frail",
                         "interval",
                         "intervention",
                         "label",
                         "leverage",
                         "lower",
                         "m_prop",
                         "Order",
                         "poor",
                         "risk",
                         "sd_value",
                         "signed_dev",
                         "signed_dev_o",
                         "stat_sign",
                         "study",
                         "studlab",
                         "treatment",
                         "treatment1",
                         "treat1",
                         "treat2",
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

Sys.setenv(`_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_` = "false")
Sys.setenv(JAGS_HOME = "C:/Program Files/JAGS/JAGS-4.2.0")
Sys.setenv(JAGS_ROOT = "C:/Program Files/JAGS/JAGS-4.2.0")
