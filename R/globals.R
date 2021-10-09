#' @import
#'   ggplot2
#'   R2jags
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
                         "EM.ref00",
                         "evidence",
                         "frail",
                         "interval",
                         "Intervention",
                         "label",
                         "leverage",
                         "lower",
                         "m.prop",
                         "Order",
                         "poor",
                         "risk",
                         "sd.value",
                         "signed.dev",
                         "signed.dev.o",
                         "stat.sign",
                         "study",
                         "studlab",
                         "treatment",
                         "treatment1",
                         "treat1",
                         "treat2",
                         "upper",
                         "value",
                         "value2",
                         "value.cum",
                         "value.rank",
                         "value.SUCRA",
                         "Var1",
                         "Var2",
                         "x",
                         "y"))

Sys.setenv(`_R_S3_METHOD_REGISTRATION_NOTE_OVERWRITES_` = "false")
suppressMessages({"Registered S3 method overwritten by 'mcmcplots':
                     method        from
                     as.mcmc.rjags R2jags"})


