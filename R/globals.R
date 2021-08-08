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
#' @importFrom netmeta pairwise
#' @importFrom pcnetmeta nma.networkplot
#' @importFrom reshape2 melt
#' @importFrom R2jags jags
#' @importFrom scales rescale percent
#' @importFrom writexl write_xlsx
#' @importFrom utils combn


utils::globalVariables(c("active",
                         "analysis",
                         "DIC",
                         "comp",
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
                         "treatment",
                         "treatment1",
                         "upper",
                         "value",
                         "value2",
                         "value.SUCRA",
                         "value.rank",
                         "value.cum",
                         "Var1",
                         "Var2",
                         "x",
                         "y"))
