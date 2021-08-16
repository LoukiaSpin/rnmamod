## A function to describe the network

describe.network <- function(data, drug.names, measure) {

  options(warn = -1)

  dat <- data.preparation(data, measure)

  two.arm.ns <- length(which(dat$na < 3))

  multi.arm.ns <- length(which(dat$na > 2))

  total.rand.partic.interv <- aggregate(unlist(dat$N), by = list(unlist(dat$t)), sum)[, 2]

  total.obs.partic.interv <- aggregate(unlist(dat$N - dat$m), by = list(unlist(dat$t)), sum)[, 2]

  pair.mod <- pairwise(as.list(dat$t), event = as.list(dat$m), n = as.list(dat$N), data = cbind(dat$t, dat$m, dat$N), studlab = 1:dat$ns)[, c(3:6, 8, 7, 9)]
  colnames(pair.mod) <- c("study", "t1", "t2", "m1", "m2", "n1", "n2")

  pair.mod[, 2] <- drug.names[pair.mod$t1]

  pair.mod[, 3] <- drug.names[pair.mod$t2]

  comp <- paste(pair.mod[, "t2"], "vs", pair.mod[, "t1"])

  direct.comp <- length(unique(comp))

  total.rand.network <- sum(apply(pair.mod[, c("n1", "n2")], 1, sum))

  total.rand.partic.comp <- aggregate(apply(pair.mod[, c("n1", "n2")], 1, sum), by = list(comp), sum)[, 2]

  total.obs.network <- sum(apply(pair.mod[, c("n1", "n2")] - pair.mod[, c("m1", "m2")], 1, sum))

  total.obs.partic.comp <- aggregate(apply(pair.mod[, c("n1", "n2")] - pair.mod[, c("m1", "m2")], 1, sum), by = list(comp), sum)[, 2]

  if (measure == "OR") {
    arm.risk <- dat$r/(dat$N - dat$m)
    rule <- apply(ifelse(arm.risk == 0.0, 1, 0), 1, sum, na.rm = T)
    trial.zero.event <- ifelse(length(which(rule > 0)) == 0, 0, which(rule > 0))
    trial.all.zero.event <- ifelse(length(which(is.element(rule, dat$na) == T)) == 0, 0, which(is.element(rule, dat$na) == T))

    total.event.network <- sum(unlist(dat$r), na.rm = T)

    total.event.interv <- aggregate(unlist(dat$r), by = list(unlist(dat$t)), sum)[, 2]

    total.risk.interv <- round(total.event.interv/total.obs.partic.interv, 2)

    min.risk.interv <- round(aggregate(unlist(dat$r)/unlist(dat$N - dat$m), by = list(unlist(dat$t)), min)[, 2], 2)

    median.risk.interv <- round(aggregate(unlist(dat$r)/unlist(dat$N - dat$m), by = list(unlist(dat$t)), median)[, 2], 2)

    max.risk.interv <- round(aggregate(unlist(dat$r)/unlist(dat$N - dat$m), by = list(unlist(dat$t)), max)[, 2], 2)

    table.interv.bin <- data.frame(drug.names, as.data.frame(table(unlist(dat$t)))[, 2], total.rand.partic.interv, total.obs.partic.interv, total.event.interv,
                                   total.risk.interv, min.risk.interv, median.risk.interv, max.risk.interv)
    colnames(table.interv.bin) <- c("Interventions", "Total trials", "Total randomised", "Total completers", "Total events", "Total events (%)",
                                    "Minimum % events", "Median % events", "Maximum % events")


    pair.bin <- pairwise(as.list(dat$t), event = as.list(dat$r), n = as.list(dat$N), data = cbind(dat$t, dat$r, dat$N), studlab = 1:dat$ns)[, c(3:6, 8, 7, 9)]
    colnames(pair.bin) <- c("study", "t1", "t2", "r1", "r2", "n1", "n2")

    total.event.comp <- aggregate(apply(pair.bin[, c("r1", "r2")], 1, sum), by = list(comp), sum)[, 2]

    total.risk.comp <- round(total.event.comp/total.obs.partic.comp, 2)

    min.risk.comp <- round(aggregate(apply(pair.bin[, c("r1", "r2")], 1, sum)/(apply(pair.bin[, c("n1", "n2")], 1, sum) - apply(pair.mod[, c("m1", "m2")], 1, sum)), by = list(comp), min)[, 2], 2)

    median.risk.comp <- round(aggregate(apply(pair.bin[, c("r1", "r2")], 1, sum)/(apply(pair.bin[, c("n1", "n2")], 1, sum) - apply(pair.mod[, c("m1", "m2")], 1, sum)), by = list(comp), median)[, 2], 2)

    max.risk.comp <- round(aggregate(apply(pair.bin[, c("r1", "r2")], 1, sum)/(apply(pair.bin[, c("n1", "n2")], 1, sum) - apply(pair.mod[, c("m1", "m2")], 1, sum)), by = list(comp), max)[, 2], 2)

    table.comp.bin <- data.frame(as.data.frame(table(comp))[, 1], as.data.frame(table(comp))[, 2], total.rand.partic.comp, total.obs.partic.comp, total.event.comp,
                                 total.risk.comp, min.risk.comp, median.risk.comp, max.risk.comp)
    colnames(table.comp.bin) <- c("Comparisons", "Total trials", "Total randomised", "Total completers", "Total events", "Total events (%)",
                                  "Minimum % events", "Median % events", "Maximum % events")

    ## Write the tables as .xlsx
    writexl::write_xlsx(table.interv.bin, paste0("Table.interventions.xlsx"))
    writexl::write_xlsx(table.comp.bin, paste0("Table.comparisons.xlsx"))
  } else {

    min.t.interv <- round(aggregate(unlist(dat$y0)/unlist(dat$se0), by = list(unlist(dat$t)), min)[, 2], 2)

    median.t.interv <- round(aggregate(unlist(dat$y0)/unlist(dat$se0), by = list(unlist(dat$t)), median)[, 2], 2)

    max.t.interv <- round(aggregate(unlist(dat$y0)/unlist(dat$se0), by = list(unlist(dat$t)), max)[, 2], 2)

    table.interv.con <- data.frame(drug.names, as.data.frame(table(unlist(dat$t)))[, 2], total.rand.partic.interv, total.obs.partic.interv,
                                   min.t.interv, median.t.interv, max.t.interv)
    colnames(table.interv.con) <- c("Interventions", "Total trials", "Total randomised", "Total completers", "Minimum y/se", "Median y/se", "Maximum y/se")


    pair.con <- pairwise(as.list(dat$t),  n = as.list(dat$N), mean = as.list(dat$y0), sd = as.list(dat$se0), data = cbind(dat$t, dat$N, dat$y0, dat$se0), studlab = 1:dat$ns)[, c(3:5, 7, 10, 8, 11, 6, 9)]
    colnames(pair.con) <- c("study", "t1", "t2", "y1", "y2", "se1", "se2","n1", "n2")

    min.t.comp <- round(aggregate((pair.con$y2 - pair.con$y1)/(pair.con$se2 + pair.con$se1), by = list(comp), min)[, 2], 2)

    median.t.comp <- round(aggregate((pair.con$y2 - pair.con$y1)/(pair.con$se2 + pair.con$se1), by = list(comp), median)[, 2], 2)

    max.t.comp <- round(aggregate((pair.con$y2 - pair.con$y1)/(pair.con$se2 + pair.con$se1), by = list(comp), max)[, 2], 2)

    table.comp.con <- data.frame(as.data.frame(table(comp))[, 1], as.data.frame(table(comp))[, 2], total.rand.partic.comp, total.obs.partic.comp,
                                 min.t.comp, median.t.comp, max.t.comp)
    colnames(table.comp.con) <- c("Comparisons", "Total trials", "Total randomised", "Total completers", "Minimum MD/se", "Median MD/se", "Maximum MD/se")

    ## Write the tables as .xlsx
    writexl::write_xlsx(table.interv.con, paste0("Table.interventions.xlsx"))
    writexl::write_xlsx(table.comp.con, paste0("Table.comparisons.xlsx"))

  }


  results <- list(direct.comp = direct.comp,
                  two.arm.ns = two.arm.ns,
                  multi.arm.ns = multi.arm.ns,
                  total.rand.network = total.rand.network,
                  total.obs.network = total.obs.network)

  results <- if (measure == "OR") {
    append(results, list(total.event.network = total.event.network,
                         trial.zero.event = trial.zero.event,
                         trial.all.zero.event = trial.all.zero.event,
                         Table.interventions = knitr::kable(table.interv.bin),
                         Table.comparisons = knitr::kable(table.comp.bin)))
  } else {
    append(results, list(Table.interventions = knitr::kable(table.interv.con),
                         Table.comparisons = knitr::kable(table.comp.con)))
  }

  return(results)

}
