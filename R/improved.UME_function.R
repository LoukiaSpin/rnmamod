## Detect the non-baseline comparisons in multi-arm trials that are not informed by two-arm trials, and hence,
## the posterior distribution coincides with the prior distribution whichi means implausible posterior SD.


improved.UME <- function(t, m, N, ns, na){

  ## Turn into contrast-level data: one row per possible comparison in each trial ('netmeta')
  (wide.format <- pairwise(as.list(t), event = as.list(m), n = as.list(N), data = cbind(t, m, N), studlab = 1:ns)[, c(3:6, 8, 7, 9)])
  colnames(wide.format) <- c("study", "t1", "t2", "m1", "m2", "n1", "n2")
  wide.format$comp <- paste0(wide.format$t1, "vs", wide.format$t2)


  ## Repeat 'na' as many times as the length of the corresponding trial
  arms0 <- list()
  for(i in 1:ns) {
    arms0[[i]] <- rep(na[i], dim(combn(na[i], 2))[2] )
  }


  ## Indicate whether a trial is two-arm or multi-arm
  (wide.format$arms <- unlist(arms0))
  wide.format[wide.format$arms > 2, 9] <- "multi-arm"
  wide.format[wide.format$arms < 3, 9] <- "two-arm"


  ## The frequency of each observed comparisons in two-arm and multi-arm trials
  (tab.comp.arms0 <- xtabs(~ comp + arms, data = wide.format))

  # Turn into data-frame into order to select the comparisons
  (tab.comp.arms <- data.frame(names(tab.comp.arms0[, 1]), tab.comp.arms0[, 1], tab.comp.arms0[, 2]))
  colnames(tab.comp.arms) <- c("comp", "multi", "two")
  rownames(tab.comp.arms) <- NULL


  ## Keep only those comparisons not studied in two-arm trials
  tab.comp.arms$select <- ifelse(tab.comp.arms$two == 0 & tab.comp.arms$multi != 0, T, F)
  subs <- subset(tab.comp.arms, select == T, select = comp)


  ## Match the selected comparisons with the study id
  pairwise.n0 <- list()
  for(i in 1:length(subs[, 1])) {
    pairwise.n0[[i]] <- wide.format[which(wide.format$comp == unlist(subs)[i]), 1:3]
  }
  pairwise.n <- do.call(rbind, pairwise.n0)


  ## Sort by the study id in increasing order
  (final0 <- pairwise.n[order(pairwise.n$study),])


  ## Finally, reduce to comparisons between non-baseline interventions
  (final <- final0[!is.element(final0$t1, unique(t[unique(final0$study), 1])),])


  ## Add also the baseline treatment for each selected trial
  base <- rep(NA, length(final[, 1]))
  for(i in 1:length(final[, 1])){
    final$base[i] <- unique(t[final$study[i], 1])
  }
  nbase.multi <- length(final[, 1])


  return(list(nbase.multi = nbase.multi, t1.bn = final$t1, t2.bn = final$t2, base = final$base, obs.comp = tab.comp.arms))


}

