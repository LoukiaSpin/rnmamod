possible.observed.frail.comparisons <- function(drug.names, obs.comp) {


  ## Obtain all unique pairwise comparisons using the 'combn' functions
  nt <- length(drug.names)
  poss.pair.comp <- data.frame(t(combn(1:nt, 2))[, 2], t(combn(1:nt, 2))[, 1])
  poss.pair.comp$comp <- paste0(poss.pair.comp[, 1], "vs", poss.pair.comp[, 2])
  colnames(poss.pair.comp) <- c("treat1", "treat2", "comp")


  ## Restrict to the observed comparisons in the network
  observed.comp <- poss.pair.comp[is.element(poss.pair.comp$comp, obs.comp), ]


  ## Replace intervention id with their original name
  # All possible comparisons - Treat1 (non-baseline arm)
  for(i in sort(unique(unlist(poss.pair.comp[, 1])))) {
    poss.pair.comp[poss.pair.comp$treat1 == i, 1] <- drug.names[i]
  }
  # Observed comparisons - Treat2 (baseline arm)
  for(i in sort(unique(unlist(poss.pair.comp[, 2])))) {
    poss.pair.comp[poss.pair.comp$treat2 == i, 2] <- drug.names[i]
  }
  poss.pair.comp$comp.name <- paste(poss.pair.comp[, 1], "vs", poss.pair.comp[, 2])


  # Observed comparisons - Treat1 (non-baseline arm)
  for(i in sort(unique(unlist(observed.comp[, 1])))) {
    observed.comp[observed.comp$treat1 == i, 1] <- drug.names[i]
  }
  # Observed comparisons - Treat2 (baseline arm)
  for(i in sort(unique(unlist(observed.comp[, 2])))) {
    observed.comp[observed.comp$treat2 == i, 2] <- drug.names[i]
  }
  observed.comp$comp.name <- paste(observed.comp[, 1], "vs", observed.comp[, 2])


  return(list(poss.comp = data.frame(ID = 1:length(poss.pair.comp$comp), poss.pair.comp),
              obs.comp = observed.comp))

}
