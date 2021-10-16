possible_observed_comparisons <- function(drug_names, obs_comp) {

  if (length(drug_names) < 3) {
    stop("This function is *not* relevant for a pairwise meta-analysis",
         call. = F)
  }

  # Obtain all unique pairwise comparisons using the 'combn' functions
  nt <- length(drug_names)
  poss_pair_comp <- data.frame(t(combn(1:nt, 2))[, 2], t(combn(1:nt, 2))[, 1])
  poss_pair_comp$comp <- paste0(poss_pair_comp[, 1], "vs", poss_pair_comp[, 2])
  colnames(poss_pair_comp) <- c("treat1", "treat2", "comp")

  # Restrict to the observed comparisons in the network
  observed_comp <- poss_pair_comp[is.element(poss_pair_comp$comp, obs_comp), ]

  # Replace intervention id with their original name
  # All possible comparisons - Treat1 (non-baseline arm)
  for (i in sort(unique(unlist(poss_pair_comp[, 1])))) {
    poss_pair_comp[poss_pair_comp$treat1 == i, 1] <- drug_names[i]
  }
  # Observed comparisons - Treat2 (baseline arm)
  for (i in sort(unique(unlist(poss_pair_comp[, 2])))) {
    poss_pair_comp[poss_pair_comp$treat2 == i, 2] <- drug_names[i]
  }
  poss_pair_comp$comp_name <-
    paste(poss_pair_comp[, 1], "vs", poss_pair_comp[, 2])

  # Observed comparisons - Treat1 (non-baseline arm)
  for (i in sort(unique(unlist(observed_comp[, 1])))) {
    observed_comp[observed_comp$treat1 == i, 1] <- drug_names[i]
  }
  # Observed comparisons - Treat2 (baseline arm)
  for (i in sort(unique(unlist(observed_comp[, 2])))) {
    observed_comp[observed_comp$treat2 == i, 2] <- drug_names[i]
  }
  observed_comp$comp_name <- paste(observed_comp[, 1], "vs", observed_comp[, 2])

  return(list(poss_comp = data.frame(ID = seq_len(length(poss_pair_comp$comp)),
                                     poss_pair_comp),
              obs_comp = observed_comp))

}
