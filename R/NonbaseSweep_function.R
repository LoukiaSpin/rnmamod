#############################################################################
#  LOAD DATA MANIPULATING FUNCTIONS
#  PREPARATION  FOR DIAS NODE-SPLITTING APPROACH #
#  Author: Dias et al., 2010 [doi: 10.1002/sim.3767]
#############################################################################

NonbaseSweep <- function(index, na) {

  # gives na-1 indexes to sweep non-baseline arms only
  N <- NROW(na)
  C <- max(na)
  out <- matrix(nrow = N, ncol = C)
  for (i in 1:N) {
    for (k in 2:na[i]) {
      # Refinement when the third arm is the baseline
      out[i, k] <- ifelse(k == index[i, "b"], 1,
                         ifelse(k < index[i, "b"], k + 1, k))
    }
  }
  out
}
