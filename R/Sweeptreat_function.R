##############################################################################
#  LOAD DATA MANIPULATING FUNCTIONS - P
#  REPARATION  FOR DIAS NODE-SPLITTING APPROACH
#  Author: Dias et al., 2010 [doi: 10.1002/sim.3767]
##############################################################################
#
Sweeptreat <- function(treat, m) {
  # Builds matrix with non-baseline treatments
  N <- NROW(treat)
  C <- NCOL(m)
  out <- matrix(nrow = N, ncol = C)
  for (i in 1:N) {
    for (k in 2:C) {
      out[i, k] <- treat[i, m[i, k]]
    }
  }
  out
}
