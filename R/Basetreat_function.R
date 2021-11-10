##############################################################################
#  LOAD DATA MANIPULATING FUNCTIONS -
#  PREPARATION  FOR DIAS NODE-SPLITTING APPROACH #
#  Author: Dias et al., 2010 [doi: 10.1002/sim.3767]
##############################################################################
#
Basetreat <- function(treat, b) {
  # Builds vector with baseline treatments
  N <- nrow(treat)
  out <- rep(0, N)
  for (i in 1:N) {
    out[i] <- treat[i, b[i]]
  }
  out
}
