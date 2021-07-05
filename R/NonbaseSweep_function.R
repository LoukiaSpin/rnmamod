#####################################################################################
#  LOAD DATA MANIPULATING FUNCTIONS - PREPARATION  FOR DIAS NODE-SPLITTING APPROACH #
#  Author: Dias et al., 2010 [doi: 10.1002/sim.3767]                                #
#####################################################################################
#
NonbaseSweep <- function(index, na)
  # gives na-1 indexes to sweep non-baseline arms only
{
  N <- NROW(na)
  C <- max(na)
  out <- matrix(nrow=N, ncol=C)
  for (i in 1:N) {
    for (k in 2:na[i]) {
      out[i,k] <- k - (index[i,"b"] >= k)
    }
  }
  out
}
