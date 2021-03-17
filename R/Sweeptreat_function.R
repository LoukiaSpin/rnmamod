#####################################################################################
#  LOAD DATA MANIPULATING FUNCTIONS - PREPARATION  FOR DIAS NODE-SPLITTING APPROACH #
#####################################################################################
#
Sweeptreat <- function(treat, m)
  # Builds matrix with non-baseline treatments
{
  N <- NROW(treat)
  C <- NCOL(m)
  out <- matrix(nrow=N, ncol=C)
  for (i in 1:N) {
    for (k in 2:C) {
      out[i,k] <- treat[i,m[i,k]]
    }
  }
  out
}
