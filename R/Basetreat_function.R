#####################################################################################
#  LOAD DATA MANIPULATING FUNCTIONS - PREPARATION  FOR DIAS NODE-SPLITTING APPROACH #
#####################################################################################
#
Basetreat <- function(treat, b)
  # Builds vector with baseline treatments
{
  N <- nrow(treat)
  out <- rep(0,N)
  for (i in 1:N) {
    out[i] <- treat[i,b[i]]
  }
  out
}