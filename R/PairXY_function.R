############################################################################
#  LOAD DATA MANIPULATING FUNCTIONS
#  PREPARATION  FOR DIAS NODE-SPLITTING APPROACH #
#  Author: Dias et al., 2010 [doi: 10.1002/sim.3767]
############################################################################
#
PairXY <- function(treat, pair) {
  # Check if pair(X,Y) in row i of data
  # and give baseline for data row i

  N <- nrow(treat)
  out <- cbind(split = rep(0, N), b = rep(0, N))
  for (i in 1:N) {
    # returns positions of matches to elements of pair in t[i,]
    # or zero if not present
    pos <- match(pair, treat[i,], nomatch = 0)   # lenght = length(pair) = 2
    out[i, 1] <- ifelse(prod(pos) > 0, 1, 0)      # 1 if pair in line i, 0 o.w.
    out[i, 2] <- ifelse(prod(pos) == 0, 1, pos[1])
  }
  out
}
