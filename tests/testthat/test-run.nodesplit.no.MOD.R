test_that("The summary log ORs and between-trial sd agree (no missing outcome data)", {

  set.seed(123)

  ## Example taken from the 'mtc.network' function of the 'gemtc' R-package
  data <- read.table(textConnection('
  study  treatment  responders  sampleSize
  01     A          2           100
  01     B          5           100
  02     B          6           110
  02     C          1           110
  03     A          3           60
  03     C          4           80
  03     B          7           80'), header = TRUE)

  network <- gemtc::mtc.network(data)

  ## Fixed-effect NMA
  #node.fe <- gemtc::mtc.nodesplit(network, linearModel = "fixed")
  #gemtc.fe <- rbind(summary(node.fe)$dir.effect[3:5], summary(node.fe)$ind.effect[3:5])
  #colnames(gemtc.fe) <- rownames(gemtc.fe) <- NULL

  ## Random-effects NMA
  node.re <- gemtc::mtc.nodesplit(network, linearModel = "random", hy.prior =  gemtc::mtc.hy.prior("std.dev", "dhnorm", 0, 1))
  gemtc.re <- rbind(summary(node.re)$dir.effect[3:5], summary(node.re)$ind.effect[3:5])
  colnames(gemtc.re) <- rownames(gemtc.re) <- NULL


  ## Adjust the example above to run the 'run_model' function of the 'rnmamod' R-package
  data.rnmamod <- read.table(textConnection('
  t1  t2  t3  r1  r2  r3  n1  n2  n3
  1   2   NA  2   5   NA  100 100 NA
  2   3   NA  6   1   NA  110 110 NA
  1   2   3   3   7   4   60  80  80'), header=TRUE)

  ## Fixed-effect NMA
  #res.fe <- run_model(data = data.rnmamod, measure = "OR", model = "FE", D = 0)
  #node.fe2 <- run_nodesplit(res.fe)
  #resnode.fe <- rbind(node.fe2$direct[c(3, 5:6)], node.fe2$indirect[c(3, 5:6)])
  #colnames(resnode.fe) <- rownames(resnode.fe) <- NULL

  ## Random-effects NMA
  res.re <- run_model(data = data.rnmamod, measure = "OR", model = "RE", heter_prior = list("halfnormal", 0, 1), D = 0, ref = 1)
  node.re2 <- run_nodesplit(res.re)
  resnode.re <- rbind(node.re2$direct[c(3, 5:6)], node.re2$indirect[c(3, 5:6)])
  colnames(resnode.re) <- rownames(resnode.re) <- NULL

  #expect_equal(resnode.fe, gemtc.fe, tolerance = 4e-2)
  #expect_equal(resnode.re, gemtc.re, tolerance = 4e-2)
})

