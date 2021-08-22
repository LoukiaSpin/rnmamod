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
  model.fe <- gemtc::mtc.model(network, linearModel = "fixed")
  gemtc.fe <- gemtc::mtc.run(model.fe)

  ## Random-effects NMA
  model.re <- gemtc::mtc.model(network, linearModel = "random", hy.prior =  gemtc::mtc.hy.prior("std.dev", "dhnorm", 0, 1))
  gemtc.re <- gemtc::mtc.run(model.re)


  ## Adjust the example above to run the 'run.model' function of the 'rnmamod' R-package
  data.rnmamod <- read.table(textConnection('
  t1  t2  t3  r1  r2  r3  n1  n2  n3
  1   2   NA  2   5   NA  100 100 NA
  2   3   NA  6   1   NA  110 110 NA
  1   2   3   3   7   4   60  80  80'), header=TRUE)

  res.fe <- run.model(data = data.rnmamod, measure = "OR", model = "FE", D = 0)
  res.re <- run.model(data = data.rnmamod, measure = "OR", model = "RE", heter.prior = list("halfnormal", 0, 1), D = 0)

  expect_equal( as.vector(res.fe$EM[1:2, 1:2]), as.vector(summary(gemtc.fe)$summaries$statistics[1:2, 1:2]), tolerance = 4e-2)
  expect_equal( as.vector(res.re$EM[1:2, 1:2]), as.vector(summary(gemtc.re)$summaries$statistics[1:2, 1:2]), tolerance = 4e-2)
  expect_equal( as.vector(res.re$tau[1:2]), as.vector(summary(gemtc.re)$summaries$statistics[3, 1:2]), tolerance = 4e-2)
})

