

test_that("The argument 'mean.misspar' is correctly specified", {
  expect_equal( missingness.param.prior(assumption = "IDE-ARM", mean.misspar = c(0, 1)), c(0.0001, 1) )
  expect_equal( missingness.param.prior(assumption = "IDE-TRIAL", mean.misspar = 0), 0 )
  expect_equal( missingness.param.prior(assumption = "IDE-COMMON", mean.misspar = 0), 0 )
  expect_equal( missingness.param.prior(assumption = "HIE-ARM", mean.misspar = c(0, 0)), c(0.0001, 0.0001) )
  expect_equal( missingness.param.prior(assumption = "HIE-TRIAL", mean.misspar = 0), 0 )
  expect_equal( missingness.param.prior(assumption = "HIE-COMMON", mean.misspar = 0), 0 )
  expect_equal( missingness.param.prior(assumption = "IND-CORR", mean.misspar = 0), 0 )
  expect_equal( missingness.param.prior(assumption = "IND-UNCORR", mean.misspar = 0), 0 )
  expect_equal( missingness.param.prior(assumption = "IND-UNCORR"), 0.0001 )
  expect_equal( missingness.param.prior(assumption = "IDE-ARM"), c(0.0001, 0.0001) )
  expect_equal( missingness.param.prior(assumption = "HIE-ARM"), c(0.0001, 0.0001) )
})



test_that("It gives the correct value of 'mean.misspar' when this argument is missing", {
  expect_equal( missingness.param.prior(assumption = "IND-UNCORR"), 0.0001 )
  expect_equal( missingness.param.prior(assumption = "IDE-ARM"), c(0.0001, 0.0001) )
  expect_equal( missingness.param.prior(assumption = "HIE-ARM"), c(0.0001, 0.0001) )
})



test_that("The argument 'mean.misspar' is erroneously specified", {
  expect_error( missingness.param.prior(assumption = "IDE-ARM", mean.misspar = c(0, 0, 0)) )
  expect_error( missingness.param.prior(assumption = "IDE-TRIAL", mean.misspar = c(0, 0, 0)) )
  expect_error( missingness.param.prior(assumption = "IDE-COMMON", mean.misspar = c(0, 0, 0)) )
  expect_error( missingness.param.prior(assumption = "HIE-ARM", mean.misspar = c(0, 0, 0)) )
  expect_error( missingness.param.prior(assumption = "HIE-TRIAL", mean.misspar = c(0, 0, 0)) )
  expect_error( missingness.param.prior(assumption = "HIE-COMMON", mean.misspar = c(0, 0, 0)) )
  expect_error( missingness.param.prior(assumption = "IND-CORR", mean.misspar = c(0, 0, 0)) )
  expect_error( missingness.param.prior(assumption = "IND-UNCORR", mean.misspar = c(0, 0, 0)) )
  expect_error( missingness.param.prior(assumption = "IND-UNCOR", mean.misspar = c(0, 0, 0)) )
})



test_that("The argument 'assumption' is erroneously specified", {
  expect_error( missingness.param.prior(assumption = "IND-UNCOR", mean.misspar = c(0, 0, 0)) )
})

