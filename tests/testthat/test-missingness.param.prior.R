

test_that("The argument 'mean_misspar' is correctly specified", {
  expect_equal( missingness_param_prior(assumption = "IDE-ARM", mean_misspar = c(0, 1)), c(0.0001, 1) )
  expect_equal( missingness_param_prior(assumption = "IDE-TRIAL", mean_misspar = 0), 0 )
  expect_equal( missingness_param_prior(assumption = "IDE-COMMON", mean_misspar = 0), 0 )
  expect_equal( missingness_param_prior(assumption = "HIE-ARM", mean_misspar = c(0, 0)), c(0.0001, 0.0001) )
  expect_equal( missingness_param_prior(assumption = "HIE-TRIAL", mean_misspar = 0), 0 )
  expect_equal( missingness_param_prior(assumption = "HIE-COMMON", mean_misspar = 0), 0 )
  expect_equal( missingness_param_prior(assumption = "IND-CORR", mean_misspar = 0), 0 )
  expect_equal( missingness_param_prior(assumption = "IND-UNCORR", mean_misspar = 0), 0 )
  expect_equal( missingness_param_prior(assumption = "IND-UNCORR"), 0.0001 )
  expect_equal( missingness_param_prior(assumption = "IDE-ARM"), c(0.0001, 0.0001) )
  expect_equal( missingness_param_prior(assumption = "HIE-ARM"), c(0.0001, 0.0001) )
})



test_that("It gives the correct value of 'mean_misspar' when this argument is missing", {
  expect_equal( missingness_param_prior(assumption = "IND-UNCORR"), 0.0001 )
  expect_equal( missingness_param_prior(assumption = "IDE-ARM"), c(0.0001, 0.0001) )
  expect_equal( missingness_param_prior(assumption = "HIE-ARM"), c(0.0001, 0.0001) )
})



test_that("The argument 'mean_misspar' is erroneously specified", {
  expect_error( missingness_param_prior(assumption = "IDE-ARM", mean_misspar = c(0, 0, 0)) )
  expect_error( missingness_param_prior(assumption = "IDE-TRIAL", mean_misspar = c(0, 0, 0)) )
  expect_error( missingness_param_prior(assumption = "IDE-COMMON", mean_misspar = c(0, 0, 0)) )
  expect_error( missingness_param_prior(assumption = "HIE-ARM", mean_misspar = c(0, 0, 0)) )
  expect_error( missingness_param_prior(assumption = "HIE-TRIAL", mean_misspar = c(0, 0, 0)) )
  expect_error( missingness_param_prior(assumption = "HIE-COMMON", mean_misspar = c(0, 0, 0)) )
  expect_error( missingness_param_prior(assumption = "IND-CORR", mean_misspar = c(0, 0, 0)) )
  expect_error( missingness_param_prior(assumption = "IND-UNCORR", mean_misspar = c(0, 0, 0)) )
  expect_error( missingness_param_prior(assumption = "IND-UNCOR", mean_misspar = c(0, 0, 0)) )
})



test_that("The argument 'assumption' is erroneously specified", {
  expect_error( missingness_param_prior(assumption = "IND-UNCOR", mean_misspar = c(0, 0, 0)) )
})

