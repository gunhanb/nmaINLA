
context("Checking analysis example: TB dataset")

### load TB dataset
data(TBdat, package="nmainla")

### Create the dataset suitable for INLA
TBdatINLA <- create_INLA_dat_pair(TBdat$TRT, TBdat$CON, TBdat$TRTTB, TBdat$CONTB)


test_that("results are correct for the data preparation function.", {

  ### compare with results on page 408
  expect_equivalent(round(TBdatINLA$data.cont$Y[3], 3), -1.386)
  expect_equivalent(TBdatINLA$data.arm$Y[3], 29)

})
