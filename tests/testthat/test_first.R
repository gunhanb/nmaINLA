
context("Checking meta-analysis example: TB dataset")

### load TB dataset
data(TBdat, package="nmainla")

### Create the dataset suitable for INLA
TBdatINLA <- create_INLA_dat_pair(TBdat$TRT, TBdat$CON, TBdat$TRTTB, TBdat$CONTB)


test_that("results are correct for the data preparation (pairiwise meta-analysis) function.", {

  ### compare with results
  expect_equivalent(round(TBdatINLA$data.cont$Y[3], 3), -1.386)
  expect_equivalent(TBdatINLA$data.arm$Y[3], 29)

})

## Fitting a pairwise random effects meta-analysis model
if(requireNamespace('INLA', quietly = TRUE)){
       require('INLA', quietly = TRUE)
       fit.TB.RE.INLA <- meta_inla(TBdatINLA, type = 'RE',
                                   approach = 'arm-level',
                                   tau.prior = 'uniform', tau.par = c(0, 5))
       }


test_that("results are correct for fitting (pairiwise meta-analysis) function.", {

  ### compare with results
  expect_equivalent(round(fit.TB.RE.INLA$nu[1], 3), -0.761)
  expect_equivalent(round(fit.TB.RE.INLA$nu[2], 3), 0.218)

})

context("Checking a NMA example: Smoking dataset")

### load TB dataset
data("Smokdat", package="nmainla")

### Create the dataset suitable for INLA
SmokdatINLA <- create_INLA_dat(dat = Smokdat,
                               armVars = c('treatment' = 't', 'responders' = 'r',
                                           'sampleSize' = 'n'),
                               nArmsVar = 'na')

test_that("results are correct for the data preparation (network meta-analysis) function.", {

  ### compare with results
  expect_equivalent(SmokdatINLA$responders[1:3], c(9, 23, 10))
  expect_equivalent(SmokdatINLA$baseline[47:50], c(3, 3, 3, 3))

})

## Fitting a pairwise random effects meta-analysis model
if(requireNamespace('INLA', quietly = TRUE)){
  require('INLA', quietly = TRUE)
  fit.Smok.cons.INLA <- nma_inla(SmokdatINLA, likelihood = 'binomial', type = 'consistency',
                                 tau.prior = 'uniform', tau.par = c(0, 5))
}


test_that("results are correct for fitting (network meta-analysis) function.", {

  ### compare with results
  expect_equivalent(round(fit.Smok.cons.INLA$d_params[, 1], 2), c(0.49, 0.84, 1.10))
  expect_equivalent(round(fit.Smok.cons.INLA$hyperpar[1:2], 2), c(0.84, 0.18))

})
