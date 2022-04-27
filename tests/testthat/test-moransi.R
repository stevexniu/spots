# Test for Moran's I calculation
library(stats)
library(ape)
set.seed(42)
test.X <- matrix(runif(300), 30, 10)
test.W <- matrix(runif(900), 30, 30)
test.r <- apply(test.X, 2, ape::Moran.I, weight = test.W)
test.moransi <- unlist(lapply(test.r, "[[", "observed"))
test.mean <- unlist(lapply(test.r, "[[", "expected"))
test.sd <- unlist(lapply(test.r, "[[", "sd"))
test.pval <- unlist(lapply(test.r, "[[", "p.value"))
test.z <- (test.moransi - test.mean) / test.sd

test_that("Univariate Moran's I calculated correctly", {
  test.univar <- UnivariateMoransI(test.X, test.W)
  expect_equal(test.univar$Morans.I, test.moransi)
  expect_equal(test.univar$Expected.I, unique(test.mean))
  expect_equal(test.univar$SD.I, test.sd)
  expect_equal(test.univar$Z.I, test.z)
  expect_equal(test.univar$p.val, test.pval)
})

test_that("Permutation Moran's I calculated correctly", {
  test.perm <- apply(test.X, 2, PermutationMoransI, W = test.W)
  expect_equal(unlist(lapply(test.perm, "[[", "Morans.I")), test.moransi)
})

test_that("Bivariate Moran's I calculated correctly", {
  test.bivar <- BivariateMoransI(test.X, test.W)
  expect_equal(diag(test.bivar$Expected.I), test.mean)
  expect_true(cor.test(diag(test.bivar$Morans.I), test.moransi)$p.value < 0.001)
  expect_true(cor.test(diag(test.bivar$SD.I), test.sd)$p.value < 0.001)
  expect_true(cor.test(diag(test.bivar$p.val), test.pval)$p.value < 0.001)
  expect_true(cor.test(diag(test.bivar$Z.I), test.z)$p.value < 0.001)
})

test_that("OLS Moran's I calculated correctly", {
  test.OLS <- OLSMoransI(test.X, test.W)
  expect_equal(diag(test.OLS$Morans.I), test.moransi)
  expect_equal(diag(test.OLS$Expected.I), test.mean)
  expect_equal(diag(test.OLS$SD.I), test.sd)
  expect_equal(diag(test.OLS$Z.I), test.z)
  expect_equal(diag(test.OLS$p.val), test.pval)
})
