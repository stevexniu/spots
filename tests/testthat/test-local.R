# Test for local Moran's I
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

test_that("Local Moran's I calculated correctly", {
  test.locali<- LocalMoransI(test.X, test.W)
  expect_equal(colSums(test.locali$Local.Morans.I), test.moransi)
})

test_that("Permutation local Moran's I calculated correctly", {
  test.permi <- apply(test.X, 2, PermutationLocalI, W = test.W)
  expect_equal(unlist(lapply(test.permi, function(x) sum(x[['Local.Morans.I']]))), test.moransi)
})
