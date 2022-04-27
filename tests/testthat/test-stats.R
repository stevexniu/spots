# Test for statistics calculation
library(stats)
set.seed(42)
test.x <- runif(1000)
test.y <- runif(1000)
test.q <- matrix(seq(0.05, 0.95, 0.05))
test.norm <- qnorm(test.q)

test_that("Permutation correlation calculated correctly", {
  test.stats <- PermutationCorr(test.x, test.y)
  test.r <- cor.test(test.x, test.y)
  expect_equal(test.stats$correlation, test.r$estimate, ignore_attr = TRUE)
  expect_equal(test.stats$p.val > 0.05, test.r$p.value > 0.05)
})

test_that("Z score p-values calculated correctly", {
  test.stats <- ZPvalue(Z = test.norm, alternative = "less")
  expect_equal(test.stats$p.val, test.q)
  test.stats <- ZPvalue(Z = test.norm, alternative = "greater")
  expect_equal(test.stats$p.val, 1-test.q)
  test.stats <- ZPvalue(Z = test.norm, alternative = "two.sided", p.adjust.method = "bonferroni")
  expect_equal(test.stats$p.val, ifelse(test.q > 0.5, 2*(1-test.q), 2*test.q))
  expect_equal(test.stats$p.adj, p.adjust(test.stats$p.val, method = "bonferroni"), ignore_attr = TRUE)
})
