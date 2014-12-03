context("Computing probability distribution")

test_that("Partial binomial probability is computed correctly",{
  p <- 0.4
  m <- 100
  n <- 20:40
  probs_biogram <- biogram:::calc_partial_binomial_log_probability(p, m, n)
  probs_dbinom <- log(dbinom(n, m, p)/(1-p)^(m-n))
  expect_equal(sum(abs(probs_biogram-probs_dbinom)), 0, tolerance = 1e-12, scale = 1)
  
  p <- 0.001
  m <- 1000
  n <- 500
  expect_equal(dbinom(n, m, p, log=T)- log(1-p)*(m-n), 
               biogram:::calc_partial_binomial_log_probability(p, m, n), 
               tolerance = 1e-12)
})