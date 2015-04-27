context("n-gram validation")

test_that("Real n-grams", {
  
  expect_true(is_ngram("1_1.1.1_0.0"))

})