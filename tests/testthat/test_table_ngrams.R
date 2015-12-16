context("n-gram tabulation")

test_that("tabulating unpositioned n-grams",{
  
  data(human_cleave)
  
  expect_identical(table_ngrams(human_cleave[, -10], c("15.15_0", "11.11_0"), human_cleave[, 10])[["target1"]], 
                   c(5, 0))
})

test_that("tabulating positioned n-grams",{
  
  data(human_cleave)
  
  expect_identical(table_ngrams(human_cleave[, -10], c("1_15.15_0", "1_11.11_0"), human_cleave[, 10])[["target1"]], 
                   c(2, 0))
})
