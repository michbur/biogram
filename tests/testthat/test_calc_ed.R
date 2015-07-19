context("Distance between encodings")

test_that("Distance 0",{
  
  l1 <- list('1' = c("a", "b", "c"),
             '2' = c("d", "e"))
  l2 <- list('1' = c("a", "b", "c"),
             '2' = c("d", "e"))
  
  expect_equal(calc_ed(l1, l2), 0)
  
  l3 <- list('1' = c("a", "b", "c"),
             '2' = c("d", "e"),
             '3' = c(),
             '4' = c())
  
  expect_equal(calc_ed(l3, l1), 0)
})

test_that("Distance 1, including empty groups",{
  
  l1 <- list('1' = c("a", "b", "c"),
             '2' = c("d", "e"))
  l2 <- list('1' = c("a", "b"),
             '2' = c("d", "e", "c"))
  
  expect_equal(calc_ed(l1, l2), 1)
  
#   l3 <- list('1' = c("a", "b", "c", "d"),
#              '2' = c("e"))
#   l4 <- list('1' = c("a", "b", "d", "e", "c"))
#   expect_equal(calc_ed(l3, l4), 1)
#   
#   l5 <- list('1' = c("a", "b", "d", "e", "c"),
#              '2' = c())
#   expect_equal(calc_ed(l3, l5), 1)
})

