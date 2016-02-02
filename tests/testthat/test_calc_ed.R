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
  
  l3 <- list('1' = c("a", "b", "c", "d"),
             '2' = c("e"))
  l4 <- list('1' = c("a", "b", "d", "e", "c"))
  expect_equal(calc_ed(l3, l4), 1)
  
  l5 <- list('1' = c("a", "b", "d", "e", "c"),
             '2' = c())
  expect_equal(calc_ed(l3, l5), 1)
})


test_that("Distance 2, including empty groups",{
  
  l1 <- list('1' = c("a", "b", "c"),
             '2' = c("d", "e"))
  l2 <- list('1' = c("a", "e"),
             '2' = c("d", "c", "b"))

  expect_equal(calc_ed(l1, l2), 2)
  
  l3 <- list('1' = c("a", "b", "c"),
             '2' = c("d", "e"),
             '3' = c())
  expect_equal(calc_ed(l3, l2), 2)
})

test_that("Distance 2, single groups",{
  
  l1 <- list('1' = "a",
             '2' = c("d", "e", "f"),
             '3' = "c",
             '4' = "b")
  
  l2 <- list('1' = c("a", "b", "c"),
             '2' = c("d", "e", "f"))
  
  expect_equal(calc_ed(l1, l2), 2)
  

  l3 <- list('1' = c("a", "b"),
             '2' = c("d", "e"),
             '3' = "c",
             '4' = "f")
  
  expect_equal(calc_ed(l3, l2), 2)
})


test_that("Distance 3",{
  l1 <- list('1' = c("a", "b"),
             '2' = c("d", "e"),
             '3' = c("g", "h"),
             '4' = "i",
             '5' = "f",
             '6' = "c")
  
  l2 <- list('1' = c("a", "b", "c"),
             '2' = c("d", "e", "f"),
             '3' = c("g", "h", "i"))
  
  expect_equal(calc_ed(l1, l2), 3)
 })


test_that("Long identical groups",{
  
  aa1 = list(`1` = c("g", "a", "p", "v", "h", "l", "i"), 
             `2` = c("k", "m"), 
             `3` = c("d", "r"), 
             `4` = c("f", "e", "w", "y", "s", "t", "c", "n", "q"))
  
  aa2 = list(`1` = c("g", "a", "p", "v", "h", "l", "i"), 
             `2` = c("k"), 
             `3` = c("d", "r", "m"), 
             `4` = c("f", "e", "w", "y", "s", "t", "c", "n", "q"))
  expect_equal(calc_ed(aa1, aa2), 1)  
})

test_that("Half-half case",{
  aa1 = list(`1` = c("g", "a", "p", "v", "h", "l", "i","k", "m"), 
             `2` = c("d", "r", "f", "e", "w", "y", "s", "t", "c", "n", "q"))
  
  aa2 = list(`1` = c("g", "a", "p", "v", "m", "l", "q"), 
             `2` = c("k", "h", "d", "e", "i"), 
             `3` = c("f", "r", "w", "y", "s", "t", "c", "n"))
  expect_equal(calc_ed(aa2, aa1), 8)
})
