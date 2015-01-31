context("Converting n-grams to positions")

test_that("Extract ngrams for different distances",{
  seqs <- structure(c(4L, 2L, 1L, 1L, 4L, 3L, 1L, 3L, 1L, 1L), .Dim = c(2L, 
                                                                        5L))
  expect_equal(position_ngrams("2_1.1.2_0.1"), 
               structure(list(`2` = 1, `3` = 1, `5` = 2), 
                         .Names = c("2", "3", "5")))
  
  #multiple n-grams
  expect_equal(position_ngrams(c("2_1.1.2_0.1", "3_1.1.2_0.0", "3_2.2.2_0.0")),
               structure(list(`2` = 1, `3` = c(1, 1, 2), `4` = c(1, 2), 
                              `5` = c(2, 2, 2)), .Names = c("2", "3", "4", "5")))
  #no overlap in positions
  expect_equal(position_ngrams(c("2_1.1.2_0.1", "7_1.1.2_0.0", "10_2.2.2_0.0")),
               structure(list(`2` = 1, `3` = 1, `5` = 2, `7` = 1, `8` = 1, `9` = 2, 
                              `10` = 2, `11` = 2, `12` = 2), 
                         .Names = c("2", "3", "5", "7", "8", "9", "10", "11", "12")))
})
