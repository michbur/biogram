context("Counting n-grams")

test_that("Count ngrams for different distances",{
  sample_seq <- c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 
                  2L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 2L)
  
  len1 <- sum(count_ngrams(sample_seq, 3, 1L:4, d = 0))
  len2 <- sum(count_ngrams(sample_seq, 3, 1L:4, d = 1))
  len3 <- sum(count_ngrams(sample_seq, 3, 1L:4, d = list(c(2, 1))))
  len4 <- sum(count_ngrams(sample_seq, 3, 1L:4, d = list(c(2, 2))))
  
  expect_equal(len1, 28)
  expect_equal(len2, 26)
  expect_equal(len3, 25)
  expect_equal(len4, 24)
  
  expect_equal(count_ngrams(sample_seq, 1, 1L:2, d = 0), 
               structure(c(13, 17), .Dim = 1:2, .Dimnames = list(NULL, c("1_0", "2_0"))))
})


test_that("Count ngrams for different positions",{
  sample_seq <- c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L)
  
  all_pos <- count_ngrams(sample_seq, 1, 1L:2, pos = TRUE)
  expect_equal(all_pos[, all_pos != 0], 
               structure(c(1, 1, 1, 1, 1, 1, 1, 1), 
                         .Names = c("1_1_0", "2_1_0", "3_1_0", "8_1_0", 
                                    "4_2_0", "5_2_0", "6_2_0", "7_2_0"))) 
})