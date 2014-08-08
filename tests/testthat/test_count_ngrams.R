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
  
#   expect_equal(count_ngrams(sample_seq, 1, 1L:2, d = 0), c("1_0" = 13, "2_0" = 17))
})

# 
# test_that("Count ngrams for different positions",{
#   sample_seq <- c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L)
#   
#   all_pos <- count_ngrams(sample_seq, 1, 1L:2, pos = TRUE)
#   expect_equal(all_pos[all_pos != 0], structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), 
#                                                 .Names = c("1_1", "2_1", "3_1", "8_1", "4_2", 
#                                                            "5_2", "6_2", "7_2")))
#   
# })