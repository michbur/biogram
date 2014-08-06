context("Counting n-grams")

test_that("Count ngrams for different distances",{
  sample_seq <- c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 
                  2L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 2L)
  ngrams1 <- c("1.1.1", "2.1.1", "1.2.1", "2.2.1", "1.1.2", "2.1.2", "1.2.2", 
               "2.2.2")
  count_ngrams(sample_seq, ngrams1, 3, d = 1)
  
  len1 <- sum(count_ngrams(sample_seq, ngrams1, 3, d = 0))
  len2 <- sum(count_ngrams(sample_seq, ngrams1, 3, d = 1))
  len3 <- sum(count_ngrams(sample_seq, ngrams1, 3, d = c(2, 1)))
  len4 <- sum(count_ngrams(sample_seq, ngrams1, 3, d = c(2, 2)))
  
  expect_equal(len1, 28)
  expect_equal(len2, 26)
  expect_equal(len3, 25)
  expect_equal(len4, 24)
  
  expect_equal(count_ngrams(sample_seq, c("1", "2"), 1, d = 0), c("1" = 13, "2" = 17))
})


test_that("Count ngrams for different positions",{
  sample_seq <- c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L)
  
  all_pos <- count_ngrams(sample_seq, 
                          create_ngrams(1, 1:4, possible_grams = 8), 
                          1, 
                          pos = TRUE)
  expect_equal(all_pos[all_pos != 0], structure(c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L), 
                                                .Names = c("1_1", "2_1", "3_1", "8_1", "4_2", 
                                                           "5_2", "6_2", "7_2")))
  
})