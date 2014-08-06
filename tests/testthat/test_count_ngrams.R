context("Counting n-grams")

test_that("Count unigrams for different distances",{
  sample_seq <- c(1L, 1L, 1L, 2L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 1L, 
                  2L, 1L, 2L, 1L, 1L, 2L, 2L, 1L, 2L, 2L, 2L, 1L, 2L, 1L, 2L)
  ngrams1 <- c("1.1.1", "2.1.1", "1.2.1", "2.2.1", "1.1.2", "2.1.2", "1.2.2", 
               "2.2.2")
  count_ngrams(sample_seq, ngrams1, 3, d = 1)
  
  len1 <- sum(count_ngrams(sample_seq, ngrams1, 3, d = 0))
  len2 <- sum(count_ngrams(sample_seq, ngrams1, 3, d = c(1, 1)))
  len3 <- sum(count_ngrams(sample_seq, ngrams1, 3, d = c(2, 1)))
  len4 <- sum(count_ngrams(sample_seq, ngrams1, 3, d = c(2, 2)))
  expect_equal(len1, 28)
  expect_equal(len2, 26)
  expect_equal(len3, 25)
  expect_equal(len4, 24)
})