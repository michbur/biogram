context("Converting sequence to n-grams")

test_that("Extract ngrams for different distances",{
  seqs <- structure(c(4L, 2L, 1L, 1L, 4L, 3L, 1L, 3L, 1L, 1L), .Dim = c(2L, 
                                                                        5L))
  expect_equal(seq2ngrams(seqs, 3, 1L:4), 
               structure(c("4.1.4_0_0", "2.1.3_0_0", "1.4.1_0_0", "1.3.3_0_0", 
                           "4.1.1_0_0", "3.3.1_0_0"), .Dim = 2:3))
  
  #the same with distance
  expect_equal(seq2ngrams(seqs, 3, 1L:4, d = c(0, 1)), 
               structure(c("4.1.1_0_1", "2.1.3_0_1", "1.4.1_0_1", "1.3.1_0_1"), 
                         .Dim = c(2L, 2L)))
  
  #only one n-gram possible
  expect_equal(seq2ngrams(seqs, 3, 1L:4, d = 1), 
               structure(c("4.4.1_1_1", "2.3.1_1_1"), .Dim = c(2L, 1L)))
  
})