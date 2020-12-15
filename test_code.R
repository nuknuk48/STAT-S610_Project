context("Check ABC Functions")
source("Raw_code.R")
## This Test_code.R file works with Raw_code.R file together

library(testthat)

## simulate a count matrix using original data H3N2_78
H3N2_78<- matrix(data = c(66, 87, 25, 22, 4, 13, 14, 15,
                          9, 4, 0, 4, 4, 9, 1, 0, 0, 4, 3, 1,
                          0, 0, 0, 1, 1,
                          0,0,0,0,0), ncol = 5, byrow = TRUE)
H3N2_81<- matrix(data = c(44, 10, 0, 0, 0,0,
                          62, 13,9, 0,0,0,
                          47,8,2,3,0,0,
                          38,11,7,5,1,0,
                          9,5,3,1,0,1), ncol = 5)

##----------------------------------------------------------------------
## Testing Data_generator_helper Function, which returns a count matrix
generated = Data_generator_helper(runif(4),5,2,H3N2_78)

test_that("Testing the class of generated data structure", {
  expect_equal(class(generated),c("matrix", "array"))
})

test_that("Testing the dimension of generated the same as original data", {
  expect_equal(dim(generated), dim(H3N2_78))
})

test_that("Testing the places in first column of generated where zeros should be", {
  expect_equal(generated[(3:6),1],rep(0,4))
})

test_that("Testing the places in second column of generated where zeros should be", {
  expect_equal(generated[(4:6),2],rep(0,3))
})

test_that("Testing whether the column sums of generated data is the same as the original data",{
  expect_true(sum(colSums(generated)==colSums(H3N2_78))==5)
})

##------------------------------------------------------------------------
## Testing Data_generator Function, which returns a list for two datasets
data_generated <- Data_generator(runif(4), 5, 5, H3N2_78, H3N2_81)

test_that("data_generated is a list structure", {
  expect_equal(class(data_generated),"list")
})

test_that("Testing the dimension of first matrix in data_generated the same as original data", {
  expect_equal(dim(data_generated$c1), dim(H3N2_78))
})

test_that("Testing the dimension of second matrix in data_generated the same as original data", {
  expect_equal(dim(data_generated$c2), dim(H3N2_81))
})
##------------------------------------------------------------------------
## Testing Distance Function
dist <- Distance(H3N2_78,H3N2_81,data_generated)

test_that("dist is a single numeric", {
  expect_equal(class(dist),"numeric")
})

##-------------------------------------------------------------------------
## Testing ABC Function

m1 = ABC(50, 1000, H3N2_78, H3N2_81)
## m1_matrix contains the probabilities of qh1, qc1, qh2, qc2
m1_matrix = matrix(unlist(m1), ncol = 4, byrow = TRUE)

test_that("Testing m1 is a list structure", {
  expect_equal(class(m1),"list")
})

test_that("Testing each in m1_matrix is a probability", {
  expect_true(all(m1_matrix>=0 & m1_matrix<=1))
})









