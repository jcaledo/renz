library(renz)
context("Linear Transformations")


## ----------------------------------------------- ##
#         Testing the function lb()                 #
## ----------------------------------------------- ##
test_that('lb() works properly with and withouth weighting', {

  data <- ONPG
  data <- data[, c(1,3)]
  data$v2[8] <- NA
  a <- lb(data)
  b <- lb(data, weighting = TRUE)

  expect_is(a, 'list')
  expect_equal(length(a), 5)
  expect_equivalent(a$fitted_parameters[1], -2.81)
  expect_equivalent(a$fitted_parameters[2], -75.12)

  expect_is(b, 'list')
  expect_equivalent(length(b), 5)
  expect_equivalent(b$fitted_parameters[1], 3.12)
  expect_equivalent(b$fitted_parameters[2], 129.82)
})


## ----------------------------------------------- ##
#         Testing the function hw()                 #
## ----------------------------------------------- ##
test_that('hw() function works properly', {

  data <- ONPG
  data <- data[, c(1,3)]
  data$v2[8] <- NA
  a <- hw(data)

  expect_is(a, 'list')
  expect_equal(length(a), 5)
  expect_equivalent(a$fitted_parameters[1], 3.87)
  expect_equivalent(a$fitted_parameters[2], 134.94)
})

## ----------------------------------------------- ##
#         Testing the function eh()                 #
## ----------------------------------------------- ##
test_that('eh() function works properly', {

  data <- ONPG
  data <- data[, c(1,3)]
  data$v2[8] <- NA
  a <- eh(data)

  expect_is(a, 'list')
  expect_equal(length(a), 5)
  expect_equivalent(a$fitted_parameters[1], 2.29)
  expect_equivalent(a$fitted_parameters[2], 104.53)
})


## ----------------------------------------------- ##
#         Testing the function ecb()                #
## ----------------------------------------------- ##
test_that('ecb() works properly', {

  data <- ONPG
  data <- data[, c(1,3)]
  data$v2[8] <- NA
  a <- ecb(data)

  expect_is(a, 'list')
  expect_equal(length(a), 3)
  expect_equivalent(a$fitted_parameters[1], 2.373)
  expect_equivalent(a$fitted_parameters[2], 168.65)
})
