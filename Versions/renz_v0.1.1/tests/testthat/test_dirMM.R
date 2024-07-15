library(renz)
context("Direct fitting")


## ----------------------------------------------- ##
#       Testing the function dir.MM()               #
## ----------------------------------------------- ##
test_that('dir.MM() works properly', {

  data <- ONPG
  data <- data[, c(1,2)]
  a <- dir.MM(data)

  expect_is(a, 'list')
  expect_equal(length(a), 2)
  expect_equivalent(a$parameters[1], 3.355)
  expect_equivalent(a$parameters[2], 156.013)
  expect_is(a$data, 'data.frame')
  expect_equal(dim(a$data), c(10,3))
})

test_that('dir.MM() works properly when NA are present', {

  data <- ONPG
  data <- data[, c(1,3)]
  data[8,2] <- NA
  a <- dir.MM(data)

  expect_is(a, 'list')
  expect_equal(length(a), 2)
  expect_equivalent(a$parameters[1], 4.162)
  expect_equivalent(a$parameters[2], 135.114)
  expect_is(a$data, 'data.frame')
  expect_equal(dim(a$data), c(9,3))
})
