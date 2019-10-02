library(renz)
context("Simulated Enzyme-catalyzed Reaction Progress Curves")

## ----------------------------------------------- ##
#       Testing the function sEprogress             #
## ----------------------------------------------- ##
a <- sEprogress(So = 3, time = 10, Km = 1, Vm = 1, replicates = 3, error = 'a')
b <- sEprogress(So = 3, time = 10, Km = 1, Vm = 1, replicates = 3, error = 'r')

test_that('sEprogress works properly', {
  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 64)
  expect_equal(ncol(a), 7)
  expect_equal(max(a), 6.3)
})

test_that('sEprogress works properly', {
  expect_is(a, 'data.frame')
  expect_equal(nrow(b), 64)
  expect_equal(ncol(b), 7)
  expect_equal(max(b), 6.3)
})

a <- sEprogress(So = 15, time = 10, Km = 5, Vm = 50, replicates = 3, error = 'a')
b<- sEprogress(So = 15, time = 10, Km = 5, Vm = 50, replicates = 3, error = 'r')

test_that('sEprogress works properly', {
  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 8)
  expect_equal(ncol(a), 7)
  expect_equal(max(a), 15)
})

test_that('sEprogress works properly', {
  expect_is(a, 'data.frame')
  expect_equal(nrow(b), 8)
  expect_equal(ncol(b), 7)
  expect_equal(max(b), 15)
})
