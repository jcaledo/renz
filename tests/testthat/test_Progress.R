library(renz)
context("Simulated Enzyme-catalyzed Reaction Progress Curves")

## ----------------------------------------------- ##
#       Testing the function sE.progress()          #
## ----------------------------------------------- ##
test_that('sE.progress() works properly', {

  a <- sE.progress(So = 3, time = 10, Km = 1, Vm = 1, replicates = 3, error = 'a')
  b <- sE.progress(So = 3, time = 10, Km = 1, Vm = 1, replicates = 3, error = 'r')

  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 87)
  expect_equal(ncol(a), 7)
  expect_lt(max(a), 10)

  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 87)
  expect_equal(ncol(b), 7)
  expect_lt(max(b), 10)
})


test_that('sE.progress() works properly', {

  a <- sE.progress(So = 15, time = 10, Km = 5, Vm = 50, replicates = 3, error = 'a')
  b<- sE.progress(So = 15, time = 10, Km = 5, Vm = 50, replicates = 3, error = 'r')

  expect_is(a, 'data.frame')
  expect_gt(nrow(a), 5)
  expect_equal(ncol(a), 7)
  expect_equal(max(a), 15)

  expect_is(b, 'data.frame')
  expect_gt(nrow(b), 5)
  expect_equal(ncol(b), 7)
  expect_equal(max(b), 15)
})

