library(renz)
context("Multisubstrate")


## ----------------------------------------------- ##
#       Testing the function bibi()                 #
## ----------------------------------------------- ##
test_that('bibi() works properly', {

  ab <- bibi(hk)
  ba <- bibi(hk, vice_versa = TRUE)
  tc <- bibi("./TcTS.txt")

  expect_is(ab, 'list')
  expect_equal(length(ab), 3)
  expect_is(ab[[1]], 'character')
  expect_is(ab[[2]], 'numeric')
  expect_is(ab[[3]], 'numeric')

  expect_is(ba, 'list')
  expect_equal(length(ab), 3)
  expect_is(ba[[1]], 'character')
  expect_is(ba[[2]], 'numeric')
  expect_is(ba[[3]], 'numeric')

  expect_equal(ab[[1]][1], ba[[1]][1])
  expect_equivalent(ab[[1]][3], ba[[1]][4])

  expect_is(tc, 'list')
  expect_equal(length(ab), 3)
  expect_is(ab[[1]], 'character')
  expect_is(ab[[2]], 'numeric')
  expect_is(ab[[3]], 'numeric')
})
