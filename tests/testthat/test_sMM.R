library(renz)
context("Simulated Kinetics")

## ----------------------------------------------- ##
#         Testing the function sMM                  #
## ----------------------------------------------- ##
a <- sMM(Km = 1, Vm = 1, replicates = 3, error = 'a', plot = F)

test_that('sMM works properly', {
  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 7)
  expect_equal(ncol(a), 7)
  expect_equal(max(a), 10)
  expect_lt(max(a[,-1]), 1)
  expect_gt(sum(a$v_sd), 0.35)
})

b <- sMM(Km = 1, Vm = 1, replicates = 3, error = 'f', plot = F)

test_that('sMM works properly', {
  expect_is(b, 'data.frame')
  expect_equal(nrow(b), 7)
  expect_equal(ncol(b), 7)
  expect(max(b), 10)
  expect_lt(max(b[,-1]), 1)
  expect_equal(sum(b$v_sd), 0.35)
})

c <- sMM(Km = 1, Vm = 1, replicates = 2, error = 'r', plot = F)

test_that('sMM works properly', {
  expect_is(c, 'data.frame')
  expect_equal(nrow(c), 7)
  expect_equal(ncol(c), 6)
  expect(max(c), 10)
  expect_lt(max(c[,-1]), 1)
  expect_lt(sum(c$v_sd), 0.35)
})

d <- sMM(Km = 1, Vm = 1, replicates = 1, error = 'a', plot = F)

test_that('sMM works properly', {
  expect_is(d, 'data.frame')
  expect_equal(nrow(d), 7)
  expect_equal(ncol(d), 5)
  expect_equal(sum(d$v_sd), 0)
})

## ----------------------------------------------- ##
#         Testing the function rMM                  #
## ----------------------------------------------- ##

a <- rMM(replicates = 3, error = 'a', plot = F)

test_that('rMM works properly', {
  expect_is(a, 'data.frame')
  expect_equal(nrow(a), 7)
  expect_equal(ncol(a), 7)
  expect_gt(max(a), 83)
  expect_equal(max(a[,1]), 73)
  expect_gt(sum(a$v_sd), 50)
  expect_equal(attributes(a)$Km, 7.3)
  expect_equal(attributes(a)$Vm, 91)
})

