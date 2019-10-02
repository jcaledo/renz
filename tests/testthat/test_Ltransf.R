library(renz)
context("Linear Transformations")

## ----------------------------------------------- ##
#         Testing the function double.rec           #
## ----------------------------------------------- ##
a <- sMM(Km = 2, Vm = 100, error = 'f')[,c(1,3:5)]
b <- double.rec(data = a)

test_that('double.rec works properly', {
  expect_is(b, 'list')
  expect_equal(length(b), 4)
  expect_equal(length(b$Kms), 3)
  expect_equal(b$Kms[2], 1.05)
  expect_equal(length(b$Vms), 3)
  expect_equal(b$Vms[2], 87.44)
  expect_is(b$inverse_data, 'data.frame')
  expect_equal(ncol(b$inverse_data), 6)
  expect_equal(nrow(b$inverse_data),7)
})

a <- sMM(Km = 2, Vm = 100, replicates = 1, error = 'r', sd = 0.1)[,c(1,3)]
b <- double.rec(data = a)

test_that('double.rec works properly', {
  expect_is(b, 'list')
  expect_equal(length(b), 4)
  expect_equal(length(b$Kms), 1)
  expect_equal(b$Kms, 2.51)
  expect_equal(length(b$Vms), 1)
  expect_equal(b$Vms, 116.32)
  expect_is(b$inverse_data, 'data.frame')
  expect_equal(ncol(b$inverse_data), 2)
  expect_equal(nrow(b$inverse_data),7)
})

## ----------------------------------------------- ##
#         Testing the function HW                   #
## ----------------------------------------------- ##
a <- sMM(Km = 2, Vm = 100, error = 'f')[,c(1,3:5)]
b <- HW(data = a)

test_that('HW function works properly', {
  expect_is(b, 'list')
  expect_equal(length(b), 4)
  expect_equal(length(b$Kms), 3)
  expect_equal(b$Kms[2], 1.66)
  expect_equal(length(b$Vms), 3)
  expect_equal(b$Vms[2], 103.78)
  expect_is(b$transformed_data, 'data.frame')
  expect_equal(ncol(b$transformed_data), 6)
  expect_equal(nrow(b$transformed_data),7)
})

a <- sMM(Km = 2, Vm = 100, replicates = 1, error = 'r', sd = 0.1)[,c(1,3)]
b <- HW(data = a)

test_that('HW function works properly', {
  expect_is(b, 'list')
  expect_equal(length(b), 4)
  expect_equal(length(b$Kms), 1)
  expect_equal(b$Kms, 2.08)
  expect_equal(length(b$Vms), 1)
  expect_equal(b$Vms, 107.51)
  expect_is(b$transformed_data, 'data.frame')
  expect_equal(ncol(b$transformed_data), 2)
  expect_equal(nrow(b$transformed_data),7)
})

## ----------------------------------------------- ##
#         Testing the function EH                   #
## ----------------------------------------------- ##
a <- sMM(Km = 2, Vm = 100, error = 'a', plot = F)[,c(1,3:5)]
b <- EH(data = a)

test_that('EH function works properly', {
  expect_is(b, 'list')
  expect_equal(length(b), 4)
  expect_equal(length(b$Kms), 3)
  expect_equal(b$Kms[2], 2.20)
  expect_equal(length(b$Vms), 3)
  expect_equal(b$Vms[2], 102.57)
  expect_is(b$transformed_data, 'data.frame')
  expect_equal(ncol(b$transformed_data), 10)
  expect_equal(nrow(b$transformed_data),7)
})

a <- sMM(Km = 2, Vm = 100, replicates = 1, error = 'r', sd = 0.1)[,c(1,3)]
b <- EH(data = a)

test_that('HW function works properly', {
  expect_is(b, 'list')
  expect_equal(length(b), 4)
  expect_equal(length(b$Kms), 1)
  expect_equal(b$Kms, 2.31)
  expect_equal(length(b$Vms), 1)
  expect_equal(b$Vms, 111.72)
  expect_is(b$transformed_data, 'data.frame')
  expect_equal(ncol(b$transformed_data), 2)
  expect_equal(nrow(b$transformed_data),7)
})

## ----------------------------------------------- ##
#         Testing the function ECB                  #
## ----------------------------------------------- ##
a <- sMM(Km = 1, Vm = 1, replicates = 3, error = 'r', plot = F)[,c(1,3:5)]
b <- ECB(data = a)

test_that('ECB works properly', {
  expect_is(b, 'list')
  expect_equal(length(b), 2)
  expect_true(grepl("Km:  0.953", b[[1]]))
  expect_true(grepl("Vm:  1.348", b[[2]]))
})

a <- sMM(Km = 1, Vm = 1, replicates = 1, error = 'f', plot = F)[, c(1,3)]
b <- ECB(data = a)

test_that('ECB works properly', {
  expect_is(b, 'list')
  expect_equal(length(b), 2)
  expect_true(grepl("Km:  1", b[[1]]))
  expect_true(grepl("Vm:  1.296", b[[2]]))
})
