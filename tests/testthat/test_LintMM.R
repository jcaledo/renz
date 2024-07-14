library(renz)
context("Linearized integrated MM equation")


## ----------------------------------------------- ##
#       Testing the function int.MM()               #
## ----------------------------------------------- ##
test_that('int.MM() works properly', {
  a <-int.MM(data = sE.progress(So = 10, time = 5,
                           Km = 4, Vm = 50)[, c(1,3)])
  expect_is(a, 'list')
  expect_equal(length(a), 2)
  expect_equivalent(a$parameters[1], 3.686)
  expect_equivalent(a$parameters[2], 48.303)
  expect_is(a$data, 'data.frame')
  expect_equal(dim(a$data), c(15,4))
})
