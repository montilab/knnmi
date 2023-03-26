test_that("mutual_inf_cc_2d returns expected result", {

  set.seed(21)
  x <- runif(10, min=0, max=1)

  M = cbind(x + rnorm(10, 0, 0.1), x + rnorm(10, 0, 0.1),
                   x + rnorm(10, 0, 0.1))


  result <- mutual_inf_cc_2d(x, M, k=3)

  expect_length(result, 3L)
  expect_type(result, "double")
  expect_equal(round(result,5), c(0.73397, 0.72397, 0.79563))

})

test_that("mutual_inf_cc_2d issues error messages when needed", {

  set.seed(21)
  x <- runif(10, min=0, max=1)

  M = cbind(rnorm(9, 0, 1), rnorm(9, 0, 1),
            rnorm(9, 0, 1))


  expect_error(  mutual_inf_cc_2d(x, M, k=3) )
})
