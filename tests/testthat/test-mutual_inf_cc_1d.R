test_that("mutual_inf_cc_1d returns expected result", {

  x <- c(0.91 , 0.08, 0.13, 0.01, 0.08,
  0.25, 0.86, 0.13 , 0.26, 0.7, 0.18, 0.7, 0.26, 0.27 , 0.46)

  y <- c(1.75,  1.14,  0.99,  0.96,  1.08,
         1.18,  1.63,  1.03,  0.95, -0.31, 0.89,  1.45,  1.02,  0.97,  1.25)


  result <- mutual_inf_cc_1d(x, y, k=3)
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(round(result,5), 0.35588)
})


test_that("mutual_inf_cc_1d issues error messages when needed", {

  x <- c(0.91 , 0.08, 0.13, 0.01, 0.08,
         0.25, 0.86, 0.13 , 0.26, 0.7, 0.18, 0.7, 0.26, 0.27 , 0.46)

  y <- c(1.75,  1.14,  0.99,  0.96,  1.08)


  expect_error(  mutual_inf_cc_1d(x, y, k=3) )
})
