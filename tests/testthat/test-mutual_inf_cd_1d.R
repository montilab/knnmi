test_that("mutual_inf_cd_1d returns expected result", {


  x <- c(0.91 , 0.08, 0.13, 0.01, 0.08,
         0.25, 0.86, 0.13 , 0.26, 0.7)

  y <- c(1, 0, 0, 0, 0, 1, 0, 0, 0, 1)


  result <- mutual_inf_cd_1d(x, y, k=3)
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(round(result,5), 0.36385)
})


test_that("mutual_inf_cd_1d issues error messages when needed", {

  x <- c(0.91 , 0.08, 0.13, 0.01, 0.08,
         0.25, 0.86, 0.13 , 0.26, 0.7)

  y <- c(1, 0, 0, 0, 0, 1, 0, 0, 0, 1, 2)

  expect_error(  mutual_inf_cd_1d(x, y, k=3) )
})
