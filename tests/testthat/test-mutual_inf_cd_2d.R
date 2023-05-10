test_that("mutual_inf_cd_2d returns expected result", {

  data(mutual_info_df)

  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  result <- mutual_inf_cd_2d(mutual_info_df$Zc_XcYc, M, k=3)

  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(round(result,5), c(0., 0.))

})

test_that("mutual_inf_cd_2d issues error messages when needed", {

  data(mutual_info_df)

  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  expect_error( mutual_inf_cd_2d(mutual_info_df$Zd_XdYd, M, k=3) )

  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  expect_error( mutual_inf_cd_2d(mutual_info_df$Zc_XcYc[-1], M, k=3) )
})
