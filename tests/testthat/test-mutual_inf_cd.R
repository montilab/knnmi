test_that("mutual_inf_cd returns expected result", {
  set.seed(654321)
  #' data(mutual_info_df)
  #'
  #' M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  #' mutual_inf_cd(mutual_info_df$Zc_XdYdWd, M)
  #' ## 0.1070804 0.1041177
  
  
  data(mutual_info_df)
  
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  result <- mutual_inf_cd(mutual_info_df$Zc_XdYdWd, M)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.1070804, 0.1041177), tolerance = 0.00001)
  
})

test_that("mutual_inf_cd issues error messages when needed", {
  
  data(mutual_info_df)
  
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  expect_error( mutual_inf_cd(mutual_info_df$Zd_XdYd, M, k=3) )
  
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  expect_error( mutual_inf_cd(mutual_info_df$Zc_XcYc[-1], M, k=3) )
})

test_that("mutual_inf_cd returns expected result", {
  set.seed(654321)
  
  data(mutual_info_df)
  
  result <- mutual_inf_cd(mutual_info_df$Zc_XdYd, mutual_info_df$Xd)
  
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.128029, tolerance = 0.00001)
})


test_that("mutual_inf_cd issues error messages when vector sizes are different", {
  
  data(mutual_info_df)
  
  expect_error(  mutual_inf_cd(mutual_info_df$Zc_XdYd[-1],
                                   mutual_info_df$Xd) )
})


test_that("mutual_inf_cd issues error messages when the value of k is too large", {
  
  data(mutual_info_df)
  
  
  expect_error(  mutual_inf_cd(mutual_info_df$Zc_XdYd,
                                   mutual_info_df$Xd, k=101) )
})
