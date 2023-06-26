test_that("cond_mutual_inf_cdd returns expected result", {
  set.seed(654321)
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  ZM <- cbind(mutual_info_df$Yd, mutual_info_df$Wd)
  result <- cond_mutual_inf_cdd(mutual_info_df$Zc_XdYdWd, M=M, Z=ZM)
  
  expect_length(result, 2L)
  expect_type(result, "double")
  expect_equal(result, c(0.1757598, 0.1086227), tolerance = 0.000001)
  
})

test_that("cond_mutual_inf_cdd issues error messages when vector and matrix have different sizes", {
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  ZM <- cbind(mutual_info_df$Yd, mutual_info_df$Wd)
  expect_error( cond_mutual_inf_cdd(mutual_info_df$Zc_XdYdWd[-1],
                                    M,
                                    ZM))
  expect_error( cond_mutual_inf_ccc(mutual_info_df$Zc_XdYdWd,
                                    M, ZM[-1,]))
  expect_error( cond_mutual_inf_ccc(mutual_info_df$Zc_XdYdWd,
                                    M, ZM[, 1]))
  expect_error( cond_mutual_inf_ccc(mutual_info_df$Zc_XdYdWd,
                                    M[-1,], ZM))
  expect_error( cond_mutual_inf_ccc(mutual_info_df$Zc_XdYdWd,
                                    M[, 1], ZM))
})


test_that("cond_mutual_inf_cdd issues error messages when the value of k is too large", {
  
  data(mutual_info_df)
  M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
  ZM <- cbind(mutual_info_df$Yd, mutual_info_df$Wd)
  expect_error( cond_mutual_inf_cdd(mutual_info_df$Zc_XdYdWd[-1],
                                    M,
                                    ZM, k=101))
})

test_that("cond_mutual_inf_cdd returns expected result", {
  set.seed(654321)
  
  data(mutual_info_df)
  
  result <- cond_mutual_inf_cdd(mutual_info_df$Zc_XdYd, mutual_info_df$Xd,
                                mutual_info_df$Yd)
  
  expect_length(result, 1L)
  expect_type(result, "double")
  expect_equal(result, 0.1338664, tolerance = 0.00001)
})


test_that("cond_mutual_inf_cdd issues error messages when vector sizes are different", {
  
  data(mutual_info_df)
  
  expect_error( cond_mutual_inf_cdd(mutual_info_df$Zc_XdYd[-1],
                                    mutual_info_df$Xd,
                                    mutual_info_df$Yd))
  expect_error( cond_mutual_inf_cdd(mutual_info_df$Zc_XdYd,
                                    mutual_info_df$Xd[-1],
                                    mutual_info_df$Yd))
  expect_error( cond_mutual_inf_cdd(mutual_info_df$Zc_XdYd,
                                    mutual_info_df$Xd,
                                    mutual_info_df$Yd[-1]))
})

 
