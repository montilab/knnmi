#' @title 
#' Mutual information estimation 
#'
#' @description
#' Estimate the mutual information MI(X;Y) of the target \code{X} and features \code{Y}  
#' where \code{X} and \code{Y} are both continuous using k-nearest neighbor distances.
#'
#' @references
#' Alexander Kraskov, Harald Stögbauer, and Peter Grassberger. Phys. Rev. E **69**, 066138 (2004).
#' \doi{10.1103/PhysRevE.69.066138}
#'
#' @param target input vector.
#' @param features input vector or matrix. 
#' @param k Integer number of nearest neighbors. The default value is 3.
#' 
#' @details 
#' The features argument is a vector of the same size as the target vector, or a matrix whose column dimension matches
#' the size of the target vector. 
#' 
#' @useDynLib knnmi _mutual_inf_cc
#'
#' @return Returns the estimated mutual information. The return value is a vector of size 1 if the features argument 
#' is a vector. If the features argument is a matrix then the return value is a vector whose size matches the 
#' number of rows in the matrix.
#' 
#' @examples
#'
#' data(mutual_info_df)
#' set.seed(654321)
#' mutual_inf_cc(mutual_info_df$Yc, t(mutual_info_df$Zc_XcYc))
#' mutual_inf_cc(mutual_info_df$Xc, t(mutual_info_df$Zc_XcYc), k=5)
#' 
#' 
#'
#' @export
mutual_inf_cc <- function(target, features, k=3L){
  if (is.vector(features)) {
    stopifnot( "target and features vectors must have the same length"=length(target) == length(features) )
  } else {
    stopifnot( "number of columns in features matrix should be equal to length of the target vector"=
                 length(target) == ncol(features) )
  }
  stopifnot("k must be less than the length of target"=k < length(target))
  
  # Make sure the inputs are double precision floats
  if (!is.double(target)) {
    storage.mode(target) <- "double"
  }
  if (!is.double(features)) {
    storage.mode(features) <- "double"
  }
  
  res <- .Call('_mutual_inf_cc', target, features, as.integer(k))
  res
}


#' @title 
#' Mutual information estimation 
#'
#' @description
#' Estimate the mutual information MI(X;Y) of the target \code{X} and features \code{Y}  
#' where \code{X} is continuous or discrete and \code{Y} is discrete using k-nearest neighbor distances.
#'
#' @references
#' Ross BC (2014) Mutual Information between Discrete and Continuous Data Sets. PLoS ONE 9(2): e87357. 
#' \doi{10.1371/journal.pone.0087357}
#' 
#' @param target input vector.
#' @param features input vector or matrix. 
#' @param k Integer number of nearest neighbors. The default value is 3.
#' 
#' @details
#' The features argument is a vector of the same size as the target vector, or a matrix whose column dimension matches
#' the size of the target vector. Discrete values for the features or targets must be numeric or integer types.
#' 
#' @useDynLib knnmi _mutual_inf_cd
#'
#' @return Returns the estimated mutual information. The return value is a vector of size 1 if the features argument is 
#' a vector. If the features argument is a matrix then the return value is a vector whose size matches the 
#' number of rows in the matrix.
#'
#' @examples
#'
#' data(mutual_info_df)
#' set.seed(654321)
#' mutual_inf_cd(mutual_info_df$Zc_XdYd, t(mutual_info_df$Xd))
#' 
#' M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
#' mutual_inf_cd(mutual_info_df$Zc_XdYdWd, t(M))
#' 
#' 
#' @export
mutual_inf_cd <- function(target, features, k=3L){
  if (is.vector(features)) {
    stopifnot( "target and features vectors must have the same length"=length(target) == length(features) )
  } else {
    stopifnot( "number of columns in features matrix should be equal to length of target"=
                 length(target) == ncol(features) )
  }
  stopifnot("k must be less than the length of target"=k < length(target))
  
  # Make sure the target is double precision floats
  # This will convert integers to doubles if needed.
  if (!is.double(target)) {
    storage.mode(target) <- "double"
  }
  # Make sure features is an integer array - a numeric array containing
  # integer values isn't quite sufficient, they need to be actual
  # integers so the C/C++ code works properly.
  if (!is.integer(features)) {
    storage.mode(features) <- "integer"
  }
  res <- .Call('_mutual_inf_cd', target, features, as.integer(k))
  res
}


#' @title Conditional mutual information estimation
#' 
#' @description
#' Conditional mutual information estimation CMI(X;Y|Z) where X is a continuous vector.
#' The input Y and conditional input Z can be vectors or matrices. If Y and Z
#' are discrete then they must be numeric or integer valued.
#'
#' @references
#' Alkiviadis Tsimpiris, Ioannis Vlachos, Dimitris Kugiumtzis,
#' Nearest neighbor estimate of conditional mutual information in feature selection,
#' Expert Systems with Applications,
#' Volume 39, Issue 16, 2012, Pages 12697-12708
#' \doi{10.1016/j.eswa.2012.05.014}
#' 
#' @param X input vector.
#' @param Y input vector or matrix.
#' @param Z conditional input vector or matrix.
#' @param k number of nearest neighbors.
#' 
#' 
#' @details
#' Argument Y is a vector of the same size as vector X, or a matrix whose column dimension matches
#' the size of X.  Argument Z is also a vector of the same size as vector X, or a matrix whose column dimension matches
#' the size of X. If Y and Z are both matrices they must additionally have the same number of rows.  If Y and/or Z are
#' discrete values they must have a numeric or integer type. 
#' 
#' @useDynLib knnmi _cond_mutual_inf
#'
#' @return Returns the estimated conditional mutual information. The return value is a vector of size 1 if both Y and Z 
#' are vectors. If either Y or Z are matrices the return value is a vector whose size is the number of rows in the 
#' matrix.
#'
#' @examples
#' data(mutual_info_df)
#' set.seed(654321)
#' cond_mutual_inf(mutual_info_df$Zc_XcYc,
#'                        mutual_info_df$Xc, t(mutual_info_df$Yc))
#' 
#' M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
#' ZM <- cbind(mutual_info_df$Yc, mutual_info_df$Wc)
#' cond_mutual_inf(mutual_info_df$Zc_XcYcWc, t(M), t(ZM))
#'
#'
#' @export
cond_mutual_inf <- function(X, Y, Z, k=3L){
  # conditional mutual information:  CMI(X, Y|Z)
  # X: vector of length N
  # Y: vector of length N or matrix of size MxN
  # Z: vector of length N or matrix of size MxN
  
  # Check that sizes match before continuing. 
  stopifnot("k must be less than the length of the input vector X"=k < length(X))
  
  # The "case" value is used in the C code to let it pick the right code
  # path without having to re-check for vector vs matrix in C.
  case <- 0L
  if (is.vector(Y) && is.vector(Z)) {
    # When they're both vectors make sure their length is the same as X.  
    stopifnot( "X and Y must have the same length"=length(X) == length(Y) )
    stopifnot( "X and Z must have the same length"=length(X) == length(Z) )
  } else if (is.vector(Y) && is.matrix(Z)) {
    case <- 1L
    # mixed vector & matrix
    stopifnot( "X and Y must have the same length"=length(X) == length(Y) )
    stopifnot( "Number of Z columns must be the same as the length of X"=length(X) == ncol(Z) )
  } else if (is.vector(Z) && is.matrix(Y)) {
    # mixed vector & matrix
    case <- 2L
    stopifnot( "X and Z must have the same length"=length(X) == length(Z) )
    stopifnot( "Number of Y columns must be the same as the length of x"=length(X) == ncol(Y) )    
  }  else if (is.matrix(Z) && is.matrix(Y)) {
    # Both Y and Z are matrices.
    case <- 3L
    stopifnot( "X and Y must have the same length"=length(X) == ncol(Y) )
    stopifnot( "X and Z must have the same length"=length(X) == ncol(Z) )
    stopifnot( "Y and Z must be the same size"=dim(Y) == dim(Z) )
  } else {
    # Incorrect arguments...
    stop("Y and Z must be vectors or matrices of the correct dimensionality.")
  }
  
  if (!is.double(X)) {
    storage.mode(X) <- "double"
  }
  if (!is.double(Y)) {
    storage.mode(Y) <- "double"
  }
  if (!is.double(Z)) {
    storage.mode(Z) <- "double"
  }
  
  res <- .Call('_cond_mutual_inf', X, Y, Z, as.integer(k), case )
  res
}


