#'
#' Mutual information estimation when the target and features are continuous.
#'
#' Compute mutual information of \code{target} and \code{features}
#' where \code{target} and \code{features} are both continuous
#' @param target input vector of length N.
#' @param features input vector of length N or a matrix of size NxM.
#' @param k number of nearest neighbors.
#' @useDynLib knnmi _mutual_inf_cc
#'
#' @return a double-precision value - mutual information estimation for
#' vectors \code{target} and \code{features}.
#' @examples
#'
#' data(mutual_info_df)
#' set.seed(654321)
#' mutual_inf_cc(mutual_info_df$Xc, mutual_info_df$Zc_XcYc)
#' ## 0
#'
#' mutual_inf_cc(mutual_info_df$Yc, mutual_info_df$Zc_XcYc)
#' ## 0.2738658
#'
#' @export
mutual_inf_cc <- function(target, features, k=3L){
  if (is.vector(features)) {
    stopifnot( "target and features vectors must have the same length"=length(target) == length(features) )
  } else {
    stopifnot( "number of rows in features matrix should be equal to length of target"=
                 length(target) == nrow(features) )
  }
  stopifnot("k must be less than the length of target"=k < length(target))
  
  res <- .Call('_mutual_inf_cc', target, features, as.integer(k))
  res
}


#'
#' Mutual information estimation when the target is continuous and
#' the features are discrete.
#'
#' Compute mutual information of \code{target} and \code{y}
#' where the \code{target} is continous and \code{features} are discrete.
#' @param target input vector of length N.
#' @param features input vector of length N or a matrix of size NxM.
#' @param k number of nearest neighbors.
#' @useDynLib knnmi _mutual_inf_cd
#'
#' @return a double-precision vector - mutual information estimation for
#' vectors \code{target} and \code{features}.
#' @examples
#'
#' data(mutual_info_df)
#' set.seed(654321)
#' mutual_inf_cd(mutual_info_df$Zc_XdYd, mutual_info_df$Xd)
#' ## 0.128029
#'
#' M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
#' mutual_inf_cd(mutual_info_df$Zc_XdYdWd, M)
#' ## 0.1070804 0.1041177
#' 
#' @export
mutual_inf_cd <- function(target, features, k=3L){
  if (is.vector(features)) {
    stopifnot( "target and features vectors must have the same length"=length(target) == length(features) )
  } else {
    stopifnot( "number of rows in features matrix should be equal to length of target"=
                 length(target) == nrow(features) )
  }
  stopifnot("k must be less than the length of target"=k < length(target))
  if (!is.integer(features)) {
    storage.mode(features) <- "integer"
  }
  res <- .Call('_mutual_inf_cd', target, features, as.integer(k))
  res
}



#'
#' Conditional mutual information estimation for continuous case of a vector and matrix, given matrix.
#'
#' Compute conditional mutual information of vector\code{x}, matrix \code{M}
#' given matrix\code{Z}
#' where \code{x}, \code{M} and \code{Z} are all continuous
#' @param x input vector of size N.
#' @param M input vector of length N or a matrix of size NxM.
#' @param Z conditional input vector of length N or a matrix of size NxM.
#' @param k number of nearest neighbors.
#' @useDynLib knnmi _cond_mutual_inf_ccc 
#'
#' @return a double-precision vector
#'
#' @examples
#' data(mutual_info_df)
#' set.seed(654321)
#' cond_mutual_inf_ccc(mutual_info_df$Zc_XcYc,
#'                        mutual_info_df$Xc, mutual_info_df$Yc)
#' ## 0.2936858
#' 
#' M <- cbind(mutual_info_df$Xc, mutual_info_df$Yc)
#' ZM <- cbind(mutual_info_df$Yc, mutual_info_df$Wc)
#' cond_mutual_inf_ccc(mutual_info_df$Zc_XcYcWc, M, ZM)
#' ## 0.1171533 0.2192397
#'
#' @export
cond_mutual_inf_ccc <- function(x, M, Z, k=3L){
  # TODO:  We want to make sure that x is a vector,
  # and the M & Z are both vectors or are both matrices.
  # Once those types are checked verify correct dimensions.
  # Best/fastest way in R?!?
  
  # Quit if the M & Z types are not the same. 
  if (xor(is.vector(M),is.vector(Z))) {
    stop("M and Z must have the same type - numeric vectors or matrices")
  }
  
  # When they're both vectors make sure their length is the same as x.  
  if (is.vector(M) && is.vector(Z)) {
    stopifnot( "x and M must have the same length"=length(x) == length(M) )
    stopifnot( "x and Z must have the same length"=length(x) == length(Z) )
    stopifnot("k must be less than the length of input vectors"=k < length(x))
  } else { # Both are matrices.
    stopifnot( "M must be a matrix"= (class(M)[1]=="matrix"))
    stopifnot( "Z must be a matrix"= (class(Z)[1]=="matrix"))
    stopifnot( "x and M must have the same length"=length(x) == nrow(M) )
    stopifnot( "x and Z must have the same length"=length(x) == nrow(Z) )
    stopifnot( "M and Z must be the same size"=dim(M) == dim(Z) )
    stopifnot("k must be less than the length of x"=k < length(x))
  }
  res <- .Call('_cond_mutual_inf_ccc', x, M, Z, as.integer(k))
  res
}


#'
#' Conditional mutual information estimation for a continuous vector
#' and a discrete matrix, given another discrete matrix.
#'
#' Compute conditional mutual information of \code{x},\code{M} given \code{Z}
#' where \code{x} is continuous, \code{M} and \code{Z} are discrete
#' @param x input vector of size N.
#' @param M input integer vector of length N or a matrix of size NxM.
#' @param Z conditional input integer vector of length N or a matrix of size NxM.
#' @param k number of nearest neighbors.
#' @useDynLib knnmi _cond_mutual_inf_cdd
#'
#'
#' @return a double-precision vector - mutual information estimation for
#' vector \code{x} and matrix \code{M}, given matrix \code{Z}.
#'
#' @examples
#' data(mutual_info_df)
#' set.seed(654321)
#' cond_mutual_inf_cdd(mutual_info_df$Zc_XdYd, mutual_info_df$Xd,
#'                   mutual_info_df$Yd)
#' ## 0.1338664
#' 
#' M <- cbind(mutual_info_df$Xd, mutual_info_df$Yd)
#' ZM <- cbind(mutual_info_df$Yd, mutual_info_df$Wd)
#' cond_mutual_inf_cdd(mutual_info_df$Zc_XdYdWd, M, ZM)
#' ## 0.1757598 0.1086227
#'
#' @export
cond_mutual_inf_cdd <- function(x, M, Z, k=3L){

  # TODO:  We want to make sure that x is a vector,
  # and the M & Z are both vectors or are both matrices.
  # Once those types are checked verify correct dimensions.
  # Best/fastest way in R?!?
  
  # Quit if the M & Z types are not the same. 
  if (xor(is.vector(M),is.vector(Z))) {
    stop("M and Z must have the same type - numeric vectors or matrices")
  }
  
  # When they're both vectors make sure their length is the same as x.  
  if (is.vector(M) && is.vector(Z)) {
    stopifnot( "x and M must have the same length"=length(x) == length(M) )
    stopifnot( "x and Z must have the same length"=length(x) == length(Z) )
    stopifnot("k must be less than the length of input vectors"=k < length(x))
  } else { # Both are matrices.
    stopifnot( "M must be a matrix"= (class(M)[1]=="matrix"))
    stopifnot( "Z must be a matrix"= (class(Z)[1]=="matrix"))
    stopifnot( "x and M must have the same length"=length(x) == nrow(M) )
    stopifnot( "x and Z must have the same length"=length(x) == nrow(Z) )
    stopifnot( "M and Z must be the same size"=dim(M) == dim(Z) )
    stopifnot("k must be less than the length of x"=k < length(x))
  }

  # TODO: i the R data type of M&Z is numeric, this should
  # just return the call to cond_mutual_inf_ccc. 
  # If they're both integers it should proceed and call the C interface
  # function _cond_mutual_inf_cdd which handles the integer->double
  # conversion on a column-by-column basis.
  if (!is.integer(M)) {
    #M <- matrix(as.integer(M), nrow=nrow(M))
    # storage.mode appears to be 5x faster...
    storage.mode(M) <- "integer"
  }
  if (!is.integer(Z)) {
   # Z <- matrix(as.integer(Z), nrow=nrow(Z))
    storage.mode(Z) <- "integer"
  }

  res <- .Call('_cond_mutual_inf_cdd', x, M, Z, as.integer(k))
  res
}


