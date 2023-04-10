#'
#' Mutual information estimation for continuous case of 2 vectors.
#'
#' Compute mutual information of \code{x} and \code{y}
#' where \code{x} and \code{y} are both continuous
#' @param x input vector.
#' @param y input vector of the same length as x.
#' @param k number of nearest neighbors.
#' @useDynLib knnmi _mutual_inf_cc_1d
#'
#' @return a double-precision value - mutual information estimation for
#' vectors \code{x} and \code{y}.
#' @examples
#'
#' x <- c(0.9065555 , 0.08110139, 0.13248763, 0.01142746, 0.07952208,
#' 0.24536721, 0.86055117, 0.1314828 , 0.26211815, 0.69562545,
#' 0.18115096, 0.69699724, 0.25942256, 0.2728219 , 0.45985588)
#'
#' y <- c(1.75049046,  1.14456434,  0.99217277,  0.96046432,  1.08373406,
#'        1.17513354,  1.63132858,  1.02685152,  0.94713564, -0.31240944,
#'        0.88935441,  1.44502339,  1.02039948,  0.97471144,  1.25480345)
#'
#'
#' mutual_inf_cc_1d(x, y, k=3)
#'
#' @export
mutual_inf_cc_1d <- function(x, y, k=3L){
  stopifnot( "x and y must have the same length"=length(x) == length(y) )
  res <- .Call('_mutual_inf_cc_1d', x, y, as.integer(k))
  res

}



#'
#' Mutual information estimation for continuous case of an NxM matrix.
#'
#' Compute conditional mutual information of \code{x} and matrix \code{M}
#' where \code{x} and \code{M} are both continuous
#' @param x input vector.
#' @param M input matrix. The number of rows should be equal to the
#' length of vector x
#' @param k number of nearest neighbors.
#' @useDynLib knnmi _mutual_inf_cc_2d
#'
#' @return a vector of length equal to the number of columns of matrix M.
#' @examples
#'
#' set.seed(21)
#' x <- runif(10, min=0, max=1)
#'
#' M = cbind(x + rnorm(10, 0, 0.1), x + rnorm(10, 0, 0.1),
#'                  x + rnorm(10, 0, 0.1))
#'
#'
#' mutual_inf_cc_2d(x, M, k=3)
#'
#' @export
mutual_inf_cc_2d <- function(x, M, k=3L){

  stopifnot( "number of rows in matrix M should be equal to length of x"=
               length(x) == nrow(M) )
  res <- .Call('_mutual_inf_cc_2d', x, M, as.integer(k))
  res

}



#'
#' Mutual information estimation for continuous/discrete case.
#'
#' Compute mutual information of \code{x} and \code{y}
#' where \code{x} is continuous vector and \code{y} are discrete
#' @param x input (continuous) vector.
#' @param y input (discrete) vector. It should have the same length as x.
#' @param k number of nearest neighbors.
#' @param use_cc (logical) if TRUE the algorithm falls into continuous/continuous case
#' @useDynLib knnmi _mutual_inf_cd_1d
#'
#' @return a double-precision value - mutual information estimation for
#' vectors \code{x} and \code{y}.
#' @export
mutual_inf_cd_1d <- function(x, y, k=3L, use_cc=FALSE){
  stopifnot( "x and y must have the same length"=length(x) == length(y) )
  if (!is.integer(y)) {
    y <- as.integer(y)
  }
  res <- .Call('_mutual_inf_cd_1d', x, as.integer(y),
               as.integer(k), as.logical(use_cc))
  res
}


#'
#' Mutual information estimation for continuous/discrete case.
#'
#' Compute mutual information of \code{x} and \code{y}
#' where \code{x} is continuous vector and \code{y} are discrete
#' @param x input (continuous) vector.
#' @param M input (discrete) matrix. It should have the same number of rows as \code{length(x)}.
#' @param k number of nearest neighbors.
#' @param use_cc (logical) if TRUE the algorithm falls into continuous/continuous case
#' @useDynLib knnmi _mutual_inf_cd_2d
#'
#' @return vector of length m, where m is the number of columns in M
#' vectors \code{x} and \code{y}.
#' @export
mutual_inf_cd_2d <- function(x, M, k=3L, use_cc=FALSE){
  stopifnot( "x and M must have the same length"=length(x) == nrow(M) )
  if (!is.integer(M)) {
    M <- matrix(as.integer(M), nrow=nrow(M))
  }
  res <- .Call('_mutual_inf_cd_2d', x, M,
               as.integer(k), as.logical(use_cc))
  res
}


#'
#' Conditional mutual information estimation for continuous case of 3 vectors.
#'
#' Compute conditional mutual information of \code{x},\code{y} given \code{z}
#' where \code{x}, \code{y} and \code{z} are all continuous
#' @param x input vector.
#' @param y input vector of the same length as x.
#' @param z conditional input vector of the same length as x.
#' @param k number of nearest neighbors.
#' @useDynLib knnmi _cond_mutual_inf_ccc_1d
#'
#' @return a double-precision value - mutual information estimation for
#' vectors \code{x} and \code{y}.
#'
#' @export
cond_mutual_inf_ccc_1d <- function(x, y, z, k=3L){

  stopifnot( "x and y must have the same length"=length(x) == length(y) )
  res <- .Call('_cond_mutual_inf_ccc_1d', x, y, z, as.integer(k))
  res

}


#'
#' Conditional mutual information estimation for continuous/discrete case of 3 vectors.
#'
#' Compute conditional mutual information of \code{x},\code{y} given \code{z}
#' where \code{x} is continuous, \code{y} and \code{z} are discrete
#' @param x input vector.
#' @param y input vector of the same length as x.
#' @param z conditional input vector of the same length as x.
#' @param k number of nearest neighbors.
#' @useDynLib knnmi _cond_mutual_inf_cdd_1d
#'
#' @return a double-precision value - mutual information estimation for
#' vectors \code{x} and \code{y}, given \code{z}.
#'
#' @export
cond_mutual_inf_cdd_1d <- function(x, y, z, k=3L){

  stopifnot( "x and y must have the same length"=length(x) == length(y) )
  res <- .Call('_cond_mutual_inf_cdd_1d', x, y, z, as.integer(k))
  res

}



#'
#' Conditional mutual information estimation for continuous case of a vector and matrix, given matrix.
#'
#' Compute conditional mutual information of vector\code{x}, matrix \code{M}
#' given matrix\code{Z}
#' where \code{x}, \code{M} and \code{Z} are all continuous
#' @param x input vector.
#' @param M input matrix of the same length as x.
#' @param Z conditional input matrix of the same length as x.
#' @param k number of nearest neighbors.
#' @useDynLib knnmi _cond_mutual_inf_ccc_2d
#'
#' @return a double-precision vector
#'
#' @export
cond_mutual_inf_ccc_2d <- function(x, M, Z, k=3L){

  res <- .Call('_cond_mutual_inf_ccc_2d', x, M, Z, as.integer(k))
  res

}


#'
#' Conditional mutual information estimation for a continuous vector
#' and a discrete matrix, given another discrete matrix.
#'
#' Compute conditional mutual information of \code{x},\code{M} given \code{Z}
#' where \code{x} is continuous, \code{M} and \code{Z} are discrete
#' @param x input vector.
#' @param M input matrix of the same length as x.
#' @param Z conditional input matrix of the same length as x.
#' @param k number of nearest neighbors.
#' @useDynLib knnmi _cond_mutual_inf_cdd_2d
#'
#'
#' @return a double-precision vector - mutual information estimation for
#' vector \code{x} and matrix \code{M}, given matrix \code{Z}.
#'
#' @export
cond_mutual_inf_cdd_2d <- function(x, M, Z, k=3L){
  stopifnot( "x and M must have the same length"=length(x) == length(M) )
  if (!is.integer(M)) {
    M <- matrix(as.integer(M), nrow=nrow(M))
  }
  if (!is.integer(Z)) {
    Z <- matrix(as.integer(Z), nrow=nrow(Z))
  }
  res <- .Call('_cond_mutual_inf_cdd_2d', x, M, Z, as.integer(k))
  res
}


