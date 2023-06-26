#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rinternals.h>
#include <Rmath.h>

extern int mutual_inf_cc_vec( double *input_x,  double *input_y,  int n_elems, int k, double *mi) ;
extern int mutual_inf_cd_vec( double *input_x,  int *input_y,  int n_elems,  int k, double *mi ) ;
extern int cond_mutual_inf_vec( double *input_x,   double *input_y,  double *input_z,  int n_elems,  int k, double *mi) ;

/* Functions defined here:
 *  _mutual_inf_cc - MI(x;y) c-c case where y is an NxM matrix. Returns an MI vector of length M.
 *  _mutual_inf_cd  - mutual information MI(x;y) c-d case where y is an NxM matrix.
 *  _cond_mutual_inf_ccc  - CMI(x;y|z), where x,y,z are continuous values. x is of length N, y and z are vectors
 *    of size N or matrices of size NxM.
 *  _cond_mutual_inf_cdd  - CMI(x;y|z), where x is continuous and both y and z are discrete. x is of length N,
 *      y and z are integer vectors of size N or integer matrices of size NxM.
 *
 *  The 1d functions all return a single numeric value. The 2d functions return a numeric vector of length M.
 *
 *  The R code that calls these functions will implement checks to make sure that vector/matrix sizes are
 *  correct.
 *
 */


SEXP _mutual_inf_cc( SEXP r_input_x, SEXP r_input_y, SEXP k) {
  /* R C wrapper for:
   * int mutual_inf_cc_vec( double *input_x,  double *input_y,  int n_elems, int k, double *mi)
   * r_input_x is the target vector. It can be an R vector of length N, or matrix of size Nx1.
   * r_input_y is the features matrix. It can be an R vector or length N or matrix of size NxM.
   * Assume that the size double-checking happened in R.
   * k is the number of nearest neighbors.
   * The mutual information is returned as a numeric vector of length 1.
   */
  
  /* declare the output mutual information */
  SEXP mi;
  double *p_mi, *p_y, *p_x ;
  int n_rows, n_cols,   k_value, i  ;
  
  if (!isMatrix(r_input_x)) {
    /* Assume it's a vector */
    n_rows = LENGTH(r_input_x);
  } else {
    /* It's a matrix, get the number of rows */
    n_rows = Rf_nrows(r_input_x);
  }
  /* In the R wrapper make sure this is passed as an integer. */
  k_value = INTEGER(k)[0] ;
  
  /* Check if r_input_y is a vector or a matrix. */
  n_cols = 1 ;
  if (isMatrix(r_input_y)) {
    /* It's a matrix, get the number of columns */
    n_cols = Rf_ncols(r_input_y);
  } 

  /* R memory allocation */
  mi = PROTECT(allocVector(REALSXP, n_cols));
  
  /* get pointers */
  p_mi = REAL(mi);
  p_y = REAL(r_input_y) ;
  p_x = REAL(r_input_x) ;
  for (i = 0 ; i < n_cols ; i++) {
   mutual_inf_cc_vec(p_x, &p_y[n_rows * i], n_rows, k_value, &p_mi[i]);
  }
  /* Lift R garbage protection */
  UNPROTECT(1);  
  /* Return the mi vector */
  return( mi );
}

SEXP _mutual_inf_cd( SEXP r_input_x, SEXP r_input_y, SEXP k) {
  /* R C wrapper for:
   * int mutual_inf_cd(double *input_x, int x_elems, int *input_y, int y_nrows, int y_ncols, int k, double *mi) ;
   *  r_input_x - the vector of continuous data. size N.
   *  r_input_y - input matrix. size NxM
   *  k - number of nearest neighbors.
   *  The mutual information is returned as a vector of length M.
   */
  
  /* declare the output mutual information */
  SEXP mi;
  double *p_mi;
  int n_rows, n_cols, k_value ;
  int *p_y ;
  double *p_x ;
  int prot_ctr = 0 ;
  int i ;
  
  /* In the R wrapper make sure this is passed as an integer. */
  k_value = INTEGER(k)[0] ;
  
  /* Get the dimensions of r_input_y */
  n_rows = Rf_nrows(r_input_y);
  /* Check if r_input_y is a vector or a matrix. */
  n_cols = 1 ;
  if (isMatrix(r_input_y)) {
    /* It's a matrix, get the number of columns */
    n_cols = Rf_ncols(r_input_y);
  }
  /* R memory allocation */
  mi = PROTECT(allocVector(REALSXP, n_cols));
  prot_ctr++ ;
  
  p_mi = REAL(mi);
  p_x = REAL(r_input_x);
  p_y = INTEGER(r_input_y) ;
  
  /* As with the 1D case try using the c-d algorithm unless a negative MI is returned.
   * If it is...re-do the entire calculation with the c-c algorithm, doing type conversions
   * of the integer columns as we go. */
  for (i = 0 ; i < n_cols ; i++) {
    mutual_inf_cd_vec(p_x, &p_y[n_rows * i], n_rows, k_value, &p_mi[i]);
  }
  
  /* Lift R garbage protection */
  UNPROTECT(prot_ctr);
  
  return( mi );
}

SEXP _cond_mutual_inf_ccc( SEXP r_input_x, SEXP r_input_y, SEXP r_input_z, SEXP k) {
  /* R C wrapper for:
   *  int cond_mutual_inf( double *input_x,  int x_elems,  double *input_y,  double *input_z,  int nrows,  int ncols,  int k, double *mi) ;
   *  r_input_x - the vector of continuous target data. size N.
   *  r_input_y - input vector or matrix. size N or size NxM
   *  r_input_z - input vector or matrix. size N or size NxM
   *  k - number of nearest neighbors.
   *
   *  The mutual information is returned as a vector of length M.
   *
   *  in CaDrA  terms  input_score = x, input_mat=y, expression_score=z
   */
  
  /* declare the output mutual information */
  SEXP mi;
  double *p_mi;
  int  yz_nrows, yz_ncols, k_value ;
  int j ;
  double *p_x, *p_y, *p_z ;
  
  /* In the R wrapper make sure this is passed as an integer. */
  k_value = INTEGER(k)[0] ;
  
  /* Get the dimensions of r_input_y */
  /* r_input_z must have the same dimensions, that should all be checked
   * on the R side. */
  yz_nrows = Rf_nrows(r_input_y);
  /* Check if r_input_y is a vector or a matrix. */
  yz_ncols = 1 ;
  if (isMatrix(r_input_y)) {
    /* It's a matrix, get the number of columns */
    yz_ncols = Rf_ncols(r_input_y);
  }
  
  
  /* R memory allocation */
  mi = PROTECT(allocVector(REALSXP, yz_ncols));
  p_mi = REAL(mi);
  
  p_x = REAL(r_input_x) ;
  p_y = REAL(r_input_y) ;
  p_z = REAL(r_input_z) ;
  for (j = 0 ; j < yz_ncols ; j++) {
    cond_mutual_inf_vec(p_x, &p_y[j * yz_nrows], &p_z[j * yz_nrows], yz_nrows, k_value, &p_mi[j]);
  }
  
  /* Lift R garbage protection */
  UNPROTECT(1);
  
  return( mi );
}



SEXP _cond_mutual_inf_cdd( SEXP r_input_x, SEXP r_input_y, SEXP r_input_z, SEXP k) {
  /* R C wrapper for:
   * int cond_mutual_inf( double *input_x,  int x_elems,  double *input_y,  double *input_z,  int nrows,  int ncols,  int k, double *mi) ;
   *  r_input_x - the vector of continuous data. size N.
   *  r_input_y - integer vector or matrix. size N or size NxM
   *  r_input_z - integer vector or matrix. size N or size NxM
   *  k - number of nearest neighbors.
   *
   *  The mutual information is returned as a vector of length M.
   *
   *  in CaDrA  terms  input_score = x, input_mat=y, expression_score=z
   */
  
  /* declare the output mutual information */
  SEXP mi;
  double *p_mi;
  int yz_nrows, yz_ncols, k_value ;
  double *p_dbl_y, *p_dbl_z ;
  SEXP dbl_y, dbl_z ;
  int i, j ;
  int *p_y, *p_z ;
  double *p_x ;
  
  /* In the R wrapper make sure this is passed as an integer. */
  k_value = INTEGER(k)[0] ;
  
  /* Get the dimensions of r_input_y */
  /* r_input_z must have the same dimensions, that should all be checked
   * on the R side. */
  yz_nrows = Rf_nrows(r_input_y);
  /* Check if r_input_y is a vector or a matrix. */
  yz_ncols = 1 ;
  if (isMatrix(r_input_y)) {
    /* It's a matrix, get the number of columns */
    yz_ncols = Rf_ncols(r_input_y);
  }
  
  /* R memory allocation */
  mi = PROTECT(allocVector(REALSXP, yz_ncols));
  p_mi = REAL(mi);
  
  /* To avoid excess memory usage, loop here and call cond_mutual_inf_vec while
   * doing integer -> double conversions as the conditional mutual inf. function
   * is only defined for all-double precision inputs.*/
  dbl_y = PROTECT(allocVector(REALSXP, yz_nrows));
  dbl_z = PROTECT(allocVector(REALSXP, yz_nrows));
  p_dbl_y = REAL(dbl_y) ;
  p_dbl_z = REAL(dbl_z) ;
  p_y = INTEGER(r_input_y) ;
  p_z = INTEGER(r_input_z) ;
  p_x = REAL(r_input_x) ;
  for (j = 0 ; j < yz_ncols ; j++) {
    /* Convert the jth column of the input y & z to double */
    for (i = 0 ; i < yz_nrows ; i++) {
      p_dbl_y[i] = (double) p_y[i + j * yz_nrows] ;
      p_dbl_z[i] = (double) p_z[i + j * yz_nrows] ;
    }
    cond_mutual_inf_vec(p_x, p_dbl_y, p_dbl_z, yz_nrows, k_value, &p_mi[j]);
  }
  
  /* Lift R garbage protection */
  UNPROTECT(3);
  
  return( mi );
} 

