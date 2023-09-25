# Introduction

knnmi is a mutual information library for R based on the k-nearest neighbor algorithm. There are three functions:
* mutual_inf_cc - Compute MI(X;Y) where X is a vector and Y is a vector or matrix. X and Y are both continuous.
* mutual_inf_cd  - Compute MI(X;Y) where X is a vector and Y is a vector or matrix. X is continuous and Y contains discrete numeric values.
* cond_mutual_inf - Compute the condtional mutual information CMI(X;Y|Z) where X is a continous vector. Y and Z are continous and be any mix of vectors and matrices.

# Implementation
The k-nearest neighbor (KNN) approach to mutual information (MI) is used in popular libraries like 
[scikit-learn](https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.mutual_info_regression.html#sklearn.feature_selection.mutual_info_regression) 
for Python. 
Here, the KNN computations are handled by the high performance [nanoflann](https://github.com/jlblancoc/nanoflann) C++ library with the addition of a custom Chebyshev distance metric. The Eigen C++ library is used to handle vectors and an effort was made to minimize the number of copies of the incoming data from R to preserve memory. 
A newer Eigen version (3.4.0) was used for some useful features and a minimal copy of that library is part of this repository. The newest Eigen version built 
into RcppEigen is 3.3.4. 

The publications used as the basis for the knnmi library are:
* mutual_inf_cc - [A. Kraskov, H. Stögbauer, and P. Grassberger, Phys. Rev. E 69, 066138 (2004)](https://doi.org/10.1103/PhysRevE.69.066138)
* mutual_inf_cd - [Ross BC (2014) Mutual Information between Discrete and Continuous Data Sets. PLoS ONE 9(2): e87357](https://doi.org/10.1371/journal.pone.0087357)
* cond_mutual_inf - [A. Tsimpiris, I. Vlachos, and D. Kugiumtzis. Expert Systems with Applications, Volume 39, Issue 16, (2012)](https://doi.org/10.1016/j.eswa.2012.05.014)

# Calling knnmi functions 
`mutual_inf_cc <- function(target, features, k=3L)` The target and features inputs should both be continuous. If the length of the target vector is (N) the dimensions of the features must be N or if it's a matrix MxN. The size of the neighborhood (k) can be changed but should probably be in the range of 2-5 (for more on this range see Kraskov, et.al. 2004).  The target and features vectors are scaled by their standard deviation and a small amount of random noise is added (with a magnitude of <1e-12) if they are continuous values.  If they are not actually continuous but are actually discrete numeric values the scaling and added noise are skipped. Discrete features should be used with the `mutual_inf_cd` function as it is approximately 2x faster. If the features variable is a matrix of size MxN the code will loop over the rows and the output MI will be a vector of size M.

`mutual_inf_cd <- function(target, features, k=3L)` The target should both be continuous and the features should be discrete numeric values. If the length of the target vector is (N) the dimensions of the features must be N or if it's a matrix MxN. The neighborhood size k is the same as in `mutual_inf_cc`. The target vector is scaled as in `mutual_inf_cc`, unless it actually contains discrete numeric values. The features input is assumed to contain discrete numeric values and is converted to an integer type. If the features variable is a matrix of size MxN the code will loop over the rows and the output MI will be a vector of size M.

`cond_mutual_inf <- function(X, Y, Z, k=3L)` X, Y, and Z should be continous values. If the length of the X vector is (N) the dimensions of Y and Z must be vectors of length N or matrices of size MxN. The inputs are treated as they are in `mutual_inf_cc`. If the Y is a matrix of size MxN and Z is a vector of size N then the return value will be a vector of length M. It is the same for the case where Y is a vector and Z is a matrix. For the case where both Y and Z are matrices of size MxN the return value will be a vector of size M where each row of Y and Z are used together to compute the CMI. 

The random noise that is added to continous inputs in all three functions is computed by calling the R random number generator. For repeatability, a seed value should be set in R before calling the knnmi functions:
```
library(knnmi)
set.seed(654321)
mi <- mutual_inf_cc(X,Y)
```
There is no parallelism implemented in this library. Parallel computations can be implemented in R code calling the knnmi library functions. 

# TODO
Add some examples, and some comparison calculations (including timing info as this is faster) with scikit-learn. 
