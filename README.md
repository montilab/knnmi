# Installation

``` r
install.packages("knnmi")
# install the latest development version
devtools::install_github("montilab/knnmi")
```


# Introduction

knnmi is a mutual information (MI) library for R based on the k-nearest neighbor (KNN) algorithm. There are three functions:

* `mutual_inf_cc` - Compute MI(X;Y) where X is a vector and Y is a vector or matrix. X and Y are both continuous.
* `mutual_inf_cd`  - Compute MI(X;Y) where X is a vector and Y is a vector or matrix. X is continuous and Y contains discrete numeric values.
* `cond_mutual_inf` - Compute the conditional mutual information CMI(X;Y|Z) where X is a continuous vector. Y and Z are continuous or discrete and be any combination of vectors and matrices.

The KNN approach to mutual information has been compared to approaches based on kernel density estimation (KDE). The KNN approach was determined to perform better for data sets of 1000 sample points or greater at all noise levels, and for small data sets (≤ 100 points) at low noise levels. The KDE MI algorithm was the better choice for small data sets with high noise levels. For more details, please refer to this paper: [S. Khan, S. Bandyopadhyay, A. R. Ganguly, S. Saigal, D. J. Erickson, III, V. Protopopescu, and G. Ostrouchov
Phys. Rev. E 76, 026209 (2007)](https://link.aps.org/doi/10.1103/PhysRevE.76.026209)


# Implementation
The k-nearest neighbor (KNN) approach to mutual information (MI) is a well estabished technique, and this library uses a C++ implementation for performance. The KNN computations are handled by the high performance [nanoflann](https://github.com/jlblancoc/nanoflann) C++ library with the addition of a custom Chebyshev distance metric. The [Eigen C++ library](https://eigen.tuxfamily.org/index.php?title=Main_Page) is used to handle vectors and an effort was made to minimize the number of copies of the incoming data from R to preserve memory. 

A newer Eigen version (3.4.0) was used for some useful features and a minimal copy of that library is part of this repository. This approach was taken as these features are not available in the 3.3.4 version of Eigen which is built into the  version of RcppEigen that was available when this library was developed. 

The publications used as the basis for the knnmi library are:
* `mutual_inf_cc` - [A. Kraskov, H. Stögbauer, and P. Grassberger, Phys. Rev. E 69, 066138 (2004)](https://doi.org/10.1103/PhysRevE.69.066138)
* `mutual_inf_cd` - [Ross BC (2014) Mutual Information between Discrete and Continuous Data Sets. PLoS ONE 9(2): e87357](https://doi.org/10.1371/journal.pone.0087357)
* `cond_mutual_inf` - [A. Tsimpiris, I. Vlachos, and D. Kugiumtzis. Expert Systems with Applications, Volume 39, Issue 16, (2012)](https://doi.org/10.1016/j.eswa.2012.05.014)

The knnmi code matches the results of the sample code provided with those publications. The algorithms for `mutual_inf_cc` and `mutual_inf_cd` are also implemented in the popular [`scikit-learn`](https://scikit-learn.org/) library for the Python language. There are small deviations between the results calculated by this library and the scikit-learn implementation, due to a larger magnitude of random noise added in the scikit-learn functions.  

# Calling knnmi functions 
### mutual_inf_cc

`mutual_inf_cc <- function(target, features, k=3L)` The target and features inputs should both be continuous. If the length of the target vector is N the the features argument must be a vector of length N or a matrix of size MxN. 

The size of the neighborhood (k) can be changed but should probably be in the range of 2-5 (Kraskov, et.al. 2004). The target and features vectors are scaled by their standard deviation and a small amount of uniform random noise is added (with a magnitude of 1e-12) if they are continuous values. This functions treats the features input as continuous but can detect the case where features contains only discrete values stored as the R `numeric` datatype. If the features input are discrete the scaling and added noise are not applied to that argument. Discrete features should be used with the `mutual_inf_cd` function as it is approximately 2x faster than `mutual_inf_cc`. 

If the features variable is a matrix of size MxN the code will loop over the rows and the output MI will be a vector of size M.

### mutual_inf_cd
`mutual_inf_cd <- function(target, features, k=3L)` The target should be continuous and the features should be discrete numeric values. If the length of the target vector is N the the features argument must be a vector of length N or a matrix of size MxN. 

The neighborhood size k is the same as in `mutual_inf_cc`. 

The target vector is scaled as in `mutual_inf_cc`, unless it actually contains discrete numeric values. The features input is assumed to contain discrete numeric values and is converted to an R integer type. 

If the features variable is a matrix of size MxN the code will loop over the rows and the output MI will be a vector of size M.

### cond_mutual_inf
`cond_mutual_inf <- function(X, Y, Z, k=3L)` X is a continuous value. Y and Z can be continuous or discrete. If the length of the X vector is N the dimensions of Y and Z must be vectors of length N or matrices of size MxN. The inputs are treated for scaling and added noise as they are in `mutual_inf_cc`. 

If the Y is a matrix of size MxN and Z is a vector of size N then the return value will be a vector of length M, likewise if Y is a vector and Z is a matrix. For the case where both Y and Z are matrices of size MxN the return value will be a vector of size M where each row of Y and Z are paired to compute the CMI. 

### Random Noise

The random noise that is added to continuous inputs in all three functions is computed by calling the R random number generator. The magnitude of the added noise is 1e-12 (Kraskov, et.al. 2004). For reproducibility, a seed value can be set in R before calling the knnmi functions:
```
library(knnmi)
set.seed(654321)
x <- c(0.1, 0.3, 0.6, 0.2, 0.3, 0.4)
y <- c(0,0,1,1,0,1)
mi <- mutual_inf_cd(x,y)
print(mi)
# [1] 0.4638889
```


# Example

### mutual_inf_cc
```
target <- c(0.91 , 0.08, 0.13, 0.01, 0.08,
       0.25, 0.86, 0.13 , 0.26, 0.7, 0.18, 0.7, 0.26, 0.27 , 0.46)

features <- c(1.75,  1.14,  0.99,  0.96,  1.08,
       1.18,  1.63,  1.03,  0.95, -0.31, 0.89,  1.45,  1.02,  0.97,  1.25)       
set.seed(654321)
# the default for k is 3
result <- mutual_inf_cc(target, features)
print(result)
# [1] 0.18977
# Try a larger value for k
result_4 <- mutual_inf_cc(target, features, k=4)
print(result_4)
# [1] 0.205987
```


```
target <- c(0.4,0.5,0.6,0.1,0.05,0.02)
# discrete values for the features argument *can* be used but this 
# function runs 2x slower than mutual_inf_cd, so this is not 
# recommended.
features <- c(1,1,1,0,0,0)
set.seed(1234)
result <- mutual_inf_cc(target,features)
print(result)
# [1] 0.225
```

### mutual_inf_cd

```
target <- c(0.75, 1.07, 2.25,0.77, 1.94, 1.63, -0.28, 0.78, 1.55, 0.35)
features <- c(0, 0, 0, 1, 0, 1, 0, 0, 1, 0)
set.seed(1234)
result <- mutual_inf_cd(target, features)
print(result)
# [1] 0.1006349
```

### cond_mutual_inf

```
x <- c(0.9279281, 1.3391192, 1.0468714, 0.2661599, -1.2366708, 1.1937904, -1.1125911, -3.5276204, 0.5838281, 0.4193065)
y <- c(-0.56047565, -0.23017749, 1.55870831, 0.07050839, 0.12928774, 1.71506499, 0.46091621, -1.26506123, -0.68685285, -0.44566197)
z <- c(-0.71040656, 0.25688371, -0.24669188, -0.34754260, -0.95161857, -0.04502772, -0.78490447, -1.66794194, -0.38022652, 0.91899661)
set.seed(1234)
result <- cond_mutual_inf(x, y, z)
print(result)
# [1] 0.05071429
```
