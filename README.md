# Introduction

knnmi is a mutual information library for R based on the k-nearest neighbor algorithm. There are three functions:
* mutual_inf_cc - Compute MI(X;Y) where X is a vector and Y is a vector or matrix. X and Y are both continuous.
* mutual_inf_cd  - Compute MI(X;Y) where X is a vector and Y is a vector or matrix. X is continuous and Y contains discrete numeric values.
* cond_mutual_inf - Compute the condtional mutual information CMI(X;Y|Z) where X is a continous vector. Y and Z are continous and be any mix of vectors and matrices.

# Implementation
The k-nearest neighbor (KNN) approach to mutual information (MI) is used in popular libraries like 
[scikit-learn](https://scikit-learn.org/stable/modules/generated/sklearn.feature_selection.mutual_info_regression.html#sklearn.feature_selection.mutual_info_regression) 
for Python. 
Here, the KNN computations are handled by the high performance [nanoflann](https://github.com/jlblancoc/nanoflann) C++ library with the addition of a custom Chebyshev distance metric. 
The Eigen C++ library is used to handle vectors and an effort was made to minimize the number of copies of the incoming data from R to preserve memory. 
A newer Eigen version (3.4.0) was used for some useful features and a minimal copy of that library is part of this repository. The newest Eigen version built 
into RcppEigen is 3.3.4. 

The publications used as the basis for the knnmi library are:
* mutual_inf_cc - [A. Kraskov, H. Stögbauer, and P. Grassberger, Phys. Rev. E 69, 066138 (2004)](https://doi.org/10.1103/PhysRevE.69.066138)
* mutual_inf_cd - [Ross BC (2014) Mutual Information between Discrete and Continuous Data Sets. PLoS ONE 9(2): e87357](https://doi.org/10.1371/journal.pone.0087357)
* cond_mutual_inf - [A. Tsimpiris, I. Vlachos, and D. Kugiumtzis. Expert Systems with Applications, Volume 39, Issue 16, (2012)](https://doi.org/10.1016/j.eswa.2012.05.014)

# Calling Conventions 
