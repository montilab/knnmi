//
// Created by bgregor on 11/1/22.
//

/* Copyright (C) 2022 - Trustees of Boston University
 
 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.
 
 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.
 
 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <https://www.gnu.org/licenses/>. */


#include "MutualInformation.h"

#include <cmath>
#include <array>
#include <numeric>
#include <functional>
#include <list>
#include <set>
#include <map>
#include <functional> 
#include <algorithm>

#include <Rmath.h>
#include "R.h"
#include "Rinternals.h"
#include <R_ext/Random.h>

#include "nanoflann.hpp"

namespace CaDrA {

MutualInformation::MutualInformation(const int mK) : MutualInformationBase(mK) {}

MutualInformation::~MutualInformation() {}

 
 
double MutualInformation::compute(const ArrayXd &x, const ArrayXd& y) {
  // compute mutual information based on k-nearest neighbors
  // This implements the algorithm described in: https://doi.org/10.1103/PhysRevE.69.066138
  // Alexander Kraskov, Harald Stögbauer, and Peter Grassberger
  // Phys. Rev. E 69, 066138 (2004).
  
  // Create a KD-tree using nanoflann and its Eigen matrix adapter
  // Store the size of the vectors for easy reference
  
  auto N = x.size() ;
  
  Array2col  tmp_mat(N, 2) ;
  // If not an integer valued array in double representation
  // scale it. Either way add a bit of noise to avoid duplicate
  // values.
  tmp_mat.col(0) = scale(x, !check_if_int(x)) ;
  tmp_mat.col(1) = scale(y, !check_if_int(y)) ;
  
  // Map the double array pointer to an Eigen vector without a copy.
  MapArrayConst x_scale(tmp_mat.col(0).data(), N) ;
  MapArrayConst y_scale(tmp_mat.col(1).data(), N) ;
  
  vector<double> dists  = calc_distances2d(N, tmp_mat) ;
 
  double x_digamma_sum = sum_digamma_from_neighbors(x_scale,dists) ;
  double y_digamma_sum = sum_digamma_from_neighbors(y_scale,dists) ;
  
  // mutual info computation
  double Nd = N ;
  double mi = digamma(Nd) + digamma(m_k) - (x_digamma_sum + y_digamma_sum) / Nd;
  
  // Can't return less than 0.
  return std::max(0.0,mi) ;
}

vector<double> MutualInformation::calc_distances2d(const long N, const Array2col &tmp_mat) const {
  // Calculate the Chebyshev distances and numbers of neighbors for the 2D array tmp_mat.
  // Adjust the max size of leaves in the kd-tree if you want, 10 is the default.
  int leaf_max_size = 10 ;
  kd_tree_2d mat_index(tmp_mat.cols(),std::cref(tmp_mat),leaf_max_size) ;
  
  // We want N neighbors in addition to the point itself so
  // add 1 to the # of neighbors.
  int real_k = m_k + 1 ;
  
  // Chebyshev distance
  vector<double> dists(N) ;
  // Number of neighbors
  vector<long> neighbors(N) ;
  
  // a query point.
  array<double,2> query_pt ;
  for (long i = 0 ; i < N ; ++i) {
    // store indexes and distances
    vector<Eigen::Index> ret_indexes(real_k, 0.0);
    vector<double> out_dists(real_k,0.0);
    
    query_pt[0] = tmp_mat(i, 0);
    query_pt[1] = tmp_mat(i, 1);
    
    neighbors[i] = mat_index.index->knnSearch(&query_pt[0], real_k,
                                              &ret_indexes[0], &out_dists[0]) ;
    out_dists.resize(neighbors[i]) ;
    dists[i] = std::nextafter(*max_element(std::begin(out_dists), std::end(out_dists)),0.0)  ;  
  }
  return dists ;
}




} // CaDrA
