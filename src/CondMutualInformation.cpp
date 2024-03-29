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
  


#include "CondMutualInformation.h"

#include <cmath>
#include <array>
#include <numeric>
#include <functional>
#include <list>
#include <set>
#include <map>

#include <Rmath.h>
#include "R.h"
#include "Rinternals.h"
#include <R_ext/Random.h>

#include "nanoflann.hpp"

namespace CaDrA {

CondMutualInformation::CondMutualInformation(const int mK) : MutualInformationBase(mK) {}

CondMutualInformation::~CondMutualInformation() {}

double CondMutualInformation::compute(const ArrayXd &x, const ArrayXd& y, const ArrayXd& z) {
  // This implements the CMI algorithm described in: https://doi.org/10.1016/j.eswa.2012.05.014
  // 
  // Alkiviadis Tsimpiris, Ioannis Vlachos, Dimitris Kugiumtzis,
  // Nearest neighbor estimate of conditional mutual information in feature selection,
  // Expert Systems with Applications,
  // Volume 39, Issue 16, 2012, Pages 12697-12708
  //
  // Compute the conditional mutual information.
  // CMI = psi(m_k) - mean( psi(n_xz[i] + psy(n_yz[i])-psi(n_z[i]) ) for all elements i
  auto N = x.size() ;
  
  // Get the 3-dimensional distance, and then re-use that for the xz, yz, and calculations, so it
  // proceeds pretty much exactly as before.
  
  Array3col tmp_mat(N, 3) ;
  
  // If not an integer valued array in double representation
  // scale it. Either way add a bit of noise to avoid duplicate
  // values.
  tmp_mat.col(0) = scale(x, !check_if_int(x)) ;
  tmp_mat.col(1) = scale(y, !check_if_int(y)) ;
  tmp_mat.col(2) = scale(z, !check_if_int(z)) ;
  
  // Calculating the distances also calculates the number of neighbors, so
  // calc_distances2 returns both.
  vector<double> dists = calc_distances3d(N, tmp_mat) ;
  
  // Get the digamma_f values...These take single vector arguments
  // so map the 3D temp array without a copy.
  MapArrayConst x_scale(tmp_mat.col(0).data(), N) ;
  MapArrayConst y_scale(tmp_mat.col(1).data(), N) ;
  MapArrayConst z_scale(tmp_mat.col(2).data(), N) ;

  // The z-projection is done using the same approach used in MutualInformation
  double z_digamma_sum = MutualInformationBase::sum_digamma_from_neighbors(z_scale, dists) ;

  // and for the 2D neighbors use the local function.
  double xz_digamma_sum = sum_digamma_from_neighbors(x_scale, z_scale, dists) ;
  double yz_digamma_sum = sum_digamma_from_neighbors(y_scale, z_scale, dists) ;


  // mutual info computation
  double mi = digamma(m_k) - (xz_digamma_sum + yz_digamma_sum - z_digamma_sum) / N;
  
  // Can't return less than 0.
  return std::max(0.0,mi) ; 
}

double CondMutualInformation::sum_digamma_from_neighbors(MapArrayConst &vec1, MapArrayConst &vec2, const vector<double> &dists) {
  // Sum of digamma_f functions over neighbor counts for 2D.
  long N = dists.size() ;
  
  // KD-Tree for this vector
  Array2col tmp_mat(N, 2) ;
  tmp_mat.col(0) = vec1 ;
  tmp_mat.col(1) = vec2 ;
  // 
  nanoflann::KDTreeEigenMatrixAdaptor<Array2col ,-1, metric_Chebyshev> vec_tree(2, tmp_mat) ;
  
  std::vector<std::pair<Eigen::Index, double>> ret_matches;
  ret_matches.reserve(dists.size()) ;
  array<double,2> pt ;
  double sum = 0.0 ;
  for (long i = 0 ; i < N ; ++i) {
    pt[0] = tmp_mat(i,0) ;
    pt[1] = tmp_mat(i,1) ;
    double radius_count = vec_tree.index->radiusSearch(pt.data(), dists[i] , ret_matches, nanoflann::SearchParams());
    sum += digamma(radius_count) ;
    ret_matches.clear() ;
  }
  return sum ;
}


vector<double> CondMutualInformation::calc_distances3d(const long N, const Array3col &tmp_mat) const {
  // Calculate the Chebyshev distances and numbers of neighbors for the 3D array tmp_mat.
  kd_tree_3d mat_index(tmp_mat.cols(),std::cref(tmp_mat),10) ;
  // We want N neighbors in addition to the point itself so
  // add 1 to the # of neighbors.
  int real_k = m_k + 1  ;
  
  // Chebyshev distance
  vector<double> dists(N) ;
  // Number of neighbors
  vector<long> neighbors(N) ;
  
  // a query point.
  array<double,3> query_pt ;
  for (long i = 0 ; i < N ; ++i) {
    // store indexes and distances
    vector<Eigen::Index> ret_indexes(real_k, 0.0);
    vector<double> out_dists(real_k,0.0);
    
    query_pt[0] = tmp_mat(i, 0);
    query_pt[1] = tmp_mat(i, 1);
    query_pt[2] = tmp_mat(i, 2);

    neighbors[i] = mat_index.index->knnSearch(&query_pt[0], real_k,
                                              &ret_indexes[0], &out_dists[0]) ;
    out_dists.resize(neighbors[i]) ;
    dists[i] = std::nextafter(*max_element(std::begin(out_dists), std::end(out_dists)),0.0) ;
  }
  return  dists  ;
}

} // CaDrA
