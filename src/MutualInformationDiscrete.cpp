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


#include "MutualInformationDiscrete.h"

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

MutualInformationDiscrete::MutualInformationDiscrete(const int mK) : MutualInformationBase(mK) {}

MutualInformationDiscrete::~MutualInformationDiscrete() {}

double MutualInformationDiscrete::compute(const ArrayXd &x, const ArrayXi &y) {
  // Compute the mutual information for a continuous vector x and a
  // discrete vector y.
  // This implements the algorithm described in: https://doi.org/10.1371/journal.pone.0087357
  // Ross BC (2014) Mutual Information between Discrete and Continuous Data Sets. PLoS ONE 9(2): e87357. 
  auto N = x.size() ;
  ArrayXd x_scale =  scale(x, !check_if_int(x)) ;
 
  // Make a kdtree for x_scale, it'll be needed later.
  kd_tree_1d xscale_index_tree(1, x_scale, 10);
  
  // Chebyshev distance
  ArrayXd dists(N) ;
  // Number of neighbors
  ArrayXd neighbors(N) ;
  
  // Get the unique labels in y by creating an STL set
  std::set<int> unique_labels{y.data(), y.data() + y.size()};

  // For each unique label store its count.
  vector<double> label_counts ;
  
  // Store the k's used for each label
  vector<vector<double>> k_all ;
  
  // Store the neighbor counts m for each label
  vector<double> m_all ;
  // master vector of all indices that are not of unique labels.
  vector<int> all_indices ;
  
  for (const auto label : unique_labels) {
    int count = (y == label).count();
    if (count > 1) {
      // Only process the non-unique labels
      label_counts.push_back(count) ;
      // Adjust real_k as necessary.
      auto real_k = min(m_k, count) + 1 ;
      vector<double> tmp = {static_cast<double>(count),static_cast<double>(real_k)};
      k_all.push_back(tmp) ;
      // Store the indices of y where they match
      // the current label.
      vector<int> label_indices(count);
      int j = 0 ;
      // Eigen doesn't seem to have a convenient way to do this
      // so just use a for loop.
      for (int i = 0; i < N; ++i) {
        if (y[i] == label) {
          label_indices[j] = i ; ++j ;
          all_indices.push_back(i) ;
        }
      }
      // Use Eigen 3.4's method of providing a vector of indices
      // to produce a sub-vector.
      ArrayXd masked_x = x_scale(label_indices) ;
      // Make a lookup tree for the points of this label.
      kd_tree_1d label_index_tree(1, masked_x, 10);
      // Get all of the distances for each point for this label.
      for (int i = 0; i < count; ++i) {
        vector<Eigen::Index> ret_indexes(real_k, 0.0);
        vector<double> out_dists(real_k,0.0);
        
        double query_pt[1] = {masked_x[i]};
        // out_dists stores the distances for this label. Get the max one.
        auto neighbors = label_index_tree.index->knnSearch(query_pt, 
                                                           real_k,
                                                           &ret_indexes[0], 
                                                           &out_dists[0]);
        out_dists.resize(neighbors) ;
        // The last one is out_dists is the furthest distance.
        auto max_dist = out_dists.back() ;
        
        std::vector<std::pair<Eigen::Index, double>> ret_matches;
        m_all.push_back(xscale_index_tree.index->radiusSearch(query_pt, 
                                                              max_dist, 
                                                              ret_matches , 
                                                              nanoflann::SearchParams(10))) ;
      }
    }
  }
  double N_mod = all_indices.size() ;
  double digamma_N = digamma(N) ;
  // Calculate the mean of digammas over the count of samples for each label.
  double digamma_labels = std::accumulate(label_counts.begin(),label_counts.end(),0.0,
                                          [&](double m, double n){ return m + 
                                            digamma(n) * n / N_mod ; } ) ;
  // Get the same for the k's used at each label.
  double digamma_k = std::accumulate(k_all.begin(),k_all.end(),0.0,
                                     [&](double m, vector<double> &n){ return m +
                                       digamma(std::max(n[1] - 1.0, 1.0)) * n[0] / N_mod; } ) ;
  double digamma_m = std::accumulate(m_all.begin(),m_all.end(), 0.0,  
                                     [&](double m, double n){return m + digamma(n); });
  digamma_m = digamma_m / N_mod ;
  
  // mutual info computation
  double mi = digamma_N - digamma_labels + digamma_k - digamma_m ;
  // Can't return less than 0.
  return std::max(0.0,mi) ;
}


 
} // CaDrA
