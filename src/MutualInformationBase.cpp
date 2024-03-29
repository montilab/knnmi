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


#include "MutualInformationBase.h"

#include <cmath>
#include <array>
#include <numeric>
#include <limits> 
#include <cstdint> 

#include <Rmath.h>
#include "R.h"
#include "Rinternals.h"
#include <R_ext/Random.h>

#include "nanoflann.hpp"

namespace CaDrA {

MutualInformationBase::MutualInformationBase(const int mK) : m_k(mK) {
  // https://cran.r-project.org/doc/manuals/R-exts.html#Random-numbers
  GetRNGstate();
}

MutualInformationBase::~MutualInformationBase() {
  PutRNGstate();
}

int MutualInformationBase::get_k() const {
  return m_k;
}


void MutualInformationBase::set_k(int mK) {
  m_k = mK;
}
 

ArrayXd MutualInformationBase::scale(const ArrayXd &x, const bool apply_scale, const bool add_noise) const {
  // Scale the vector x by its standard deviation. Add a bit of random 
  // noise.
  ArrayXd x_scale = x ;
  double mean_x ; 
  auto size_x = x.size() ;
  if (apply_scale) {
    double std_dev = std::sqrt(x_scale.pow(2).sum() / (size_x - 1)) ;
    x_scale = x_scale  / std_dev ; 
  }  
  // add a wee bit of noise as suggested in
  // Kraskov et. al. 2004. For the CaDrA package the suggested value
  // of 1e-10 seems too high and leads to too much variance in the MI
  // calculation. The point of the noise is to avoid points too close
  // together which may lead to accidentally double-counting them.
  // Trying 1e-12 which seems ok.
  if (add_noise) {
    // Get a uniform value from the R RNG
    // https://cran.r-project.org/doc/manuals/R-exts.html#Random-numbers
    mean_x = x_scale.mean() ;
    for (auto i = 0 ; i < size_x ; i++) {
      x_scale[i] += 1e-12 * mean_x * unif_rand();
    }
  }
  return x_scale;
}

bool MutualInformationBase::check_if_int(const ArrayXd &x) {
  // check if an array is made up of integers masquerading as doubles.
  // This is used as a decision as to whether or not to center & scale
  // the array.
  
  // Cast each element of the array to a long integer. Subtract 
  // from x. If at any point the difference is greater than the 
  // double-prec eps return false. Otherwise return true. This
  // will quit immediately if a non-integer value is found.
  auto eps = std::numeric_limits<double>::epsilon() ;
  for (long i = 0 ; i < x.size() ; ++i) {
    auto val = static_cast<std::int64_t>(x[i]) ; 
    if ((x[i] - val) > eps) {
      return false ;
    }
  }
  return true ;
}

 

double MutualInformationBase::sum_digamma_from_neighbors(MapArrayConst &vec, const vector<double> &dists) {
  // This one is called from mutual_information_cc and cond_mutual_information for the neighbors
  // for a single vector.
  long N = dists.size() ;
  double sum = 0.0 ;
  
  // KD-Tree for this vector
  nanoflann::KDTreeEigenMatrixAdaptor<MapArrayConst,-1,nanoflann::metric_L1> vec_tree(1, vec) ;
  
  std::vector<std::pair<Eigen::Index, double>> ret_matches;
  for (long i = 0 ; i < N ; ++i) {
    double pt = vec(i) ; // avoids type issues with the compiler and the radiusSearch.
    double tmp = vec_tree.index->radiusSearch(&pt, dists[i] , ret_matches , nanoflann::SearchParams(10));
    sum += digamma(tmp) ;
    ret_matches.clear() ;
  }
  return sum ;
}

 
 
// placeholder
double MutualInformationBase::compute(){ return 0.0 ;}

} // CaDrA
