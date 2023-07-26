//
// Created by bgregor on 11/1/22.
//
#include "MutualInformationBase.h"

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
 

ArrayXd MutualInformationBase::add_noise(const ArrayXd &x) const {
  // Scale the vector x by its standard deviation. Add a bit of random 
  // noise.
  /*auto size_v = x.size() ;
  double mean_x = x.mean() ;
  double sum = 0.0 ;
  for (auto i = 0 ; i < size_v ; i++) {
    double tmp = x[i] - mean_x ;
    sum += tmp * tmp ;
  }
  double std_dev = std::sqrt(sum / (size_v-1)) ;
  ArrayXd x_scale = (x - mean_x)  / std_dev ;*/
  
  // add a wee bit of noise as suggested in
  // Kraskov et. al. 2004.
  //if (add_noise) {
  ArrayXd x_noise = x ;
  double mean_xs = x.mean() ;
  // Get a uniform value from the R RNG
  // https://cran.r-project.org/doc/manuals/R-exts.html#Random-numbers
  for (auto i = 0 ; i < x.size() ; i++) {
    x_noise[i] += 1e-10 * mean_xs * unif_rand();
  }
  //}
  return x_noise;
}

vector<double> MutualInformationBase::count_neighbors(MapArrayConst &vec, const vector<double> &dists) {
  // This one is called from mutual_information_cc and cond_mutual_information for the neighbors
  // for a single vector.
  long N = dists.size() ;
  vector<double> neighbors(N,0) ;

  // KD-Tree for this vector. L1 distance is identical to the Chebyshev distance in 1D. 
  nanoflann::KDTreeEigenMatrixAdaptor<MapArrayConst,-1,nanoflann::metric_L1> vec_tree(1, vec) ;
  
  std::vector<std::pair<Eigen::Index, double>> ret_matches;
  ret_matches.reserve(dists.size());
  
  for (long i = 0 ; i < N ; ++i) {
    double pt = vec(i) ; // avoids type issues with the compiler and the radiusSearch.
    double radius_count = vec_tree.index->radiusSearch(&pt, dists[i] , ret_matches, nanoflann::SearchParams());
    neighbors[i] = radius_count ;
    ret_matches.clear() ;
  }
  return neighbors ;
}


// Calls R's digamma function for all elements of a vector.
vector<double> MutualInformationBase::digamma_vec(vector<double> counts) const {
    vector<double> result(counts.size()) ;
    std::transform (counts.begin(), counts.end(), result.begin(), digamma);
    return result ;
}
 
 
 
// placeholder
double MutualInformationBase::compute(){ return 0.0 ;}

} // CaDrA
