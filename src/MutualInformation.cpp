//
// Created by bgregor on 11/1/22.
//
#include "MutualInformation.h"

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

#include "nanoflann.hpp"

namespace CaDrA {

MutualInformation::MutualInformation(const int mK, const int seed) : m_k(mK), m_seed(seed) {
  // Construct this object with the neighborhood size k and a 
  // seed. If the seed is <= 0 then use a random seed.
  
  if (seed > 0) {
    // Reverse the bits in seed.
    // https://www.geeksforgeeks.org/write-an-efficient-c-program-to-reverse-bits-of-a-number/
    unsigned int useed = seed ;
    unsigned int count = sizeof(useed) * 8 - 1;
    unsigned int reverse_useed = useed;
    
    useed >>= 1;
    while (useed) {
      reverse_useed <<= 1;
      reverse_useed |= useed & 1;
      useed >>= 1;
      count--;
    }
    reverse_useed <<= count;
    // use seed and its reversed bits representation to help
    // build a seed sequence.
    std::seed_seq seq = {seed, ~seed, static_cast<int>(reverse_useed), 
                         ~static_cast<int>(reverse_useed)};
    m_rng = new pcg64(seq);
  } else {
    pcg_extras::seed_seq_from<std::random_device> ran_seed_source;
    m_rng = new pcg64(ran_seed_source);
  }
}

MutualInformation::~MutualInformation() {
  if (m_rng != NULL) {
    delete m_rng ;
  }
}

int MutualInformation::get_k() const {
  return m_k;
}

int MutualInformation::get_seed() const {
  return m_seed;
}

void MutualInformation::set_k(int mK) {
  m_k = mK;
}

double MutualInformation::mutual_information_cc(const ArrayXd &x, const ArrayXd& y) {
  // compute mutual information based on k-nearest neighbors
  // This implements the algorithm described in: https://doi.org/10.1103/PhysRevE.69.066138
  // Alexander Kraskov, Harald Stogbauer, and Peter Grassberger
  // Phys. Rev. E 69, 066138 ?? Published 23 June 2004; Erratum Phys. Rev. E 83, 019903 (2011)
  
  // Creat a KD-tree using nanoflann and its Eigen matrix adapter
  // Store the size of the vectors for easy reference
  
  auto N = x.size() ;
  
  Array2col  tmp_mat(N, 2) ;
  tmp_mat.col(0) = scale(x) ;
  tmp_mat.col(1) = scale(y) ;
  // Map the double array pointer to an Eigen vector without a copy.
  MapArrayConst x_scale(tmp_mat.col(0).data(), N) ;
  MapArrayConst y_scale(tmp_mat.col(1).data(), N) ;
  
  vector<double> dists  = calc_distances2d(N, tmp_mat).first;
  
  double x_digamma_sum = sum_digamma_from_neighbors(x_scale, dists) ;
  double y_digamma_sum = sum_digamma_from_neighbors(y_scale, dists) ;
  
  // mutual info computation
  double mi = MutualInformation::digamma_f(N)
    + MutualInformation::digamma_f(m_k)
    - (x_digamma_sum + y_digamma_sum) / N;
    
    return std::max(0.0,std::min(mi,1.0)) ;
}

double MutualInformation::mutual_information_cd(const ArrayXd &x, const ArrayXi &y) {
  // Compute the mutual information for a continuous vector x and a
  // discrete vector y.
  // This implements the algorithm described in: https://doi.org/10.1371/journal.pone.0087357
  // Ross BC (2014) Mutual Information between Discrete and Continuous Data Sets. PLoS ONE 9(2): e87357. 
  
  auto N = x.size() ;
  ArrayXd x_scale = scale(x) ;
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
    if (count > 0) {
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
  double digamma_N = digamma_f(N) ;
  // Calculate the mean of digammas over the count of samples for each label.
  double digamma_labels = std::accumulate(label_counts.begin(),label_counts.end(),0.0,
                                          [N_mod](double m, double n){ return m +
                                            MutualInformation::digamma_f(n) * n / N_mod ; } ) ;
  // Get the same for the k's used at each label.
  double digamma_k = std::accumulate(k_all.begin(),k_all.end(),0.0,
                                     [N_mod](double m, vector<double> &n){ return m +
                                       MutualInformation::digamma_f(std::max(n[1] - 1.0, 1.0)) * n[0] / N_mod; } ) ;
  double digamma_m = std::accumulate(m_all.begin(),m_all.end(),
                                     0.0,  [](double m, double n){
                                       return m + MutualInformation::digamma_f(n); });
  digamma_m = digamma_m / N_mod ;
  
  // mutual info computation
  // Matlab:   4.602 - 4.3965 + 0.9228 - 1.0961
  double mi = digamma_N - digamma_labels + digamma_k - digamma_m ;
  return std::max(0.0,std::min(mi,1.0)) ;
}

double MutualInformation::cond_mutual_information_ccc(const ArrayXd &x, const ArrayXd& y, const ArrayXd& z) {
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
  
  tmp_mat.col(0) = scale(x) ;
  tmp_mat.col(1) = scale(y) ;
  tmp_mat.col(2) = scale(z) ;
  
  // Map the double array pointer to an Eigen vector without a copy.
  MapArrayConst x_scale(tmp_mat.col(0).data(), N) ;
  MapArrayConst y_scale(tmp_mat.col(1).data(), N) ;
  MapArrayConst z_scale(tmp_mat.col(2).data(), N) ;
  
  // Calculating the distances also calculates the number of neighbors, so
  // calc_distances2 returns both.
  vector<double> dists = calc_distances3d(N, tmp_mat).first ;
  
  // Get the digamma_f values...
  double xz_digamma_sum = sum_digamma_from_neighbors(x_scale, z_scale, dists) ;
  double yz_digamma_sum = sum_digamma_from_neighbors(y_scale, z_scale, dists) ;
  double z_digamma_sum = sum_digamma_from_neighbors(z_scale, dists) ;
  
  // mutual info computation
  double mi = MutualInformation::digamma_f(m_k)
    - (xz_digamma_sum + yz_digamma_sum - z_digamma_sum) / N;
  
  return std::max(0.0,std::min(mi,1.0)) ;
}

double MutualInformation::cond_mutual_information_cdd(const ArrayXd &x, const ArrayXi& y, const ArrayXi& z) {
  // Implement conditional mutual information between continuous x and discrete y & z.
  // Convert the y & z arrays over to double precision, call the
  // cond_mutual_information_ccc function.
  
  return cond_mutual_information_ccc(x, y.cast<double>(), z.cast<double>()) ;
}

double MutualInformation::sum_digamma_from_neighbors(MapArrayConst &vec, const vector<double> &dists) {
  // This one is called from mutual_information_cc and cond_mutual_information for the neighbors
  // for a single vector.
  long N = dists.size() ;
  double sum = 0.0 ;
  
  // KD-Tree for this vector
  nanoflann::KDTreeEigenMatrixAdaptor<MapArrayConst,-1,metric_Chebyshev> vec_tree(1, vec, 10) ;
  
  std::vector<std::pair<Eigen::Index, double>> ret_matches;
  for (long i = 0 ; i < N ; ++i) {
    double pt = vec(i) ; // avoids type issues with the compiler and the radiusSearch.
    double tmp = vec_tree.index->radiusSearch(&pt, dists[i] , ret_matches , nanoflann::SearchParams(10));
    sum += MutualInformation::digamma_f(tmp) ;
    ret_matches.clear() ;
  }
  return sum ;
}

double MutualInformation::sum_digamma_from_neighbors(MapArrayConst &vec1, MapArrayConst &vec2, const vector<double> &dists) {
  // Sum of digamma_f functions over neighbor counts for 2D.
  long N = dists.size() ;
  double sum = 0.0 ;
  
  // KD-Tree for this vector
  Array2col tmp_mat(N, 2) ;
  tmp_mat.col(0) = vec1 ;
  tmp_mat.col(1) = vec2 ;
  nanoflann::KDTreeEigenMatrixAdaptor<Array2col ,-1,metric_Chebyshev> vec_tree(2, tmp_mat, 10) ;
  
  std::vector<std::pair<Eigen::Index, double>> ret_matches;
  array<double,2> pt ;
  for (long i = 0 ; i < N ; ++i) {
    pt[0] = tmp_mat(i,0) ;
    pt[1] = tmp_mat(i,1) ;
    double tmp = vec_tree.index->radiusSearch(pt.data(), dists[i] , ret_matches , nanoflann::SearchParams(10));
    sum += MutualInformation::digamma_f(tmp) ;
    ret_matches.clear() ;
  }
  return sum ;
}



ArrayXd MutualInformation::scale(const ArrayXd &x, bool add_noise) const {
  // Scale the vector x by its standard deviation
  
  auto size_v = x.size() ;
  double mean_x = x.mean() ;
  double sum = 0.0 ;
  for (auto i = 0 ; i < size_v ; i++) {
    double tmp = x[i] - mean_x ;
    sum += tmp * tmp ;
  }
  double std_dev = std::sqrt(sum / (size_v-1)) ;
  ArrayXd x_scale = x  / std_dev ;
  
  
  // add a wee bit of noise as suggested in
  // Kraskov et. al. 2004.
  if (add_noise) {
    double mean_xs = x_scale.mean() ;
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    for (auto i = 0 ; i < x_scale.size() ; i++) {
      x_scale[i] += 1e-10 * mean_xs * dist(*m_rng)  ;
    }
  }
  return x_scale;
}


pair<vector<double>,vector<long>>  MutualInformation::calc_distances3d(const long N, 
                                                                       const Array3col &tmp_mat) const {
  // Calculate the Chebyshev distances and numbers of neighbors for the 3D array tmp_mat.
  kd_tree_3d mat_index(tmp_mat.cols(),std::cref(tmp_mat),20) ;
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
    vector<double> out_dists_sqr(real_k,0.0);
    
    for (long j = 0 ; j < 3 ; ++j)
      query_pt[j] = tmp_mat(i, j);
    
    neighbors[i] = mat_index.index->knnSearch(&query_pt[0], real_k,
                                              &ret_indexes[0], &out_dists_sqr[0]) ;
    out_dists_sqr.resize(neighbors[i]) ;
    dists[i] = out_dists_sqr.back() ;
  }
  return make_pair(dists,neighbors) ;
}



pair<vector<double>,vector<long>> MutualInformation::calc_distances2d(const long N, 
                                                                      const Array2col &tmp_mat) const {
  // Calculate the Chebyshev distances and numbers of neighbors for the 2D array tmp_mat.
  kd_tree_2d mat_index(tmp_mat.cols(),std::cref(tmp_mat),20) ;
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
    vector<double> out_dists_sqr(real_k,0.0);
    
    query_pt[0] = tmp_mat(i, 0);
    query_pt[1] = tmp_mat(i, 1);
    
    neighbors[i] = mat_index.index->knnSearch(&query_pt[0], real_k,
                                              &ret_indexes[0], &out_dists_sqr[0]) ;
    out_dists_sqr.resize(neighbors[i]) ;
    auto max_dist= std::nextafter(*max_element(std::begin(out_dists_sqr), std::end(out_dists_sqr)),0.0)  ; // following sklearn
    dists[i] = max_dist ;
  }
  return make_pair(dists,neighbors) ;
}

// For development outside of R use the Boost library's digamma_f funtion. 
// For R integration this can is switched to call the digamma_f
// function that comes with R.
inline double MutualInformation::digamma_f(const double x) {
#ifndef BOOST_DIGAMMA
  return digamma(x) ;
#else
  return boost::math::digamma(x) ;
#endif
}

} // CaDrA
