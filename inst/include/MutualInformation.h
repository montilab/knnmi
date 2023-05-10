//
// Created by bgregor on 11/1/22.
//

#ifndef REVEALER_MUTUALINFORMATION_H
#define REVEALER_MUTUALINFORMATION_H

#include <eigen3/Eigen/Core>
#include <random>
#include <vector>

#include "nanoflann.hpp"
#include "ChebyshevMetric.h"
#include "pcg_random.hpp"



using namespace Eigen;
using namespace std ;

namespace CaDrA {

    // A const map that allows a double* array to be read as
    // an Eigen ArrayXd.
    typedef Map<const ArrayXd> MapArrayConst;
    typedef Map<const ArrayXi> MapArrayIConst;
    
    typedef Array<double, Dynamic, 2> Array2col  ;
    typedef Array<double, Dynamic, 3> Array3col  ;
    typedef nanoflann::KDTreeEigenMatrixAdaptor<ArrayXd, -1, metric_Chebyshev> kd_tree_1d ;
    typedef nanoflann::KDTreeEigenMatrixAdaptor<Array2col, -1, metric_Chebyshev> kd_tree_2d ;
    typedef nanoflann::KDTreeEigenMatrixAdaptor<Array3col, -1, metric_Chebyshev> kd_tree_3d ;

    class MutualInformation {
        // Methods to compute mutual information for
        // continuous and discrete variables.
    public:
        // Constructor
        // Initialize with a neighbor size and a seed. A seed value
        // of <= 0 results in a random seed being used.
        MutualInformation(const int k, const int seed) ;
        MutualInformation() = delete ;
        
        // Destructor
        virtual ~MutualInformation() ;
        
        // Compute conditional mutual information of 3 continuous variables
        double cond_mutual_information_ccc(const ArrayXd &x, const ArrayXd& y, const ArrayXd& z) ;

        // Compute conditional mutual information of a continuous variable and two discrete variables
        double cond_mutual_information_cdd(const ArrayXd &x, const ArrayXi& y, const ArrayXi& z) ;

        //  Computes the Mutual Information of 2 continuous variables.
        double mutual_information_cc(const ArrayXd &x, const ArrayXd& y) ;

        //  Computes the Mutual Information of a continuous and a discrete variable
        double mutual_information_cd(const ArrayXd &x, const ArrayXi& y) ;

        // getter/setter for m_k.
        int get_k() const;
        void set_k(int mK);
        // getter for the seed value.
        int get_seed() const ;
    protected:
        int m_k ;
        int m_seed ;
        // The PCG64 RNG is dynamically allocated.
        pcg64 *m_rng ;

        ArrayXd scale(const ArrayXd &x, bool add_noise=true) const;

        // Compute the sum of digamma_f functions as nearest neighbors are calculated.  1D and 2D versions.
        double sum_digamma_from_neighbors(MapArrayConst &vec, const vector<double> &dists) ;
        double sum_digamma_from_neighbors(MapArrayConst &vec1, MapArrayConst &vec2, const vector<double> &dists) ;

        // Calculate distances and nearest neighbors in 2D and 3D.
        pair<vector<double>,vector<long>>  calc_distances2d(const long N, const Array2col &tmp_mat) const;
        pair<vector<double>,vector<long>>  calc_distances3d(const long N, const Array3col &tmp_mat) const;

        static double digamma_f(const double x)  ;
    };

} // CaDrA

#endif //REVEALER_MUTUALINFORMATION_H
