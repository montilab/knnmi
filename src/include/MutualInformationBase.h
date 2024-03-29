//
// Created by bgregor on 11/1/22.
//

#ifndef REVEALER_MUTUALINFORMATIONBASE_H
#define REVEALER_MUTUALINFORMATIONBASE_H

#include <eigen3/Eigen/Core>
#include <vector>

#include "nanoflann.hpp"
#include "ChebyshevMetric.h"

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

    class MutualInformationBase {
        // Methods to compute mutual information for
        // continuous and discrete variables.
    public:
        // Constructor
        // Initialize with a neighbor size.
        MutualInformationBase(const int k) ;
        MutualInformationBase() = delete ;
        
        // Destructor
        virtual ~MutualInformationBase() ;

        // getter/setter for m_k.
        int get_k() const;
        void set_k(int mK);

        // Computation functions. Tag with _c and _d for continuous
        // and discrete forms...the compiler is giving complaints about
        // overloads and this is the easy way out.
        virtual double compute() ;

    protected:
        int m_k ;
        
        // center & scale & add a bit of noise
        virtual ArrayXd scale(const ArrayXd &x, const bool apply_scale=true, const bool add_noise=true) const ;
        
        // check if an array is made up of integers masquerading as doubles.
        virtual bool check_if_int(const ArrayXd &x) ;
        
        // digamma function
        vector<double> digamma_vec(vector<double> counts) const ;
        
        // Get the sum of digamma functions for the number of nearest neighbors for all points.
        // Used in both the MutualInformation and CondMutualInformation subclasses.
        virtual double sum_digamma_from_neighbors(MapArrayConst &vec, const vector<double> &dists) ;
    };

} // CaDrA

#endif //REVEALER_MUTUALINFORMATIONBASE_H
