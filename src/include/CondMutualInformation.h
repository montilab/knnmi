//
// Created by bgregor on 11/1/22.
//

#ifndef REVEALER_CONDMUTUALINFORMATION_H
#define REVEALER_CONDMUTUALINFORMATION_H

#include <eigen3/Eigen/Core>
#include <vector>

#include "MutualInformationBase.h"
#include "nanoflann.hpp"
#include "ChebyshevMetric.h"

using namespace Eigen;
using namespace std ;

namespace CaDrA {
 

    class CondMutualInformation : public MutualInformationBase {
        // Methods to compute mutual information for
        // continuous and discrete variables.
    public:
        // Constructor
        // Initialize with a neighbor size.
        CondMutualInformation(const int k) ;
        CondMutualInformation() = delete ;
        
        // Destructor
        virtual ~CondMutualInformation() ;
        
        // Compute conditional mutual information of 3 continuous variables
        virtual double compute(const ArrayXd &x, const ArrayXd& y, const ArrayXd& z) ;

    protected:

        // Compute the sum of digamma_f functions as nearest neighbors are calculated. 
        virtual double sum_digamma_from_neighbors(MapArrayConst &vec1, MapArrayConst &vec2, const vector<double> &dists) ;

        // Calculate distances and nearest neighbors in 3D.
        vector<double> calc_distances3d(const long N, const Array3col &tmp_mat) const;
    };

} // CaDrA

#endif //REVEALER_CONDMUTUALINFORMATION_H
