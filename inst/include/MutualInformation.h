//
// Created by bgregor on 11/1/22.
//

#ifndef REVEALER_MUTUALINFORMATION_H
#define REVEALER_MUTUALINFORMATION_H

#include <eigen3/Eigen/Core>
#include <vector>

#include "MutualInformationBase.h"
#include "nanoflann.hpp"
#include "ChebyshevMetric.h"

using namespace Eigen;
using namespace std ;

namespace CaDrA {
 

class MutualInformation : public MutualInformationBase {
  // Methods to compute mutual information for
  // continuous and discrete variables.
public:
  // Constructor
  // Initialize with a neighbor size.
  MutualInformation(const int k) ;
  MutualInformation() = delete ;
  
  // Destructor
  virtual ~MutualInformation() ;
  
  //  Computes the Mutual Information of 2 continuous variables.
  virtual double compute(const ArrayXd &x, const ArrayXd& y) ;
  
protected:
  // Calculate distances in 2D.  
  vector<double>  calc_distances2d(const long N, const Array2col &tmp_mat) const;
};

} // CaDrA

#endif //REVEALER_MUTUALINFORMATION_H
