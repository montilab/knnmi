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
