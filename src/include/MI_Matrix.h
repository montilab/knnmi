//
// Created by bgregor on 9/19/22.
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

#ifndef REVEALER_MI_MATRIX_H
#define REVEALER_MI_MATRIX_H

// This provides C-callable wrappers for the mutual information implementation
// which is in C++. The functions here are quite simple, they wrap the incoming
// arrays as Eigen Array types to call the

extern "C" {
        // and for when the input score is just 1 vector
        int mutual_inf_cc_vec(const double *input_x, const double *input_y, 
                              const int n_elems,  const int k, double *mi) ;

        // case where input_y is discrete integers
        int mutual_inf_cd_vec(const double *input_x, const int *input_y, const int n_elems, 
                              const int k, double *mi ) ;

        //  But - there's only 1 algorithm for the conditional mutual information, so there's
        // only a need for 1 function for this.
        // Conditional mutual information.
       int cond_mutual_inf_vec(const double *input_x,  const double *input_y, const double *input_z, 
                               const int n_elems, const int k, double *mi) ;
}
#endif //REVEALER_MI_MATRIX_H
