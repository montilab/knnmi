//
// Created by bgregor on 10/3/22.
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


#ifndef REVEALER_CHEBYSHEVMETRIC_H
#define REVEALER_CHEBYSHEVMETRIC_H

#include "nanoflann.hpp"
#include <eigen3/Eigen/Core>

#include <vector>

namespace CaDrA {
/** Chebyshev distance metric
 *
 *  Following the example to implement a custom metric for nanoflann from: 
 *  
 *  https://github.com/jlblancoc/nanoflann/blob/v1.4.2/examples/pointcloud_custom_metric.cpp
 *
 *
 *  T Type of the elements (e.g. double, float, uint8_t)
 *  DataSource Source of the data, i.e. where the vectors are stored
 *  _DistanceType Type of distance variables (must be signed)
 *  AccessorType Type of the arguments with which the data can be
 * accessed (e.g. float, double, int64_t, T*)
 */
    template<
            class T, class DataSource, typename _DistanceType = T,
            typename AccessorType = uint32_t>
    struct Chebyshev_Adaptor {
        // This gets marked as unused here by the IDE but it gets used when the
        // template is instantiated by nanoflann. Don't remove it.
        using ElementType = T;
        using DistanceType = _DistanceType;

        const DataSource &data_source;

        Chebyshev_Adaptor(const DataSource &_data_source): data_source(_data_source) {}

        // Note the size argument has type Eigen::Index. This avoids template errors between
        // different compilers and operating systems.
        inline DistanceType evalMetric(const T *a, const AccessorType b_idx, Eigen::Index size) const {
            // Chebyshev distance: max difference between points in a and b
            DistanceType max_diff = std::numeric_limits<DistanceType>::min();
            for (auto i = 0; i < size; ++i) {
                DistanceType diff = std::fabs(a[i] - data_source.kdtree_get_pt(b_idx, i));
                if (diff > max_diff) {
                    max_diff = diff ;
                }
            }
           return max_diff ;
        }

        template<typename U, typename V>
        inline DistanceType accum_dist(const U a, const V b, const size_t) const {
            return std::fabs(a - b);
        }
    };


    struct metric_Chebyshev : public nanoflann::Metric
    {
        template <class T, class DataSource, typename AccessorType = uint32_t>
        struct traits
        {
            using distance_t = Chebyshev_Adaptor<T, DataSource, T, AccessorType>;
        };
    };


}

#endif //REVEALER_CHEBYSHEVMETRIC_H

 
