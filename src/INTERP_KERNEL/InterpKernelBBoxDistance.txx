// Copyright (C) 2026  CEA, EDF
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//

#pragma once

#include <array>
#include <algorithm>

// fmt: off
namespace INTERP_KERNEL
{
inline void
DistInternal(double x, double boxLower, double boxUpper, double &xmin, double &xmax)
{
    if (x < boxLower)
    {
        xmin = (boxLower - x) * (boxLower - x);
        xmax = (boxUpper - x) * (boxUpper - x);
        return;
    }
    if (x > boxUpper)
    {
        xmin = (x - boxUpper) * (x - boxUpper);
        xmax = (x - boxLower) * (x - boxLower);
        return;
    }
    else
    {
        xmin = 0.0;
        xmax = std::max(boxUpper - x, x - boxLower) * std::max(boxUpper - x, x - boxLower);
        return;
    }
}

template <int dim>
inline void
HighLevelDistInternal(const double pt[dim], const double boundary[2 * dim], double &xmin, double &xmax)
{
    double xminTmp, xmaxTmp;
    DistInternal(pt[0], boundary[0], boundary[1], xminTmp, xmaxTmp);
    double boundary2[2 * (dim - 1)];
    double pt2[dim - 1];
    for (int i = 0; i < dim - 1; ++i)
    {
        pt2[i] = pt[i + 1];
        boundary2[2 * i] = boundary[2 * (i + 1)];
        boundary2[2 * i + 1] = boundary[2 * (i + 1) + 1];
    }
    double xminTmp2, xmaxTmp2;
    HighLevelDistInternal<dim - 1>(pt2, boundary2, xminTmp2, xmaxTmp2);
    xmin = xminTmp + xminTmp2;
    xmax = xmaxTmp + xmaxTmp2;
}

template <>
inline void
HighLevelDistInternal<1>(const double pt[1], const double boundary[2], double &xmin, double &xmax)
{
    DistInternal(pt[0], boundary[0], boundary[1], xmin, xmax);
}

inline void
BBDistInternal(double boxLowerA, double boxUpperA, double boxLowerB, double boxUpperB, double &xmin, double &xmax)
{
    xmin = std::max(0.0, std::max(boxLowerB - boxUpperA, boxLowerA - boxUpperB));
    xmax = std::max(std::abs(boxUpperB - boxLowerA), std::abs(boxUpperA - boxLowerB));
    xmin = xmin * xmin;
    xmax = xmax * xmax;
}

template <int dim>
inline void
HighLevelBBDistInternal(
    const std::array<double, 2 * dim> &bba, const std::array<double, 2 * dim> &bbb, double &xmin, double &xmax
)
{
    double xminTmp, xmaxTmp;
    BBDistInternal(bba[0], bba[1], bbb[0], bbb[1], xminTmp, xmaxTmp);
    std::array<double, 2 * (dim - 1)> bba2, bbb2;
    for (int i = 0; i < dim - 1; ++i)
    {
        bba2[2 * i] = bba[2 * (i + 1)];
        bba2[2 * i + 1] = bba[2 * (i + 1) + 1];
        bbb2[2 * i] = bbb[2 * (i + 1)];
        bbb2[2 * i + 1] = bbb[2 * (i + 1) + 1];
    }
    double xminTmp2, xmaxTmp2;
    HighLevelBBDistInternal<dim - 1>(bba2, bbb2, xminTmp2, xmaxTmp2);
    xmin = xminTmp + xminTmp2;
    xmax = xmaxTmp + xmaxTmp2;
}

template <>
inline void
HighLevelBBDistInternal<1>(
    const std::array<double, 2> &bba, const std::array<double, 2> &bbb, double &xmin, double &xmax
)
{
    BBDistInternal(bba[0], bba[1], bbb[0], bbb[1], xmin, xmax);
}
}  // namespace INTERP_KERNEL
// fmt: on
