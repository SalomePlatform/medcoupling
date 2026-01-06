// Copyright (C) 2007-2026  CEA, EDF
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

#include "INTERPKERNELDefines.hxx"
#include <iostream>

namespace INTERP_KERNEL
{

/**
 * \brief Class representing the bounding box of a number of points.
 *
 */
template <int SPACEDIM>
class INTERPKERNEL_EXPORT BoundingBoxT
{
   public:
    /// Enumeration representing the six coordinates that define the bounding box
    enum BoxCoord
    {
        XMIN = 0,
        YMIN = 1,
        ZMIN = 2,
        XMAX = 3,
        YMAX = 4,
        ZMAX = 5
    };

    BoundingBoxT() = default;

    BoundingBoxT(const double **pts, const unsigned numPts);

    BoundingBoxT(const BoundingBoxT<SPACEDIM> &box1, const BoundingBoxT<SPACEDIM> &box2);

    ~BoundingBoxT() = default;

    void fillInXMinXmaxYminYmaxZminZmaxFormat(double data[2 * SPACEDIM]) const;

    void initializeWith(const double **pts, const unsigned numPts);

    bool isDisjointWith(const BoundingBoxT<SPACEDIM> &box) const;

    inline void setCoordinate(int coord, double value);

    inline double getCoordinate(int coord) const;

    void updateWithPoint(const double *pt);

    inline void dumpCoords() const;

    void toCompactData(double data[2 * SPACEDIM]) const;

    BoundingBoxT<SPACEDIM> &operator=(const BoundingBoxT<SPACEDIM> &box) = delete;

   private:
    bool isValid() const;

    /// disallow copying
    BoundingBoxT(const BoundingBoxT &box);

    /// Vector containing the coordinates of the box
    /// interlaced in the order XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX
    double _coords[2 * SPACEDIM];
};

/**
 * Sets a coordinate of the box to a given value.
 *
 * @param coord coordinate to set
 * @param value new value for coordinate
 *
 */
template <int SPACEDIM>
inline void
BoundingBoxT<SPACEDIM>::setCoordinate(int coord, double value)
{
    _coords[coord] = value;
}

/**
 * Gets a coordinate of the box
 *
 * @param coord coordinate to get
 * @return value of coordinate
 *
 */
template <int SPACEDIM>
inline double
BoundingBoxT<SPACEDIM>::getCoordinate(int coord) const
{
    return _coords[coord];
}

/**
 * Prints the coordinates of the box to std::cout
 *
 */
template <int SPACEDIM>
inline void
BoundingBoxT<SPACEDIM>::dumpCoords() const
{
    if (SPACEDIM >= 1)
    {
        std::cout << "[xmin, xmax] = [" << _coords[XMIN] << ", " << _coords[XMIN + SPACEDIM] << "]" << " | ";
    }
    if (SPACEDIM >= 2)
    {
        std::cout << "[ymin, ymax] = [" << _coords[YMIN] << ", " << _coords[YMIN + SPACEDIM] << "]" << " | ";
    }
    if (SPACEDIM >= 3)
    {
        std::cout << "[zmin, zmax] = [" << _coords[ZMIN] << ", " << _coords[ZMIN + SPACEDIM] << "]";
    }
    std::cout << std::endl;
}

using BoundingBox = BoundingBoxT<3>;

}  // namespace INTERP_KERNEL
