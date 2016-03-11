// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

#ifndef __BOUNDINGBOX_HXX__
#define __BOUNDINGBOX_HXX__

#include "INTERPKERNELDefines.hxx"
#include <iostream>

namespace INTERP_KERNEL
{

  /**
   * \brief Class representing the bounding box of a number of points.
   *
   */
  class INTERPKERNEL_EXPORT BoundingBox 
  {
  public:

    /// Enumeration representing the six coordinates that define the bounding box
    enum BoxCoord { XMIN = 0, YMIN = 1, ZMIN = 2, XMAX = 3, YMAX = 4, ZMAX = 5 };
        
    BoundingBox(const double** pts, const unsigned numPts);

    BoundingBox(const BoundingBox& box1, const BoundingBox& box2);

    ~BoundingBox();

    bool isDisjointWith(const BoundingBox& box) const;
    
    inline void setCoordinate(const BoxCoord coord, double value);

    inline double getCoordinate(const BoxCoord coord) const;

    void updateWithPoint(const double* pt);

    inline void dumpCoords() const;

  private:
    
    bool isValid() const;

    /// disallow copying
    BoundingBox(const BoundingBox& box);
    
    /// disallow assignment
    BoundingBox& operator=(const BoundingBox& box);
    
    /// Vector containing the coordinates of the box
    /// interlaced in the order XMIN, YMIN, ZMIN, XMAX, YMAX, ZMAX
    double* _coords;

  };

  /**
   * Sets a coordinate of the box to a given value.
   * 
   * @param coord coordinate to set
   * @param value new value for coordinate
   *
   */
  inline void BoundingBox::setCoordinate(const BoxCoord coord, double value)
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
  inline double BoundingBox::getCoordinate(const BoxCoord coord) const
  {
    return _coords[coord];
  }

  /**
   * Prints the coordinates of the box to std::cout
   *
   */
  inline void BoundingBox::dumpCoords() const
  {
    std::cout << "[xmin, xmax] = [" << _coords[XMIN] << ", " << _coords[XMAX] << "]" << " | ";
    std::cout << "[ymin, ymax] = [" << _coords[YMIN] << ", " << _coords[YMAX] << "]" << " | ";
    std::cout << "[zmin, zmax] = [" << _coords[ZMIN] << ", " << _coords[ZMAX] << "]";
    std::cout << std::endl;
  }

}

#endif
