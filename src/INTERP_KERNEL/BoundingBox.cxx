// Copyright (C) 2007-2024  CEA, EDF
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

#include "BoundingBox.hxx"

#include <iostream>
#include <algorithm>
#include <cassert>

namespace INTERP_KERNEL
{
  
  /**
   * Constructor creating box from an array of the points corresponding
   * to the vertices of the element.
   * Each point is represented by an array of three doubles.
   *
   * @param pts     array of points 
   * @param numPts  number of vertices
   *
   */
  BoundingBox::BoundingBox(const double** pts, const unsigned numPts)
  {
    initializeWith(pts,numPts);
  }

  void BoundingBox::fillInXMinXmaxYminYmaxZminZmaxFormat(double data[6]) const
  {
    data[0] = this->getCoordinate(BoundingBox::XMIN);
    data[1] = this->getCoordinate(BoundingBox::XMAX);
    data[2] = this->getCoordinate(BoundingBox::YMIN);
    data[3] = this->getCoordinate(BoundingBox::YMAX);
    data[4] = this->getCoordinate(BoundingBox::ZMIN);
    data[5] = this->getCoordinate(BoundingBox::ZMAX);
  }

  /**
   * Constructor creating box from an array of the points corresponding
   * to the vertices of the element.
   * Each point is represented by an array of three doubles.
   *
   * @param pts     array of points 
   * @param numPts  number of vertices
   *
   */
  void BoundingBox::initializeWith(const double** pts, const unsigned numPts)
  {
    // initialize with first two points
    const double *pt0(pts[0]);

    for(BoxCoord c = XMIN ; c <= ZMIN ; c = BoxCoord(c + 1))
      {
        _coords[c] = pt0[c];
        _coords[c + 3] = pt0[c];
      }

    for(unsigned i = 1 ; i < numPts ; ++i)
      {
        updateWithPoint(pts[i]);
      }
  
    assert(isValid());
  }

  /**
   * Constructor creating box from union of two boxes, resulting in a box that encloses both of them
   *
   * @param  box1  the first box
   * @param  box2  the second box
   */
  BoundingBox::BoundingBox(const BoundingBox& box1, const BoundingBox& box2)
  {
    for(BoxCoord c = XMIN ; c <= ZMIN ; c = BoxCoord(c + 1))
      {
        _coords[c] = std::min(box1._coords[c], box2._coords[c]);
        _coords[c + 3] = std::max(box1._coords[c + 3], box2._coords[c + 3]);
      }
    
    assert(isValid());
  }

  /**
   * Determines if the intersection with a given box is empty
   * 
   * @param    box   BoundingBox with which intersection is tested
   * @return  true if intersection between boxes is empty, false if not
   */
  bool BoundingBox::isDisjointWith(const BoundingBox& box) const
  {
    for(BoxCoord c = XMIN ; c <= ZMIN ; c = BoxCoord(c + 1))
      {
        const double otherMinCoord = box.getCoordinate(c);
        const double otherMaxCoord = box.getCoordinate(BoxCoord(c + 3));
       
        // boxes are disjoint if there exists a direction in which the 
        // minimum coordinate of one is greater than the maximum coordinate of the other

        // more stable version ?
        // const double tol = 1.0e-2*_coords[c];
        // if(_coords[c] > otherMaxCoord + tol 
        //   || _coords[c + 3] < otherMinCoord - tol)
       
       
        if(_coords[c] > otherMaxCoord 
           || _coords[c + 3] < otherMinCoord)
       
          {
            return true;
          }
       
      }
    return false;
  }
    
  

  /**
   * Updates the bounding box to include a given point
   * 
   * @param pt    point to be included
   *
   */
  void BoundingBox::updateWithPoint(const double* pt)
  {
    for(BoxCoord c = XMIN ; c <= ZMIN ; c = BoxCoord(c + 1))
      {
        const double ptVal = pt[c];

        // update min and max coordinates
        _coords[c] = std::min(_coords[c], ptVal);
        _coords[c + 3] = std::max(_coords[c + 3], ptVal);

      }
  }
  
  /**
   * Checks if the box is valid, which it is if its minimum coordinates are
   * smaller than its maximum coordinates in all directions.
   *
   * @return  true if the box is valid, false if not
   */
  bool BoundingBox::isValid() const
  {
    bool valid = true;
    for(BoxCoord c = XMIN ; c < ZMIN ; c = BoxCoord(c + 1))
      {
        if(_coords[c] > _coords[c + 3])
          {
            std::cout << "+++ Error in  BoundingBox |: coordinate " << c << " is invalid : "
                      <<_coords[c] << " > " << _coords[c+3] << std::endl;
            valid = false;
          }
      }
    return valid;
  }

  void BoundingBox::toCompactData(double data[6]) const
  {
    data[0]=_coords[XMIN];
    data[1]=_coords[XMAX];
    data[2]=_coords[YMIN];
    data[3]=_coords[YMAX];
    data[4]=_coords[ZMIN];
    data[5]=_coords[ZMAX];
  }

}
