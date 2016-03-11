// Copyright (C) 2009-2016  OPEN CASCADE
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

#ifndef __DIRECTEDBOUNDINGBOX_HXX__
#define __DIRECTEDBOUNDINGBOX_HXX__

#include "INTERPKERNELDefines.hxx"

#include <vector>

namespace INTERP_KERNEL
{

  /**
   * \brief Class representing the bounding box of a number of points
   *  with box axes parallel to principal axes of inertia of points
   */
  class DirectedBoundingBox
  {
  public:

    INTERPKERNEL_EXPORT DirectedBoundingBox();

    INTERPKERNEL_EXPORT DirectedBoundingBox(const double* pts, const unsigned numPts, const unsigned dim);

    INTERPKERNEL_EXPORT DirectedBoundingBox(const double** pts, const unsigned numPts, const unsigned dim);

    //~DirectedBoundingBox();

    INTERPKERNEL_EXPORT void enlarge(const double tol);
    
    INTERPKERNEL_EXPORT bool isDisjointWith(const DirectedBoundingBox& box) const;

    INTERPKERNEL_EXPORT bool isDisjointWith(const double* box) const;

    INTERPKERNEL_EXPORT bool isOut(const double* point) const;


    // return internal data
    INTERPKERNEL_EXPORT std::vector<double> getData() const;

    // initialize with data returned by getData()
    INTERPKERNEL_EXPORT void setData(const double* data);

    // return size of internal data
    INTERPKERNEL_EXPORT static int dataSize(int dim);

  private:

    //void computeAxes3D(const std::vector<double>& tensor);

    //void computeAxes2D(const std::vector<double>& tensor);

    inline void addPointToBox(const double* coord);

    void toLocalCS(const double* p, double* pLoc) const;

    void fromLocalCS(const double* p, double* pGlob) const;

    inline bool isLocalOut(const double* pLoc) const;

    void getCorners(std::vector<double>& corners, const double* minmax) const;

    unsigned _dim;

    std::vector<double> _axes; //!< principal axes of inertia in full interlace
    std::vector<double> _minmax; //!< pairs of min an max coordinates along the axes

  };

  //================================================================================
  /*!
   * \brief Test point in local CS against box extremities
   * 
   */
  //================================================================================

  inline bool DirectedBoundingBox::isLocalOut(const double* pLoc) const
    {
      for ( int i = 0; i < (int)_dim; ++i )
        if ( pLoc[i] < _minmax[i*2] || pLoc[i] > _minmax[i*2+1] )
          return true;
      return false;
    }

  //================================================================================
  /*!
   * \brief Update box extremities
   */
  //================================================================================

  inline void DirectedBoundingBox::addPointToBox(const double* coord)
  {
    for ( int i = 0; i < (int)_dim; ++i )
      {
        double c = 0;
        for ( int j = 0; j < (int)_dim; ++j ) c += coord[j]*_axes[i*_dim+j];
        if ( c < _minmax[i*2] )   _minmax[i*2] = c;
        if ( c > _minmax[i*2+1] ) _minmax[i*2+1] = c;
      }
  }
}
#endif
