// Copyright (C) 2007-2025  CEA, EDF
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
// Author : Anthony Geay (CEA/DEN)

#ifndef __TARGETINTERSECTOR__TXX__
#define __TARGETINTERSECTOR__TXX__

#include "TargetIntersector.hxx"
#include <limits>

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  void TargetIntersector<MyMeshType,MyMatrix>::adjustBoundingBoxes(std::vector<double>& bbox, double adjustmentEps, double adjustmentEpsAbs)
  {
    this->adjustBoundingBoxes(bbox.data(),bbox.size(),adjustmentEps,adjustmentEpsAbs);
  }

  /*! Readjusts a set of bounding boxes so that they are extended
    in all dimensions for avoiding missing interesting intersections

    @param bbox vector containing the bounding boxes
    @param adjustmentEps relative adjustment value (a percentage of the maximal BBox dimension)
    @param adjustmentEpsAbs absolute adjustment value (added on each side of the BBox in each dimension)
  */
  template<class MyMeshType, class MyMatrix>
  void TargetIntersector<MyMeshType,MyMatrix>::adjustBoundingBoxes(double *bbox, std::size_t sz, double adjustmentEps, double adjustmentEpsAbs)
  {
    std::size_t size = sz/(2*SPACEDIM);
    for (std::size_t i=0; i<size; i++)
      {
        double max=- std::numeric_limits<double>::max();
        for(int idim=0; idim<SPACEDIM; idim++)
          {
            double Dx=bbox[i*2*SPACEDIM+1+2*idim]-bbox[i*2*SPACEDIM+2*idim];
            max=(max<Dx)?Dx:max;
          }
        for(int idim=0; idim<SPACEDIM; idim++)
          {
            bbox[i*2*SPACEDIM+2*idim  ] -= adjustmentEps*max+adjustmentEpsAbs;
            bbox[i*2*SPACEDIM+2*idim+1] += adjustmentEps*max+adjustmentEpsAbs;
          }
      }
  }

}

#endif
