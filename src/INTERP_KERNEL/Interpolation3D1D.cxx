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
// Author : Adrien Bruneton (CEA/DEN)

#include "Interpolation3D1D.hxx"
#include "Interpolation3D1D.txx"

namespace INTERP_KERNEL
{
  Interpolation3D1D::Interpolation3D1D() { }

  Interpolation3D1D::Interpolation3D1D(const InterpolationOptions& io):Interpolation<Interpolation3D1D>(io) { }

  void Interpolation3D1D::adjustBoundingBoxes(double *bbox, std::size_t sz)
  {
    const int SPACE_DIM = 3;
    const double adj = getBoundingBoxAdjustmentAbs();
    const double adjRel = getBoundingBoxAdjustment();

    std::size_t size = sz/(2*SPACE_DIM);
    for (std::size_t i=0; i<size; i++)
      {
        double max=- std::numeric_limits<double>::max();
        for(int idim=0; idim<SPACE_DIM; idim++)
          {
            double Dx=bbox[i*2*SPACE_DIM+1+2*idim]-bbox[i*2*SPACE_DIM+2*idim];
            max=(max<Dx)?Dx:max;
          }
        for(int idim=0; idim<SPACE_DIM; idim++)
          {
            bbox[i*2*SPACE_DIM+2*idim  ] -= adjRel*max+adj;
            bbox[i*2*SPACE_DIM+2*idim+1] += adjRel*max+adj;
          }
      }
  }

  /**
   * Inspired from PlanarIntersector<MyMeshType,MyMatrix>::adjustBoundingBoxes
   */
  void Interpolation3D1D::adjustBoundingBoxes(std::vector<double>& bbox)
  {
    adjustBoundingBoxes(bbox.data(),bbox.size());
  }
}
