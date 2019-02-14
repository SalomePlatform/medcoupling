// Copyright (C) 2018-2019  CEA/DEN, EDF R&D
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
// Author : Anthony GEAY (EDF R&D)

#include "Interpolation1D0D.hxx"
#include "Interpolation1D0D.txx"

namespace INTERP_KERNEL
{
  /**
   * \class Interpolation1D0D
   * \brief Class used to calculate the interpolation between a 1D mesh and 0D mesh (in 3D space)
   * Can be seen as a specialization of Interpolation3D, and allows notably the adjustment of bounind boxes.
   * 
   */

  Interpolation1D0D::Interpolation1D0D()
  {}

  Interpolation1D0D::Interpolation1D0D(const InterpolationOptions& io):Interpolation<Interpolation1D0D>(io)
  {}

  /**
   * Inspired from PlanarIntersector<MyMeshType,MyMatrix>::adjustBoundingBoxes
   */
  void Interpolation1D0D::adjustBoundingBoxes(std::vector<double>& bbox)
  {
    const int SPACE_DIM = 3;
    const double adj(getPrecision());// here precision is used instead of getBoundingBoxAdjustment and getBoundingBoxAdjustmentAbs because in the context only precision is relevant

    long size = bbox.size()/(2*SPACE_DIM);
    for (int i=0; i<size; i++)
      {
        for(int idim=0; idim<SPACE_DIM; idim++)
          {
            bbox[i*2*SPACE_DIM+2*idim  ] -= adj;
            bbox[i*2*SPACE_DIM+2*idim+1] += adj;
          }
      }
  }
}
