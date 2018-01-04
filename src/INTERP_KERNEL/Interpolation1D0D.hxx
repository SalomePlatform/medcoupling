// Copyright (C) 2018  CEA/DEN, EDF R&D
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
// Author : A Geay (EDF R&D)

#ifndef __INTERPOLATION1D0D_HXX__
#define __INTERPOLATION1D0D_HXX__

#include "INTERPKERNELDefines.hxx"
#include "Interpolation.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpolationOptions.hxx"

#include <vector>

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT Interpolation1D0D : public Interpolation<Interpolation1D0D>
  {
  public:
    Interpolation1D0D();
    Interpolation1D0D(const InterpolationOptions& io);
    template<class MyMeshType, class MatrixType>
    int interpolateMeshes(const MyMeshType& srcMesh, const MyMeshType& targetMesh, MatrixType& result, const std::string& method);
  private:
    void adjustBoundingBoxes(std::vector<double>& bbox);
  };
}

#endif
