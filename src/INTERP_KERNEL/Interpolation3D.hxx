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
// Author : Anthony Geay (CEA/DEN)

#ifndef __INTERPOLATION3D_HXX__
#define __INTERPOLATION3D_HXX__

#include "INTERPKERNELDefines.hxx"
#include "Interpolation.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpolationOptions.hxx"

namespace INTERP_KERNEL
{
  /**
   * \class Interpolation3D
   * \brief Class used to calculate the volumes of intersection between the elements of two 3D meshes.
   */

  class INTERPKERNEL_EXPORT Interpolation3D : public Interpolation<Interpolation3D>
  {
  public:
    Interpolation3D();
    Interpolation3D(const InterpolationOptions& io);
    template<class MyMeshType, class MatrixType>
    typename MyMeshType::MyConnType interpolateMeshes(const MyMeshType& srcMesh, const MyMeshType& targetMesh, MatrixType& result, const std::string& method);
  };
}

#endif
