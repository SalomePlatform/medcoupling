// Copyright (C) 2009-2024  OPEN CASCADE
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

// File      : InterpolationCU.hxx
// Created   : Mon Dec 14 16:52:53 2009
// Author    : Edward AGAPOV (eap)
//
#ifndef __InterpolationCU_HXX__
#define __InterpolationCU_HXX__

#include "Interpolation.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelUtilities.hxx"

namespace INTERP_KERNEL
{
  class InterpolationCU : public Interpolation< InterpolationCU >
  {
  public:
    InterpolationCU();
    InterpolationCU(const InterpolationOptions & io);

    template<class MyCMeshType, class MyUMeshType, class MatrixType>
    typename MyCMeshType::MyConnType interpolateMeshes(const MyCMeshType& meshS, const MyUMeshType& meshT, MatrixType& result, const char *method);

    template<class MyUMeshType, class MyCMeshType, class MatrixType>
    typename MyUMeshType::MyConnType interpolateMeshesRev(const MyUMeshType& meshS, const MyCMeshType& meshT, MatrixType& result, const char *method);

  };
}

#endif
