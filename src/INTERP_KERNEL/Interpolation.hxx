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
// Author : Anthony Geay (CEA/DEN)

#ifndef __INTERPOLATION_HXX__
#define __INTERPOLATION_HXX__

#include "INTERPKERNELDefines.hxx"
#include "InterpolationOptions.hxx"
#include "InterpKernelException.hxx"

#include <string>

namespace INTERP_KERNEL
{
  template<class TrueMainInterpolator>
  class Interpolation : public InterpolationOptions
  {
  public:
    Interpolation() { }
    Interpolation(const InterpolationOptions& io) :InterpolationOptions(io){}
    //interpolation of two triangular meshes.
    template<class MatrixType, class MyMeshType>
    int interpolateMeshes(const MyMeshType& meshS, const MyMeshType& meshT, MatrixType& result)
    { return asLeaf().interpolateMeshes(meshS,meshT,result); }
    template<class MyMeshType, class MatrixType>
    int fromIntegralUniform(const MyMeshType& meshT, MatrixType& result, const std::string& method) { return fromToIntegralUniform(false,meshT,result,method); }
    template<class MyMeshType, class MatrixType>
    int toIntegralUniform(const MyMeshType& meshS, MatrixType& result, const std::string& method) { return fromToIntegralUniform(true,meshS,result,method); }
    template<class MyMeshType>
    static double CalculateCharacteristicSizeOfMeshes(const MyMeshType& myMeshS, const MyMeshType& myMeshT, const int printLevel);
  protected:
    template<class MyMeshType, class MatrixType>
    int fromToIntegralUniform(bool fromTo, const MyMeshType& mesh, MatrixType& result, const std::string& method);
  protected:
    TrueMainInterpolator& asLeaf() { return static_cast<TrueMainInterpolator&>(*this); }
  };
}

#endif
