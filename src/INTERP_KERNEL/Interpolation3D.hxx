//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D, OPEN CASCADE
//
//  Copyright (C) 2003-2007  OPEN CASCADE, EADS/CCR, LIP6, CEA/DEN,
//  CEDRAT, EDF R&D, LEG, PRINCIPIA R&D, BUREAU VERITAS
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef __INTERPOLATION3D_HXX__
#define __INTERPOLATION3D_HXX__

#include "Interpolation.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "IntersectorHexa.hxx"
#include "InterpolationOptions.hxx"
namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT Interpolation3D : public Interpolation<Interpolation3D>
  {
  public:
    Interpolation3D();
    Interpolation3D(const InterpolationOptions& io);
    template<class MatrixType, class MyMeshType>
    void interpolateMeshes(const MyMeshType& srcMesh, const MyMeshType& targetMesh, MatrixType& result);
  private:
    SplittingPolicy _splitting_policy;
  };
}

#endif
