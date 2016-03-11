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
#ifndef __INTERSECTOR3DP1P0_TXX__
#define __INTERSECTOR3DP1P0_TXX__

#include "Intersector3DP1P0.hxx"
#include "Intersector3D.txx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  Intersector3DP1P0<MyMeshType,MyMatrix>::Intersector3DP1P0(const MyMeshType& targetMesh, const MyMeshType& srcMesh):Intersector3D<MyMeshType,MyMatrix>(targetMesh,srcMesh)
  {
  }

  template<class MyMeshType, class MyMatrix>
  int Intersector3DP1P0<MyMeshType,MyMatrix>::getNumberOfRowsOfResMatrix() const
  {
    return Intersector3D<MyMeshType,MyMatrix>::_target_mesh.getNumberOfElements();
  }
  
  template<class MyMeshType, class MyMatrix>
  int Intersector3DP1P0<MyMeshType,MyMatrix>::getNumberOfColsOfResMatrix() const
  {
    return Intersector3D<MyMeshType,MyMatrix>::_src_mesh.getNumberOfNodes();
  }
}

#endif
