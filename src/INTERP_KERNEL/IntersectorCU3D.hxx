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

// File      : IntersectorCU3D.hxx
// Created   : Thu Dec 17 14:10:00 2009
// Author    : Edward AGAPOV (eap)
//
#ifndef __IntersectorCU3D_HXX__
#define __IntersectorCU3D_HXX__

#include "IntersectorCU.hxx"
#include "SplitterTetra.hxx"

namespace INTERP_KERNEL
{
  class _Cartesian3D2UnstructHexMesh;

  template<class MyCMeshType, class MyUMeshType, class MyMatrix >
  class IntersectorCU3D : public IntersectorCU<MyCMeshType,MyUMeshType,MyMatrix,IntersectorCU3D<MyCMeshType,MyUMeshType,MyMatrix> >
  {
  public:
    typedef typename MyUMeshType::MyConnType UConnType;
    typedef typename MyCMeshType::MyConnType CConnType;
  public:
    IntersectorCU3D(const MyCMeshType& meshS, const MyUMeshType& meshT, SplittingPolicy splitting_policy);
    ~IntersectorCU3D();
    double intersectGeometry(UConnType icellT, const std::vector<CConnType>& icellC);

  private:

    typedef SplitterTetra2<MyUMeshType, _Cartesian3D2UnstructHexMesh > TSplitter;
    typedef SplitterTetra <_Cartesian3D2UnstructHexMesh >              TTetra;
    _Cartesian3D2UnstructHexMesh* _uHexMesh;
    TSplitter*                    _split;
  };
}


#endif
