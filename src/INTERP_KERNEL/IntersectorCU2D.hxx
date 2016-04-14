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

// File      : IntersectorCU2D.hxx
// Created   : Thu Dec 17 14:10:00 2009
// Author    : Edward AGAPOV (eap)
//
#ifndef __IntersectorCU2D_HXX__
#define __IntersectorCU2D_HXX__

#include "IntersectorCU.hxx"

namespace INTERP_KERNEL
{
  template<class MyCMeshType, class MyUMeshType, class MyMatrix >
  class IntersectorCU2D : public IntersectorCU<MyCMeshType,MyUMeshType,MyMatrix,IntersectorCU2D<MyCMeshType,MyUMeshType,MyMatrix> >
  {
  public:
    typedef typename MyUMeshType::MyConnType UConnType;
    typedef typename MyCMeshType::MyConnType CConnType;
  public:
    IntersectorCU2D(const MyCMeshType& meshS, const MyUMeshType& meshT);
    double intersectGeometry(UConnType icellT, const std::vector<CConnType>& icellC);

  private:
    TriangulationIntersector<MyUMeshType,MyMatrix,PlanarIntersectorP0P0> _intersector;
  };
}


#endif
