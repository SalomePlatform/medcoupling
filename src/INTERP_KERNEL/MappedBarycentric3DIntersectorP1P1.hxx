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
// Author : Adrien Bruneton (CEA/DEN)

#ifndef __MappedBarycentric3DIntersectorP1P1_HXX__
#define __MappedBarycentric3DIntersectorP1P1_HXX__

#include "Intersector3DP1P1.hxx"
#include "NormalizedUnstructuredMesh.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  class MappedBarycentric3DIntersectorP1P1 : public Intersector3DP1P1<MyMeshType,MyMatrix>
  { 
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
  public:
    MappedBarycentric3DIntersectorP1P1(const MyMeshType& targetMesh, const MyMeshType& srcMesh, double precision);
    ~MappedBarycentric3DIntersectorP1P1();
    void intersectCells(ConnType targetCell, const std::vector<ConnType>& srcCells, MyMatrix& res);
  protected:
    double _precision;
  };
}

#endif
