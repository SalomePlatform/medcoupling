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

#ifndef __INTERSECTOR3D_HXX__
#define __INTERSECTOR3D_HXX__

#include "TargetIntersector.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  class Intersector3D : public TargetIntersector<MyMeshType,MyMatrix>
  {
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
  public:
    Intersector3D(const MyMeshType& targetMesh, const MyMeshType& srcMesh);
    void getRealTargetCoordinates(ConnType icellT, std::vector<double>& coordsT) const;
    void getRealSourceCoordinates(ConnType icellT, std::vector<double>& coordsT) const;
    const ConnType *getStartConnOfTargetCell(ConnType icellT) const;
    const ConnType *getStartConnOfSourceCell(ConnType icellS) const;
    void getConnOfSourceCell(ConnType icellS, typename std::vector<ConnType>& res) const;
  protected:
    const MyMeshType& _target_mesh;
    const MyMeshType& _src_mesh;
  };
}

#endif
