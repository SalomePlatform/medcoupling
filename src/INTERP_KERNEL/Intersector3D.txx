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
#ifndef __INTERSECTOR3D_TXX__
#define __INTERSECTOR3D_TXX__

#include "Intersector3D.hxx"

#include <algorithm>

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix>
  Intersector3D<MyMeshType,MyMatrix>::Intersector3D(const MyMeshType& targetMesh, const MyMeshType& srcMesh):_target_mesh(targetMesh),_src_mesh(srcMesh)
  {
  }

  /*!
   * @param icellT in format of MyMeshType.
   */
  template<class MyMeshType, class MyMatrix>
  void Intersector3D<MyMeshType,MyMatrix>::getRealTargetCoordinates(ConnType icellT, std::vector<double>& coordsT) const
  {
    int nbNodesT=_target_mesh.getNumberOfNodesOfElement(icellT);
    coordsT.resize(SPACEDIM*nbNodesT);
    std::vector<double>::iterator iter=coordsT.begin();
    for (ConnType iT=0; iT<nbNodesT; iT++)
      {
        const double *coordsCur=getCoordsOfNode(iT,icellT,_target_mesh);
        iter=std::copy(coordsCur,coordsCur+SPACEDIM,iter);
      }
  }

  /*!
   * @param icellT in format of MyMeshType.
   */
  template<class MyMeshType, class MyMatrix>
  void Intersector3D<MyMeshType,MyMatrix>::getRealSourceCoordinates(ConnType icellS, std::vector<double>& coordsS) const
  {
    int nbNodesS=_src_mesh.getNumberOfNodesOfElement(icellS);
    coordsS.resize(SPACEDIM*nbNodesS);
    std::vector<double>::iterator iter=coordsS.begin();
    for (ConnType iS=0; iS<nbNodesS; iS++)
      {
        const double *coordsCur=getCoordsOfNode(iS,icellS,_src_mesh);
        iter=std::copy(coordsCur,coordsCur+SPACEDIM,iter);
      }
  }

  /*!
   * @param icellT in C format.
   * @return is in format of MyMeshType
   */
  template<class MyMeshType, class MyMatrix>
  const typename MyMeshType::MyConnType *Intersector3D<MyMeshType,MyMatrix>::getStartConnOfTargetCell(ConnType icellT) const
  {
    const ConnType *myConectT=_target_mesh.getConnectivityPtr();
    const ConnType *myConIndexT=_target_mesh.getConnectivityIndexPtr();
    return myConectT+OTT<ConnType,numPol>::conn2C(myConIndexT[icellT]);
  }

  /*!
   * @param icellT in C format.
   * @return is in format of MyMeshType
   */
  template<class MyMeshType, class MyMatrix>
  const typename MyMeshType::MyConnType *Intersector3D<MyMeshType,MyMatrix>::getStartConnOfSourceCell(ConnType icellS) const
  {
    const ConnType *myConectS=_src_mesh.getConnectivityPtr();
    const ConnType *myConIndexS=_src_mesh.getConnectivityIndexPtr();
    return myConectS+OTT<ConnType,numPol>::conn2C(myConIndexS[icellS]);
  }

  /*!
   * @param icellS in format of MyMeshType.
   * @param res ; out param in format of MyMeshType.
   */
  template<class MyMeshType, class MyMatrix>
  void Intersector3D<MyMeshType,MyMatrix>::getConnOfSourceCell(ConnType icellS, typename std::vector<ConnType>& res) const
  {
    const ConnType *myConectS=_src_mesh.getConnectivityPtr();
    const ConnType *myConIndexS=_src_mesh.getConnectivityIndexPtr();
    ConnType start=myConIndexS[OTT<ConnType,numPol>::ind2C(icellS)];
    ConnType end=myConIndexS[OTT<ConnType,numPol>::ind2C(icellS)+1];
    int nbNodesS=end-start;
    res.resize(nbNodesS);
    std::copy(myConectS+OTT<ConnType,numPol>::conn2C(start),myConectS+OTT<ConnType,numPol>::conn2C(end),res.begin());
  }
}

#endif
