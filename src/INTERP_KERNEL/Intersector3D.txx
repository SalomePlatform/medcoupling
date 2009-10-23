//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
    const ConnType *myConectT=_target_mesh.getConnectivityPtr();
    const ConnType *myConIndexT=_target_mesh.getConnectivityIndexPtr();
    const double *myCoordsT=_target_mesh.getCoordinatesPtr();
    int nbNodesT=myConIndexT[OTT<ConnType,numPol>::ind2C(icellT)+1]-myConIndexT[OTT<ConnType,numPol>::ind2C(icellT)];
    coordsT.resize(SPACEDIM*nbNodesT);
    for (ConnType iT=0; iT<nbNodesT; iT++)
      for(int idim=0; idim<SPACEDIM; idim++)
        coordsT[SPACEDIM*iT+idim]=myCoordsT[SPACEDIM*OTT<ConnType,numPol>::coo2C(myConectT[OTT<ConnType,numPol>::conn2C(myConIndexT[OTT<ConnType,numPol>::ind2C(icellT)]+iT)])+idim];
  }

  /*!
   * @param icellT in format of MyMeshType.
   */
  template<class MyMeshType, class MyMatrix>
  void Intersector3D<MyMeshType,MyMatrix>::getRealSourceCoordinates(ConnType icellS, std::vector<double>& coordsS) const
  {
    const ConnType *myConectS=_src_mesh.getConnectivityPtr();
    const ConnType *myConIndexS=_src_mesh.getConnectivityIndexPtr();
    const double *myCoordsS=_src_mesh.getCoordinatesPtr();
    int nbNodesS=myConIndexS[OTT<ConnType,numPol>::ind2C(icellS)+1]-myConIndexS[OTT<ConnType,numPol>::ind2C(icellS)];
    coordsS.resize(SPACEDIM*nbNodesS);
    for (ConnType iS=0; iS<nbNodesS; iS++)
      for(int idim=0; idim<SPACEDIM; idim++)
        coordsS[SPACEDIM*iS+idim]=myCoordsS[SPACEDIM*OTT<ConnType,numPol>::coo2C(myConectS[OTT<ConnType,numPol>::conn2C(myConIndexS[OTT<ConnType,numPol>::ind2C(icellS)]+iS)])+idim];
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
