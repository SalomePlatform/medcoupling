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
#ifndef __GEOMETRIC2DINTERSECTOR_TXX__
#define __GEOMETRIC2DINTERSECTOR_TXX__

#include "Geometric2DIntersector.hxx"
#include "QuadraticPolygon.hxx"
#include "EdgeArcCircle.hxx"
#include "EdgeLin.hxx"
#include "Node.hxx"

namespace INTERP_KERNEL
{
  /*namespace QUADRATIC_PLANAR
  {
    extern double _precision;
  }*/

  template<class MyMeshType>
  Geometric2DIntersector<MyMeshType>::Geometric2DIntersector(const MyMeshType& mesh_A, const MyMeshType& mesh_B,
                                                             double dimCaracteristic, double precision):
    PlanarIntersector<MyMeshType>(dimCaracteristic, precision, 0., false, 0),_meshA(mesh_A),_meshB(mesh_B)
  {
    _connectA= mesh_A.getConnectivityPtr();
    _connectB= mesh_B.getConnectivityPtr();
    _connIndexA= mesh_A.getConnectivityIndexPtr();
    _connIndexB= mesh_B.getConnectivityIndexPtr();
    _coordsA = mesh_A.getCoordinatesPtr();
    _coordsB = mesh_B.getCoordinatesPtr();
    QUADRATIC_PLANAR::_precision=dimCaracteristic*precision;
  }
  
  template<class MyMeshType>
  double Geometric2DIntersector<MyMeshType>::intersectCells(ConnType icell_A, ConnType icell_B, 
                                                            int nb_NodesA, int nb_NodesB)
  {
    NormalizedCellType tA=_meshA.getTypeOfElement(icell_A);
    NormalizedCellType tB=_meshA.getTypeOfElement(icell_B);
    QuadraticPolygon *p1=buildPolygonAFrom(icell_A,nb_NodesA,tA);
    QuadraticPolygon *p2=buildPolygonBFrom(icell_B,nb_NodesB,tB);
    double ret=p1->intersectWith(*p2);
    delete p1; delete p2;
    return ret;
  }

  template<class MyMeshType>
  QuadraticPolygon *Geometric2DIntersector<MyMeshType>::buildPolygonAFrom(ConnType cell, int nbOfPoints, NormalizedCellType type)
  {
    const ConnType *startOfCellNodeConn=_connectA+OTT<ConnType,numPol>::conn2C(_connIndexA[OTT<ConnType,numPol>::ind2C(cell)]);
    std::vector<Node *> nodes(nbOfPoints);
    for(int i=0;i<nbOfPoints;i++)
      nodes[i]=new Node(_coordsA+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[i])*SPACEDIM);
    if(type!=NORM_TRI6 && type!=NORM_QUAD8)
      return QuadraticPolygon::buildLinearPolygon(nodes);
    else
      return QuadraticPolygon::buildArcCirclePolygon(nodes);
  }

  template<class MyMeshType>
  QuadraticPolygon *Geometric2DIntersector<MyMeshType>::buildPolygonBFrom(ConnType cell, int nbOfPoints, NormalizedCellType type)
  {
    const ConnType *startOfCellNodeConn=_connectB+OTT<ConnType,numPol>::conn2C(_connIndexB[OTT<ConnType,numPol>::ind2C(cell)]);
    std::vector<Node *> nodes(nbOfPoints);
    for(int i=0;i<nbOfPoints;i++)
      nodes[i]=new Node(_coordsB+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[i])*SPACEDIM);
    if(type!=NORM_TRI6 && type!=NORM_QUAD8)
      return QuadraticPolygon::buildLinearPolygon(nodes);
    else
      return QuadraticPolygon::buildArcCirclePolygon(nodes);
  }
}

#endif
