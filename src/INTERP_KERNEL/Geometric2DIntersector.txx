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
#include "PlanarIntersectorP0P0.txx"
#include "PlanarIntersectorP0P1.txx"
#include "PlanarIntersectorP1P0.txx"
#include "PlanarIntersectorP1P1.txx"
#include "CellModel.hxx"

#include "QuadraticPolygon.hxx"
#include "EdgeArcCircle.hxx"
#include "EdgeLin.hxx"
#include "Node.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  Geometric2DIntersector<MyMeshType,MyMatrix,InterpType>::Geometric2DIntersector(const MyMeshType& meshT, const MyMeshType& meshS,
                                                                          double dimCaracteristic, double medianPlane, double precision, int orientation):
    InterpType<MyMeshType,MyMatrix,Geometric2DIntersector<MyMeshType,MyMatrix,InterpType> >(meshT,meshS,dimCaracteristic, precision, medianPlane, true, orientation, 0)
  {
    QUADRATIC_PLANAR::_precision=dimCaracteristic*precision;
  }
  
  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  double Geometric2DIntersector<MyMeshType,MyMatrix,InterpType>::intersectGeometry(ConnType icellT, ConnType icellS, ConnType nbNodesT, ConnType nbNodesS)
  {
    int orientation = 1;
    std::vector<double> CoordsT;
    std::vector<double> CoordsS;
    PlanarIntersector<MyMeshType,MyMatrix>::getRealCoordinates(icellT,icellS,nbNodesT,nbNodesS,CoordsT,CoordsS,orientation);
    NormalizedCellType tT=PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getTypeOfElement(icellT);
    NormalizedCellType tS=PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getTypeOfElement(icellS);
    QuadraticPolygon *p1=buildPolygonFrom(CoordsT,tT);
    QuadraticPolygon *p2=buildPolygonFrom(CoordsS,tS);
    double ret=p1->intersectWith(*p2);
    delete p1; delete p2;
    return ret;
  }

  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  double Geometric2DIntersector<MyMeshType,MyMatrix,InterpType>::intersectGeometryWithQuadrangle(const double *quadrangle, const std::vector<double>& sourceCoords, bool isSourceQuad)
  {
    std::vector<Node *> nodes(4);
    nodes[0]=new Node(quadrangle[0],quadrangle[1]);
    nodes[1]=new Node(quadrangle[SPACEDIM],quadrangle[SPACEDIM+1]);
    nodes[2]=new Node(quadrangle[2*SPACEDIM],quadrangle[2*SPACEDIM+1]);
    nodes[3]=new Node(quadrangle[3*SPACEDIM],quadrangle[3*SPACEDIM+1]);
    int nbOfSourceNodes=sourceCoords.size()/SPACEDIM;
    std::vector<Node *> nodes2(nbOfSourceNodes);
    for(int i=0;i<nbOfSourceNodes;i++)
      nodes2[i]=new Node(sourceCoords[i*SPACEDIM],sourceCoords[i*SPACEDIM+1]);
    QuadraticPolygon *p1=QuadraticPolygon::buildLinearPolygon(nodes);
    QuadraticPolygon *p2;
    if(!isSourceQuad)
      p2=QuadraticPolygon::buildLinearPolygon(nodes2);
    else
      p2=QuadraticPolygon::buildArcCirclePolygon(nodes2);
    double ret=p1->intersectWith(*p2);
    delete p1; delete p2;
    return ret;
  }

  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  double Geometric2DIntersector<MyMeshType,MyMatrix,InterpType>::intersectGeometryGeneral(const std::vector<double>& targetCoords, const std::vector<double>& sourceCoords)
  {
    int nbOfTargetNodes=targetCoords.size()/SPACEDIM;
    std::vector<Node *> nodes(nbOfTargetNodes);
    for(int i=0;i<nbOfTargetNodes;i++)
      nodes[i]=new Node(targetCoords[i*SPACEDIM],targetCoords[i*SPACEDIM+1]);
    int nbOfSourceNodes=sourceCoords.size()/SPACEDIM;
    std::vector<Node *> nodes2(nbOfSourceNodes);
    for(int i=0;i<nbOfSourceNodes;i++)
      nodes2[i]=new Node(sourceCoords[i*SPACEDIM],sourceCoords[i*SPACEDIM+1]);
    QuadraticPolygon *p1=QuadraticPolygon::buildLinearPolygon(nodes);
    QuadraticPolygon *p2=QuadraticPolygon::buildLinearPolygon(nodes2);
    double ret=p1->intersectWith(*p2);
    delete p1; delete p2;
    return ret;
  }

  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  QuadraticPolygon *Geometric2DIntersector<MyMeshType,MyMatrix,InterpType>::buildPolygonFrom(const std::vector<double>& coords, NormalizedCellType type)
  {
    int nbNodes=coords.size()/SPACEDIM;
    std::vector<Node *> nodes(nbNodes);
    for(int i=0;i<nbNodes;i++)
      nodes[i]=new Node(coords[i*SPACEDIM],coords[i*SPACEDIM+1]);
    if(!CellModel::getCellModel(type).isQuadratic())
      return QuadraticPolygon::buildLinearPolygon(nodes);
    else
      return QuadraticPolygon::buildArcCirclePolygon(nodes);
  }

  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  QuadraticPolygon *Geometric2DIntersector<MyMeshType,MyMatrix,InterpType>::buildPolygonAFrom(ConnType cell, int nbOfPoints, NormalizedCellType type)
  {
    const ConnType *startOfCellNodeConn=PlanarIntersector<MyMeshType,MyMatrix>::_connectT+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexT[OTT<ConnType,numPol>::ind2C(cell)]);
    std::vector<Node *> nodes(nbOfPoints);
    for(int i=0;i<nbOfPoints;i++)
      nodes[i]=new Node(PlanarIntersector<MyMeshType,MyMatrix>::_coordsT+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[i])*SPACEDIM);
    if(CellModel::getCellModel(type).isQuadratic())
      return QuadraticPolygon::buildLinearPolygon(nodes);
    else
      return QuadraticPolygon::buildArcCirclePolygon(nodes);
  }

  template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>
  QuadraticPolygon *Geometric2DIntersector<MyMeshType,MyMatrix,InterpType>::buildPolygonBFrom(ConnType cell, int nbOfPoints, NormalizedCellType type)
  {
    const ConnType *startOfCellNodeConn=PlanarIntersector<MyMeshType,MyMatrix>::_connectS+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexS[OTT<ConnType,numPol>::ind2C(cell)]);
    std::vector<Node *> nodes(nbOfPoints);
    for(int i=0;i<nbOfPoints;i++)
      nodes[i]=new Node(PlanarIntersector<MyMeshType,MyMatrix>::_coordsS+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[i])*SPACEDIM);
    if(type!=NORM_TRI6 && type!=NORM_QUAD8)
      return QuadraticPolygon::buildLinearPolygon(nodes);
    else
      return QuadraticPolygon::buildArcCirclePolygon(nodes);
  }
}

#endif
