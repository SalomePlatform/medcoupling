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
#ifndef __POINTLOCATORINTERSECTOR_TXX__
#define __POINTLOCATORINTERSECTOR_TXX__

#include "PointLocator2DIntersector.hxx"
#include "PlanarIntersectorP0P0.txx"
#include "PlanarIntersectorP0P1.txx"
#include "PlanarIntersectorP1P0.txx"
#include "PlanarIntersectorP1P1.txx"
#include "PlanarIntersectorP1P0Bary.txx"
#include "CellModel.hxx"

#include "InterpKernelGeo2DQuadraticPolygon.hxx"
#include "PointLocatorAlgos.txx"

#define PTLOC2D_INTERSECTOR PointLocator2DIntersector<MyMeshType,MyMatrix,InterpType>
#define INTERSECTOR_TEMPLATE template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>

namespace INTERP_KERNEL
{
  INTERSECTOR_TEMPLATE
  PTLOC2D_INTERSECTOR::PointLocator2DIntersector(const MyMeshType& meshT, const MyMeshType& meshS,
                                               double dimCaracteristic, double md3DSurf, double minDot3DSurf, double medianPlane,
                                               double precision, int orientation):
    InterpType<MyMeshType,MyMatrix,PTLOC2D_INTERSECTOR >(meshT,meshS,dimCaracteristic, precision, md3DSurf, minDot3DSurf, medianPlane, true, orientation, 0)
  {
  }
  
  INTERSECTOR_TEMPLATE
  double PTLOC2D_INTERSECTOR::intersectGeometry(ConnType icellT,   ConnType icellS,
                                                ConnType nbNodesT, ConnType nbNodesS)
  {
    int orientation = 1;
    std::vector<double> CoordsT;
    std::vector<double> CoordsS;
    PlanarIntersector<MyMeshType,MyMatrix>::getRealCoordinates(icellT,icellS,nbNodesT,nbNodesS,CoordsT,CoordsS,orientation);
    NormalizedCellType tT=PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getTypeOfElement(icellT);
    QuadraticPolygon *pT=buildPolygonFrom(CoordsT,tT);
    double baryT[SPACEDIM];
    pT->getBarycenterGeneral(baryT);
    delete pT;
    if(PointLocatorAlgos<MyMeshType>::isElementContainsPointAlg2D(baryT,&CoordsS[0],nbNodesS,InterpType<MyMeshType,MyMatrix,PTLOC2D_INTERSECTOR >::_precision))
      return 1.;
    return 0.;
  }

  INTERSECTOR_TEMPLATE
  double PTLOC2D_INTERSECTOR::intersectGeometryWithQuadrangle(const double             * quadrangle,
                                                              const std::vector<double>& sourceCoords,
                                                              bool                       isSourceQuad)
  {
    int nbOfSourceNodes=sourceCoords.size()/SPACEDIM;
    std::vector<Node *> nodes2(nbOfSourceNodes);
    for(int i=0;i<nbOfSourceNodes;i++)
      nodes2[i]=new Node(sourceCoords[i*SPACEDIM],sourceCoords[i*SPACEDIM+1]);
    QuadraticPolygon *p2;
    if(!isSourceQuad)
      p2=QuadraticPolygon::BuildLinearPolygon(nodes2);
    else
      p2=QuadraticPolygon::BuildArcCirclePolygon(nodes2);
    double bary[SPACEDIM];
    p2->getBarycenter(bary);
    delete p2;
    if( PointLocatorAlgos<MyMeshType>::isElementContainsPointAlg2D(bary,quadrangle,4) )
      return 1.;
    return 0.;
  }

  INTERSECTOR_TEMPLATE
  double PTLOC2D_INTERSECTOR::intersectGeometryGeneral(const std::vector<double>& targetCoords,
                                                       const std::vector<double>& sourceCoords)
  {
    int nbOfTargetNodes=targetCoords.size()/SPACEDIM;
    int nbOfSourceNodes=sourceCoords.size()/SPACEDIM;
    std::vector<Node *> nodes2(nbOfSourceNodes);
    for(int i=0;i<nbOfSourceNodes;i++)
      nodes2[i]=new Node(sourceCoords[i*SPACEDIM],sourceCoords[i*SPACEDIM+1]);
    QuadraticPolygon *p=QuadraticPolygon::BuildLinearPolygon(nodes2);
    double bary[SPACEDIM];
    p->getBarycenterGeneral(bary);
    delete p;
    if( PointLocatorAlgos<MyMeshType>::isElementContainsPointAlg2D(bary,&targetCoords[0],nbOfTargetNodes) )
      return 1.;
    return 0.;
  }

  //================================================================================
  /*!
   * \brief Intersect a triangle and a polygon for P1P0 barycentric algorithm
   *  \param targetCell - list of coordinates of target polygon in full interlace
   *  \param targetCellQuadratic - specifies if target polygon is quadratic or not
   *  \param sourceTria - list of coordinates of source triangle
   *  \param res - coefficients a,b and c associated to nodes of sourceTria
   */
  //================================================================================

  INTERSECTOR_TEMPLATE
  double PTLOC2D_INTERSECTOR::intersectGeoBary(const std::vector<double>& targetCell,
                                               bool                       targetCellQuadratic,
                                               const double *             sourceTria,
                                               std::vector<double>&       res)
  {
    throw INTERP_KERNEL::Exception("intersectGeoBary incompatible with PointLocator. Desactivate P1P0Bary to avoid the problem");
    return 0.;
  }

  INTERSECTOR_TEMPLATE
  QuadraticPolygon *PTLOC2D_INTERSECTOR::buildPolygonFrom(const std::vector<double>& coords, NormalizedCellType type)
  {
    int nbNodes=coords.size()/SPACEDIM;
    std::vector<Node *> nodes(nbNodes);
    for(int i=0;i<nbNodes;i++)
      nodes[i]=new Node(coords[i*SPACEDIM],coords[i*SPACEDIM+1]);
    if(!CellModel::GetCellModel(type).isQuadratic())
      return QuadraticPolygon::BuildLinearPolygon(nodes);
    else
      return QuadraticPolygon::BuildArcCirclePolygon(nodes);
  }

  INTERSECTOR_TEMPLATE
  QuadraticPolygon *PTLOC2D_INTERSECTOR::buildPolygonAFrom(ConnType cell, int nbOfPoints, NormalizedCellType type)
  {
    const ConnType *startOfCellNodeConn=PlanarIntersector<MyMeshType,MyMatrix>::_connectT+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexT[OTT<ConnType,numPol>::ind2C(cell)]);
    std::vector<Node *> nodes(nbOfPoints);
    for(int i=0;i<nbOfPoints;i++)
      nodes[i]=new Node(PlanarIntersector<MyMeshType,MyMatrix>::_coordsT+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[i])*SPACEDIM);
    if(CellModel::GetCellModel(type).isQuadratic())
      return QuadraticPolygon::BuildLinearPolygon(nodes);
    else
      return QuadraticPolygon::BuildArcCirclePolygon(nodes);
  }

  INTERSECTOR_TEMPLATE
  QuadraticPolygon *PTLOC2D_INTERSECTOR::buildPolygonBFrom(ConnType cell, int nbOfPoints, NormalizedCellType type)
  {
    const ConnType *startOfCellNodeConn=PlanarIntersector<MyMeshType,MyMatrix>::_connectS+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexS[OTT<ConnType,numPol>::ind2C(cell)]);
    std::vector<Node *> nodes(nbOfPoints);
    for(int i=0;i<nbOfPoints;i++)
      nodes[i]=new Node(PlanarIntersector<MyMeshType,MyMatrix>::_coordsS+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[i])*SPACEDIM);
    if(type!=NORM_TRI6 && type!=NORM_QUAD8)
      return QuadraticPolygon::BuildLinearPolygon(nodes);
    else
      return QuadraticPolygon::BuildArcCirclePolygon(nodes);
  }
}

#endif
