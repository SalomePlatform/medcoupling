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
#ifndef __GEOMETRIC2DINTERSECTOR_TXX__
#define __GEOMETRIC2DINTERSECTOR_TXX__

#include "Geometric2DIntersector.hxx"
#include "PlanarIntersectorP0P0.txx"
#include "Planar2D1DIntersectorP0P0.txx"
#include "PlanarIntersectorP0P1.txx"
#include "PlanarIntersectorP1P0.txx"
#include "PlanarIntersectorP1P1.txx"
#include "PlanarIntersectorP1P0Bary.txx"
#include "PlanarIntersectorP0P1Bary.txx"
#include "CellModel.hxx"

#include "InterpKernelGeo2DQuadraticPolygon.hxx"
#include "InterpKernelGeo2DEdgeArcCircle.hxx"
#include "InterpKernelGeo2DEdgeLin.hxx"
#include "InterpKernelGeo2DNode.hxx"

#define GEO2D_INTERSECTOR    Geometric2DIntersector<MyMeshType,MyMatrix,InterpType>
#define INTERSECTOR_TEMPLATE template<class MyMeshType, class MyMatrix, template <class MeshType, class TheMatrix, class ThisIntersector> class InterpType>

namespace INTERP_KERNEL
{
  INTERSECTOR_TEMPLATE
  GEO2D_INTERSECTOR::Geometric2DIntersector(const MyMeshType& meshT, const MyMeshType& meshS,
                                            double dimCaracteristic, double md3DSurf, double minDot3DSurf, double medianPlane,
                                            double precision, int orientation):
    InterpType<MyMeshType,MyMatrix,GEO2D_INTERSECTOR >(meshT,meshS,dimCaracteristic, precision, md3DSurf, minDot3DSurf, medianPlane, true, orientation, 0),
    _precision(precision)
  {
  }
  
  INTERSECTOR_TEMPLATE
  double GEO2D_INTERSECTOR::intersectGeometry(ConnType icellT,   ConnType icellS,
                                              ConnType nbNodesT, ConnType nbNodesS)
  {
    int orientation = 1;
    std::vector<double> CoordsT;
    std::vector<double> CoordsS;
    PlanarIntersector<MyMeshType,MyMatrix>::getRealCoordinates(icellT,icellS,nbNodesT,nbNodesS,CoordsT,CoordsS,orientation);
    NormalizedCellType tT=PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getTypeOfElement(icellT);
    NormalizedCellType tS=PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getTypeOfElement(icellS);
    QuadraticPolygon *p1=buildPolygonFrom(CoordsT,tT);
    QuadraticPolygon *p2=buildPolygonFrom(CoordsS,tS);
    double ret=p1->intersectWithAbs(*p2);
    delete p1; delete p2;
    return orientation*ret;
  }

  INTERSECTOR_TEMPLATE
  double GEO2D_INTERSECTOR::intersectGeometry1D(ConnType icellT,   ConnType icellS,
                                                ConnType nbNodesT, ConnType nbNodesS,
                                                bool& isColinear)
  {
    int orientation = 1;
    std::vector<double> CoordsT;
    std::vector<double> CoordsS;
    PlanarIntersector<MyMeshType,MyMatrix>::getRealCoordinates(icellT,icellS,nbNodesT,nbNodesS,CoordsT,CoordsS,orientation);
    NormalizedCellType tT=PlanarIntersector<MyMeshType,MyMatrix>::_meshT.getTypeOfElement(icellT);
    NormalizedCellType tS=PlanarIntersector<MyMeshType,MyMatrix>::_meshS.getTypeOfElement(icellS);
    QuadraticPolygon *p1=buildPolygonFrom(CoordsT,tT);
    QuadraticPolygon *p2=buildPolygonOfOneEdgeFrom(CoordsS,tS);
    double ret=p1->intersectWithAbs1D(*p2, isColinear);
    delete p1; delete p2;
    return orientation*ret;
  }

  INTERSECTOR_TEMPLATE
  double GEO2D_INTERSECTOR::intersectGeometryWithQuadrangle(const double             * quadrangle,
                                                            const std::vector<double>& sourceCoords,
                                                            bool                       isSourceQuad)
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
    QuadraticPolygon *p1=QuadraticPolygon::BuildLinearPolygon(nodes);
    QuadraticPolygon *p2;
    if(!isSourceQuad)
      p2=QuadraticPolygon::BuildLinearPolygon(nodes2);
    else
      p2=QuadraticPolygon::BuildArcCirclePolygon(nodes2);
    double ret=p1->intersectWithAbs(*p2);
    delete p1; delete p2;
    return ret;
  }

  INTERSECTOR_TEMPLATE
  double GEO2D_INTERSECTOR::intersectGeometryGeneral(const std::vector<double>& targetCoords,
                                                     const std::vector<double>& sourceCoords)
  {
    int nbOfTargetNodes=targetCoords.size()/SPACEDIM;
    std::vector<Node *> nodes(nbOfTargetNodes);
    for(int i=0;i<nbOfTargetNodes;i++)
      nodes[i]=new Node(targetCoords[i*SPACEDIM],targetCoords[i*SPACEDIM+1]);
    int nbOfSourceNodes=sourceCoords.size()/SPACEDIM;
    std::vector<Node *> nodes2(nbOfSourceNodes);
    for(int i=0;i<nbOfSourceNodes;i++)
      nodes2[i]=new Node(sourceCoords[i*SPACEDIM],sourceCoords[i*SPACEDIM+1]);
    QuadraticPolygon *p1=QuadraticPolygon::BuildLinearPolygon(nodes);
    QuadraticPolygon *p2=QuadraticPolygon::BuildLinearPolygon(nodes2);
    double ret=p1->intersectWithAbs(*p2);
    delete p1; delete p2;
    return ret;
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
  double GEO2D_INTERSECTOR::intersectGeoBary(const std::vector<double>& targetCell,
                                             bool                       targetCellQuadratic,
                                             const double *             sourceTria,
                                             std::vector<double>&       res)
  {
    std::vector<Node *> nodes(3);
    nodes[0]=new Node(sourceTria[0*SPACEDIM],sourceTria[0*SPACEDIM+1]);
    nodes[1]=new Node(sourceTria[1*SPACEDIM],sourceTria[1*SPACEDIM+1]);
    nodes[2]=new Node(sourceTria[2*SPACEDIM],sourceTria[2*SPACEDIM+1]);
    int nbOfTargetNodes=targetCell.size()/SPACEDIM;
    std::vector<Node *> nodes2(nbOfTargetNodes);
    for(int i=0;i<nbOfTargetNodes;i++)
      nodes2[i]=new Node(targetCell[i*SPACEDIM],targetCell[i*SPACEDIM+1]);
    QuadraticPolygon *p1=QuadraticPolygon::BuildLinearPolygon(nodes);
    QuadraticPolygon *p2;
    if(!targetCellQuadratic)
      p2=QuadraticPolygon::BuildLinearPolygon(nodes2);
    else
      p2=QuadraticPolygon::BuildArcCirclePolygon(nodes2);
    double barycenter[2];
    double ret=p1->intersectWithAbs(*p2,barycenter);
    delete p1; delete p2;
    if ( ret > std::numeric_limits<double>::min() )
    {
      std::vector<const double* > sourceCell(3);
      sourceCell[0] = &sourceTria[0];
      sourceCell[1] = &sourceTria[SPACEDIM];
      sourceCell[2] = &sourceTria[SPACEDIM*2];
      res.resize(3);
      barycentric_coords( sourceCell, barycenter, &res[0]);
      res[0] *= ret;
      res[1] *= ret;
      res[2] *= ret;
    }
    else
    {
      ret = 0;
    }
    return ret;
  }

  INTERSECTOR_TEMPLATE
  QuadraticPolygon *GEO2D_INTERSECTOR::buildPolygonFrom(const std::vector<double>& coords, NormalizedCellType type)
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
  QuadraticPolygon *GEO2D_INTERSECTOR::buildPolygonOfOneEdgeFrom(const std::vector<double>& coords, NormalizedCellType type)
  {
    if(type==NORM_SEG2)
      {
        Node *node0=new Node(coords[0],coords[1]);
        Node *node1=new Node(coords[SPACEDIM],coords[SPACEDIM+1]);
        QuadraticPolygon *ret=new QuadraticPolygon;
        ret->pushBack(new EdgeLin(node0,node1));
        node0->decrRef(); node1->decrRef();
        return ret;
      }
    else if(type==NORM_SEG3)
      {
        Node *nodeBg=new Node(coords[0],coords[1]);
        Node *nodeEnd=new Node(coords[SPACEDIM],coords[SPACEDIM+1]);
        Node *nodeMiddle=new Node(coords[2*SPACEDIM],coords[2*SPACEDIM+1]);
        QuadraticPolygon *ret=new QuadraticPolygon;
        ret->pushBack(new EdgeArcCircle(nodeBg,nodeMiddle,nodeEnd));
        nodeBg->decrRef(); nodeEnd->decrRef(); nodeMiddle->decrRef();
        return ret;
      }
    else
      throw INTERP_KERNEL::Exception("buildPolygonOfOneEdgeFrom : trying to build such non close QuadraticPolygon with 1D type !");
  }

  INTERSECTOR_TEMPLATE
  QuadraticPolygon *GEO2D_INTERSECTOR::buildPolygonAFrom(ConnType cell, int nbOfPoints, NormalizedCellType type)
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
  QuadraticPolygon *GEO2D_INTERSECTOR::buildPolygonBFrom(ConnType cell, int nbOfPoints, NormalizedCellType type)
  {
    const ConnType *startOfCellNodeConn=PlanarIntersector<MyMeshType,MyMatrix>::_connectS+OTT<ConnType,numPol>::conn2C(PlanarIntersector<MyMeshType,MyMatrix>::_connIndexS[OTT<ConnType,numPol>::ind2C(cell)]);
    std::vector<Node *> nodes(nbOfPoints);
    for(int i=0;i<nbOfPoints;i++)
      nodes[i]=new Node(PlanarIntersector<MyMeshType,MyMatrix>::_coordsS+OTT<ConnType,numPol>::coo2C(startOfCellNodeConn[i])*SPACEDIM);
    const CellModel& cm=CellModel::GetCellModel(type);
    if(!cm.isQuadratic())
      return QuadraticPolygon::BuildLinearPolygon(nodes);
    else
      return QuadraticPolygon::BuildArcCirclePolygon(nodes);
  }
}

#endif
