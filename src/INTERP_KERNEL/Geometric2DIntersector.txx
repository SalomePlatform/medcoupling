#ifndef __GEOMETRIC2DINTERSECTOR_TXX__
#define __GEOMETRIC2DINTERSECTOR_TXX__

#include "Geometric2DIntersector.hxx"
#include "QuadraticPolygon.hxx"
#include "EdgeArcCircle.hxx"
#include "EdgeLin.hxx"
#include "Node.hxx"

namespace INTERP_KERNEL
{
  namespace QUADRATIC_PLANAR
  {
    extern double _precision;
  }

  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  Geometric2DIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::Geometric2DIntersector(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh_A,
                                                                                              const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh_B,
                                                                                              double dimCaracteristic, double precision):
    PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>(dimCaracteristic, precision, 0., false, 0),_meshA(mesh_A),_meshB(mesh_B)
  {
    _connectA= mesh_A.getConnectivityPtr();
    _connectB= mesh_B.getConnectivityPtr();
    _connIndexA= mesh_A.getConnectivityIndexPtr();
    _connIndexB= mesh_B.getConnectivityIndexPtr();
    _coordsA = mesh_A.getCoordinatesPtr();
    _coordsB = mesh_B.getCoordinatesPtr();
    QUADRATIC_PLANAR::_precision=dimCaracteristic*precision;
  }
  
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  double Geometric2DIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::intersectCells(ConnType icell_A, ConnType icell_B, 
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

  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  QuadraticPolygon *Geometric2DIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::buildPolygonAFrom(ConnType cell, int nbOfPoints, NormalizedCellType type)
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

  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  QuadraticPolygon *Geometric2DIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>::buildPolygonBFrom(ConnType cell, int nbOfPoints, NormalizedCellType type)
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
