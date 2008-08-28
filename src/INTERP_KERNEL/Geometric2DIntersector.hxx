#ifndef __GEOMETRIC2DINTERSECTOR_HXX__
#define __GEOMETRIC2DINTERSECTOR_HXX__

#include "PlanarIntersector.hxx"

namespace INTERP_KERNEL
{
  class QuadraticPolygon;

  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  class INTERPKERNEL_EXPORT Geometric2DIntersector : public PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>
  {
  public:
    Geometric2DIntersector(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh_A,
                           const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh_B,
                           double dimCaracteristic, double precision);
    double intersectCells(ConnType icell_A, ConnType icell_B, int nb_NodesA, int nb_NodesB);
  private:
    QuadraticPolygon *buildPolygonAFrom(ConnType cell, int nbOfPoints, NormalizedCellType type);
    QuadraticPolygon *buildPolygonBFrom(ConnType cell, int nbOfPoints, NormalizedCellType type);
  private:
    const ConnType *_connectA;
    const ConnType *_connectB;
    const double *_coordsA;
    const double *_coordsB;
    const ConnType *_connIndexA;
    const ConnType *_connIndexB;
    const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& _meshA;
    const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& _meshB;
  };
}

#endif
