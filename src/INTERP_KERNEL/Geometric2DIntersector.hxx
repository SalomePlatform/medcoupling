#ifndef __GEOMETRIC2DINTERSECTOR_HXX__
#define __GEOMETRIC2DINTERSECTOR_HXX__

#include "PlanarIntersector.hxx"

namespace INTERP_KERNEL
{
  class QuadraticPolygon;

  template<class MyMeshType>
  class INTERPKERNEL_EXPORT Geometric2DIntersector : public PlanarIntersector<MyMeshType>
  {
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
  public:
    Geometric2DIntersector(const MyMeshType& mesh_A, const MyMeshType& mesh_B,
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
    const MyMeshType& _meshA;
    const MyMeshType& _meshB;
  };
}

#endif
