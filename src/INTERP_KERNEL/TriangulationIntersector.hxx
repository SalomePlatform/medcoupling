#ifndef _TRIANGULATION_INTERSECTOR_HXX_
#define _TRIANGULATION_INTERSECTOR_HXX_

#include "PlanarIntersector.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType>
  class INTERPKERNEL_EXPORT TriangulationIntersector : public PlanarIntersector<MyMeshType>
  {
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
  public:
    TriangulationIntersector(const MyMeshType& mesh_A, const MyMeshType& mesh_B,
                             double dimCaracteristic, double precision, double medianPlane, int printLevel);
    double intersectCells(ConnType icell_A, ConnType icell_B, int nb_NodesA, int nb_NodesB);
  private :
    const ConnType *_connectA;
    const ConnType *_connectB;
    const double *_coordsA;
    const double *_coordsB;
    const ConnType *_connIndexA;
    const ConnType *_connIndexB;
  };
}

#endif
