#ifndef _TRIANGULATION_INTERSECTOR_HXX_
#define _TRIANGULATION_INTERSECTOR_HXX_

#include "PlanarIntersector.hxx"

namespace INTERP_KERNEL
{
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  class INTERPKERNEL_EXPORT TriangulationIntersector : public PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>
  {
  public:
    TriangulationIntersector(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh_A,
                             const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh_B,
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
