#ifndef __CONVEXINTERSECTOR_HXX__
#define __CONVEXINTERSECTOR_HXX__

#include "PlanarIntersector.hxx"
#include "InterpolationUtils.hxx"

namespace INTERP_KERNEL
{
  template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
  class INTERPKERNEL_EXPORT ConvexIntersector: public PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>
  {
  public:
    ConvexIntersector(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh_A,
                      const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh_B, 
                      double dimCaracteristic, double precision, double medianPlane,
                      bool doRotate, int printLevel);
    double intersectCells(ConnType icell_A, ConnType icell_B, int nb_NodesA, int nb_NodesB);
  private :
    double _epsilon;
    const ConnType *_connectA;
    const ConnType *_connectB;
    const ConnType *_connIndexA;
    const ConnType *_connIndexB;
    const double *_coordsA;
    const double *_coordsB;
  };
}

#endif
