#ifndef __CONVEXINTERSECTOR_HXX__
#define __CONVEXINTERSECTOR_HXX__

#include "PlanarIntersector.hxx"
#include "InterpolationUtils.hxx"

namespace INTERP_KERNEL
{
  template<class MyMeshType>
  class INTERPKERNEL_EXPORT ConvexIntersector : public PlanarIntersector<MyMeshType>
  {
  public:
    static const int SPACEDIM=MyMeshType::MY_SPACEDIM;
    static const int MESHDIM=MyMeshType::MY_MESHDIM;
    typedef typename MyMeshType::MyConnType ConnType;
    static const NumberingPolicy numPol=MyMeshType::My_numPol;
  public:
    ConvexIntersector(const MyMeshType& mesh_A, const MyMeshType& mesh_B, 
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
