#ifndef _INTERPOLATION2D_HXX_
#define _INTERPOLATION2D_HXX_

#include "InterpolationPlanar.hxx"

namespace INTERP_KERNEL
{
  class Interpolation2D : public InterpolationPlanar<Interpolation2D>
  {
  public:
    Interpolation2D() { }
  public:
    bool doRotate() const { return false; }
    double medianPlane() const { return 0.; }
    template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
    void performAdjustmentOfBB(PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>* intersector, std::vector<double>& bbox) const
    {}
  };
}

#endif
