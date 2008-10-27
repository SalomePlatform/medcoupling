#ifndef _INTERPOLATION2D_HXX_
#define _INTERPOLATION2D_HXX_

#include "InterpolationPlanar.hxx"

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT Interpolation2D : public InterpolationPlanar<Interpolation2D>
  {
  public:
    Interpolation2D() { }
		Interpolation2D(const InterpolationOptions& io):InterpolationPlanar<Interpolation2D>(io){};
  public:
    bool doRotate() const { return false; }
    double medianPlane() const { return 0.; }
    template<class MyMeshType>
    void performAdjustmentOfBB(PlanarIntersector<MyMeshType>* intersector, std::vector<double>& bbox) const
    {}
  };
}

#endif
