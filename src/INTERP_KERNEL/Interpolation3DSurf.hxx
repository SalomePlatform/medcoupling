#ifndef __INTERPOLATION3DSURF_HXX__
#define __INTERPOLATION3DSURF_HXX__

#include "InterpolationPlanar.hxx"
#include "InterpolationOptions.hxx"

namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT Interpolation3DSurf : public InterpolationPlanar<Interpolation3DSurf>
  {
  public:
    Interpolation3DSurf();
		Interpolation3DSurf(const InterpolationOptions& io);
    void setOptions(double precision, int printLevel, double medianPlane, 
                    IntersectionType intersectionType, bool doRotate, int orientation=0);

  public:
    bool doRotate() const { return _doRotate; }
    double medianPlane() const { return _medianPlane; }
    template<class MyMeshType>
    void performAdjustmentOfBB(PlanarIntersector<MyMeshType>* intersector, std::vector<double>& bbox) const
    { intersector->adjustBoundingBoxes(bbox,_surf3DAdjustmentEps); }
  protected:
    bool _doRotate;
    double _medianPlane;
    double _surf3DAdjustmentEps;
    static const double DFT_MEDIAN_PLANE;
    static const double DFT_SURF3D_ADJ_EPS;
  };
}

#endif
