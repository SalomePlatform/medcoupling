#ifndef __INTERPOLATION3DSURF_HXX__
#define __INTERPOLATION3DSURF_HXX__

#include "InterpolationPlanar.hxx"

namespace INTERP_KERNEL
{
  class Interpolation3DSurf : public InterpolationPlanar<Interpolation3DSurf>
  {
  public:
    Interpolation3DSurf();
    void setOptions(double precision, int printLevel, double medianPlane, 
                    IntersectionType intersectionType, bool doRotate);
  public:
    bool doRotate() const { return _doRotate; }
    double medianPlane() const { return _medianPlane; }
    template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
    void performAdjustmentOfBB(PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>* intersector, std::vector<double>& bbox) const
    { intersector->adjustBoundingBoxes(bbox,_surf3DAdjustmentEps); }
  protected:
    bool _doRotate;
    double _medianPlane;
    double _surf3DAdjustmentEps;
    static const double DFT_MEDIAN_PLANE=0.5;
    static const double DFT_SURF3D_ADJ_EPS=1e-4;
  };
}

#endif
