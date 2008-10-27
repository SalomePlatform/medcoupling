#ifndef _INTERPOLATIONPLANAR_HXX_
#define _INTERPOLATIONPLANAR_HXX_

#include "Interpolation.hxx"
#include "PlanarIntersector.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpolationOptions.hxx"

namespace INTERP_KERNEL
{
  template<class RealPlanar>
  class INTERPKERNEL_EXPORT InterpolationPlanar : public Interpolation< InterpolationPlanar<RealPlanar> >
  {
  private:
    double _dimCaracteristic;
    static const double DEFAULT_PRECISION;

  public:
    InterpolationPlanar();
    InterpolationPlanar(const InterpolationOptions & io);

    // geometric precision, debug print level, coice of the median plane, intersection etc ...
    void setOptions(double precision, int printLevel,
                    IntersectionType intersectionType, int orientation=0);

    // Main function to interpolate triangular and quadratic meshes
    template<class MatrixType, class MyMeshType>
    void interpolateMeshes(const MyMeshType& mesh1, const MyMeshType& mesh2, MatrixType& result);

  public:
    bool doRotate() const { return asLeafInterpPlanar().doRotate(); }
    double medianPlane() const { return asLeafInterpPlanar().medianPlane(); }
    template<class MyMeshType>
    void performAdjustmentOfBB(PlanarIntersector<MyMeshType>* intersector, std::vector<double>& bbox) const
    { return asLeafInterpPlanar().performAdjustmentOfBB(intersector,bbox); }

  protected:
    RealPlanar& asLeafInterpPlanar() { return static_cast<RealPlanar&>(*this); }
    const RealPlanar& asLeafInterpPlanar() const { return static_cast< const RealPlanar& >(*this); }
  };
}

#endif
