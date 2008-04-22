#ifndef _INTERPOLATIONPLANAR_HXX_
#define _INTERPOLATIONPLANAR_HXX_

#include "Interpolation.hxx"
#include "PlanarIntersector.hxx"
#include "NormalizedUnstructuredMesh.hxx"

namespace INTERP_KERNEL
{
  typedef enum {Triangulation, Convex, Geometric2D, Generic} IntersectionType;

  template<class RealPlanar>
  class InterpolationPlanar : public Interpolation< InterpolationPlanar<RealPlanar> >
  {
  private: 
    int  _printLevel;
    double _precision;
    double _dimCaracteristic;
    IntersectionType _intersectionType;
    static const double DEFAULT_PRECISION=1.e-12;
  public:
    InterpolationPlanar();
  
    // geometric precision, debug print level, coice of the median plane, intersection etc ...
    void setOptions(double precision, int printLevel,
                    IntersectionType intersectionType);
  
    // Main function to interpolate triangular and quadratic meshes
    template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MatrixType, class MyMeshType>
    void interpolateMeshes(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh1, 
                           const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh2,
                           MatrixType& result);
  public:
    bool doRotate() const { return asLeafInterpPlanar().doRotate(); }
    double medianPlane() const { return asLeafInterpPlanar().medianPlane(); }
    template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MyMeshType>
    void performAdjustmentOfBB(PlanarIntersector<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>* intersector, std::vector<double>& bbox) const
    { return asLeafInterpPlanar().performAdjustmentOfBB(intersector,bbox); }
  protected:
    RealPlanar& asLeafInterpPlanar() { return static_cast<RealPlanar&>(*this); }
    const RealPlanar& asLeafInterpPlanar() const { return static_cast< const RealPlanar& >(*this); }
  };
}

#endif
