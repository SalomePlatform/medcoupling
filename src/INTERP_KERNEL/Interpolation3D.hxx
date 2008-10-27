#ifndef __INTERPOLATION3D_HXX__
#define __INTERPOLATION3D_HXX__

#include "Interpolation.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "IntersectorHexa.hxx"
#include "InterpolationOptions.hxx"
namespace INTERP_KERNEL
{
  class INTERPKERNEL_EXPORT Interpolation3D : public Interpolation<Interpolation3D>
  {
  public:
    Interpolation3D();
    Interpolation3D(const InterpolationOptions& io);
    template<class MatrixType, class MyMeshType>
    void interpolateMeshes(const MyMeshType& srcMesh, const MyMeshType& targetMesh, MatrixType& result);
  private:
    SplittingPolicy _splitting_policy;
  };
}

#endif
