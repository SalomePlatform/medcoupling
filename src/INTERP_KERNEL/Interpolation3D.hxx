#ifndef __INTERPOLATION3D_HXX__
#define __INTERPOLATION3D_HXX__

#include "Interpolation.hxx"
#include "NormalizedUnstructuredMesh.hxx"

namespace INTERP_KERNEL
{

  class INTERPKERNEL_EXPORT Interpolation3D : public Interpolation<Interpolation3D>
  {
  public:
    Interpolation3D();
    template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MatrixType, class MyMeshType>
    void interpolateMeshes(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& srcMesh, 
                           const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& targetMesh,
                           MatrixType& result);
  };
}

#endif
