#ifndef _INTERPOLATION_HXX_
#define _INTERPOLATION_HXX_

#include "NormalizedUnstructuredMesh.hxx"

/**
 * \mainpage
 * Status : documentation of 3D - part of intersection matrix calculation more or less complete
 *
 *
 */
namespace INTERP_KERNEL
{
  template<class TrueMainInterpolator>
  class Interpolation
  {
  public:
    Interpolation() { }
  
    //interpolation of two triangular meshes.
    template<int SPACEDIM, int MESHDIM, class ConnType, NumberingPolicy numPol, class MatrixType, class MyMeshType>
    void interpolateMeshes(const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh1,
                           const NormalizedUnstructuredMesh<SPACEDIM,MESHDIM,ConnType,numPol,MyMeshType>& mesh2,
                           MatrixType& result)
    { return asLeaf().interpolateMeshes(mesh1,mesh2,result); }
  protected:
    TrueMainInterpolator& asLeaf() { return static_cast<TrueMainInterpolator&>(*this); }
  };
}

#endif
