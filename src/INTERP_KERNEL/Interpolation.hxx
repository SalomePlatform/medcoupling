#ifndef _INTERPOLATION_HXX_
#define _INTERPOLATION_HXX_

/**
 * \mainpage
 * Status : documentation of 3D - part of intersection matrix calculation more or less complete
 *
 *
 */
#include "INTERPKERNEL_defines.hxx"
#include "InterpolationOptions.hxx"

namespace INTERP_KERNEL
{
  template<class TrueMainInterpolator>
  class INTERPKERNEL_EXPORT Interpolation : public InterpolationOptions
  {
  public:
    Interpolation() { 
	// 	double InterpolationOptions::_precision=1e-12;
// 		int InterpolationOptions::_printLevel=0;
// 		InterpolationOptions::setIntersectionType(Triangulation);
// 		InterpolationOptions::setMedianPlane(0.5);
// 		InterpolationOptions::setDoRotate(true);
// 		InterpolationOptions::setBoundingBoxAdjustment(0.1);
// 		InterpolationOptions::setSplittingPolicy(GENERAL_48);
		}
		Interpolation(const InterpolationOptions& io) :InterpolationOptions(io){}
    //interpolation of two triangular meshes.
    template<class MatrixType, class MyMeshType>
    void interpolateMeshes(const MyMeshType& mesh1, const MyMeshType& mesh2, MatrixType& result)
    { return asLeaf().interpolateMeshes(mesh1,mesh2,result); }
  protected:
    TrueMainInterpolator& asLeaf() { return static_cast<TrueMainInterpolator&>(*this); }
  };
}

#endif
