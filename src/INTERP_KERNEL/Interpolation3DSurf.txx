#ifndef __INTERPOLATION3DSURF_TXX__
#define __INTERPOLATION3DSURF_TXX__

#include "Interpolation3DSurf.hxx"
#include "InterpolationPlanar.txx"

namespace INTERP_KERNEL
{
  const double Interpolation3DSurf::DFT_MEDIAN_PLANE=0.5;
  const double Interpolation3DSurf::DFT_SURF3D_ADJ_EPS=1e-4;
  
  Interpolation3DSurf::Interpolation3DSurf():_doRotate(true)
                                            ,_medianPlane(DFT_MEDIAN_PLANE)
                                            ,_surf3DAdjustmentEps(DFT_SURF3D_ADJ_EPS)
  {
  }

  /**
     \brief  Function used to set the options for the intersection calculation
     \details The following options can be modified:
     -# Intersection_type: the type of algorithm to be used in the computation of the cell-cell intersections.
     - Values: Triangle, Convex.
     - Default: Triangle.
     -# MedianPlane: Position of the median plane where both cells will be projected
     - Values: between 0 and 1.
     - Default: 0.5.
     -# DoRotate: rotate the coordinate system such that the target cell is in the Oxy plane.
     - Values: true (necessarilly if Intersection_type=Triangle), false.
     - Default: true (as default Intersection_type=Triangle)
     -# Precision: Level of precision of the computations is precision times the characteristic size of the mesh.
     - Values: positive real number.
     - Default: 1.0E-12.
     -# PrintLevel: Level of verboseness during the computations.
     - Values: interger between 0 and 3.
     - Default: 0.
  */
  void Interpolation3DSurf::setOptions(double precision, int printLevel, double medianPlane, 
                                       IntersectionType intersectionType, bool doRotate)
  {
    InterpolationPlanar<Interpolation3DSurf>::setOptions(precision,printLevel,intersectionType);
    _doRotate=doRotate;
    _medianPlane=medianPlane;
  }
}

#endif
