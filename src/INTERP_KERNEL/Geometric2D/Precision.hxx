#ifndef __PRECISION_HXX__
#define __PRECISION_HXX__

#include "Geometric2D_defines.hxx"

namespace INTERP_KERNEL
{
  class GEOMETRIC2D_EXPORT QUADRATIC_PLANAR
  {
  public:
    static double _precision;
    static double _arcDetectionPrecision;
    static void setPrecision(double precision);
    static void setArcDetectionPrecision(double precision);
  };
}

#endif
