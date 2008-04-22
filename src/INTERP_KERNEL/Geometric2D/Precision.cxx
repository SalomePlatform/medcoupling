#include "Precision.hxx"

double INTERP_KERNEL::QUADRATIC_PLANAR::_precision=1e-14;

double INTERP_KERNEL::QUADRATIC_PLANAR::_arcDetectionPrecision=1e-14;

void INTERP_KERNEL::QUADRATIC_PLANAR::setPrecision(double precision)
{ 
  INTERP_KERNEL::QUADRATIC_PLANAR::_precision=precision;
}

void INTERP_KERNEL::QUADRATIC_PLANAR::setArcDetectionPrecision(double precision)
{
  INTERP_KERNEL::QUADRATIC_PLANAR::_arcDetectionPrecision=precision;
}
