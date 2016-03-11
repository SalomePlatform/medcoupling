// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
// See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
// Author : Anthony Geay (CEA/DEN)

#include "Interpolation2DCurve.hxx"
#include "InterpolationCurve.txx"

namespace INTERP_KERNEL
{
  Interpolation2DCurve::Interpolation2DCurve()
  {
    // to have non-zero default thickness of target element
    InterpolationOptions::setBoundingBoxAdjustmentAbs( InterpolationOptions::getPrecision() );
  }

  Interpolation2DCurve::Interpolation2DCurve
    (const InterpolationOptions& io):InterpolationCurve<Interpolation2DCurve>(io)
  {
    // to have non-zero default thickness of target element
    InterpolationOptions::setBoundingBoxAdjustmentAbs( InterpolationOptions::getPrecision() );
  }

  /**
   *  \brief  Function used to set the options for the intersection calculation
   * \details The following options can be modified:
   *  -# Precision: Level of precision of the computations.
   *   - Values: positive real number.
   *   - Default: 1.0E-12.
   *  -# Tolerance: Thickness of target element.
   *   - Values: positive real number.
   *   - Default: 1.0E-12.
   *  -# Median line: Position of the median line where both segments will be projected.
   *   - Values: real number between 0.0 and 1.0.
   *   - Default: 0.5
   */
  void Interpolation2DCurve::setOptions (double precision,
                                         double tolerance,
                                         double medianLine)
  {
    InterpolationOptions::setPrecision(precision);
    InterpolationOptions::setBoundingBoxAdjustmentAbs(tolerance);
    InterpolationOptions::setMedianPlane(medianLine);
  }
}
