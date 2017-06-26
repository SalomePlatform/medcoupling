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

#ifndef __INTERPKERNELGEO2DPRECISION_HXX__
#define __INTERPKERNELGEO2DPRECISION_HXX__

#include "INTERPKERNELDefines.hxx"

namespace INTERP_KERNEL
{
  /* !!TODO: a more global review of the code should be done, so that eps is always a parameter of all methods
     instead of being stored as a static attribute.
  */

  /** Class storing the precision for the detection of colinear segments, coincident points, etc ...
   * in Geometric2D computations.
   *
   * RAII pattern allowing to temporarily override Geometric2D precision.
   * When the instance is destroyed, the previous precision is set back.
   *
   */
  class QuadraticPlanarPrecision
  {
  public:
    QuadraticPlanarPrecision(double prec);
    virtual ~QuadraticPlanarPrecision();

    static void setPrecision(double precision);
    inline static double getPrecision() { return _precision; }
  private:
    static double _precision;
    double _initial_precision;
  };

  /** Class storing the precision for the detection of degenerated arc of circles
   * in Geometric2D computations.
   *
   * RAII pattern allowing to temporarily override Geometric2D precision.
   * When the instance is destroyed, the previous precision is set back.
   */
  class QuadraticPlanarArcDetectionPrecision
  {
  public:
    QuadraticPlanarArcDetectionPrecision(double prec);
    virtual ~QuadraticPlanarArcDetectionPrecision();

    static void setArcDetectionPrecision(double precision);
    inline static double getArcDetectionPrecision() { return _arc_detection_precision; }
  private:
    static double _arc_detection_precision;
    double _initial_arc_detection_precision;

  };
}

#endif

