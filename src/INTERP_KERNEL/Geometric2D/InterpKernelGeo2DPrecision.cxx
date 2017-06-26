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

#include "InterpKernelGeo2DPrecision.hxx"

double INTERP_KERNEL::QuadraticPlanarPrecision::_precision=1e-14;
double INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::_arc_detection_precision=1e-14;

INTERP_KERNEL::QuadraticPlanarPrecision::QuadraticPlanarPrecision(double precision):
    _initial_precision(_precision)
{
  _precision=precision;
}

INTERP_KERNEL::QuadraticPlanarPrecision::~QuadraticPlanarPrecision()
{
  _precision = _initial_precision;
}

INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::QuadraticPlanarArcDetectionPrecision(double precision):
    _initial_arc_detection_precision(_arc_detection_precision)
{
  _arc_detection_precision=precision;
}

INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::~QuadraticPlanarArcDetectionPrecision()
{
  _arc_detection_precision = _initial_arc_detection_precision;
}


void INTERP_KERNEL::QuadraticPlanarPrecision::setPrecision(double precision)
{ 
  _precision=precision;
}

void INTERP_KERNEL::QuadraticPlanarArcDetectionPrecision::setArcDetectionPrecision(double precision)
{
  _arc_detection_precision=precision;
}
