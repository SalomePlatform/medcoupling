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

double INTERP_KERNEL::QUADRATIC_PLANAR::_precision=1e-14;

double INTERP_KERNEL::QUADRATIC_PLANAR::_arc_detection_precision=1e-14;

void INTERP_KERNEL::QUADRATIC_PLANAR::setPrecision(double precision)
{ 
  INTERP_KERNEL::QUADRATIC_PLANAR::_precision=precision;
}

void INTERP_KERNEL::QUADRATIC_PLANAR::setArcDetectionPrecision(double precision)
{
  INTERP_KERNEL::QUADRATIC_PLANAR::_arc_detection_precision=precision;
}
