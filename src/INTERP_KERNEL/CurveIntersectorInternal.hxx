// Copyright (C) 2024  CEA, EDF
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

#pragma once

#include "INTERPKERNELDefines.hxx"

INTERPKERNEL_EXPORT double InternalProjectionFrom3DTo2D(const double coordsT[6], const double coordsS[6], double tolerance, double coordsTOut[4], double coordsSOut[4]);

INTERPKERNEL_EXPORT bool CurveIntersectorInternalProjectionThis2D(const double *coordsT, const double *coordsS, double tolerance, double precision, double medianLine,
                                                                  double& xs0, double& xs1, double& xt0, double& xt1);
