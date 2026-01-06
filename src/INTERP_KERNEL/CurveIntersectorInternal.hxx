// Copyright (C) 2024-2026  CEA, EDF
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

#ifndef __CurveIntersectorInternal_HXX__
#define __CurveIntersectorInternal_HXX__

#include "INTERPKERNELDefines.hxx"

#include <array>

namespace INTERP_KERNEL
{
/*!
 * \class CurveIntersectorInternal
 * Class defining some functions fir internal projections
 */
class CurveIntersectorInternal
{
   public:
    CurveIntersectorInternal() {};
    ~CurveIntersectorInternal() {};
    INTERPKERNEL_EXPORT static double InternalProjectionFrom3DTo2D(
        const double coordsT[6], const double coordsS[6], double tolerance, double coordsTOut[4], double coordsSOut[4]
    );
    INTERPKERNEL_EXPORT static bool CurveIntersectorInternalProjectionThis2D(
        const double *coordsT,
        const double *coordsS,
        double tolerance,
        double precision,
        double medianLine,
        double &xs0,
        double &xs1,
        double &xt0,
        double &xt1
    );
    static double InternalProjectionFrom3DTo2DHelper(
        const std::array<double, 3> &t01, double norm_t01, const std::array<double, 3> &t0s0, double s0x
    );
};
}  // namespace INTERP_KERNEL

#endif
