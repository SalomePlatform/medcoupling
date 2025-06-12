// Copyright (C) 2007-2025  CEA, EDF
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

#include "CurveIntersectorInternal.hxx"
#include <cmath>
#include <limits>
#include <array>

#include "VectorUtils.hxx"

double
INTERP_KERNEL::CurveIntersectorInternal::InternalProjectionFrom3DTo2DHelper(
    const std::array<double, 3> &t01, double norm_t01, const std::array<double, 3> &t0s0, double s0x
)
{
    std::array<double, 3> t0s0_ortho(t01);
    INTERP_KERNEL::scaleVector(t0s0_ortho.data(), -s0x / norm_t01);
    INTERP_KERNEL::add(t0s0.data(), t0s0_ortho.data());
    double dotForZeSign = INTERP_KERNEL::dot(t01.data(), t0s0_ortho.data());
    double zeSign = std::signbit(dotForZeSign) ? -1.0 : 1.0;
    double s0y = INTERP_KERNEL::norm(t0s0_ortho.data()) * zeSign;
    return s0y;
}

double
INTERP_KERNEL::CurveIntersectorInternal::InternalProjectionFrom3DTo2D(
    const double coordsT[6], const double coordsS[6], double tolerance, double coordsTOut[4], double coordsSOut[4]
)
{  // EDF31137
    std::array<double, 3> s01 = {coordsS[3] - coordsS[0], coordsS[4] - coordsS[1], coordsS[5] - coordsS[2]};
    double norm_s01 = INTERP_KERNEL::norm(s01.data());
    std::array<double, 3> t01 = {coordsT[3] - coordsT[0], coordsT[4] - coordsT[1], coordsT[5] - coordsT[2]};
    double norm_t01 = INTERP_KERNEL::norm(t01.data());
    std::array<double, 3> s01_cross_t01;
    INTERP_KERNEL::cross(s01.data(), t01.data(), s01_cross_t01.data());
    double norm_s01_cross_t01 = INTERP_KERNEL::norm(s01_cross_t01.data());
    double areSTParallel = norm_s01_cross_t01 / (norm_s01 * norm_t01);
    std::array<double, 3> t0s0{coordsS[0] - coordsT[0], coordsS[1] - coordsT[1], coordsS[2] - coordsT[2]};
    std::array<double, 3> t0s1{coordsS[3] - coordsT[0], coordsS[4] - coordsT[1], coordsS[5] - coordsT[2]};
    double *s0x(coordsSOut), *s0y(coordsSOut + 1), *s1x(coordsSOut + 2), *s1y(coordsSOut + 3);
    coordsTOut[0] = 0.0;
    coordsTOut[1] = 0.0;
    coordsTOut[2] = norm_t01;
    coordsTOut[3] = 0.0;
    if (areSTParallel < tolerance)
    {
        *s0x = INTERP_KERNEL::dot(t01.data(), t0s0.data()) / norm_t01;
        *s1x = INTERP_KERNEL::dot(t01.data(), t0s1.data()) / norm_t01;
        *s0y = InternalProjectionFrom3DTo2DHelper(t01, norm_t01, t0s0, *s0x);
        *s1y = InternalProjectionFrom3DTo2DHelper(t01, norm_t01, t0s1, *s1x);
        return 0.0;
    }
    std::array<double, 3> w(s01_cross_t01);
    INTERP_KERNEL::scaleVector(w.data(), 1.0 / (areSTParallel * norm_s01 * norm_t01));
    std::array<double, 3> v;
    INTERP_KERNEL::cross(w.data(), t01.data(), v.data());
    INTERP_KERNEL::scaleVector(v.data(), 1.0 / norm_t01);
    std::array<double, 3> u(t01);
    INTERP_KERNEL::scaleVector(u.data(), 1.0 / norm_t01);
    double distW(INTERP_KERNEL::dot(t0s0.data(), w.data()));
    std::array<double, 3> t0s0_inplane(t0s0);
    {
        std::array<double, 3> tmp(w);
        INTERP_KERNEL::scaleVector(tmp.data(), -distW);
        INTERP_KERNEL::add(tmp.data(), t0s0_inplane.data());
    }
    *s0x = INTERP_KERNEL::dot(t0s0_inplane.data(), u.data());
    *s0y = INTERP_KERNEL::dot(t0s0_inplane.data(), v.data());
    std::array<double, 3> t0s1_inplane(t0s1);
    {
        std::array<double, 3> tmp(w);
        INTERP_KERNEL::scaleVector(tmp.data(), -distW);
        INTERP_KERNEL::add(tmp.data(), t0s1_inplane.data());
    }
    *s1x = INTERP_KERNEL::dot(t0s1_inplane.data(), u.data());
    *s1y = INTERP_KERNEL::dot(t0s1_inplane.data(), v.data());
    return fabs(distW);
}

bool
INTERP_KERNEL::CurveIntersectorInternal::CurveIntersectorInternalProjectionThis2D(
    const double *coordsT,
    const double *coordsS,
    double tolerance,
    double precision,
    double medianLine,
    double &xs0,
    double &xs1,
    double &xt0,
    double &xt1
)
{
    // Pass 2D->1D

    enum
    {
        X = 0,
        Y
    };

    // check if two segments overlap in 2D within tolerance

    const double *t0 = coordsT;
    const double *t1 = coordsT + 2;
    double t01[2] = {t1[X] - t0[X], t1[Y] - t0[Y]};          // tgt segment direction
    double tSize = sqrt(t01[X] * t01[X] + t01[Y] * t01[Y]);  // tgt segment size
    if (tSize < precision)
        return false;                  // degenerated segment
    t01[X] /= tSize, t01[Y] /= tSize;  // normalize t01

    const double *s0 = coordsS;
    const double *s1 = coordsS + 2;
    double t0s0[2] = {s0[X] - t0[X], s0[Y] - t0[Y]};
    double t0s1[2] = {s1[X] - t0[X], s1[Y] - t0[Y]};
    double nt01_x_t0s0 = t0s0[X] * t01[Y] - t0s0[Y] * t01[X];  // t0s0 dot norm of t01
    double nt01_x_t0s1 = t0s1[X] * t01[Y] - t0s1[Y] * t01[X];  // t0s1 dot norm of t01
    double dist_ts0 = fabs(nt01_x_t0s0);                       // dist from tgt seg to s0
    double dist_ts1 = fabs(nt01_x_t0s1);                       // dist from tgt seg to s1
    bool s0_out_of_tol = (dist_ts0 > tolerance);
    bool s1_out_of_tol = (dist_ts1 > tolerance);
    if (nt01_x_t0s0 * nt01_x_t0s1 > 0 && (s0_out_of_tol || s1_out_of_tol))
        return false;  // tgt segment is to far from src segment

    double S0[2] = {s0[X], s0[Y]};
    double S1[2] = {s1[X], s1[Y]};
    if (s0_out_of_tol)  // put s0 within tolerance
    {
        double t = tolerance * nt01_x_t0s0 / dist_ts0;  // signed tolerance
        double r = (nt01_x_t0s0 - t) / (nt01_x_t0s0 - nt01_x_t0s1);
        S0[X] = s0[X] * (1. - r) + s1[X] * r;
        S0[Y] = s0[Y] * (1. - r) + s1[Y] * r;
    }
    if (s1_out_of_tol)  // put s1 within tolerance
    {
        double t = tolerance * nt01_x_t0s1 / dist_ts1;  // signed tolerance
        double r = (nt01_x_t0s1 - t) / (nt01_x_t0s1 - nt01_x_t0s0);
        S1[X] = s1[X] * (1. - r) + s0[X] * r;
        S1[Y] = s1[Y] * (1. - r) + s0[Y] * r;
    }

    // project tgt and src segments to median line

    double s01[2] = {S1[X] - S0[X], S1[Y] - S0[Y]};          // src segment direction
    double sSize = sqrt(s01[X] * s01[X] + s01[Y] * s01[Y]);  // src segment size
    if (sSize < precision)
        return false;                  // degenerated segment
    s01[X] /= sSize, s01[Y] /= sSize;  // normalize s01

    // make t01 and s01 codirected
    double t01_x_s01 = t01[X] * s01[X] + t01[Y] * s01[Y];  // t01 dot s01
    if (t01_x_s01 < 0)
        s01[X] = -s01[X], s01[Y] = -s01[Y];

    double medianDir[2] = {
        t01[X] * (1. - medianLine) + s01[X] * medianLine, t01[Y] * (1. - medianLine) + s01[Y] * medianLine
    };
    double medianSize = sqrt(medianDir[X] * medianDir[X] + medianDir[Y] * medianDir[Y]);
    if (medianSize < std::numeric_limits<double>::min())
        return false;  // strange...
    medianDir[X] /= medianSize, medianDir[Y] /= medianSize;

    xt0 = t0[X] * medianDir[X] + t0[Y] * medianDir[Y];
    xt1 = t1[X] * medianDir[X] + t1[Y] * medianDir[Y];
    xs0 = S0[X] * medianDir[X] + S0[Y] * medianDir[Y];
    xs1 = S1[X] * medianDir[X] + S1[Y] * medianDir[Y];
    return true;
}
