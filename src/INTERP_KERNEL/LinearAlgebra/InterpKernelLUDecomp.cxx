// Copyright (C) 2022-2026  CEA, EDF
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

// Implementation coming from Numerical Recipes in C of 1994 (version 2.04)

#include "InterpKernelLUDecomp.hxx"
#include "InterpKernelException.hxx"

#include <cmath>
#include <sstream>

using namespace INTERP_KERNEL;

LUDecomp::LUDecomp(const INTERP_KERNEL::DenseMatrix &a) : n(a.nrows()), lu(a), aref(a), indx(n)
{
    const double TINY = 1.0e-40;
    mcIdType i, imax, j, k;
    double big, temp;
    std::vector<double> vv(n);
    d = 1.0;
    for (i = 0; i < n; i++)
    {
        big = 0.0;
        for (j = 0; j < n; j++)
            if ((temp = std::abs(lu[i][j])) > big)
                big = temp;
        if (big == 0.0)
            THROW_IK_EXCEPTION("Singular matrix in LUDecomp");
        vv[i] = 1.0 / big;
    }
    for (k = 0; k < n; k++)
    {
        big = 0.0;
        imax = k;
        for (i = k; i < n; i++)
        {
            temp = vv[i] * std::abs(lu[i][k]);
            if (temp > big)
            {
                big = temp;
                imax = i;
            }
        }
        if (k != imax)
        {
            for (j = 0; j < n; j++)
            {
                temp = lu[imax][j];
                lu[imax][j] = lu[k][j];
                lu[k][j] = temp;
            }
            d = -d;
            vv[imax] = vv[k];
        }
        indx[k] = imax;
        if (lu[k][k] == 0.0)
            lu[k][k] = TINY;
        for (i = k + 1; i < n; i++)
        {
            temp = lu[i][k] /= lu[k][k];
            for (j = k + 1; j < n; j++) lu[i][j] -= temp * lu[k][j];
        }
    }
}
void
LUDecomp::solve(const std::vector<double> &b, std::vector<double> &x)
{
    mcIdType i, ii = 0, ip, j;
    double sum;
    if (b.size() != ToSizeT(n) || x.size() != ToSizeT(n))
        THROW_IK_EXCEPTION("LUDecomp::solve bad sizes");
    for (i = 0; i < n; i++) x[i] = b[i];
    for (i = 0; i < n; i++)
    {
        ip = indx[i];
        sum = x[ip];
        x[ip] = x[i];
        if (ii != 0)
            for (j = ii - 1; j < i; j++) sum -= lu[i][j] * x[j];
        else if (sum != 0.0)
            ii = i + 1;
        x[i] = sum;
    }
    for (i = n - 1; i >= 0; i--)
    {
        sum = x[i];
        for (j = i + 1; j < n; j++) sum -= lu[i][j] * x[j];
        x[i] = sum / lu[i][i];
    }
}

void
LUDecomp::solve(const INTERP_KERNEL::DenseMatrix &b, INTERP_KERNEL::DenseMatrix &x)
{
    mcIdType i, j, m(b.ncols());
    if (b.nrows() != n || x.nrows() != n || b.ncols() != x.ncols())
        THROW_IK_EXCEPTION("LUDecomp::solve bad sizes");
    std::vector<double> xx(n);
    for (j = 0; j < m; j++)
    {
        for (i = 0; i < n; i++) xx[i] = b[i][j];
        solve(xx, xx);
        for (i = 0; i < n; i++) x[i][j] = xx[i];
    }
}

void
LUDecomp::inverse(INTERP_KERNEL::DenseMatrix &ainv)
{
    mcIdType i, j;
    ainv.resize(n, n);
    for (i = 0; i < n; i++)
    {
        for (j = 0; j < n; j++) ainv[i][j] = 0.;
        ainv[i][i] = 1.;
    }
    solve(ainv, ainv);
}

double
LUDecomp::det()
{
    double dd = d;
    for (mcIdType i = 0; i < n; i++) dd *= lu[i][i];
    return dd;
}

void
LUDecomp::mprove(const std::vector<double> &b, std::vector<double> &x)
{
    mcIdType i, j;
    std::vector<double> r(n);
    for (i = 0; i < n; i++)
    {
        long double sdp = -b[i];
        for (j = 0; j < n; j++) sdp += (long double)aref[i][j] * (long double)x[j];
        r[i] = double(sdp);
    }
    solve(r, r);
    for (i = 0; i < n; i++) x[i] -= r[i];
}
