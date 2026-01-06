// Copyright (C) 2007-2026  CEA, EDF
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

#ifndef __INTERPKERNELMATRIXTOOLS_HXX__
#define __INTERPKERNELMATRIXTOOLS_HXX__

#include "INTERPKERNELDefines.hxx"
#include "MCIdType.hxx"

namespace INTERP_KERNEL
{
void INTERPKERNEL_EXPORT
matrixProduct(const double *A, mcIdType n1, mcIdType p1, const double *B, mcIdType n2, mcIdType p2, double *C);
void INTERPKERNEL_EXPORT
inverseMatrix(const double *A, mcIdType n, double *iA);
void INTERPKERNEL_EXPORT
daxpy(mcIdType n, double da, const double *dx, mcIdType incx, double *dy, mcIdType incy);
}  // namespace INTERP_KERNEL

#endif
