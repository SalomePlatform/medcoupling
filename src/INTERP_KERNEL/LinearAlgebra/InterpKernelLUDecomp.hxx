// Copyright (C) 2022  CEA/DEN, EDF R&D
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

#include "InterpKernelDenseMatrix.hxx"
#include <vector>

namespace INTERP_KERNEL
{
  class LUDecomp
  {
  private:
    mcIdType n;
    INTERP_KERNEL::DenseMatrix lu;
    std::vector<mcIdType> indx;
    double d;
    const INTERP_KERNEL::DenseMatrix& aref;
  public:
    LUDecomp (const INTERP_KERNEL::DenseMatrix& a);
    void solve(const std::vector<double>& b, std::vector<double> &x);
    void solve(const INTERP_KERNEL::DenseMatrix& b, INTERP_KERNEL::DenseMatrix &x);
    void inverse(INTERP_KERNEL::DenseMatrix &ainv);
    double det();
    void mprove(const std::vector<double>& b, std::vector<double> &x);
  };
}
