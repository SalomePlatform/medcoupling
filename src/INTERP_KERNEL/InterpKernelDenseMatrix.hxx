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

#include "INTERPKERNELDefines.hxx"
#include "MCIdType.hxx"

namespace INTERP_KERNEL
{
  template <class T>
  class INTERPKERNEL_EXPORT DenseMatrixT
  {
  private:
    mcIdType nn;
    mcIdType mm;
    T **v;
  public:
    DenseMatrixT();
    DenseMatrixT(mcIdType n, mcIdType m);
    DenseMatrixT(mcIdType n, mcIdType m, const T &a);	
    DenseMatrixT(mcIdType n, mcIdType m, const T *a);
    DenseMatrixT(const DenseMatrixT &rhs);
    DenseMatrixT & operator=(const DenseMatrixT &rhs);
    using value_type = T;
    //! subscripting: pointer to row i
    T* operator[](const mcIdType i) { return v[i]; }
    const T* operator[](const mcIdType i) const { return v[i]; }
    mcIdType nrows() const { return nn; }
    mcIdType ncols() const { return mm; }
    void resize(mcIdType newn, mcIdType newm);
    void assign(mcIdType newn, mcIdType newm, const T &a);
    ~DenseMatrixT();
    T determinant() const;
    T toJacobian() const;
  };

  using DenseMatrix = DenseMatrixT<double>;
}
