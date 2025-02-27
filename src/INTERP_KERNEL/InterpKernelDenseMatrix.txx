// Copyright (C) 2022-2025  CEA, EDF
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
#include "InterpKernelException.hxx"
#include "VectorUtils.hxx"

#include <cmath>

namespace INTERP_KERNEL
{
  template <class T>
  DenseMatrixT<T>::DenseMatrixT() : nn(0), mm(0), v(nullptr) {}

  template <class T>
  DenseMatrixT<T>::DenseMatrixT(mcIdType n, mcIdType m) : nn(n), mm(m), v(n>0 ? new T*[n] : nullptr)
  {
    mcIdType i,nel=m*n;
    if (v) v[0] = nel>0 ? new T[nel] : nullptr;
    for (i=1;i<n;i++) v[i] = v[i-1] + m;
  }

  template <class T>
  DenseMatrixT<T>::DenseMatrixT(mcIdType n, mcIdType m, const T &a) : nn(n), mm(m), v(n>0 ? new T*[n] : nullptr)
  {
    mcIdType i,j,nel=m*n;
    if (v) v[0] = nel>0 ? new T[nel] : nullptr;
    for (i=1; i< n; i++) v[i] = v[i-1] + m;
    for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = a;
  }

  template <class T>
  DenseMatrixT<T>::DenseMatrixT(mcIdType n, mcIdType m, const T *a) : nn(n), mm(m), v(n>0 ? new T*[n] : nullptr)
  {
    mcIdType i,j,nel=m*n;
    if (v) v[0] = nel>0 ? new T[nel] : nullptr;
    for (i=1; i< n; i++) v[i] = v[i-1] + m;
    for (i=0; i< n; i++) for (j=0; j<m; j++) v[i][j] = *a++;
  }

  template <class T>
  DenseMatrixT<T>::DenseMatrixT(const DenseMatrixT &rhs) : nn(rhs.nn), mm(rhs.mm), v(nn>0 ? new T*[nn] : nullptr)
  {
    mcIdType i,j,nel=mm*nn;
    if (v) v[0] = nel>0 ? new T[nel] : nullptr;
    for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
    for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
  }

  /*!
  * postcondition: normal assignment via copying has been performed
  * if matrix and rhs were different sizes, matrix
  * has been resized to match the size of rhs
  */
  template <class T>
  DenseMatrixT<T> & DenseMatrixT<T>::operator=(const DenseMatrixT<T> &rhs)
  {
    if (this != &rhs) {
      mcIdType i,j,nel;
      if (nn != rhs.nn || mm != rhs.mm) {
        if ( v ) {
          delete[] (v[0]);
          delete[] (v);
        }
        nn=rhs.nn;
        mm=rhs.mm;
        v = nn>0 ? new T*[nn] : nullptr;
        nel = mm*nn;
        if (v) v[0] = nel>0 ? new T[nel] : nullptr;
        for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
      }
      for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = rhs[i][j];
    }
    return *this;
  }

  template <class T>
  void DenseMatrixT<T>::resize(mcIdType newn, mcIdType newm)
  {
    mcIdType i,nel;
    if (newn != nn || newm != mm)
    {
      if ( v )
      {
        delete[] (v[0]);
        delete[] (v);
      }
      nn = newn;
      mm = newm;
      v = nn>0 ? new T*[nn] : nullptr;
      nel = mm*nn;
      if (v) v[0] = nel>0 ? new T[nel] : nullptr;
      for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
    }
  }

  template <class T>
  void DenseMatrixT<T>::assign(mcIdType newn, mcIdType newm, const T& a)
  {
    mcIdType i,j,nel;
    if (newn != nn || newm != mm)
    {
      if ( v )
      {
        delete[] (v[0]);
        delete[] (v);
      }
      nn = newn;
      mm = newm;
      v = nn>0 ? new T*[nn] : nullptr;
      nel = mm*nn;
      if (v) v[0] = nel>0 ? new T[nel] : nullptr;
      for (i=1; i< nn; i++) v[i] = v[i-1] + mm;
    }
    for (i=0; i< nn; i++) for (j=0; j<mm; j++) v[i][j] = a;
  }

  template <class T>
  DenseMatrixT<T>::~DenseMatrixT()
  {
    if ( v )
    {
      delete[] (v[0]);
      delete[] (v);
    }
  }

  template<class T>
  T Determinant22(const T *m)
  { return m[0]*m[3]-m[1]*m[2]; }

  template<class T>
  T Determinant33(const T *m)
  { return m[0]*(m[4]*m[8]-m[7]*m[5])-m[1]*(m[3]*m[8]-m[6]*m[5])+m[2]*(m[3]*m[7]-m[6]*m[4]);}

  template <class T>
  T DenseMatrixT<T>::determinant() const
  {
    if(nn==1 && mm==1)
      return v[0][0];
    if(nn==2 && mm==2)
      return Determinant22(v[0]);
    if(nn==3 && mm==3)
      return Determinant33(v[0]);
    THROW_IK_EXCEPTION("DenseMatrixT::determinant : only 1x1, 2x2 and 3x3 implemented !");
  }

  template <class T>
  T DenseMatrixT<T>::toJacobian() const
  {
    if(nn == mm)
      return determinant();
    const T *vPtr(this->v[0]);
    if(nn==3 && mm==1)
      return norm(vPtr);
    if(nn==2 && mm==1)
      return std::sqrt(vPtr[0]*vPtr[0] + vPtr[1]*vPtr[1]);
    if(nn==3 && mm==2)
    {
      T tmp[3];
      T VA[3] = {v[0][0],v[1][0],v[2][0]};
      T VB[3] = {v[0][1],v[1][1],v[2][1]};
      cross(VA,VB,tmp);
      return norm(tmp);
    }
    THROW_IK_EXCEPTION("DenseMatrixT::toJacobian : only 1x1, 2x1, 3x1, 3x2, 2x2 and 3x3 implemented !");
  }
}
