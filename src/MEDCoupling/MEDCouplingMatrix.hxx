// Copyright (C) 2007-2024  CEA, EDF
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
// Author : Anthony Geay

#ifndef __PARAMEDMEM_MEDCOUPLINGMATRIX_HXX__
#define __PARAMEDMEM_MEDCOUPLINGMATRIX_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingTimeLabel.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MCAuto.hxx"

#include "InterpKernelException.hxx"

namespace MEDCoupling
{
  /*!
   * The aim of this class is \b NOT to reimplement all linear algebra but only to store a dense matrix.
   * It only provides basic set/get and basic operations and bindings to linear algebra libraries (numpy/scipy) and a compatible format to Petsc.
   */
  class DenseMatrix : public RefCountObject, public TimeLabel
  {
  public:
    MEDCOUPLING_EXPORT static DenseMatrix *New(mcIdType nbRows, mcIdType nbCols);
    MEDCOUPLING_EXPORT static DenseMatrix *New(DataArrayDouble *array, mcIdType nbRows, mcIdType nbCols);
    MEDCOUPLING_EXPORT std::string getClassName() const override { return std::string("DenseMatrix"); }
    MEDCOUPLING_EXPORT DenseMatrix *deepCopy() const;
    MEDCOUPLING_EXPORT DenseMatrix *shallowCpy() const;
    MEDCOUPLING_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDCOUPLING_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDCOUPLING_EXPORT void updateTime() const;
    //
    MEDCOUPLING_EXPORT mcIdType getNumberOfRows() const { return _nb_rows; }
    MEDCOUPLING_EXPORT mcIdType getNumberOfCols() const { return _nb_cols; }
    MEDCOUPLING_EXPORT mcIdType getNbOfElems() const { return _nb_rows*_nb_cols; }
    MEDCOUPLING_EXPORT void reBuild(DataArrayDouble *array, mcIdType nbRows=-1, mcIdType nbCols=-1);
    MEDCOUPLING_EXPORT void reShape(mcIdType nbRows, mcIdType nbCols);
    MEDCOUPLING_EXPORT void transpose();
    //
    MEDCOUPLING_EXPORT bool isEqual(const DenseMatrix& other, double eps) const;
    MEDCOUPLING_EXPORT bool isEqualIfNotWhy(const DenseMatrix& other, double eps, std::string& reason) const;
    MEDCOUPLING_EXPORT DataArrayDouble *matVecMult(const DataArrayDouble *vec) const;
    MEDCOUPLING_EXPORT static DataArrayDouble *MatVecMult(const DenseMatrix *mat, const DataArrayDouble *vec);
    MEDCOUPLING_EXPORT static DenseMatrix *Add(const DenseMatrix *a1, const DenseMatrix *a2);
    MEDCOUPLING_EXPORT void addEqual(const DenseMatrix *other);
    MEDCOUPLING_EXPORT static DenseMatrix *Substract(const DenseMatrix *a1, const DenseMatrix *a2);
    MEDCOUPLING_EXPORT void substractEqual(const DenseMatrix *other);
    MEDCOUPLING_EXPORT static DenseMatrix *Multiply(const DenseMatrix *a1, const DenseMatrix *a2);
    MEDCOUPLING_EXPORT static DenseMatrix *Multiply(const DenseMatrix *a1, const DataArrayDouble *a2);
    //
    MEDCOUPLING_EXPORT const DataArrayDouble *getData() const { return _data; }
    MEDCOUPLING_EXPORT DataArrayDouble *getData() { return _data; }
  private:
    ~DenseMatrix();
    DenseMatrix(mcIdType nbRows, mcIdType nbCols);
    DenseMatrix(DataArrayDouble *array, mcIdType nbRows, mcIdType nbCols);
    mcIdType getNumberOfRowsExt(mcIdType nbRows) const;
    mcIdType getNumberOfColsExt(mcIdType nbCols) const;
    void checkValidData() const;
    static void CheckArraySizes(DataArrayDouble *array, mcIdType nbRows, mcIdType nbCols);
    static void CheckSameSize(const DenseMatrix *a1, const DenseMatrix *a2);
    static void CheckCompatibleSizeForMul(const DenseMatrix *a1, const DenseMatrix *a2);
  private:
    mcIdType _nb_rows;
    mcIdType _nb_cols;
    MCAuto<DataArrayDouble> _data;
  };
}

#endif
