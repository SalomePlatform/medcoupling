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
// Author : Anthony Geay

#include "MEDCouplingMatrix.hxx"

#include "InterpKernelMatrixTools.hxx"

#include <sstream>

using namespace MEDCoupling;

DenseMatrix *DenseMatrix::New(int nbRows, int nbCols)
{
  return new DenseMatrix(nbRows,nbCols);
}

DenseMatrix *DenseMatrix::New(DataArrayDouble *array, int nbRows, int nbCols)
{
  return new DenseMatrix(array,nbRows,nbCols);
}

DenseMatrix *DenseMatrix::deepCopy() const
{
  MCAuto<DataArrayDouble> arr(getData()->deepCopy());
  MCAuto<DenseMatrix> ret(DenseMatrix::New(arr,getNumberOfRows(),getNumberOfCols()));
  return ret.retn();
}

DenseMatrix *DenseMatrix::shallowCpy() const
{
  MCAuto<DenseMatrix> ret(DenseMatrix::New(const_cast<DataArrayDouble *>(getData()),getNumberOfRows(),getNumberOfCols()));
  return ret.retn();
}

std::size_t DenseMatrix::getHeapMemorySizeWithoutChildren() const
{
  return sizeof(DenseMatrix);
}

std::vector<const BigMemoryObject *> DenseMatrix::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back((const DataArrayDouble *)_data);
  return ret;
}

void DenseMatrix::updateTime() const
{
  const DataArrayDouble *pt(_data);
  if(pt)
    updateTimeWith(*pt);
}

/*!
 * This method scratch \a this to use a new input. The shape of \a this can be modified freely without any constraints.
 *
 * \param [in] array - The array containing data that is expected to be taken as new data.
 * \param [in] nbRows - The new number of rows (>0 or -1). If -1, the current number of rows will be taken.
 * \param [in] nbCols - The new number of columns (>0 or -1). If -1, the current number of cols will be taken.
 *
 * \sa reShape
 */
void DenseMatrix::reBuild(DataArrayDouble *array, int nbRows, int nbCols)
{
  int nbr(getNumberOfRowsExt(nbRows)),nbc(getNumberOfColsExt(nbCols));
  CheckArraySizes(array,nbr,nbc);
  DataArrayDouble *data(_data);
  if(data!=array)
    {
      _data=array; _data->incrRef();
      declareAsNew();
    }
  if(nbr!=_nb_rows)
    {
      _nb_rows=nbr;
      declareAsNew();
    }
  if(nbc!=_nb_cols)
    {
      _nb_cols=nbc;
      declareAsNew();
    }
}

/*!
 * This method does \b not change the content of the data in \a this. It only changes the shape (with a same number of elements in the matrix).
 * If the number of elements needs to be changed call reBuild method instead.
 *
 * \param [in] nbRows - The new number of rows (>0)
 * \param [in] nbCols - The new number of columns (>0)
 * \throw if the \c nbRows*nbCols is not equal to \c this->getNbOfElems()
 * \sa reBuild
 */
void DenseMatrix::reShape(int nbRows, int nbCols)
{
  if(nbRows<0 || nbCols<0)
    throw INTERP_KERNEL::Exception("DenseMatrix::reShape : number of rows and number of cols must be > 0 both !");
  if(nbRows*nbCols!=getNbOfElems())
    throw INTERP_KERNEL::Exception("DenseMatrix::reShape : This method is designed to change only the shape ! Number of elements must remain the same !");
  if(_nb_rows!=nbRows)
    {
      _nb_rows=nbRows;
      declareAsNew();
    }
  if(_nb_cols!=nbCols)
    {
      _nb_cols=nbCols;
      declareAsNew();
    }
}

void DenseMatrix::transpose()
{
  const MemArray<double>& mem(getData()->accessToMemArray());
  double *pt(mem.toNoInterlace(getNumberOfCols()));
  std::copy(pt,pt+getNbOfElems(),getData()->getPointer());//declareAsNew done here automatically by getPointer
  free(pt);
  std::swap(_nb_rows,_nb_cols);
  updateTime();
}

bool DenseMatrix::isEqual(const DenseMatrix& other, double eps) const
{
  std::string tmp;
  return isEqualIfNotWhy(other,eps,tmp);
}

bool DenseMatrix::isEqualIfNotWhy(const DenseMatrix& other, double eps, std::string& reason) const
{
  if(_nb_rows!=other._nb_rows)
    {
      std::ostringstream oss; oss << "Number of rows differs (" << _nb_rows << "!=" << other._nb_rows << ") !";
      reason+=oss.str();
      return false;
    }
  if(_nb_cols!=other._nb_cols)
      {
        std::ostringstream oss; oss << "Number of cols differs (" << _nb_cols << "!=" << other._nb_cols << ") !";
        reason+=oss.str();
        return false;
      }
  std::string tmp1;
  if(!_data->isEqualIfNotWhy(*other._data,eps,tmp1))
    {
      reason+="Data differs : "+tmp1;
      return false;
    }
  return true;
}

DataArrayDouble *DenseMatrix::matVecMult(const DataArrayDouble *vec) const
{
  return MatVecMult(this,vec);
}

DataArrayDouble *DenseMatrix::MatVecMult(const DenseMatrix *mat, const DataArrayDouble *vec)
{
  if(!mat || !vec)
    throw INTERP_KERNEL::Exception("DenseMatrix::MatVecMult : input matrix or vec is NULL !");
  vec->checkAllocated();
  if(vec->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DenseMatrix::MatVecMult : input vector must have only one component !");
  if(vec->getNumberOfTuples()!=mat->getNumberOfCols())
    throw INTERP_KERNEL::Exception("DenseMatrix::MatVecMult : Number of columns of this must be equal to number of tuples of vec !");
  MCAuto<DataArrayDouble> ret(DataArrayDouble::New()); ret->alloc(mat->getNumberOfRows(),1);
  INTERP_KERNEL::matrixProduct(mat->getData()->begin(),mat->getNumberOfRows(),mat->getNumberOfCols(),vec->begin(),vec->getNumberOfTuples(),1,ret->getPointer());
  return ret.retn();
}

DenseMatrix *DenseMatrix::Add(const DenseMatrix *a1, const DenseMatrix *a2)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DenseMatrix::Add : input matrices must be not NULL !");
  CheckSameSize(a1,a2);
  MCAuto<DataArrayDouble> data(DataArrayDouble::Add(a1->getData(),a2->getData()));
  MCAuto<DenseMatrix> ret(DenseMatrix::New(data,a1->getNumberOfRows(),a1->getNumberOfCols()));
  return ret.retn();
}

void DenseMatrix::addEqual(const DenseMatrix *other)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DenseMatrix::addEqual : other must be not NULL !");
  CheckSameSize(this,other);
  getData()->addEqual(other->getData());
}

DenseMatrix *DenseMatrix::Substract(const DenseMatrix *a1, const DenseMatrix *a2)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DenseMatrix::Substract : input matrices must be not NULL !");
  CheckSameSize(a1,a2);
  MCAuto<DataArrayDouble> data(DataArrayDouble::Substract(a1->getData(),a2->getData()));
  MCAuto<DenseMatrix> ret(DenseMatrix::New(data,a1->getNumberOfRows(),a1->getNumberOfCols()));
  return ret.retn();
}

void DenseMatrix::substractEqual(const DenseMatrix *other)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DenseMatrix::substractEqual : other must be not NULL !");
  CheckSameSize(this,other);
  getData()->substractEqual(other->getData());
}

DenseMatrix *DenseMatrix::Multiply(const DenseMatrix *a1, const DenseMatrix *a2)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DenseMatrix::Multiply : input matrices must be not NULL !");
  CheckCompatibleSizeForMul(a1,a2);
  int nbr(a1->getNumberOfRows()),nbc(a2->getNumberOfCols());
  MCAuto<DataArrayDouble> data(DataArrayDouble::New()); data->alloc(nbr*nbc,1);
  MCAuto<DenseMatrix> ret(DenseMatrix::New(data,a1->getNumberOfRows(),a2->getNumberOfCols()));
  INTERP_KERNEL::matrixProduct(a1->getData()->begin(),a1->getNumberOfRows(),a1->getNumberOfCols(),a2->getData()->begin(),a2->getNumberOfRows(),a2->getNumberOfCols(),data->getPointer());
  return ret.retn();
}

DenseMatrix *DenseMatrix::Multiply(const DenseMatrix *a1, const DataArrayDouble *a2)
{
  if(!a1 || !a2 || !a2->isAllocated())
    throw INTERP_KERNEL::Exception("DenseMatrix::Multiply #2 : input matrices must be not NULL and a2 allocated !");
  if(a2->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DenseMatrix::Multiply #2 : The 2nd member must have exactly one component !");
  MCAuto<DenseMatrix> a2Bis(DenseMatrix::New(const_cast<DataArrayDouble *>(a2),a2->getNumberOfTuples(),1));
  return DenseMatrix::Multiply(a1,a2Bis);
}

DenseMatrix::~DenseMatrix()
{
}

DenseMatrix::DenseMatrix(int nbRows, int nbCols):_nb_rows(nbRows),_nb_cols(nbCols),_data(DataArrayDouble::New())
{
  if(_nb_rows<0 || _nb_cols<0)
    throw INTERP_KERNEL::Exception("constructor of DenseMatrix : number of rows and number of cols must be > 0 both !");
  int nbOfTuples(_nb_rows*_nb_cols);
  _data->alloc(nbOfTuples,1);
}

DenseMatrix::DenseMatrix(DataArrayDouble *array, int nbRows, int nbCols):_nb_rows(nbRows),_nb_cols(nbCols)
{
  CheckArraySizes(array,_nb_rows,_nb_cols);
  _data=array; _data->incrRef();
}

int DenseMatrix::getNumberOfRowsExt(int nbRows) const
{
  if(nbRows<-1)
    throw INTERP_KERNEL::Exception("DenseMatrix::getNumberOfRowsExt : invalid input must be >= -1 !");
  if(nbRows==-1)
    return _nb_rows;
  else
    return nbRows;
}

int DenseMatrix::getNumberOfColsExt(int nbCols) const
{
  if(nbCols<-1)
    throw INTERP_KERNEL::Exception("DenseMatrix::getNumberOfColsExt : invalid input must be >= -1 !");
  if(nbCols==-1)
    return _nb_cols;
  else
    return nbCols;
}

void DenseMatrix::checkValidData() const
{
  if(!getData())
    throw INTERP_KERNEL::Exception("DenseMatrix::checkValidData : data is NULL !");
  if(!getData()->isAllocated())
    throw INTERP_KERNEL::Exception("DenseMatrix::checkValidData : data is not allocated !");
  if(getData()->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DenseMatrix::checkValidData : data has not 1 component !");
}

void DenseMatrix::CheckArraySizes(DataArrayDouble *array, int nbRows, int nbCols)
{
  if(nbRows<0 || nbCols<0)
    throw INTERP_KERNEL::Exception("constructor #2 of DenseMatrix : number of rows and number of cols must be > 0 both !");
  if(!array || !array->isAllocated())
    throw INTERP_KERNEL::Exception("constructor #2 of DenseMatrix : input array is empty or not allocated !");
  if(array->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("constructor #2 of DenseMatrix : input array must have exactly one component !");
  std::size_t nbr((std::size_t)nbRows),nbc((std::size_t)nbCols);
  if(nbr*nbc!=array->getNbOfElems())
    throw INTERP_KERNEL::Exception("constructor #2 of DenseMatrix : the number of elems in input array is not equal to the product of nbRows and nbCols !");
}

void DenseMatrix::CheckSameSize(const DenseMatrix *a1, const DenseMatrix *a2)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DenseMatrix::CheckSameSize : a1 or a2 is NULL !");
  a1->checkValidData(); a2->checkValidData();
  if(a1->getNumberOfRows()!=a2->getNumberOfRows())
    throw INTERP_KERNEL::Exception("DenseMatrix::CheckSameSize : number of rows mismatches !");
  if(a1->getNumberOfCols()!=a2->getNumberOfCols())
    throw INTERP_KERNEL::Exception("DenseMatrix::CheckSameSize : number of columns mismatches !");

}

void DenseMatrix::CheckCompatibleSizeForMul(const DenseMatrix *a1, const DenseMatrix *a2)
{
  if(!a1 || !a2)
      throw INTERP_KERNEL::Exception("DenseMatrix::CheckCompatibleSizeForMul : a1 or a2 is NULL !");
  a1->checkValidData(); a2->checkValidData();
  if(a1->getNumberOfCols()!=a2->getNumberOfRows())
    throw INTERP_KERNEL::Exception("DenseMatrix::CheckCompatibleSizeForMul : number of cols of a1 must be equal to number of rows of a2 !");
}

