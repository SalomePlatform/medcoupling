// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#include "MEDCouplingMemArray.txx"
#include "MEDCouplingAutoRefCountObjectPtr.hxx"

#include <set>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <functional>

using namespace ParaMEDMEM;

/*!
 * Checks if raw data is allocated. Read more on the raw data
 * in \ref MEDCouplingArrayBasicsTuplesAndCompo "DataArrays infos" for more information.
 *  \return bool - \a true if the raw data is allocated, \a false else.
 */
bool DataArrayChar::isAllocated() const
{
  return getConstPointer()!=0;
}

/*!
 * Checks if raw data is allocated and throws an exception if it is not the case.
 *  \throw If the raw data is not allocated.
 */
void DataArrayChar::checkAllocated() const
{
  if(!isAllocated())
    throw INTERP_KERNEL::Exception("DataArrayChar::checkAllocated : Array is defined but not allocated ! Call alloc or setValues method first !");
}

/*!
 * This method desallocated \a this without modification of informations relative to the components.
 * After call of this method, DataArrayChar::isAllocated will return false.
 * If \a this is already not allocated, \a this is let unchanged.
 */
void DataArrayChar::desallocate()
{
  _mem.destroy();
}

std::size_t DataArrayChar::getHeapMemorySizeWithoutChildren() const
{
  std::size_t sz(_mem.getNbOfElemAllocated());
  return DataArray::getHeapMemorySizeWithoutChildren()+sz;
}

/*!
 * Returns an integer value characterizing \a this array, which is useful for a quick
 * comparison of many instances of DataArrayInt.
 *  \return int - the hash value.
 *  \throw If \a this is not allocated.
 */
int DataArrayChar::getHashCode() const
{
  checkAllocated();
  std::size_t nbOfElems=getNbOfElems();
  int ret=nbOfElems*65536;
  int delta=3;
  if(nbOfElems>48)
    delta=nbOfElems/8;
  int ret0=0;
  const char *pt=begin();
  for(std::size_t i=0;i<nbOfElems;i+=delta)
    ret0+=pt[i];
  return ret+ret0;
}

/*!
 * Checks the number of tuples.
 *  \return bool - \a true if getNumberOfTuples() == 0, \a false else.
 *  \throw If \a this is not allocated.
 */
bool DataArrayChar::empty() const
{
  checkAllocated();
  return getNumberOfTuples()==0;
}

/*!
 * Copies all the data from another DataArrayChar. For more info see
 * \ref MEDCouplingArrayBasicsCopyDeepAssign.
 *  \param [in] other - another instance of DataArrayChar to copy data from.
 *  \throw If the \a other is not allocated.
 */
void DataArrayChar::cpyFrom(const DataArrayChar& other)
{
  other.checkAllocated();
  int nbOfTuples=other.getNumberOfTuples();
  int nbOfComp=other.getNumberOfComponents();
  allocIfNecessary(nbOfTuples,nbOfComp);
  std::size_t nbOfElems=(std::size_t)nbOfTuples*nbOfComp;
  char *pt=getPointer();
  const char *ptI=other.getConstPointer();
  for(std::size_t i=0;i<nbOfElems;i++)
    pt[i]=ptI[i];
  copyStringInfoFrom(other);
}

/*!
 * This method reserve nbOfElems elements in memory ( nbOfElems*4 bytes ) \b without impacting the number of tuples in \a this.
 * If \a this has already been allocated, this method checks that \a this has only one component. If not an INTERP_KERNEL::Exception will be thrown.
 * If \a this has not already been allocated, number of components is set to one.
 * This method allows to reduce number of reallocations on invokation of DataArrayChar::pushBackSilent and DataArrayChar::pushBackValsSilent on \a this.
 * 
 * \sa DataArrayChar::pack, DataArrayChar::pushBackSilent, DataArrayChar::pushBackValsSilent
 */
void DataArrayChar::reserve(std::size_t nbOfElems)
{
  int nbCompo=getNumberOfComponents();
  if(nbCompo==1)
    {
      _mem.reserve(nbOfElems);
    }
  else if(nbCompo==0)
    {
      _mem.reserve(nbOfElems);
      _info_on_compo.resize(1);
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayChar::reserve : not available for DataArrayChar with number of components different than 1 !");
}

/*!
 * This method adds at the end of \a this the single value \a val. This method do \b not update its time label to avoid useless incrementation
 * of counter. So the caller is expected to call TimeLabel::declareAsNew on \a this at the end of the push session.
 *
 * \param [in] val the value to be added in \a this
 * \throw If \a this has already been allocated with number of components different from one.
 * \sa DataArrayChar::pushBackValsSilent
 */
void DataArrayChar::pushBackSilent(char val)
{
  int nbCompo=getNumberOfComponents();
  if(nbCompo==1)
    _mem.pushBack(val);
  else if(nbCompo==0)
    {
      _info_on_compo.resize(1);
      _mem.pushBack(val);
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayChar::pushBackSilent : not available for DataArrayChar with number of components different than 1 !");
}

/*!
 * This method adds at the end of \a this a serie of values [\c valsBg,\c valsEnd). This method do \b not update its time label to avoid useless incrementation
 * of counter. So the caller is expected to call TimeLabel::declareAsNew on \a this at the end of the push session.
 *
 *  \param [in] valsBg - an array of values to push at the end of \this.
 *  \param [in] valsEnd - specifies the end of the array \a valsBg, so that
 *              the last value of \a valsBg is \a valsEnd[ -1 ].
 * \throw If \a this has already been allocated with number of components different from one.
 * \sa DataArrayChar::pushBackSilent
 */
void DataArrayChar::pushBackValsSilent(const char *valsBg, const char *valsEnd)
{
  int nbCompo=getNumberOfComponents();
  if(nbCompo==1)
    _mem.insertAtTheEnd(valsBg,valsEnd);
  else if(nbCompo==0)
    {
      _info_on_compo.resize(1);
      _mem.insertAtTheEnd(valsBg,valsEnd);
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayChar::pushBackValsSilent : not available for DataArrayChar with number of components different than 1 !");
}

/*!
 * This method returns silently ( without updating time label in \a this ) the last value, if any and suppress it.
 * \throw If \a this is already empty.
 * \throw If \a this has number of components different from one.
 */
char DataArrayChar::popBackSilent()
{
  if(getNumberOfComponents()==1)
    return _mem.popBack();
  else
    throw INTERP_KERNEL::Exception("DataArrayChar::popBackSilent : not available for DataArrayChar with number of components different than 1 !");
}

/*!
 * This method \b do \b not modify content of \a this. It only modify its memory footprint if the allocated memory is to high regarding real data to store.
 *
 * \sa DataArrayChar::getHeapMemorySizeWithoutChildren, DataArrayChar::reserve
 */
void DataArrayChar::pack() const
{
  _mem.pack();
}

/*!
 * Allocates the raw data in memory. If exactly as same memory as needed already
 * allocated, it is not re-allocated.
 *  \param [in] nbOfTuple - number of tuples of data to allocate.
 *  \param [in] nbOfCompo - number of components of data to allocate.
 *  \throw If \a nbOfTuple < 0 or \a nbOfCompo < 0.
 */
void DataArrayChar::allocIfNecessary(int nbOfTuple, int nbOfCompo)
{
  if(isAllocated())
    {
      if(nbOfTuple!=getNumberOfTuples() || nbOfCompo!=getNumberOfComponents())
        alloc(nbOfTuple,nbOfCompo);
    }
  else
    alloc(nbOfTuple,nbOfCompo);
}

/*!
 * Allocates the raw data in memory. If the memory was already allocated, then it is
 * freed and re-allocated. See an example of this method use
 * \ref MEDCouplingArraySteps1WC "here".
 *  \param [in] nbOfTuple - number of tuples of data to allocate.
 *  \param [in] nbOfCompo - number of components of data to allocate.
 *  \throw If \a nbOfTuple < 0 or \a nbOfCompo < 0.
 */
void DataArrayChar::alloc(int nbOfTuple, int nbOfCompo)
{
  if(nbOfTuple<0 || nbOfCompo<0)
    throw INTERP_KERNEL::Exception("DataArrayChar::alloc : request for negative length of data !");
  _info_on_compo.resize(nbOfCompo);
  _mem.alloc(nbOfCompo*(std::size_t)nbOfTuple);
  declareAsNew();
}

/*!
 * Checks if \a this and another DataArrayChar are fully equal. For more info see
 * \ref MEDCouplingArrayBasicsCompare.
 *  \param [in] other - an instance of DataArrayChar to compare with \a this one.
 *  \return bool - \a true if the two arrays are equal, \a false else.
 */
bool DataArrayChar::isEqual(const DataArrayChar& other) const
{
  std::string tmp;
  return isEqualIfNotWhy(other,tmp);
}

/*!
 * Equivalent to DataArrayChar::isEqual except that if false the reason of
 * mismatch is given.
 * 
 * \param [in] other the instance to be compared with \a this
 * \param [out] reason In case of inequality returns the reason.
 * \sa DataArrayChar::isEqual
 */
bool DataArrayChar::isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const
{
  if(!areInfoEqualsIfNotWhy(other,reason))
    return false;
  return _mem.isEqual(other._mem,0,reason);
}

/*!
 * Checks if values of \a this and another DataArrayChar are equal. For more info see
 * \ref MEDCouplingArrayBasicsCompare.
 *  \param [in] other - an instance of DataArrayChar to compare with \a this one.
 *  \return bool - \a true if the values of two arrays are equal, \a false else.
 */
bool DataArrayChar::isEqualWithoutConsideringStr(const DataArrayChar& other) const
{
  std::string tmp;
  return _mem.isEqual(other._mem,0,tmp);
}

/*!
 * Reverse the array values.
 *  \throw If \a this->getNumberOfComponents() < 1.
 *  \throw If \a this is not allocated.
 */
void DataArrayChar::reverse()
{
  checkAllocated();
  _mem.reverse(getNumberOfComponents());
  declareAsNew();
}

/*!
 * Assign zero to all values in \a this array. To know more on filling arrays see
 * \ref MEDCouplingArrayFill.
 * \throw If \a this is not allocated.
 */
void DataArrayChar::fillWithZero()
{
  checkAllocated();
  _mem.fillWithValue(0);
  declareAsNew();
}

/*!
 * Assign \a val to all values in \a this array. To know more on filling arrays see
 * \ref MEDCouplingArrayFill.
 *  \param [in] val - the value to fill with.
 *  \throw If \a this is not allocated.
 */
void DataArrayChar::fillWithValue(char val)
{
  checkAllocated();
  _mem.fillWithValue(val);
  declareAsNew();
}

/*!
 * Returns a textual and human readable representation of \a this instance of
 * DataArrayChar. This text is shown when a DataArrayChar is printed in Python.
 *  \return std::string - text describing \a this DataArrayChar.
 */
std::string DataArrayChar::repr() const
{
  std::ostringstream ret;
  reprStream(ret);
  return ret.str();
}

std::string DataArrayChar::reprZip() const
{
  std::ostringstream ret;
  reprZipStream(ret);
  return ret.str();
}

/*!
 * Changes number of tuples in the array. If the new number of tuples is smaller
 * than the current number the array is truncated, otherwise the array is extended.
 *  \param [in] nbOfTuples - new number of tuples. 
 *  \throw If \a this is not allocated.
 *  \throw If \a nbOfTuples is negative.
 */
void DataArrayChar::reAlloc(int nbOfTuples)
{
  if(nbOfTuples<0)
    throw INTERP_KERNEL::Exception("DataArrayChar::reAlloc : input new number of tuples should be >=0 !");
  checkAllocated();
  _mem.reAlloc(getNumberOfComponents()*(std::size_t)nbOfTuples);
  declareAsNew();
}

/*!
 * Creates a new DataArrayInt and assigns all (textual and numerical) data of \a this
 * array to the new one.
 *  \return DataArrayInt * - the new instance of DataArrayChar.
 */
DataArrayInt *DataArrayChar::convertToIntArr() const
{
  checkAllocated();
  DataArrayInt *ret=DataArrayInt::New();
  ret->alloc(getNumberOfTuples(),getNumberOfComponents());
  std::size_t nbOfVals=getNbOfElems();
  const char *src=getConstPointer();
  int *dest=ret->getPointer();
  std::copy(src,src+nbOfVals,dest);
  ret->copyStringInfoFrom(*this);
  return ret;
}

/*!
 * Permutes values of \a this array as required by \a old2New array. The values are
 * permuted so that \c new[ \a old2New[ i ]] = \c old[ i ]. Number of tuples remains
 * the same as in \this one.
 * If a permutation reduction is needed, substr() or selectByTupleId() should be used.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] old2New - C array of length equal to \a this->getNumberOfTuples()
 *     giving a new position for i-th old value.
 */
void DataArrayChar::renumberInPlace(const int *old2New)
{
  checkAllocated();
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  char *tmp=new char[nbTuples*nbOfCompo];
  const char *iptr=getConstPointer();
  for(int i=0;i<nbTuples;i++)
    std::copy(iptr+nbOfCompo*i,iptr+nbOfCompo*(i+1),tmp+nbOfCompo*old2New[i]);
  std::copy(tmp,tmp+nbTuples*nbOfCompo,getPointer());
  delete [] tmp;
  declareAsNew();
}

/*!
 * Permutes values of \a this array as required by \a new2Old array. The values are
 * permuted so that \c new[ i ] = \c old[ \a new2Old[ i ]]. Number of tuples remains
 * the same as in \this one.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] new2Old - C array of length equal to \a this->getNumberOfTuples()
 *     giving a previous position of i-th new value.
 */
void DataArrayChar::renumberInPlaceR(const int *new2Old)
{
  checkAllocated();
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  char *tmp=new char[nbTuples*nbOfCompo];
  const char *iptr=getConstPointer();
  for(int i=0;i<nbTuples;i++)
    std::copy(iptr+nbOfCompo*new2Old[i],iptr+nbOfCompo*(new2Old[i]+1),tmp+nbOfCompo*i);
  std::copy(tmp,tmp+nbTuples*nbOfCompo,getPointer());
  delete [] tmp;
  declareAsNew();
}

/*!
 * Returns a copy of \a this array with values permuted as required by \a old2New array.
 * The values are permuted so that  \c new[ \a old2New[ i ]] = \c old[ i ].
 * Number of tuples in the result array remains the same as in \this one.
 * If a permutation reduction is needed, renumberAndReduce() should be used.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] old2New - C array of length equal to \a this->getNumberOfTuples()
 *          giving a new position for i-th old value.
 *  \return DataArrayChar * - the new instance of DataArrayChar that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 */
DataArrayChar *DataArrayChar::renumber(const int *old2New) const
{
  checkAllocated();
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  MEDCouplingAutoRefCountObjectPtr<DataArrayChar> ret=buildEmptySpecializedDAChar();
  ret->alloc(nbTuples,nbOfCompo);
  ret->copyStringInfoFrom(*this);
  const char *iptr=getConstPointer();
  char *optr=ret->getPointer();
  for(int i=0;i<nbTuples;i++)
    std::copy(iptr+nbOfCompo*i,iptr+nbOfCompo*(i+1),optr+nbOfCompo*old2New[i]);
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * Returns a copy of \a this array with values permuted as required by \a new2Old array.
 * The values are permuted so that  \c new[ i ] = \c old[ \a new2Old[ i ]]. Number of
 * tuples in the result array remains the same as in \this one.
 * If a permutation reduction is needed, substr() or selectByTupleId() should be used.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] new2Old - C array of length equal to \a this->getNumberOfTuples()
 *     giving a previous position of i-th new value.
 *  \return DataArrayChar * - the new instance of DataArrayChar that the caller
 *          is to delete using decrRef() as it is no more needed.
 */
DataArrayChar *DataArrayChar::renumberR(const int *new2Old) const
{
  checkAllocated();
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  MEDCouplingAutoRefCountObjectPtr<DataArrayChar> ret=buildEmptySpecializedDAChar();
  ret->alloc(nbTuples,nbOfCompo);
  ret->copyStringInfoFrom(*this);
  const char *iptr=getConstPointer();
  char *optr=ret->getPointer();
  for(int i=0;i<nbTuples;i++)
    std::copy(iptr+nbOfCompo*new2Old[i],iptr+nbOfCompo*(new2Old[i]+1),optr+nbOfCompo*i);
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * Returns a shorten and permuted copy of \a this array. The new DataArrayChar is
 * of size \a newNbOfTuple and it's values are permuted as required by \a old2New array.
 * The values are permuted so that  \c new[ \a old2New[ i ]] = \c old[ i ] for all
 * \a old2New[ i ] >= 0. In other words every i-th tuple in \a this array, for which 
 * \a old2New[ i ] is negative, is missing from the result array.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] old2New - C array of length equal to \a this->getNumberOfTuples()
 *     giving a new position for i-th old tuple and giving negative position for
 *     for i-th old tuple that should be omitted.
 *  \return DataArrayChar * - the new instance of DataArrayChar that the caller
 *          is to delete using decrRef() as it is no more needed.
 */
DataArrayChar *DataArrayChar::renumberAndReduce(const int *old2New, int newNbOfTuple) const
{
  checkAllocated();
  int nbTuples=getNumberOfTuples();
  int nbOfCompo=getNumberOfComponents();
  MEDCouplingAutoRefCountObjectPtr<DataArrayChar> ret=buildEmptySpecializedDAChar();
  ret->alloc(newNbOfTuple,nbOfCompo);
  const char *iptr=getConstPointer();
  char *optr=ret->getPointer();
  for(int i=0;i<nbTuples;i++)
    {
      int w=old2New[i];
      if(w>=0)
        std::copy(iptr+i*nbOfCompo,iptr+(i+1)*nbOfCompo,optr+w*nbOfCompo);
    }
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * Returns a shorten and permuted copy of \a this array. The new DataArrayChar is
 * of size \a new2OldEnd - \a new2OldBg and it's values are permuted as required by
 * \a new2OldBg array.
 * The values are permuted so that  \c new[ i ] = \c old[ \a new2OldBg[ i ]].
 * This method is equivalent to renumberAndReduce() except that convention in input is
 * \c new2old and \b not \c old2new.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] new2OldBg - pointer to the beginning of a permutation array that gives a
 *              tuple index in \a this array to fill the i-th tuple in the new array.
 *  \param [in] new2OldEnd - specifies the end of the permutation array that starts at
 *              \a new2OldBg, so that pointer to a tuple index (\a pi) varies as this:
 *              \a new2OldBg <= \a pi < \a new2OldEnd.
 *  \return DataArrayChar * - the new instance of DataArrayChar that the caller
 *          is to delete using decrRef() as it is no more needed.
 */
DataArrayChar *DataArrayChar::selectByTupleId(const int *new2OldBg, const int *new2OldEnd) const
{
  return selectByTupleIdSafe(new2OldBg,new2OldEnd);
}

/*!
 * Returns a shorten and permuted copy of \a this array. The new DataArrayChar is
 * of size \a new2OldEnd - \a new2OldBg and it's values are permuted as required by
 * \a new2OldBg array.
 * The values are permuted so that  \c new[ i ] = \c old[ \a new2OldBg[ i ]].
 * This method is equivalent to renumberAndReduce() except that convention in input is
 * \c new2old and \b not \c old2new.
 * This method is equivalent to selectByTupleId() except that it prevents coping data
 * from behind the end of \a this array.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] new2OldBg - pointer to the beginning of a permutation array that gives a
 *              tuple index in \a this array to fill the i-th tuple in the new array.
 *  \param [in] new2OldEnd - specifies the end of the permutation array that starts at
 *              \a new2OldBg, so that pointer to a tuple index (\a pi) varies as this:
 *              \a new2OldBg <= \a pi < \a new2OldEnd.
 *  \return DataArrayChar * - the new instance of DataArrayChar that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a new2OldEnd - \a new2OldBg > \a this->getNumberOfTuples().
 */
DataArrayChar *DataArrayChar::selectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayChar> ret=buildEmptySpecializedDAChar();
  int nbComp=getNumberOfComponents();
  int oldNbOfTuples=getNumberOfTuples();
  ret->alloc((int)std::distance(new2OldBg,new2OldEnd),nbComp);
  ret->copyStringInfoFrom(*this);
  char *pt=ret->getPointer();
  const char *srcPt=getConstPointer();
  int i=0;
  for(const int *w=new2OldBg;w!=new2OldEnd;w++,i++)
    if(*w>=0 && *w<oldNbOfTuples)
      std::copy(srcPt+(*w)*nbComp,srcPt+((*w)+1)*nbComp,pt+i*nbComp);
    else
      throw INTERP_KERNEL::Exception("DataArrayChar::selectByTupleIdSafe : some ids has been detected to be out of [0,this->getNumberOfTuples) !");
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * Returns a shorten copy of \a this array. The new DataArrayChar contains every
 * (\a bg + \c i * \a step)-th tuple of \a this array located before the \a end2-th
 * tuple. Indices of the selected tuples are the same as ones returned by the Python
 * command \c range( \a bg, \a end2, \a step ).
 * This method is equivalent to selectByTupleIdSafe() except that the input array is
 * not constructed explicitly.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] bg - index of the first tuple to copy from \a this array.
 *  \param [in] end2 - index of the tuple before which the tuples to copy are located.
 *  \param [in] step - index increment to get index of the next tuple to copy.
 *  \return DataArrayChar * - the new instance of DataArrayChar that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If (\a end2 < \a bg) or (\a step <= 0).
 *  \sa DataArrayChar::substr.
 */
DataArrayChar *DataArrayChar::selectByTupleId2(int bg, int end2, int step) const
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayChar> ret=buildEmptySpecializedDAChar();
  int nbComp=getNumberOfComponents();
  int newNbOfTuples=GetNumberOfItemGivenBESRelative(bg,end2,step,"DataArrayInt::selectByTupleId2 : ");
  ret->alloc(newNbOfTuples,nbComp);
  char *pt=ret->getPointer();
  const char *srcPt=getConstPointer()+bg*nbComp;
  for(int i=0;i<newNbOfTuples;i++,srcPt+=step*nbComp)
    std::copy(srcPt,srcPt+nbComp,pt+i*nbComp);
  ret->copyStringInfoFrom(*this);
  return ret.retn();
}

/*!
 * Checks if all values in \a this array are equal to \a val.
 *  \param [in] val - value to check equality of array values to.
 *  \return bool - \a true if all values are \a val.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1
 */
bool DataArrayChar::isUniform(char val) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::isUniform : must be applied on DataArrayChar with only one component, you can call 'rearrange' method before !");
  int nbOfTuples=getNumberOfTuples();
  const char *w=getConstPointer();
  const char *end2=w+nbOfTuples;
  for(;w!=end2;w++)
    if(*w!=val)
      return false;
  return true;
}

/*!
 * Changes the number of components within \a this array so that its raw data **does
 * not** change, instead splitting this data into tuples changes.
 *  \param [in] newNbOfComp - number of components for \a this array to have.
 *  \throw If \a this is not allocated
 *  \throw If getNbOfElems() % \a newNbOfCompo != 0.
 *  \throw If \a newNbOfCompo is lower than 1.
 *  \throw If the rearrange method would lead to a number of tuples higher than 2147483647 (maximal capacity of int32 !).
 *  \warning This method erases all (name and unit) component info set before!
 */
void DataArrayChar::rearrange(int newNbOfCompo)
{
  checkAllocated();
  if(newNbOfCompo<1)
    throw INTERP_KERNEL::Exception("DataArrayChar::rearrange : input newNbOfCompo must be > 0 !");
  std::size_t nbOfElems=getNbOfElems();
  if(nbOfElems%newNbOfCompo!=0)
    throw INTERP_KERNEL::Exception("DataArrayChar::rearrange : nbOfElems%newNbOfCompo!=0 !");
  if(nbOfElems/newNbOfCompo>(std::size_t)std::numeric_limits<int>::max())
    throw INTERP_KERNEL::Exception("DataArrayChar::rearrange : the rearrangement leads to too high number of tuples (> 2147483647) !");
  _info_on_compo.clear();
  _info_on_compo.resize(newNbOfCompo);
  declareAsNew();
}

/*!
 * Returns a shorten copy of \a this array. The new DataArrayChar contains all
 * tuples starting from the \a tupleIdBg-th tuple and including all tuples located before
 * the \a tupleIdEnd-th one. This methods has a similar behavior as std::string::substr().
 * This method is a specialization of selectByTupleId2().
 *  \param [in] tupleIdBg - index of the first tuple to copy from \a this array.
 *  \param [in] tupleIdEnd - index of the tuple before which the tuples to copy are located.
 *          If \a tupleIdEnd == -1, all the tuples till the end of \a this array are copied.
 *  \return DataArrayChar * - the new instance of DataArrayChar that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a tupleIdBg < 0.
 *  \throw If \a tupleIdBg > \a this->getNumberOfTuples().
    \throw If \a tupleIdEnd != -1 && \a tupleIdEnd < \a this->getNumberOfTuples().
 *  \sa DataArrayChar::selectByTupleId2
 */
DataArrayChar *DataArrayChar::substr(int tupleIdBg, int tupleIdEnd) const
{
  checkAllocated();
  int nbt=getNumberOfTuples();
  if(tupleIdBg<0)
    throw INTERP_KERNEL::Exception("DataArrayChar::substr : The tupleIdBg parameter must be greater than 0 !");
  if(tupleIdBg>nbt)
    throw INTERP_KERNEL::Exception("DataArrayChar::substr : The tupleIdBg parameter is greater than number of tuples !");
  int trueEnd=tupleIdEnd;
  if(tupleIdEnd!=-1)
    {
      if(tupleIdEnd>nbt)
        throw INTERP_KERNEL::Exception("DataArrayChar::substr : The tupleIdBg parameter is greater or equal than number of tuples !");
    }
  else
    trueEnd=nbt;
  int nbComp=getNumberOfComponents();
  MEDCouplingAutoRefCountObjectPtr<DataArrayChar> ret=buildEmptySpecializedDAChar();
  ret->alloc(trueEnd-tupleIdBg,nbComp);
  ret->copyStringInfoFrom(*this);
  std::copy(getConstPointer()+tupleIdBg*nbComp,getConstPointer()+trueEnd*nbComp,ret->getPointer());
  return ret.retn();
}

/*!
 * Returns a shorten or extended copy of \a this array. If \a newNbOfComp is less
 * than \a this->getNumberOfComponents() then the result array is shorten as each tuple
 * is truncated to have \a newNbOfComp components, keeping first components. If \a
 * newNbOfComp is more than \a this->getNumberOfComponents() then the result array is
 * expanded as each tuple is populated with \a dftValue to have \a newNbOfComp
 * components.  
 *  \param [in] newNbOfComp - number of components for the new array to have.
 *  \param [in] dftValue - value assigned to new values added to the new array.
 *  \return DataArrayChar * - the new instance of DataArrayChar that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 */
DataArrayChar *DataArrayChar::changeNbOfComponents(int newNbOfComp, char dftValue) const
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayChar> ret=buildEmptySpecializedDAChar();
  ret->alloc(getNumberOfTuples(),newNbOfComp);
  const char *oldc=getConstPointer();
  char *nc=ret->getPointer();
  int nbOfTuples=getNumberOfTuples();
  int oldNbOfComp=getNumberOfComponents();
  int dim=std::min(oldNbOfComp,newNbOfComp);
  for(int i=0;i<nbOfTuples;i++)
    {
      int j=0;
      for(;j<dim;j++)
        nc[newNbOfComp*i+j]=oldc[i*oldNbOfComp+j];
      for(;j<newNbOfComp;j++)
        nc[newNbOfComp*i+j]=dftValue;
    }
  ret->setName(getName());
  for(int i=0;i<dim;i++)
    ret->setInfoOnComponent(i,getInfoOnComponent(i));
  ret->setName(getName());
  return ret.retn();
}

/*!
 * Returns a copy of \a this array composed of selected components.
 * The new DataArrayChar has the same number of tuples but includes components
 * specified by \a compoIds parameter. So that getNbOfElems() of the result array
 * can be either less, same or more than \a this->getNbOfElems().
 *  \param [in] compoIds - sequence of zero based indices of components to include
 *              into the new array.
 *  \return DataArrayChar * - the new instance of DataArrayChar that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If a component index (\a i) is not valid: 
 *         \a i < 0 || \a i >= \a this->getNumberOfComponents().
 *
 *  \ref py_mcdataarrayint_keepselectedcomponents "Here is a Python example".
 */
DataArray *DataArrayChar::keepSelectedComponents(const std::vector<int>& compoIds) const
{
  checkAllocated();
  MEDCouplingAutoRefCountObjectPtr<DataArrayChar> ret(buildEmptySpecializedDAChar());
  int newNbOfCompo=(int)compoIds.size();
  int oldNbOfCompo=getNumberOfComponents();
  for(std::vector<int>::const_iterator it=compoIds.begin();it!=compoIds.end();it++)
    DataArray::CheckValueInRange(oldNbOfCompo,(*it),"keepSelectedComponents invalid requested component");
  int nbOfTuples=getNumberOfTuples();
  ret->alloc(nbOfTuples,newNbOfCompo);
  ret->copyPartOfStringInfoFrom(*this,compoIds);
  const char *oldc=getConstPointer();
  char *nc=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(int j=0;j<newNbOfCompo;j++,nc++)
      *nc=oldc[i*oldNbOfCompo+compoIds[j]];
  return ret.retn();
}

/*!
 * Appends components of another array to components of \a this one, tuple by tuple.
 * So that the number of tuples of \a this array remains the same and the number of 
 * components increases.
 *  \param [in] other - the DataArrayChar to append to \a this one.
 *  \throw If \a this is not allocated.
 *  \throw If \a this and \a other arrays have different number of tuples.
 *
 *  \ref cpp_mcdataarrayint_meldwith "Here is a C++ example".
 *
 *  \ref py_mcdataarrayint_meldwith "Here is a Python example".
 */
void DataArrayChar::meldWith(const DataArrayChar *other)
{
  if(!other)
    throw INTERP_KERNEL::Exception("DataArrayChar::meldWith : DataArrayChar pointer in input is NULL !");
  checkAllocated();
  other->checkAllocated();
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples!=other->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("DataArrayChar::meldWith : mismatch of number of tuples !");
  int nbOfComp1=getNumberOfComponents();
  int nbOfComp2=other->getNumberOfComponents();
  char *newArr=(char *)malloc(nbOfTuples*(nbOfComp1+nbOfComp2)*sizeof(char));
  char *w=newArr;
  const char *inp1=getConstPointer();
  const char *inp2=other->getConstPointer();
  for(int i=0;i<nbOfTuples;i++,inp1+=nbOfComp1,inp2+=nbOfComp2)
    {
      w=std::copy(inp1,inp1+nbOfComp1,w);
      w=std::copy(inp2,inp2+nbOfComp2,w);
    }
  useArray(newArr,true,C_DEALLOC,nbOfTuples,nbOfComp1+nbOfComp2);
  std::vector<int> compIds(nbOfComp2);
  for(int i=0;i<nbOfComp2;i++)
    compIds[i]=nbOfComp1+i;
  copyPartOfStringInfoFrom2(compIds,*other);
}

/*!
 * Copy all values from another DataArrayChar into specified tuples and components
 * of \a this array. Textual data is not copied.
 * The tree parameters defining set of indices of tuples and components are similar to
 * the tree parameters of the Python function \c range(\c start,\c stop,\c step).
 *  \param [in] a - the array to copy values from.
 *  \param [in] bgTuples - index of the first tuple of \a this array to assign values to.
 *  \param [in] endTuples - index of the tuple before which the tuples to assign to
 *              are located.
 *  \param [in] stepTuples - index increment to get index of the next tuple to assign to.
 *  \param [in] bgComp - index of the first component of \a this array to assign values to.
 *  \param [in] endComp - index of the component before which the components to assign
 *              to are located.
 *  \param [in] stepComp - index increment to get index of the next component to assign to.
 *  \param [in] strictCompoCompare - if \a true (by default), then \a a->getNumberOfComponents() 
 *              must be equal to the number of columns to assign to, else an
 *              exception is thrown; if \a false, then it is only required that \a
 *              a->getNbOfElems() equals to number of values to assign to (this condition
 *              must be respected even if \a strictCompoCompare is \a true). The number of 
 *              values to assign to is given by following Python expression:
 *              \a nbTargetValues = 
 *              \c len(\c range(\a bgTuples,\a endTuples,\a stepTuples)) *
 *              \c len(\c range(\a bgComp,\a endComp,\a stepComp)).
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a this is not allocated.
 *  \throw If parameters specifying tuples and components to assign to do not give a
 *            non-empty range of increasing indices.
 *  \throw If \a a->getNbOfElems() != \a nbTargetValues.
 *  \throw If \a strictCompoCompare == \a true && \a a->getNumberOfComponents() !=
 *            \c len(\c range(\a bgComp,\a endComp,\a stepComp)).
 *
 *  \ref py_mcdataarrayint_setpartofvalues1 "Here is a Python example".
 */
void DataArrayChar::setPartOfValues1(const DataArrayChar *a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayChar::setPartOfValues1 : DataArrayChar pointer in input is NULL !");
  const char msg[]="DataArrayChar::setPartOfValues1";
  checkAllocated();
  a->checkAllocated();
  int newNbOfTuples=DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg);
  int newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
  DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
  bool assignTech=true;
  if(a->getNbOfElems()==(std::size_t)newNbOfTuples*newNbOfComp)
    {
      if(strictCompoCompare)
        a->checkNbOfTuplesAndComp(newNbOfTuples,newNbOfComp,msg);
    }
  else
    {
      a->checkNbOfTuplesAndComp(1,newNbOfComp,msg);
      assignTech=false;
    }
  char *pt=getPointer()+bgTuples*nbComp+bgComp;
  const char *srcPt=a->getConstPointer();
  if(assignTech)
    {
      for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        for(int j=0;j<newNbOfComp;j++,srcPt++)
          pt[j*stepComp]=*srcPt;
    }
  else
    {
      for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        {
          const char *srcPt2=srcPt;
          for(int j=0;j<newNbOfComp;j++,srcPt2++)
            pt[j*stepComp]=*srcPt2;
        }
    }
}

/*!
 * Assign a given value to values at specified tuples and components of \a this array.
 * The tree parameters defining set of indices of tuples and components are similar to
 * the tree parameters of the Python function \c range(\c start,\c stop,\c step)..
 *  \param [in] a - the value to assign.
 *  \param [in] bgTuples - index of the first tuple of \a this array to assign to.
 *  \param [in] endTuples - index of the tuple before which the tuples to assign to
 *              are located.
 *  \param [in] stepTuples - index increment to get index of the next tuple to assign to.
 *  \param [in] bgComp - index of the first component of \a this array to assign to.
 *  \param [in] endComp - index of the component before which the components to assign
 *              to are located.
 *  \param [in] stepComp - index increment to get index of the next component to assign to.
 *  \throw If \a this is not allocated.
 *  \throw If parameters specifying tuples and components to assign to, do not give a
 *            non-empty range of increasing indices or indices are out of a valid range
 *            for \this array.
 *
 *  \ref py_mcdataarrayint_setpartofvaluessimple1 "Here is a Python example".
 */
void DataArrayChar::setPartOfValuesSimple1(char a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp)
{
  const char msg[]="DataArrayChar::setPartOfValuesSimple1";
  checkAllocated();
  int newNbOfTuples=DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg);
  int newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
  DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
  char *pt=getPointer()+bgTuples*nbComp+bgComp;
  for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
    for(int j=0;j<newNbOfComp;j++)
      pt[j*stepComp]=a;
}


/*!
 * Copy all values from another DataArrayChar (\a a) into specified tuples and 
 * components of \a this array. Textual data is not copied.
 * The tuples and components to assign to are defined by C arrays of indices.
 * There are two *modes of usage*:
 * - If \a a->getNbOfElems() equals to number of values to assign to, then every value
 *   of \a a is assigned to its own location within \a this array. 
 * - If \a a includes one tuple, then all values of \a a are assigned to the specified
 *   components of every specified tuple of \a this array. In this mode it is required
 *   that \a a->getNumberOfComponents() equals to the number of specified components.
 * 
 *  \param [in] a - the array to copy values from.
 *  \param [in] bgTuples - pointer to an array of tuple indices of \a this array to
 *              assign values of \a a to.
 *  \param [in] endTuples - specifies the end of the array \a bgTuples, so that
 *              pointer to a tuple index <em>(pi)</em> varies as this: 
 *              \a bgTuples <= \a pi < \a endTuples.
 *  \param [in] bgComp - pointer to an array of component indices of \a this array to
 *              assign values of \a a to.
 *  \param [in] endComp - specifies the end of the array \a bgTuples, so that
 *              pointer to a component index <em>(pi)</em> varies as this: 
 *              \a bgComp <= \a pi < \a endComp.
 *  \param [in] strictCompoCompare - this parameter is checked only if the
 *               *mode of usage* is the first; if it is \a true (default), 
 *               then \a a->getNumberOfComponents() must be equal 
 *               to the number of specified columns, else this is not required.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a this is not allocated.
 *  \throw If any index of tuple/component given by <em>bgTuples / bgComp</em> is
 *         out of a valid range for \a this array.
 *  \throw In the first *mode of usage*, if <em>strictCompoCompare == true </em> and
 *         if <em> a->getNumberOfComponents() != (endComp - bgComp) </em>.
 *  \throw In the second *mode of usage*, if \a a->getNumberOfTuples() != 1 or
 *         <em> a->getNumberOfComponents() != (endComp - bgComp)</em>.
 *
 *  \ref py_mcdataarrayint_setpartofvalues2 "Here is a Python example".
 */
void DataArrayChar::setPartOfValues2(const DataArrayChar *a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp, bool strictCompoCompare)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayChar::setPartOfValues2 : DataArrayChar pointer in input is NULL !");
  const char msg[]="DataArrayChar::setPartOfValues2";
  checkAllocated();
  a->checkAllocated();
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  for(const int *z=bgComp;z!=endComp;z++)
    DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
  int newNbOfTuples=(int)std::distance(bgTuples,endTuples);
  int newNbOfComp=(int)std::distance(bgComp,endComp);
  bool assignTech=true;
  if(a->getNbOfElems()==(std::size_t)newNbOfTuples*newNbOfComp)
    {
      if(strictCompoCompare)
        a->checkNbOfTuplesAndComp(newNbOfTuples,newNbOfComp,msg);
    }
  else
    {
      a->checkNbOfTuplesAndComp(1,newNbOfComp,msg);
      assignTech=false;
    }
  char *pt=getPointer();
  const char *srcPt=a->getConstPointer();
  if(assignTech)
    {    
      for(const int *w=bgTuples;w!=endTuples;w++)
        {
          DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
          for(const int *z=bgComp;z!=endComp;z++,srcPt++)
            {    
              pt[(std::size_t)(*w)*nbComp+(*z)]=*srcPt;
            }
        }
    }
  else
    {
      for(const int *w=bgTuples;w!=endTuples;w++)
        {
          const char *srcPt2=srcPt;
          DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
          for(const int *z=bgComp;z!=endComp;z++,srcPt2++)
            {    
              pt[(std::size_t)(*w)*nbComp+(*z)]=*srcPt2;
            }
        }
    }
}

/*!
 * Assign a given value to values at specified tuples and components of \a this array.
 * The tuples and components to assign to are defined by C arrays of indices.
 *  \param [in] a - the value to assign.
 *  \param [in] bgTuples - pointer to an array of tuple indices of \a this array to
 *              assign \a a to.
 *  \param [in] endTuples - specifies the end of the array \a bgTuples, so that
 *              pointer to a tuple index (\a pi) varies as this: 
 *              \a bgTuples <= \a pi < \a endTuples.
 *  \param [in] bgComp - pointer to an array of component indices of \a this array to
 *              assign \a a to.
 *  \param [in] endComp - specifies the end of the array \a bgTuples, so that
 *              pointer to a component index (\a pi) varies as this: 
 *              \a bgComp <= \a pi < \a endComp.
 *  \throw If \a this is not allocated.
 *  \throw If any index of tuple/component given by <em>bgTuples / bgComp</em> is
 *         out of a valid range for \a this array.
 *
 *  \ref py_mcdataarrayint_setpartofvaluessimple2 "Here is a Python example".
 */
void DataArrayChar::setPartOfValuesSimple2(char a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp)
{
  checkAllocated();
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  for(const int *z=bgComp;z!=endComp;z++)
    DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
  char *pt=getPointer();
  for(const int *w=bgTuples;w!=endTuples;w++)
    for(const int *z=bgComp;z!=endComp;z++)
      {
        DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
        pt[(std::size_t)(*w)*nbComp+(*z)]=a;
      }
}

/*!
 * Copy all values from another DataArrayChar (\a a) into specified tuples and 
 * components of \a this array. Textual data is not copied.
 * The tuples to assign to are defined by a C array of indices.
 * The components to assign to are defined by three values similar to parameters of
 * the Python function \c range(\c start,\c stop,\c step).
 * There are two *modes of usage*:
 * - If \a a->getNbOfElems() equals to number of values to assign to, then every value
 *   of \a a is assigned to its own location within \a this array. 
 * - If \a a includes one tuple, then all values of \a a are assigned to the specified
 *   components of every specified tuple of \a this array. In this mode it is required
 *   that \a a->getNumberOfComponents() equals to the number of specified components.
 *
 *  \param [in] a - the array to copy values from.
 *  \param [in] bgTuples - pointer to an array of tuple indices of \a this array to
 *              assign values of \a a to.
 *  \param [in] endTuples - specifies the end of the array \a bgTuples, so that
 *              pointer to a tuple index <em>(pi)</em> varies as this: 
 *              \a bgTuples <= \a pi < \a endTuples.
 *  \param [in] bgComp - index of the first component of \a this array to assign to.
 *  \param [in] endComp - index of the component before which the components to assign
 *              to are located.
 *  \param [in] stepComp - index increment to get index of the next component to assign to.
 *  \param [in] strictCompoCompare - this parameter is checked only in the first
 *               *mode of usage*; if \a strictCompoCompare is \a true (default), 
 *               then \a a->getNumberOfComponents() must be equal 
 *               to the number of specified columns, else this is not required.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a this is not allocated.
 *  \throw If any index of tuple given by \a bgTuples is out of a valid range for 
 *         \a this array.
 *  \throw In the first *mode of usage*, if <em>strictCompoCompare == true </em> and
 *         if <em> a->getNumberOfComponents()</em> is unequal to the number of components
 *         defined by <em>(bgComp,endComp,stepComp)</em>.
 *  \throw In the second *mode of usage*, if \a a->getNumberOfTuples() != 1 or
 *         <em> a->getNumberOfComponents()</em> is unequal to the number of components
 *         defined by <em>(bgComp,endComp,stepComp)</em>.
 *  \throw If parameters specifying components to assign to, do not give a
 *            non-empty range of increasing indices or indices are out of a valid range
 *            for \this array.
 *
 *  \ref py_mcdataarrayint_setpartofvalues3 "Here is a Python example".
 */
void DataArrayChar::setPartOfValues3(const DataArrayChar *a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayChar::setPartOfValues3 : DataArrayChar pointer in input is NULL !");
  const char msg[]="DataArrayChar::setPartOfValues3";
  checkAllocated();
  a->checkAllocated();
  int newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
  int newNbOfTuples=(int)std::distance(bgTuples,endTuples);
  bool assignTech=true;
  if(a->getNbOfElems()==(std::size_t)newNbOfTuples*newNbOfComp)
    {
      if(strictCompoCompare)
        a->checkNbOfTuplesAndComp(newNbOfTuples,newNbOfComp,msg);
    }
  else
    {
      a->checkNbOfTuplesAndComp(1,newNbOfComp,msg);
      assignTech=false;
    }
  char *pt=getPointer()+bgComp;
  const char *srcPt=a->getConstPointer();
  if(assignTech)
    {
      for(const int *w=bgTuples;w!=endTuples;w++)
        for(int j=0;j<newNbOfComp;j++,srcPt++)
          {
            DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
            pt[(std::size_t)(*w)*nbComp+j*stepComp]=*srcPt;
          }
    }
  else
    {
      for(const int *w=bgTuples;w!=endTuples;w++)
        {
          const char *srcPt2=srcPt;
          for(int j=0;j<newNbOfComp;j++,srcPt2++)
            {
              DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
              pt[(std::size_t)(*w)*nbComp+j*stepComp]=*srcPt2;
            }
        }
    }
}

/*!
 * Assign a given value to values at specified tuples and components of \a this array.
 * The tuples to assign to are defined by a C array of indices.
 * The components to assign to are defined by three values similar to parameters of
 * the Python function \c range(\c start,\c stop,\c step).
 *  \param [in] a - the value to assign.
 *  \param [in] bgTuples - pointer to an array of tuple indices of \a this array to
 *              assign \a a to.
 *  \param [in] endTuples - specifies the end of the array \a bgTuples, so that
 *              pointer to a tuple index <em>(pi)</em> varies as this: 
 *              \a bgTuples <= \a pi < \a endTuples.
 *  \param [in] bgComp - index of the first component of \a this array to assign to.
 *  \param [in] endComp - index of the component before which the components to assign
 *              to are located.
 *  \param [in] stepComp - index increment to get index of the next component to assign to.
 *  \throw If \a this is not allocated.
 *  \throw If any index of tuple given by \a bgTuples is out of a valid range for 
 *         \a this array.
 *  \throw If parameters specifying components to assign to, do not give a
 *            non-empty range of increasing indices or indices are out of a valid range
 *            for \this array.
 *
 *  \ref py_mcdataarrayint_setpartofvaluessimple3 "Here is a Python example".
 */
void DataArrayChar::setPartOfValuesSimple3(char a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp)
{
  const char msg[]="DataArrayChar::setPartOfValuesSimple3";
  checkAllocated();
  int newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
  int nbComp=getNumberOfComponents();
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
  char *pt=getPointer()+bgComp;
  for(const int *w=bgTuples;w!=endTuples;w++)
    for(int j=0;j<newNbOfComp;j++)
      {
        DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
        pt[(*w)*nbComp+j*stepComp]=a;
      }
}

void DataArrayChar::setPartOfValues4(const DataArrayChar *a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp, bool strictCompoCompare)
{
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayInt::setPartOfValues4 : input DataArrayInt is NULL !");
  const char msg[]="DataArrayInt::setPartOfValues4";
  checkAllocated();
  a->checkAllocated();
  int newNbOfTuples=DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg);
  int newNbOfComp=(int)std::distance(bgComp,endComp);
  int nbComp=getNumberOfComponents();
  for(const int *z=bgComp;z!=endComp;z++)
    DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
  bool assignTech=true;
  if(a->getNbOfElems()==(std::size_t)newNbOfTuples*newNbOfComp)
    {
      if(strictCompoCompare)
        a->checkNbOfTuplesAndComp(newNbOfTuples,newNbOfComp,msg);
    }
  else
    {
      a->checkNbOfTuplesAndComp(1,newNbOfComp,msg);
      assignTech=false;
    }
  const char *srcPt=a->getConstPointer();
  char *pt=getPointer()+bgTuples*nbComp;
  if(assignTech)
    {
      for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        for(const int *z=bgComp;z!=endComp;z++,srcPt++)
          pt[*z]=*srcPt;
    }
  else
    {
      for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        {
          const char *srcPt2=srcPt;
          for(const int *z=bgComp;z!=endComp;z++,srcPt2++)
            pt[*z]=*srcPt2;
        }
    }
}

void DataArrayChar::setPartOfValuesSimple4(char a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp)
{
  const char msg[]="DataArrayInt::setPartOfValuesSimple4";
  checkAllocated();
  int newNbOfTuples=DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg);
  int nbComp=getNumberOfComponents();
  for(const int *z=bgComp;z!=endComp;z++)
    DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
  int nbOfTuples=getNumberOfTuples();
  DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
  char *pt=getPointer()+bgTuples*nbComp;
  for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
    for(const int *z=bgComp;z!=endComp;z++)
      pt[*z]=a;
}

/*!
 * Copy some tuples from another DataArrayChar into specified tuples
 * of \a this array. Textual data is not copied. Both arrays must have equal number of
 * components.
 * Both the tuples to assign and the tuples to assign to are defined by a DataArrayChar.
 * All components of selected tuples are copied.
 *  \param [in] a - the array to copy values from.
 *  \param [in] tuplesSelec - the array specifying both source tuples of \a a and
 *              target tuples of \a this. \a tuplesSelec has two components, and the
 *              first component specifies index of the source tuple and the second
 *              one specifies index of the target tuple.
 *  \throw If \a this is not allocated.
 *  \throw If \a a is NULL.
 *  \throw If \a a is not allocated.
 *  \throw If \a tuplesSelec is NULL.
 *  \throw If \a tuplesSelec is not allocated.
 *  \throw If <em>this->getNumberOfComponents() != a->getNumberOfComponents()</em>.
 *  \throw If \a tuplesSelec->getNumberOfComponents() != 2.
 *  \throw If any tuple index given by \a tuplesSelec is out of a valid range for 
 *         the corresponding (\a this or \a a) array.
 */
void DataArrayChar::setPartOfValuesAdv(const DataArrayChar *a, const DataArrayChar *tuplesSelec)
{
  if(!a || !tuplesSelec)
    throw INTERP_KERNEL::Exception("DataArrayChar::setPartOfValuesAdv : DataArrayChar pointer in input is NULL !");
  checkAllocated();
  a->checkAllocated();
  tuplesSelec->checkAllocated();
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=a->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("DataArrayChar::setPartOfValuesAdv : This and a do not have the same number of components !");
  if(tuplesSelec->getNumberOfComponents()!=2)
    throw INTERP_KERNEL::Exception("DataArrayChar::setPartOfValuesAdv : Expecting to have a tuple selector DataArrayChar instance with exactly 2 components !");
  int thisNt=getNumberOfTuples();
  int aNt=a->getNumberOfTuples();
  char *valsToSet=getPointer();
  const char *valsSrc=a->getConstPointer();
  for(const char *tuple=tuplesSelec->begin();tuple!=tuplesSelec->end();tuple+=2)
    {
      if(tuple[1]>=0 && tuple[1]<aNt)
        {
          if(tuple[0]>=0 && tuple[0]<thisNt)
            std::copy(valsSrc+nbOfComp*tuple[1],valsSrc+nbOfComp*(tuple[1]+1),valsToSet+nbOfComp*tuple[0]);
          else
            {
              std::ostringstream oss; oss << "DataArrayChar::setPartOfValuesAdv : Tuple #" << std::distance(tuplesSelec->begin(),tuple)/2;
              oss << " of 'tuplesSelec' request of tuple id #" << tuple[0] << " in 'this' ! It should be in [0," << thisNt << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayChar::setPartOfValuesAdv : Tuple #" << std::distance(tuplesSelec->begin(),tuple)/2;
          oss << " of 'tuplesSelec' request of tuple id #" << tuple[1] << " in 'a' ! It should be in [0," << aNt << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

/*!
 * Copy some tuples from another DataArrayChar (\a aBase) into contiguous tuples
 * of \a this array. Textual data is not copied. Both arrays must have equal number of
 * components.
 * The tuples to assign to are defined by index of the first tuple, and
 * their number is defined by \a tuplesSelec->getNumberOfTuples().
 * The tuples to copy are defined by values of a DataArrayChar.
 * All components of selected tuples are copied.
 *  \param [in] tupleIdStart - index of the first tuple of \a this array to assign
 *              values to.
 *  \param [in] aBase - the array to copy values from.
 *  \param [in] tuplesSelec - the array specifying tuples of \a aBase to copy.
 *  \throw If \a this is not allocated.
 *  \throw If \a aBase is NULL.
 *  \throw If \a aBase is not allocated.
 *  \throw If \a tuplesSelec is NULL.
 *  \throw If \a tuplesSelec is not allocated.
 *  \throw If <em>this->getNumberOfComponents() != aBase->getNumberOfComponents()</em>.
 *  \throw If \a tuplesSelec->getNumberOfComponents() != 1.
 *  \throw If <em>tupleIdStart + tuplesSelec->getNumberOfTuples() > this->getNumberOfTuples().</em>
 *  \throw If any tuple index given by \a tuplesSelec is out of a valid range for 
 *         \a aBase array.
 */
void DataArrayChar::setContigPartOfSelectedValues(int tupleIdStart, const DataArray *aBase, const DataArrayInt *tuplesSelec)
{
  if(!aBase || !tuplesSelec)
    throw INTERP_KERNEL::Exception("DataArrayChar::setContigPartOfSelectedValues : input DataArray is NULL !");
  const DataArrayChar *a=dynamic_cast<const DataArrayChar *>(aBase);
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayChar::setContigPartOfSelectedValues : input DataArray aBase is not a DataArrayChar !");
  checkAllocated();
  a->checkAllocated();
  tuplesSelec->checkAllocated();
  int nbOfComp=getNumberOfComponents();
  if(nbOfComp!=a->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("DataArrayChar::setContigPartOfSelectedValues : This and a do not have the same number of components !");
  if(tuplesSelec->getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::setContigPartOfSelectedValues : Expecting to have a tuple selector DataArrayChar instance with exactly 1 component !");
  int thisNt=getNumberOfTuples();
  int aNt=a->getNumberOfTuples();
  int nbOfTupleToWrite=tuplesSelec->getNumberOfTuples();
  char *valsToSet=getPointer()+tupleIdStart*nbOfComp;
  if(tupleIdStart+nbOfTupleToWrite>thisNt)
    throw INTERP_KERNEL::Exception("DataArrayChar::setContigPartOfSelectedValues : invalid number range of values to write !");
  const char *valsSrc=a->getConstPointer();
  for(const int *tuple=tuplesSelec->begin();tuple!=tuplesSelec->end();tuple++,valsToSet+=nbOfComp)
    {
      if(*tuple>=0 && *tuple<aNt)
        {
          std::copy(valsSrc+nbOfComp*(*tuple),valsSrc+nbOfComp*(*tuple+1),valsToSet);
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayChar::setContigPartOfSelectedValues : Tuple #" << std::distance(tuplesSelec->begin(),tuple);
          oss << " of 'tuplesSelec' request of tuple id #" << *tuple << " in 'a' ! It should be in [0," << aNt << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
}

/*!
 * Copy some tuples from another DataArrayChar (\a aBase) into contiguous tuples
 * of \a this array. Textual data is not copied. Both arrays must have equal number of
 * components.
 * The tuples to copy are defined by three values similar to parameters of
 * the Python function \c range(\c start,\c stop,\c step).
 * The tuples to assign to are defined by index of the first tuple, and
 * their number is defined by number of tuples to copy.
 * All components of selected tuples are copied.
 *  \param [in] tupleIdStart - index of the first tuple of \a this array to assign
 *              values to.
 *  \param [in] aBase - the array to copy values from.
 *  \param [in] bg - index of the first tuple to copy of the array \a aBase.
 *  \param [in] end2 - index of the tuple of \a aBase before which the tuples to copy
 *              are located.
 *  \param [in] step - index increment to get index of the next tuple to copy.
 *  \throw If \a this is not allocated.
 *  \throw If \a aBase is NULL.
 *  \throw If \a aBase is not allocated.
 *  \throw If <em>this->getNumberOfComponents() != aBase->getNumberOfComponents()</em>.
 *  \throw If <em>tupleIdStart + len(range(bg,end2,step)) > this->getNumberOfTuples().</em>
 *  \throw If parameters specifying tuples to copy, do not give a
 *            non-empty range of increasing indices or indices are out of a valid range
 *            for the array \a aBase.
 */
void DataArrayChar::setContigPartOfSelectedValues2(int tupleIdStart, const DataArray *aBase, int bg, int end2, int step)
{
  if(!aBase)
    throw INTERP_KERNEL::Exception("DataArrayChar::setContigPartOfSelectedValues2 : input DataArray is NULL !");
  const DataArrayChar *a=dynamic_cast<const DataArrayChar *>(aBase);
  if(!a)
    throw INTERP_KERNEL::Exception("DataArrayChar::setContigPartOfSelectedValues2 : input DataArray aBase is not a DataArrayChar !");
  checkAllocated();
  a->checkAllocated();
  int nbOfComp=getNumberOfComponents();
  const char msg[]="DataArrayChar::setContigPartOfSelectedValues2";
  int nbOfTupleToWrite=DataArray::GetNumberOfItemGivenBES(bg,end2,step,msg);
  if(nbOfComp!=a->getNumberOfComponents())
    throw INTERP_KERNEL::Exception("DataArrayChar::setContigPartOfSelectedValues2 : This and a do not have the same number of components !");
  int thisNt=getNumberOfTuples();
  int aNt=a->getNumberOfTuples();
  char *valsToSet=getPointer()+tupleIdStart*nbOfComp;
  if(tupleIdStart+nbOfTupleToWrite>thisNt)
    throw INTERP_KERNEL::Exception("DataArrayChar::setContigPartOfSelectedValues2 : invalid number range of values to write !");
  if(end2>aNt)
    throw INTERP_KERNEL::Exception("DataArrayChar::setContigPartOfSelectedValues2 : invalid range of values to read !");
  const char *valsSrc=a->getConstPointer()+bg*nbOfComp;
  for(int i=0;i<nbOfTupleToWrite;i++,valsToSet+=nbOfComp,valsSrc+=step*nbOfComp)
    {
      std::copy(valsSrc,valsSrc+nbOfComp,valsToSet);
    }
}

/*!
 * Returns a shorten copy of \a this array. The new DataArrayChar contains ranges
 * of tuples specified by \a ranges parameter.
 * For more info on renumbering see \ref MEDCouplingArrayRenumbering.
 *  \param [in] ranges - std::vector of std::pair's each of which defines a range
 *              of tuples in [\c begin,\c end) format.
 *  \return DataArrayChar * - the new instance of DataArrayChar that the caller
 *          is to delete using decrRef() as it is no more needed.
 *  \throw If \a end < \a begin.
 *  \throw If \a end > \a this->getNumberOfTuples().
 *  \throw If \a this is not allocated.
 */
DataArray *DataArrayChar::selectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const
{
  checkAllocated();
  int nbOfComp=getNumberOfComponents();
  int nbOfTuplesThis=getNumberOfTuples();
  if(ranges.empty())
    {
      MEDCouplingAutoRefCountObjectPtr<DataArrayChar> ret=buildEmptySpecializedDAChar();
      ret->alloc(0,nbOfComp);
      ret->copyStringInfoFrom(*this);
      return ret.retn();
    }
  int ref=ranges.front().first;
  int nbOfTuples=0;
  bool isIncreasing=true;
  for(std::vector<std::pair<int,int> >::const_iterator it=ranges.begin();it!=ranges.end();it++)
    {
      if((*it).first<=(*it).second)
        {
          if((*it).first>=0 && (*it).second<=nbOfTuplesThis)
            {
              nbOfTuples+=(*it).second-(*it).first;
              if(isIncreasing)
                isIncreasing=ref<=(*it).first;
              ref=(*it).second;
            }
          else
            {
              std::ostringstream oss; oss << "DataArrayChar::selectByTupleRanges : on range #" << std::distance(ranges.begin(),it);
              oss << " (" << (*it).first << "," << (*it).second << ") is greater than number of tuples of this :" << nbOfTuples << " !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayChar::selectByTupleRanges : on range #" << std::distance(ranges.begin(),it);
          oss << " (" << (*it).first << "," << (*it).second << ") end is before begin !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  if(isIncreasing && nbOfTuplesThis==nbOfTuples)
    return deepCpy();
  MEDCouplingAutoRefCountObjectPtr<DataArrayChar> ret=buildEmptySpecializedDAChar();
  ret->alloc(nbOfTuples,nbOfComp);
  ret->copyStringInfoFrom(*this);
  const char *src=getConstPointer();
  char *work=ret->getPointer();
  for(std::vector<std::pair<int,int> >::const_iterator it=ranges.begin();it!=ranges.end();it++)
    work=std::copy(src+(*it).first*nbOfComp,src+(*it).second*nbOfComp,work);
  return ret.retn();
}

/*!
 * Returns a value located at specified tuple and component.
 * This method is equivalent to DataArrayChar::getIJ() except that validity of
 * parameters is checked. So this method is safe but expensive if used to go through
 * all values of \a this.
 *  \param [in] tupleId - index of tuple of interest.
 *  \param [in] compoId - index of component of interest.
 *  \return char - value located by \a tupleId and \a compoId.
 *  \throw If \a this is not allocated.
 *  \throw If condition <em>( 0 <= tupleId < this->getNumberOfTuples() )</em> is violated.
 *  \throw If condition <em>( 0 <= compoId < this->getNumberOfComponents() )</em> is violated.
 */
char DataArrayChar::getIJSafe(int tupleId, int compoId) const
{
  checkAllocated();
  if(tupleId<0 || tupleId>=getNumberOfTuples())
    {
      std::ostringstream oss; oss << "DataArrayChar::getIJSafe : request for tupleId " << tupleId << " should be in [0," << getNumberOfTuples() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  if(compoId<0 || compoId>=getNumberOfComponents())
    {
      std::ostringstream oss; oss << "DataArrayChar::getIJSafe : request for compoId " << compoId << " should be in [0," << getNumberOfComponents() << ") !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  return _mem[tupleId*_info_on_compo.size()+compoId];
}

/*!
 * Returns the first value of \a this. 
 *  \return char - the last value of \a this array.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this->getNumberOfTuples() < 1.
 */
char DataArrayChar::front() const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::front : number of components not equal to one !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<1)
    throw INTERP_KERNEL::Exception("DataArrayChar::front : number of tuples must be >= 1 !");
  return *(getConstPointer());
}

/*!
 * Returns the last value of \a this. 
 *  \return char - the last value of \a this array.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \throw If \a this->getNumberOfTuples() < 1.
 */
char DataArrayChar::back() const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::back : number of components not equal to one !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<1)
    throw INTERP_KERNEL::Exception("DataArrayChar::back : number of tuples must be >= 1 !");
  return *(getConstPointer()+nbOfTuples-1);
}

/*!
 * Creates a new DataArrayChar containing IDs (indices) of tuples holding value equal to a
 * given one.
 *  \param [in] val - the value to find within \a this.
 *  \return DataArrayChar * - a new instance of DataArrayChar. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 */
DataArrayInt *DataArrayChar::getIdsEqual(char val) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::getIdsEqual : the array must have only one component, you can call 'rearrange' method before !");
  const char *cptr=getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr==val)
      ret->pushBackSilent(i);
  return ret.retn();
}

/*!
 * Creates a new DataArrayChar containing IDs (indices) of tuples holding value \b not
 * equal to a given one. 
 *  \param [in] val - the value to ignore within \a this.
 *  \return DataArrayChar * - a new instance of DataArrayChar. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 */
DataArrayInt *DataArrayChar::getIdsNotEqual(char val) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::getIdsNotEqual : the array must have only one component, you can call 'rearrange' method before !");
  const char *cptr=getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr!=val)
      ret->pushBackSilent(i);
  return ret.retn();
}

/*!
 * This method searches the sequence specified in input parameter \b vals in \b this.
 * This works only for DataArrayChar having number of components equal to one (if not an INTERP_KERNEL::Exception will be thrown).
 * This method differs from DataArrayChar::locateTuple in that the position is internal raw data is not considered here contrary to DataArrayChar::locateTuple.
 * \sa DataArrayChar::locateTuple
 */
int DataArrayChar::search(const std::vector<char>& vals) const
{
  checkAllocated();
  int nbOfCompo=getNumberOfComponents();
  if(nbOfCompo!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::search : works only for DataArrayChar instance with one component !");
  const char *cptr=getConstPointer();
  std::size_t nbOfVals=getNbOfElems();
  const char *loc=std::search(cptr,cptr+nbOfVals,vals.begin(),vals.end());
  if(loc!=cptr+nbOfVals)
    return std::distance(cptr,loc);
  return -1;
}

/*!
 * This method is an extension of DataArrayChar::locateValue method because this method works for DataArrayChar with
 * any number of components excepted 0 (an INTERP_KERNEL::Exception is thrown in this case).
 * This method searches in \b this is there is a tuple that matched the input parameter \b tupl.
 * If any the tuple id is returned. If not -1 is returned.
 * 
 * This method throws an INTERP_KERNEL::Exception if the number of components in \b this mismatches with the size of
 * the input vector. An INTERP_KERNEL::Exception is thrown too if \b this is not allocated.
 *
 * \return tuple id where \b tupl is. -1 if no such tuple exists in \b this.
 * \sa DataArrayChar::search.
 */
int DataArrayChar::locateTuple(const std::vector<char>& tupl) const
{
  checkAllocated();
  int nbOfCompo=getNumberOfComponents();
  if(nbOfCompo==0)
    throw INTERP_KERNEL::Exception("DataArrayChar::locateTuple : 0 components in 'this' !");
  if(nbOfCompo!=(int)tupl.size())
    {
      std::ostringstream oss; oss << "DataArrayChar::locateTuple : 'this' contains " << nbOfCompo << " components and searching for a tuple of length " << tupl.size() << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  const char *cptr=getConstPointer();
  std::size_t nbOfVals=getNbOfElems();
  for(const char *work=cptr;work!=cptr+nbOfVals;)
    {
      work=std::search(work,cptr+nbOfVals,tupl.begin(),tupl.end());
      if(work!=cptr+nbOfVals)
        {
          if(std::distance(cptr,work)%nbOfCompo!=0)
            work++;
          else
            return std::distance(cptr,work)/nbOfCompo;
        }
    }
  return -1;
}

/*!
 * This method is an extension of DataArrayChar::presenceOfValue method because this method works for DataArrayChar with
 * any number of components excepted 0 (an INTERP_KERNEL::Exception is thrown in this case).
 * This method searches in \b this is there is a tuple that matched the input parameter \b tupl.
 * This method throws an INTERP_KERNEL::Exception if the number of components in \b this mismatches with the size of
 * the input vector. An INTERP_KERNEL::Exception is thrown too if \b this is not allocated.
 * \sa DataArrayChar::locateTuple
 */
bool DataArrayChar::presenceOfTuple(const std::vector<char>& tupl) const
{
  return locateTuple(tupl)!=-1;
}

/*!
 * Returns \a true if a given value is present within \a this one-dimensional array.
 *  \param [in] value - the value to find within \a this array.
 *  \return bool - \a true in case if \a value is present within \a this array.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \sa locateValue()
 */
bool DataArrayChar::presenceOfValue(char value) const
{
  return locateValue(value)!=-1;
}

/*!
 * This method expects to be called when number of components of this is equal to one.
 * This method returns true if it exists a tuple so that the value is contained in \b vals.
 * If not any tuple contains one of the values contained in 'vals' false is returned.
 * \sa DataArrayChar::locateValue
 */
bool DataArrayChar::presenceOfValue(const std::vector<char>& vals) const
{
  return locateValue(vals)!=-1;
}

/*!
 * This method expects to be called when number of components of this is equal to one.
 * This method returns the tuple id, if it exists, of the first tuple equal to \b value.
 * If not any tuple contains \b value -1 is returned.
 * \sa DataArrayChar::presenceOfValue
 */
int DataArrayChar::locateValue(char value) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::presenceOfValue : the array must have only one component, you can call 'rearrange' method before !");
  const char *cptr=getConstPointer();
  int nbOfTuples=getNumberOfTuples();
  const char *ret=std::find(cptr,cptr+nbOfTuples,value);
  if(ret!=cptr+nbOfTuples)
    return std::distance(cptr,ret);
  return -1;
}

/*!
 * This method expects to be called when number of components of this is equal to one.
 * This method returns the tuple id, if it exists, of the first tuple so that the value is contained in \b vals.
 * If not any tuple contains one of the values contained in 'vals' false is returned.
 * \sa DataArrayChar::presenceOfValue
 */
int DataArrayChar::locateValue(const std::vector<char>& vals) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayInt::presenceOfValue : the array must have only one component, you can call 'rearrange' method before !");
  std::set<char> vals2(vals.begin(),vals.end());
  const char *cptr=getConstPointer();
  int nbOfTuples=getNumberOfTuples();
  for(const char *w=cptr;w!=cptr+nbOfTuples;w++)
    if(vals2.find(*w)!=vals2.end())
      return std::distance(cptr,w);
  return -1;
}

/*!
 * Returns the maximal value and its location within \a this one-dimensional array.
 *  \param [out] tupleId - index of the tuple holding the maximal value.
 *  \return char - the maximal value among all values of \a this array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this->getNumberOfTuples() < 1
 */
char DataArrayChar::getMaxValue(int& tupleId) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::getMaxValue : must be applied on DataArrayChar with only one component !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayChar::getMaxValue : array exists but number of tuples must be > 0 !");
  const char *vals=getConstPointer();
  const char *loc=std::max_element(vals,vals+nbOfTuples);
  tupleId=(int)std::distance(vals,loc);
  return *loc;
}

/*!
 * Returns the maximal value within \a this array that is allowed to have more than
 *  one component.
 *  \return char - the maximal value among all values of \a this array.
 *  \throw If \a this is not allocated.
 */
char DataArrayChar::getMaxValueInArray() const
{
  checkAllocated();
  const char *loc=std::max_element(begin(),end());
  return *loc;
}

/*!
 * Returns the minimal value and its location within \a this one-dimensional array.
 *  \param [out] tupleId - index of the tuple holding the minimal value.
 *  \return char - the minimal value among all values of \a this array.
 *  \throw If \a this->getNumberOfComponents() != 1
 *  \throw If \a this->getNumberOfTuples() < 1
 */
char DataArrayChar::getMinValue(int& tupleId) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::getMaxValue : must be applied on DataArrayChar with only one component !");
  int nbOfTuples=getNumberOfTuples();
  if(nbOfTuples<=0)
    throw INTERP_KERNEL::Exception("DataArrayChar::getMaxValue : array exists but number of tuples must be > 0 !");
  const char *vals=getConstPointer();
  const char *loc=std::min_element(vals,vals+nbOfTuples);
  tupleId=(int)std::distance(vals,loc);
  return *loc;
}

/*!
 * Returns the minimal value within \a this array that is allowed to have more than
 *  one component.
 *  \return char - the minimal value among all values of \a this array.
 *  \throw If \a this is not allocated.
 */
char DataArrayChar::getMinValueInArray() const
{
  checkAllocated();
  const char *loc=std::min_element(begin(),end());
  return *loc;
}

/*!
 * This method works only on data array with one component.
 * This method returns a newly allocated array storing stored ascendantly tuple ids in \b this so that
 * this[*id] in [\b vmin,\b vmax)
 * 
 * \param [in] vmin begin of range. This value is included in range.
 * \param [in] vmax end of range. This value is \b not included in range.
 * \return a newly allocated data array that the caller should deal with.
 */
DataArrayInt *DataArrayChar::getIdsInRange(char vmin, char vmax) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::getIdsInRange : this must have exactly one component !");
  const char *cptr=getConstPointer();
  MEDCouplingAutoRefCountObjectPtr<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(0,1);
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr>=vmin && *cptr<vmax)
      ret->pushBackSilent(i);
  return ret.retn();
}

/*!
 * Returns a new DataArrayChar by concatenating two given arrays, so that (1) the number
 * of tuples in the result array is <em> a1->getNumberOfTuples() + a2->getNumberOfTuples() -
 * offsetA2</em> and (2)
 * the number of component in the result array is same as that of each of given arrays.
 * First \a offsetA2 tuples of \a a2 are skipped and thus are missing from the result array.
 * Info on components is copied from the first of the given arrays. Number of components
 * in the given arrays must be the same.
 *  \param [in] a1 - an array to include in the result array.
 *  \param [in] a2 - another array to include in the result array.
 *  \param [in] offsetA2 - number of tuples of \a a2 to skip.
 *  \return DataArrayChar * - the new instance of DataArrayChar.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If either \a a1 or \a a2 is NULL.
 *  \throw If \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents().
 */
DataArrayChar *DataArrayChar::Aggregate(const DataArrayChar *a1, const DataArrayChar *a2)
{
  if(!a1 || !a2)
    throw INTERP_KERNEL::Exception("DataArrayChar::Aggregate : input DataArrayChar instance is NULL !");
  std::vector<const DataArrayChar *> v(2); v[0]=a1; v[1]=a2;
  return Aggregate(v);
}

/*!
 * Returns a new DataArrayChar by concatenating all given arrays, so that (1) the number
 * of tuples in the result array is a sum of the number of tuples of given arrays and (2)
 * the number of component in the result array is same as that of each of given arrays.
 * Info on components is copied from the first of the given arrays. Number of components
 * in the given arrays must be  the same.
 *  \param [in] arr - a sequence of arrays to include in the result array.
 *  \return DataArrayChar * - the new instance of DataArrayChar.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If all arrays within \a arr are NULL.
 *  \throw If getNumberOfComponents() of arrays within \a arr.
 */
DataArrayChar *DataArrayChar::Aggregate(const std::vector<const DataArrayChar *>& arr)
{
  std::vector<const DataArrayChar *> a;
  for(std::vector<const DataArrayChar *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
    if(*it4)
      a.push_back(*it4);
  if(a.empty())
    throw INTERP_KERNEL::Exception("DataArrayChar::Aggregate : input list must be NON EMPTY !");
  std::vector<const DataArrayChar *>::const_iterator it=a.begin();
  int nbOfComp=(*it)->getNumberOfComponents();
  int nbt=(*it++)->getNumberOfTuples();
  for(int i=1;it!=a.end();it++,i++)
    {
      if((*it)->getNumberOfComponents()!=nbOfComp)
        throw INTERP_KERNEL::Exception("DataArrayChar::Aggregate : Nb of components mismatch for array aggregation !");
      nbt+=(*it)->getNumberOfTuples();
    }
  MEDCouplingAutoRefCountObjectPtr<DataArrayChar> ret=a[0]->buildEmptySpecializedDAChar();
  ret->alloc(nbt,nbOfComp);
  char *pt=ret->getPointer();
  for(it=a.begin();it!=a.end();it++)
    pt=std::copy((*it)->getConstPointer(),(*it)->getConstPointer()+(*it)->getNbOfElems(),pt);
  ret->copyStringInfoFrom(*(a[0]));
  return ret.retn();
}

/*!
 * Returns a new DataArrayChar by aggregating two given arrays, so that (1) the number
 * of components in the result array is a sum of the number of components of given arrays
 * and (2) the number of tuples in the result array is same as that of each of given
 * arrays. In other words the i-th tuple of result array includes all components of
 * i-th tuples of all given arrays.
 * Number of tuples in the given arrays must be the same.
 *  \param [in] a1 - an array to include in the result array.
 *  \param [in] a2 - another array to include in the result array.
 *  \return DataArrayChar * - the new instance of DataArrayChar.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If both \a a1 and \a a2 are NULL.
 *  \throw If any given array is not allocated.
 *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
 */
DataArrayChar *DataArrayChar::Meld(const DataArrayChar *a1, const DataArrayChar *a2)
{
  std::vector<const DataArrayChar *> arr(2);
  arr[0]=a1; arr[1]=a2;
  return Meld(arr);
}

/*!
 * Returns a new DataArrayChar by aggregating all given arrays, so that (1) the number
 * of components in the result array is a sum of the number of components of given arrays
 * and (2) the number of tuples in the result array is same as that of each of given
 * arrays. In other words the i-th tuple of result array includes all components of
 * i-th tuples of all given arrays.
 * Number of tuples in the given arrays must be  the same.
 *  \param [in] arr - a sequence of arrays to include in the result array.
 *  \return DataArrayChar * - the new instance of DataArrayChar.
 *          The caller is to delete this result array using decrRef() as it is no more
 *          needed.
 *  \throw If all arrays within \a arr are NULL.
 *  \throw If any given array is not allocated.
 *  \throw If getNumberOfTuples() of arrays within \a arr is different.
 */
DataArrayChar *DataArrayChar::Meld(const std::vector<const DataArrayChar *>& arr)
{
  std::vector<const DataArrayChar *> a;
  for(std::vector<const DataArrayChar *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
    if(*it4)
      a.push_back(*it4);
  if(a.empty())
    throw INTERP_KERNEL::Exception("DataArrayChar::Meld : array must be NON empty !");
  std::vector<const DataArrayChar *>::const_iterator it;
  for(it=a.begin();it!=a.end();it++)
    (*it)->checkAllocated();
  it=a.begin();
  int nbOfTuples=(*it)->getNumberOfTuples();
  std::vector<int> nbc(a.size());
  std::vector<const char *> pts(a.size());
  nbc[0]=(*it)->getNumberOfComponents();
  pts[0]=(*it++)->getConstPointer();
  for(int i=1;it!=a.end();it++,i++)
    {
      if(nbOfTuples!=(*it)->getNumberOfTuples())
        throw INTERP_KERNEL::Exception("DataArrayChar::meld : mismatch of number of tuples !");
      nbc[i]=(*it)->getNumberOfComponents();
      pts[i]=(*it)->getConstPointer();
    }
  int totalNbOfComp=std::accumulate(nbc.begin(),nbc.end(),0);
  DataArrayChar *ret=a[0]->buildEmptySpecializedDAChar();
  ret->alloc(nbOfTuples,totalNbOfComp);
  char *retPtr=ret->getPointer();
  for(int i=0;i<nbOfTuples;i++)
    for(int j=0;j<(int)a.size();j++)
      {
        retPtr=std::copy(pts[j],pts[j]+nbc[j],retPtr);
        pts[j]+=nbc[j];
      }
  int k=0;
  for(int i=0;i<(int)a.size();i++)
    for(int j=0;j<nbc[i];j++,k++)
      ret->setInfoOnComponent(k,a[i]->getInfoOnComponent(j));
  return ret;
}

/*!
 * Sets a C array to be used as raw data of \a this. The previously set info
 *  of components is retained and re-sized. 
 * For more info see \ref MEDCouplingArraySteps1.
 *  \param [in] array - the C array to be used as raw data of \a this.
 *  \param [in] ownership - if \a true, \a array will be deallocated at destruction of \a this.
 *  \param [in] type - specifies how to deallocate \a array. If \a type == ParaMEDMEM::CPP_DEALLOC,
 *                     \c delete [] \c array; will be called. If \a type == ParaMEDMEM::C_DEALLOC,
 *                     \c free(\c array ) will be called.
 *  \param [in] nbOfTuple - new number of tuples in \a this.
 *  \param [in] nbOfCompo - new number of components in \a this.
 */
void DataArrayChar::useArray(const char *array, bool ownership,  DeallocType type, int nbOfTuple, int nbOfCompo)
{
  _info_on_compo.resize(nbOfCompo);
  _mem.useArray(array,ownership,type,(std::size_t)nbOfTuple*nbOfCompo);
  declareAsNew();
}

void DataArrayChar::useExternalArrayWithRWAccess(const char *array, int nbOfTuple, int nbOfCompo)
{
  _info_on_compo.resize(nbOfCompo);
  _mem.useExternalArrayWithRWAccess(array,(std::size_t)nbOfTuple*nbOfCompo);
  declareAsNew();
}

/*!
 * Returns a new instance of DataArrayByte. The caller is to delete this array
 * using decrRef() as it is no more needed. 
 */
DataArrayByte *DataArrayByte::New()
{
  return new DataArrayByte;
}

DataArrayByteIterator *DataArrayByte::iterator()
{
  return new DataArrayByteIterator(this);
}

/*!
 * Returns a full copy of \a this. For more info on copying data arrays see
 * \ref MEDCouplingArrayBasicsCopyDeep.
 *  \return DataArrayByte * - a new instance of DataArrayByte.
 */
DataArrayByte *DataArrayByte::deepCpy() const
{
  return new DataArrayByte(*this);
}

/*!
 * Returns either a \a deep or \a shallow copy of this array. For more info see
 * \ref MEDCouplingArrayBasicsCopyDeep and \ref MEDCouplingArrayBasicsCopyShallow.
 *  \param [in] dCpy - if \a true, a deep copy is returned, else, a shallow one.
 *  \return DataArrayByte * - either a new instance of DataArrayByte (if \a dCpy
 *          == \a true) or \a this instance (if \a dCpy == \a false).
 */
DataArrayByte *DataArrayByte::performCpy(bool dCpy) const
{
  if(dCpy)
    return deepCpy();
  else
    {
      incrRef();
      return const_cast<DataArrayByte *>(this);
    }
}

/*!
 * Returns the only one value in \a this, if and only if number of elements
 * (nb of tuples * nb of components) is equal to 1, and that \a this is allocated.
 *  \return char - the sole value stored in \a this array.
 *  \throw If at least one of conditions stated above is not fulfilled.
 */
char DataArrayByte::byteValue() const
{
  if(isAllocated())
    {
      if(getNbOfElems()==1)
        {
          return *getConstPointer();
        }
      else
        throw INTERP_KERNEL::Exception("DataArrayByte::byteValue : DataArrayByte instance is allocated but number of elements is not equal to 1 !");
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayByte::byteValue : DataArrayByte instance is not allocated !");
}

DataArrayChar *DataArrayByte::buildEmptySpecializedDAChar() const
{
  return DataArrayByte::New();
}

void DataArrayByte::reprStream(std::ostream& stream) const
{
  stream << "Name of byte array : \"" << _name << "\"\n";
  reprWithoutNameStream(stream);
}

void DataArrayByte::reprZipStream(std::ostream& stream) const
{
  stream << "Name of byte array : \"" << _name << "\"\n";
  reprZipWithoutNameStream(stream);
}

void DataArrayByte::reprWithoutNameStream(std::ostream& stream) const
{
  DataArray::reprWithoutNameStream(stream);
  if(_mem.reprHeader(getNumberOfComponents(),stream))
    {
      const char *data=begin();
      int nbOfTuples=getNumberOfTuples();
      int nbCompo=getNumberOfComponents();
      for(int i=0;i<nbOfTuples;i++,data+=nbCompo)
        {
          stream << "Tuple #" << i << " : ";
          std::copy(data,data+nbCompo,std::ostream_iterator<int>(stream," "));//it is not a bug int here not char because it is not ASCII here contrary to DataArrayAsciiChar
          stream << "\n";
        }
    }
}

void DataArrayByte::reprZipWithoutNameStream(std::ostream& stream) const
{
  DataArray::reprWithoutNameStream(stream);
  _mem.reprZip(getNumberOfComponents(),stream);
}

void DataArrayByte::reprCppStream(const std::string& varName, std::ostream& stream) const
{
  int nbTuples=getNumberOfTuples(),nbComp=getNumberOfComponents();
  const char *data=getConstPointer();
  stream << "DataArrayByte *" << varName << "=DataArrayByte::New();" << std::endl;
  if(nbTuples*nbComp>=1)
    {
      stream << "const char " << varName << "Data[" << nbTuples*nbComp << "]={";
      std::copy(data,data+nbTuples*nbComp-1,std::ostream_iterator<char>(stream,","));
      stream << data[nbTuples*nbComp-1] << "};" << std::endl;
      stream << varName << "->useArray(" << varName << "Data,false,CPP_DEALLOC," << nbTuples << "," << nbComp << ");" << std::endl;
    }
  else
    stream << varName << "->alloc(" << nbTuples << "," << nbComp << ");" << std::endl;
  stream << varName << "->setName(\"" << getName() << "\");" << std::endl;
}

/*!
 * Method that gives a quick overvien of \a this for python.
 */
void DataArrayByte::reprQuickOverview(std::ostream& stream) const
{
  static const std::size_t MAX_NB_OF_BYTE_IN_REPR=300;
  stream << "DataArrayByte C++ instance at " << this << ". ";
  if(isAllocated())
    {
      int nbOfCompo=(int)_info_on_compo.size();
      if(nbOfCompo>=1)
        {
          int nbOfTuples=getNumberOfTuples();
          stream << "Number of tuples : " << nbOfTuples << ". Number of components : " << nbOfCompo << "." << std::endl;
          reprQuickOverviewData(stream,MAX_NB_OF_BYTE_IN_REPR);
        }
      else
        stream << "Number of components : 0.";
    }
  else
    stream << "*** No data allocated ****";
}

void DataArrayByte::reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const
{
  const char *data=begin();
  int nbOfTuples=getNumberOfTuples();
  int nbOfCompo=(int)_info_on_compo.size();
  std::ostringstream oss2; oss2 << "[";
  std::string oss2Str(oss2.str());
  bool isFinished=true;
  for(int i=0;i<nbOfTuples && isFinished;i++)
    {
      if(nbOfCompo>1)
        {
          oss2 << "(";
          for(int j=0;j<nbOfCompo;j++,data++)
            {
              oss2 << (int)*data;
              if(j!=nbOfCompo-1) oss2 << ", ";
            }
          oss2 << ")";
        }
      else
        { oss2 << (int)*data; data++; }
      if(i!=nbOfTuples-1) oss2 << ", ";
      std::string oss3Str(oss2.str());
      if(oss3Str.length()<maxNbOfByteInRepr)
        oss2Str=oss3Str;
      else
        isFinished=false;
    }
  stream << oss2Str;
  if(!isFinished)
    stream << "... ";
  stream << "]";
}

bool DataArrayByte::isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const
{
  const DataArrayByte *otherC=dynamic_cast<const DataArrayByte *>(&other);
  if(!otherC)
    { reason="this is of type DataArrayByte whereas other is not a DataArrayByte instance"; return false; }
  return DataArrayChar::isEqualIfNotWhy(other,reason);
}

DataArrayByteIterator::DataArrayByteIterator(DataArrayByte *da):_da(da),_pt(0),_tuple_id(0),_nb_comp(0),_nb_tuple(0)
{
  if(_da)
    {
      _da->incrRef();
      if(_da->isAllocated())
        {
          _nb_comp=da->getNumberOfComponents();
          _nb_tuple=da->getNumberOfTuples();
          _pt=da->getPointer();
        }
    }
}

DataArrayByteIterator::~DataArrayByteIterator()
{
  if(_da)
    _da->decrRef();
}

DataArrayByteTuple *DataArrayByteIterator::nextt()
{
  if(_tuple_id<_nb_tuple)
    {
      _tuple_id++;
      DataArrayByteTuple *ret=new DataArrayByteTuple(_pt,_nb_comp);
      _pt+=_nb_comp;
      return ret;
    }
  else
    return 0;
}

DataArrayByteTuple::DataArrayByteTuple(char *pt, int nbOfComp):_pt(pt),_nb_of_compo(nbOfComp)
{
}

std::string DataArrayByteTuple::repr() const
{
  std::ostringstream oss; oss << "(";
  for(int i=0;i<_nb_of_compo-1;i++)
    oss << (int)_pt[i] << ", ";
  oss << _pt[_nb_of_compo-1] << ")";
  return oss.str();
}

char DataArrayByteTuple::byteValue() const
{
  if(_nb_of_compo==1)
    return *_pt;
  throw INTERP_KERNEL::Exception("DataArrayByteTuple::byteValue : DataArrayByteTuple instance has not exactly 1 component -> Not possible to convert it into an character !");
}

/*!
 * This method returns a newly allocated instance the caller should dealed with by a ParaMEDMEM::DataArrayByte::decrRef.
 * This method performs \b no copy of data. The content is only referenced using ParaMEDMEM::DataArrayByte::useArray with ownership set to \b false.
 * This method throws an INTERP_KERNEL::Exception is it is impossible to match sizes of \b this that is too say \b nbOfCompo=this->_nb_of_elem and \bnbOfTuples==1 or
 * \b nbOfCompo=1 and \bnbOfTuples==this->_nb_of_elem.
 */
DataArrayByte *DataArrayByteTuple::buildDAByte(int nbOfTuples, int nbOfCompo) const
{
  if((_nb_of_compo==nbOfCompo && nbOfTuples==1) || (_nb_of_compo==nbOfTuples && nbOfCompo==1))
    {
      DataArrayByte *ret=DataArrayByte::New();
      ret->useExternalArrayWithRWAccess(_pt,nbOfTuples,nbOfCompo);
      return ret;
    }
  else
    {
      std::ostringstream oss; oss << "DataArrayByteTuple::buildDAByte : unable to build a requested DataArrayByte instance with nbofTuple=" << nbOfTuples << " and nbOfCompo=" << nbOfCompo;
      oss << ".\nBecause the number of elements in this is " << _nb_of_compo << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}

/*!
 * Returns a new instance of DataArrayAsciiChar. The caller is to delete this array
 * using decrRef() as it is no more needed. 
 */
DataArrayAsciiChar *DataArrayAsciiChar::New()
{
  return new DataArrayAsciiChar;
}

/*!
 * Returns a new instance of DataArrayAsciiChar. The caller is to delete this array
 * using decrRef() as it is no more needed. 
 * \param [in] st the string. This input string should have a length greater than 0. If not an excpetion will be thrown.
 */
DataArrayAsciiChar *DataArrayAsciiChar::New(const std::string& st)
{
  return new DataArrayAsciiChar(st);
}

/*!
 * \param [in] st the string. This input string should have a length greater than 0. If not an excpetion will be thrown.
 */
DataArrayAsciiChar::DataArrayAsciiChar(const std::string& st)
{
  std::size_t lgth=st.length();
  if(lgth==0)
    throw INTERP_KERNEL::Exception("DataArrayAsciiChar contructor with string ! Size of input string is null !");
  alloc(1,lgth);
  std::copy(st.begin(),st.begin()+lgth,getPointer());
}

/*!
 * Returns a new instance of DataArrayAsciiChar. The caller is to delete this array
 * using decrRef() as it is no more needed. 
 * This constructor uses \a vst input vector of strings to initialize itself. For all strings whose length is lower than max length of strings in
 * \a vst the remaining locations in memory will be set to character \a defaultChar.
 *
 * \param [in] defaultChar the default character used to fill not defined locations in \a this
 * \param [in] vst vector of strings. This input vector must be non empty. \a this will have its component size set to the max lgth of strings contained
 *             in \a vst. If all strings are empty an INTERP_KERNEL::Exception will be thrown.
 *
 * \throw If input \a vst is empty.
 * \throw If all strings in \a vst are empty.
 */
DataArrayAsciiChar *DataArrayAsciiChar::New(const std::vector<std::string>& vst, char defaultChar)
{
  return new DataArrayAsciiChar(vst,defaultChar);
}

/*!
 * This constructor uses \a vst input vector of strings to initialize itself. For all strings whose length is lower than max length of strings in
 * \a vst the remaining locations in memory will be set to character \a defaultChar.
 *
 * \param [in] defaultChar the default character used to fill not defined locations in \a this
 * \param [in] vst vector of strings. This input vector must be non empty. \a this will have its component size set to the max lgth of strings contained
 *             in \a vst. If all strings are empty an INTERP_KERNEL::Exception will be thrown.
 *
 * \throw If input \a vst is empty.
 * \throw If all strings in \a vst are empty.
 */
DataArrayAsciiChar::DataArrayAsciiChar(const std::vector<std::string>& vst, char defaultChar)
{
  if(vst.empty())
    throw INTERP_KERNEL::Exception("DataArrayAsciiChar contructor with vector of strings ! Empty array !");
  std::size_t nbCompo=0;
  for(std::vector<std::string>::const_iterator it=vst.begin();it!=vst.end();it++)
    nbCompo=std::max(nbCompo,(*it).length());
  if(nbCompo==0)
    throw INTERP_KERNEL::Exception("DataArrayAsciiChar contructor with vector of strings ! All strings in not empty vector are empty !");
  int nbTuples=(int)vst.size();
  alloc(nbTuples,(int)nbCompo);
  char *pt=getPointer();
  for(int i=0;i<nbTuples;i++,pt+=nbCompo)
    {
      const std::string& tmp=vst[i];
      std::size_t sz=tmp.length();
      std::copy(tmp.begin(),tmp.begin()+sz,pt);
      std::fill(pt+sz,pt+nbCompo,defaultChar);
    }
}

DataArrayAsciiCharIterator *DataArrayAsciiChar::iterator()
{
  return new DataArrayAsciiCharIterator(this);
}

/*!
 * Returns a full copy of \a this. For more info on copying data arrays see
 * \ref MEDCouplingArrayBasicsCopyDeep.
 *  \return DataArrayAsciiChar * - a new instance of DataArrayAsciiChar.
 */
DataArrayAsciiChar *DataArrayAsciiChar::deepCpy() const
{
  return new DataArrayAsciiChar(*this);
}

/*!
 * Returns either a \a deep or \a shallow copy of this array. For more info see
 * \ref MEDCouplingArrayBasicsCopyDeep and \ref MEDCouplingArrayBasicsCopyShallow.
 *  \param [in] dCpy - if \a true, a deep copy is returned, else, a shallow one.
 *  \return DataArrayAsciiChar * - either a new instance of DataArrayAsciiChar (if \a dCpy
 *          == \a true) or \a this instance (if \a dCpy == \a false).
 */
DataArrayAsciiChar *DataArrayAsciiChar::performCpy(bool dCpy) const
{
  if(dCpy)
    return deepCpy();
  else
    {
      incrRef();
      return const_cast<DataArrayAsciiChar *>(this);
    }
}

/*!
 * Returns the only one value in \a this, if and only if number of elements
 * (nb of tuples * nb of components) is equal to 1, and that \a this is allocated.
 *  \return char - the sole value stored in \a this array.
 *  \throw If at least one of conditions stated above is not fulfilled.
 */
char DataArrayAsciiChar::asciiCharValue() const
{
  if(isAllocated())
    {
      if(getNbOfElems()==1)
        {
          return *getConstPointer();
        }
      else
        throw INTERP_KERNEL::Exception("DataArrayAsciiChar::asciiCharValue : DataArrayAsciiChar instance is allocated but number of elements is not equal to 1 !");
    }
  else
    throw INTERP_KERNEL::Exception("DataArrayAsciiChar::asciiCharValue : DataArrayAsciiChar instance is not allocated !");
}

DataArrayChar *DataArrayAsciiChar::buildEmptySpecializedDAChar() const
{
  return DataArrayAsciiChar::New();
}

void DataArrayAsciiChar::reprStream(std::ostream& stream) const
{
  stream << "Name of ASCII char array : \"" << _name << "\"\n";
  reprWithoutNameStream(stream);
}

void DataArrayAsciiChar::reprZipStream(std::ostream& stream) const
{
  stream << "Name of ASCII char array : \"" << _name << "\"\n";
  reprZipWithoutNameStream(stream);
}

void DataArrayAsciiChar::reprWithoutNameStream(std::ostream& stream) const
{
  DataArray::reprWithoutNameStream(stream);
  if(_mem.reprHeader(getNumberOfComponents(),stream))
    {
      const char *data=begin();
      int nbOfTuples=getNumberOfTuples();
      int nbCompo=getNumberOfComponents();
      for(int i=0;i<nbOfTuples;i++,data+=nbCompo)
        {
          stream << "Tuple #" << i << " : \"";
          std::copy(data,data+nbCompo,std::ostream_iterator<char>(stream));
          stream << "\"\n";
        }
    }
}

void DataArrayAsciiChar::reprZipWithoutNameStream(std::ostream& stream) const
{
  reprWithoutNameStream(stream);
}

void DataArrayAsciiChar::reprCppStream(const std::string& varName, std::ostream& stream) const
{
  int nbTuples=getNumberOfTuples(),nbComp=getNumberOfComponents();
  const char *data=getConstPointer();
  stream << "DataArrayAsciiChar *" << varName << "=DataArrayAsciiChar::New();" << std::endl;
  if(nbTuples*nbComp>=1)
    {
      stream << "const char " << varName << "Data[" << nbTuples*nbComp << "]={";
      std::copy(data,data+nbTuples*nbComp-1,std::ostream_iterator<char>(stream,","));
      stream << data[nbTuples*nbComp-1] << "};" << std::endl;
      stream << varName << "->useArray(" << varName << "Data,false,CPP_DEALLOC," << nbTuples << "," << nbComp << ");" << std::endl;
    }
  else
    stream << varName << "->alloc(" << nbTuples << "," << nbComp << ");" << std::endl;
  stream << varName << "->setName(\"" << getName() << "\");" << std::endl;
}

/*!
 * Method that gives a quick overvien of \a this for python.
 */
void DataArrayAsciiChar::reprQuickOverview(std::ostream& stream) const
{
  static const std::size_t MAX_NB_OF_BYTE_IN_REPR=300;
  stream << "DataArrayAsciiChar C++ instance at " << this << ". ";
  if(isAllocated())
    {
      int nbOfCompo=(int)_info_on_compo.size();
      if(nbOfCompo>=1)
        {
          int nbOfTuples=getNumberOfTuples();
          stream << "Number of tuples : " << nbOfTuples << ". Number of components : " << nbOfCompo << "." << std::endl;
          reprQuickOverviewData(stream,MAX_NB_OF_BYTE_IN_REPR);
        }
      else
        stream << "Number of components : 0.";
    }
  else
    stream << "*** No data allocated ****";
}

void DataArrayAsciiChar::reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const
{
  const char *data=begin();
  int nbOfTuples=getNumberOfTuples();
  int nbOfCompo=(int)_info_on_compo.size();
  std::ostringstream oss2; oss2 << "[";
  std::string oss2Str(oss2.str());
  bool isFinished=true;
  for(int i=0;i<nbOfTuples && isFinished;i++)
    {
      bool isAscii=true;
      for(int j=0;j<nbOfCompo;j++)
        if(data[j]<32) isAscii=false;
      if(isAscii)
        {
          oss2 << "\'";
          for(int j=0;j<nbOfCompo;j++,data++)
            oss2 << *data;
          oss2 << "\'";
        }
      else
        {
          oss2 << "(";
          for(int j=0;j<nbOfCompo;j++,data++)
            {
              oss2 << (int)*data;
              if(j!=nbOfCompo-1) oss2 << ", ";
            }
          oss2 << ")";
        }
      if(i!=nbOfTuples-1) oss2 << ", ";
      std::string oss3Str(oss2.str());
      if(oss3Str.length()<maxNbOfByteInRepr)
        oss2Str=oss3Str;
      else
        isFinished=false;
    }
  stream << oss2Str;
  if(!isFinished)
    stream << "... ";
  stream << "]";
}

bool DataArrayAsciiChar::isEqualIfNotWhy(const DataArrayChar& other, std::string& reason) const
{
  const DataArrayAsciiChar *otherC=dynamic_cast<const DataArrayAsciiChar *>(&other);
  if(!otherC)
    { reason="this is of type DataArrayAsciiChar whereas other is not a DataArrayAsciiChar instance"; return false; }
  return DataArrayChar::isEqualIfNotWhy(other,reason);
}

DataArrayAsciiCharIterator::DataArrayAsciiCharIterator(DataArrayAsciiChar *da):_da(da),_pt(0),_tuple_id(0),_nb_comp(0),_nb_tuple(0)
{
  if(_da)
    {
      _da->incrRef();
      if(_da->isAllocated())
        {
          _nb_comp=da->getNumberOfComponents();
          _nb_tuple=da->getNumberOfTuples();
          _pt=da->getPointer();
        }
    }
}

DataArrayAsciiCharIterator::~DataArrayAsciiCharIterator()
{
  if(_da)
    _da->decrRef();
}

DataArrayAsciiCharTuple *DataArrayAsciiCharIterator::nextt()
{
  if(_tuple_id<_nb_tuple)
    {
      _tuple_id++;
      DataArrayAsciiCharTuple *ret=new DataArrayAsciiCharTuple(_pt,_nb_comp);
      _pt+=_nb_comp;
      return ret;
    }
  else
    return 0;
}

DataArrayAsciiCharTuple::DataArrayAsciiCharTuple(char *pt, int nbOfComp):_pt(pt),_nb_of_compo(nbOfComp)
{
}

std::string DataArrayAsciiCharTuple::repr() const
{
  std::ostringstream oss;
  std::copy(_pt,_pt+_nb_of_compo,std::ostream_iterator<char>(oss));
  return oss.str();
}

char DataArrayAsciiCharTuple::asciiCharValue() const
{
  if(_nb_of_compo==1)
    return *_pt;
  throw INTERP_KERNEL::Exception("DataArrayAsciiCharTuple::asciiCharValue : DataArrayAsciiCharTuple instance has not exactly 1 component -> Not possible to convert it into an character !");
}

/*!
 * This method returns a newly allocated instance the caller should dealed with by a ParaMEDMEM::DataArrayAsciiChar::decrRef.
 * This method performs \b no copy of data. The content is only referenced using ParaMEDMEM::DataArrayAsciiChar::useArray with ownership set to \b false.
 * This method throws an INTERP_KERNEL::Exception is it is impossible to match sizes of \b this that is too say \b nbOfCompo=this->_nb_of_elem and \bnbOfTuples==1 or
 * \b nbOfCompo=1 and \bnbOfTuples==this->_nb_of_elem.
 */
DataArrayAsciiChar *DataArrayAsciiCharTuple::buildDAAsciiChar(int nbOfTuples, int nbOfCompo) const
{
  if((_nb_of_compo==nbOfCompo && nbOfTuples==1) || (_nb_of_compo==nbOfTuples && nbOfCompo==1))
    {
      DataArrayAsciiChar *ret=DataArrayAsciiChar::New();
      ret->useExternalArrayWithRWAccess(_pt,nbOfTuples,nbOfCompo);
      return ret;
    }
  else
    {
      std::ostringstream oss; oss << "DataArrayAsciiCharTuple::buildDAAsciiChar : unable to build a requested DataArrayAsciiChar instance with nbofTuple=" << nbOfTuples << " and nbOfCompo=" << nbOfCompo;
      oss << ".\nBecause the number of elements in this is " << _nb_of_compo << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
}
