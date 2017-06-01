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
// Author : Anthony Geay (CEA/DEN)

#include "MEDCouplingMemArray.txx"
#include "MCAuto.hxx"

#include <set>
#include <cmath>
#include <limits>
#include <numeric>
#include <algorithm>
#include <functional>

using namespace MEDCoupling;

template class MEDCoupling::MemArray<char>;
template class MEDCoupling::DataArrayTemplate<char>;

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
 * Appends components of another array to components of \a this one, tuple by tuple.
 * So that the number of tuples of \a this array remains the same and the number of 
 * components increases.
 *  \param [in] other - the DataArrayChar to append to \a this one.
 *  \throw If \a this is not allocated.
 *  \throw If \a this and \a other arrays have different number of tuples.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref cpp_mcdataarrayint_meldwith "Here is a C++ example".
 *
 *  \ref py_mcdataarrayint_meldwith "Here is a Python example".
 *  \endif
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
 * Creates a new DataArrayChar containing IDs (indices) of tuples holding value equal to a
 * given one.
 *  \param [in] val - the value to find within \a this.
 *  \return DataArrayChar * - a new instance of DataArrayChar. The caller is to delete this
 *          array using decrRef() as it is no more needed.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 */
DataArrayInt *DataArrayChar::findIdsEqual(char val) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::findIdsEqual : the array must have only one component, you can call 'rearrange' method before !");
  const char *cptr=getConstPointer();
  MCAuto<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
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
DataArrayInt *DataArrayChar::findIdsNotEqual(char val) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::findIdsNotEqual : the array must have only one component, you can call 'rearrange' method before !");
  const char *cptr=getConstPointer();
  MCAuto<DataArrayInt> ret(DataArrayInt::New()); ret->alloc(0,1);
  int nbOfTuples=getNumberOfTuples();
  for(int i=0;i<nbOfTuples;i++,cptr++)
    if(*cptr!=val)
      ret->pushBackSilent(i);
  return ret.retn();
}

/*!
 * This method searches the sequence specified in input parameter \b vals in \b this.
 * This works only for DataArrayChar having number of components equal to one (if not an INTERP_KERNEL::Exception will be thrown).
 * This method differs from DataArrayChar::findIdFirstEqualTuple in that the position is internal raw data is not considered here contrary to DataArrayChar::findIdFirstEqualTuple.
 * \sa DataArrayChar::findIdFirstEqualTuple
 */
int DataArrayChar::findIdSequence(const std::vector<char>& vals) const
{
  checkAllocated();
  int nbOfCompo=getNumberOfComponents();
  if(nbOfCompo!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::findIdSequence : works only for DataArrayChar instance with one component !");
  const char *cptr=getConstPointer();
  std::size_t nbOfVals=getNbOfElems();
  const char *loc=std::search(cptr,cptr+nbOfVals,vals.begin(),vals.end());
  if(loc!=cptr+nbOfVals)
    return std::distance(cptr,loc);
  return -1;
}

/*!
 * This method is an extension of DataArrayChar::findIdFirstEqual method because this method works for DataArrayChar with
 * any number of components excepted 0 (an INTERP_KERNEL::Exception is thrown in this case).
 * This method searches in \b this is there is a tuple that matched the input parameter \b tupl.
 * If any the tuple id is returned. If not -1 is returned.
 * 
 * This method throws an INTERP_KERNEL::Exception if the number of components in \b this mismatches with the size of
 * the input vector. An INTERP_KERNEL::Exception is thrown too if \b this is not allocated.
 *
 * \return tuple id where \b tupl is. -1 if no such tuple exists in \b this.
 * \sa DataArrayChar::findIdSequence.
 */
int DataArrayChar::findIdFirstEqualTuple(const std::vector<char>& tupl) const
{
  checkAllocated();
  int nbOfCompo=getNumberOfComponents();
  if(nbOfCompo==0)
    throw INTERP_KERNEL::Exception("DataArrayChar::findIdFirstEqualTuple : 0 components in 'this' !");
  if(nbOfCompo!=(int)tupl.size())
    {
      std::ostringstream oss; oss << "DataArrayChar::findIdFirstEqualTuple : 'this' contains " << nbOfCompo << " components and searching for a tuple of length " << tupl.size() << " !";
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
 * \sa DataArrayChar::findIdFirstEqualTuple
 */
bool DataArrayChar::presenceOfTuple(const std::vector<char>& tupl) const
{
  return findIdFirstEqualTuple(tupl)!=-1;
}

/*!
 * Returns \a true if a given value is present within \a this one-dimensional array.
 *  \param [in] value - the value to find within \a this array.
 *  \return bool - \a true in case if \a value is present within \a this array.
 *  \throw If \a this is not allocated.
 *  \throw If \a this->getNumberOfComponents() != 1.
 *  \sa findIdFirstEqual()
 */
bool DataArrayChar::presenceOfValue(char value) const
{
  return findIdFirstEqual(value)!=-1;
}

/*!
 * This method expects to be called when number of components of this is equal to one.
 * This method returns true if it exists a tuple so that the value is contained in \b vals.
 * If not any tuple contains one of the values contained in 'vals' false is returned.
 * \sa DataArrayChar::findIdFirstEqual
 */
bool DataArrayChar::presenceOfValue(const std::vector<char>& vals) const
{
  return findIdFirstEqual(vals)!=-1;
}

/*!
 * This method expects to be called when number of components of this is equal to one.
 * This method returns the tuple id, if it exists, of the first tuple equal to \b value.
 * If not any tuple contains \b value -1 is returned.
 * \sa DataArrayChar::presenceOfValue
 */
int DataArrayChar::findIdFirstEqual(char value) const
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
int DataArrayChar::findIdFirstEqual(const std::vector<char>& vals) const
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
 * This method works only on data array with one component.
 * This method returns a newly allocated array storing stored ascendantly tuple ids in \b this so that
 * this[*id] in [\b vmin,\b vmax)
 * 
 * \param [in] vmin begin of range. This value is included in range.
 * \param [in] vmax end of range. This value is \b not included in range.
 * \return a newly allocated data array that the caller should deal with.
 */
DataArrayInt *DataArrayChar::findIdsInRange(char vmin, char vmax) const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayChar::findIdsInRange : this must have exactly one component !");
  const char *cptr=getConstPointer();
  MCAuto<DataArrayInt> ret=DataArrayInt::New(); ret->alloc(0,1);
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
  std::size_t nbOfComp((*it)->getNumberOfComponents());
  int nbt=(*it++)->getNumberOfTuples();
  for(int i=1;it!=a.end();it++,i++)
    {
      if((*it)->getNumberOfComponents()!=nbOfComp)
        throw INTERP_KERNEL::Exception("DataArrayChar::Aggregate : Nb of components mismatch for array aggregation !");
      nbt+=(*it)->getNumberOfTuples();
    }
  MCAuto<DataArrayChar> ret=a[0]->buildEmptySpecializedDAChar();
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
DataArrayByte *DataArrayByte::deepCopy() const
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
DataArrayByte *DataArrayByte::performCopyOrIncrRef(bool dCpy) const
{
  if(dCpy)
    return deepCopy();
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

/*!
 * This method is \b NOT wrapped into python because it can be useful only for performance reasons in C++ context.
 * \throw if \a this is not allocated.
 * \throw if \a this has not exactly one component.
 */
std::vector<bool> DataArrayByte::toVectorOfBool() const
{
  checkAllocated();
  if(getNumberOfComponents()!=1)
    throw INTERP_KERNEL::Exception("DataArrayByte::toVectorOfBool : this method can be used only if this has one component !");
  int nbt(getNumberOfTuples());
  std::vector<bool> ret(nbt,false);
  const char *pt(begin());
  for(int i=0;i<nbt;i++,pt++)
    if(*pt!=0)
      ret[i]=true;
  return ret;
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
 * This method returns a newly allocated instance the caller should dealed with by a MEDCoupling::DataArrayByte::decrRef.
 * This method performs \b no copy of data. The content is only referenced using MEDCoupling::DataArrayByte::useArray with ownership set to \b false.
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
DataArrayAsciiChar *DataArrayAsciiChar::deepCopy() const
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
DataArrayAsciiChar *DataArrayAsciiChar::performCopyOrIncrRef(bool dCpy) const
{
  if(dCpy)
    return deepCopy();
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
 * This method returns a newly allocated instance the caller should dealed with by a MEDCoupling::DataArrayAsciiChar::decrRef.
 * This method performs \b no copy of data. The content is only referenced using MEDCoupling::DataArrayAsciiChar::useArray with ownership set to \b false.
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
