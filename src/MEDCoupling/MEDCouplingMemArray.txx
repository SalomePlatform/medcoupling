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

#ifndef __PARAMEDMEM_MEDCOUPLINGMEMARRAY_TXX__
#define __PARAMEDMEM_MEDCOUPLINGMEMARRAY_TXX__

#include "MEDCouplingMemArray.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"
#include "InterpolationUtils.hxx"
#include "MCAuto.hxx"

#include <sstream>
#include <cstdlib>
#include <algorithm>

namespace MEDCoupling
{
  template<class T>
  void MEDCouplingPointer<T>::setInternal(T *pointer)
  {
    _internal=pointer;
    _external=0;
  }

  template<class T>
  void MEDCouplingPointer<T>::setExternal(const T *pointer)
  {
    _external=pointer;
    _internal=0;
  }

  template<class T>
  MemArray<T>::MemArray(const MemArray<T>& other):_nb_of_elem(0),_nb_of_elem_alloc(0),_ownership(false),_dealloc(0),_param_for_deallocator(0)
  {
    if(!other._pointer.isNull())
      {
        _nb_of_elem_alloc=other._nb_of_elem;
        T *pointer=(T*)malloc(_nb_of_elem_alloc*sizeof(T));
        std::copy(other._pointer.getConstPointer(),other._pointer.getConstPointer()+other._nb_of_elem,pointer);
        useArray(pointer,true,C_DEALLOC,other._nb_of_elem);
      }
  }

  template<class T>
  void MemArray<T>::useArray(const T *array, bool ownership, DeallocType type, std::size_t nbOfElem)
  {
    destroy();
    _nb_of_elem=nbOfElem;
    _nb_of_elem_alloc=nbOfElem;
    if(ownership)
      _pointer.setInternal(const_cast<T *>(array));
    else
      _pointer.setExternal(array);
    _ownership=ownership;
    _dealloc=BuildFromType(type);
  }

  template<class T>
  void MemArray<T>::useExternalArrayWithRWAccess(const T *array, std::size_t nbOfElem)
  {
    destroy();
    _nb_of_elem=nbOfElem;
    _nb_of_elem_alloc=nbOfElem;
    _pointer.setInternal(const_cast<T *>(array));
    _ownership=false;
    _dealloc=CPPDeallocator;
  }

  template<class T>
  void MemArray<T>::writeOnPlace(std::size_t id, T element0, const T *others, std::size_t sizeOfOthers)
  {
    if(id+sizeOfOthers>=_nb_of_elem_alloc)
      reserve(2*_nb_of_elem+sizeOfOthers+1);
    T *pointer=_pointer.getPointer();
    pointer[id]=element0;
    std::copy(others,others+sizeOfOthers,pointer+id+1);
    _nb_of_elem=std::max<std::size_t>(_nb_of_elem,id+sizeOfOthers+1);
  }

  template<class T>
  template<class InputIterator>
  void MemArray<T>::insertAtTheEnd(InputIterator first, InputIterator last)
  {
    T *pointer=_pointer.getPointer();
    while(first!=last)
      {
        if(_nb_of_elem>=_nb_of_elem_alloc)
          {
            reserve(_nb_of_elem_alloc>0?2*_nb_of_elem_alloc:1);
            pointer=_pointer.getPointer();
          }
        pointer[_nb_of_elem++]=*first++;
      }
  }

  template<class T>
  void MemArray<T>::pushBack(T elem)
  {
    if(_nb_of_elem>=_nb_of_elem_alloc)
      reserve(_nb_of_elem_alloc>0?2*_nb_of_elem_alloc:1);
    T *pt=getPointer();
    pt[_nb_of_elem++]=elem;
  }

  template<class T>
  T MemArray<T>::popBack()
  {
    if(_nb_of_elem>0)
      {
        const T *pt=getConstPointer();
        return pt[--_nb_of_elem];
      }
    throw INTERP_KERNEL::Exception("MemArray::popBack : nothing to pop in array !");
  }

  template<class T>
  void MemArray<T>::pack() const
  {
    (const_cast<MemArray<T> * >(this))->reserve(_nb_of_elem);
  }

  template<class T>
  bool MemArray<T>::isEqual(const MemArray<T>& other, T prec, std::string& reason) const
  {
    std::ostringstream oss; oss.precision(15);
    if(_nb_of_elem!=other._nb_of_elem)
      {
        oss << "Number of elements in coarse data of DataArray mismatch : this=" << _nb_of_elem << " other=" << other._nb_of_elem;
        reason=oss.str();
        return false;
      }
    const T *pt1=_pointer.getConstPointer();
    const T *pt2=other._pointer.getConstPointer();
    if(pt1==0 && pt2==0)
      return true;
    if(pt1==0 || pt2==0)
      {
        oss << "coarse data pointer is defined for only one DataArray instance !";
        reason=oss.str();
        return false;
      }
    if(pt1==pt2)
      return true;
    for(std::size_t i=0;i<_nb_of_elem;i++)
      if(pt1[i]-pt2[i]<-prec || (pt1[i]-pt2[i])>prec)
        {
          oss << "The content of data differs at pos #" << i << " of coarse data ! this[i]=" << pt1[i] << " other[i]=" << pt2[i];
          reason=oss.str();
          return false;
        }
    return true;
  }

  /*!
   * \param [in] sl is typically the number of components
   * \return True if a not null pointer is present, False if not.
   */
  template<class T>
  bool MemArray<T>::reprHeader(int sl, std::ostream& stream) const
  {
    stream << "Number of tuples : ";
    if(!_pointer.isNull())
      {
        if(sl!=0)
          stream << _nb_of_elem/sl << std::endl << "Internal memory facts : " << _nb_of_elem << "/" << _nb_of_elem_alloc;
        else
          stream << "Empty Data";
      }
    else
      stream << "No data";
    stream << "\n";
    stream << "Data content :\n";
    bool ret=!_pointer.isNull();
    if(!ret)
      stream << "No data !\n";
    return ret;
  }

  /*!
   * \param [in] sl is typically the number of components
   */
  template<class T>
  void MemArray<T>::repr(int sl, std::ostream& stream) const
  {
    if(reprHeader(sl,stream))
      {
        const T *data=getConstPointer();
        if(_nb_of_elem!=0 && sl!=0)
          {
            std::size_t nbOfTuples=_nb_of_elem/std::abs(sl);
            for(std::size_t i=0;i<nbOfTuples;i++)
              {
                stream << "Tuple #" << i << " : ";
                std::copy(data,data+sl,std::ostream_iterator<T>(stream," "));
                stream << "\n";
                data+=sl;
              }
          }
        else
          stream << "Empty Data\n";
      }
  }

  /*!
   * \param [in] sl is typically the number of components
   */
  template<class T>
  void MemArray<T>::reprZip(int sl, std::ostream& stream) const
  {
    stream << "Number of tuples : ";
    if(!_pointer.isNull())
      {
        if(sl!=0)
          stream << _nb_of_elem/sl;
        else
          stream << "Empty Data";
      }
    else
      stream << "No data";
    stream << "\n";
    stream << "Data content : ";
    const T *data=getConstPointer();
    if(!_pointer.isNull())
      {
        if(_nb_of_elem!=0 && sl!=0)
          {
            std::size_t nbOfTuples=_nb_of_elem/std::abs(sl);
            for(std::size_t i=0;i<nbOfTuples;i++)
              {
                stream << "|";
                std::copy(data,data+sl,std::ostream_iterator<T>(stream," "));
                stream << "| ";
                data+=sl;
              }
            stream << "\n";
          }
        else
          stream << "Empty Data\n";
      }
    else
      stream << "No data !\n";
  }

  /*!
   * \param [in] sl is typically the number of components
   */
  template<class T>
  void MemArray<T>::reprNotTooLong(int sl, std::ostream& stream) const
  {
    if(reprHeader(sl,stream))
      {
        const T *data=getConstPointer();
        if(_nb_of_elem!=0 && sl!=0)
          {
            std::size_t nbOfTuples=_nb_of_elem/std::abs(sl);
            if(nbOfTuples<=1000)
              {
                for(std::size_t i=0;i<nbOfTuples;i++)
                  {
                    stream << "Tuple #" << i << " : ";
                    std::copy(data,data+sl,std::ostream_iterator<T>(stream," "));
                    stream << "\n";
                    data+=sl;
                  }
              }
            else
              {// too much tuples -> print the 3 first tuples and 3 last.
                stream << "Tuple #0 : ";
                std::copy(data,data+sl,std::ostream_iterator<T>(stream," ")); stream << "\n";
                stream << "Tuple #1 : ";
                std::copy(data+sl,data+2*sl,std::ostream_iterator<T>(stream," ")); stream << "\n";
                stream << "Tuple #2 : ";
                std::copy(data+2*sl,data+3*sl,std::ostream_iterator<T>(stream," ")); stream << "\n";
                stream << "...\n";
                stream << "Tuple #" << nbOfTuples-3 << " : ";
                std::copy(data+(nbOfTuples-3)*sl,data+(nbOfTuples-2)*sl,std::ostream_iterator<T>(stream," ")); stream << "\n";
                stream << "Tuple #" << nbOfTuples-2 << " : ";
                std::copy(data+(nbOfTuples-2)*sl,data+(nbOfTuples-1)*sl,std::ostream_iterator<T>(stream," ")); stream << "\n";
                stream << "Tuple #" << nbOfTuples-1 << " : ";
                std::copy(data+(nbOfTuples-1)*sl,data+nbOfTuples*sl,std::ostream_iterator<T>(stream," ")); stream << "\n";
              }
          }
        else
          stream << "Empty Data\n";
      }
  }

  template<class T>
  void MemArray<T>::fillWithValue(const T& val)
  {
    T *pt=_pointer.getPointer();
    std::fill(pt,pt+_nb_of_elem,val);
  }

  template<class T>
  T *MemArray<T>::fromNoInterlace(int nbOfComp) const
  {
    if(nbOfComp<1)
      throw INTERP_KERNEL::Exception("MemArray<T>::fromNoInterlace : number of components must be > 0 !");
    const T *pt=_pointer.getConstPointer();
    std::size_t nbOfTuples=_nb_of_elem/nbOfComp;
    T *ret=(T*)malloc(_nb_of_elem*sizeof(T));
    T *w=ret;
    for(std::size_t i=0;i<nbOfTuples;i++)
      for(int j=0;j<nbOfComp;j++,w++)
        *w=pt[j*nbOfTuples+i];
    return ret;
  }

  template<class T>
  T *MemArray<T>::toNoInterlace(int nbOfComp) const
  {
    if(nbOfComp<1)
      throw INTERP_KERNEL::Exception("MemArray<T>::toNoInterlace : number of components must be > 0 !");
    const T *pt=_pointer.getConstPointer();
    std::size_t nbOfTuples=_nb_of_elem/nbOfComp;
    T *ret=(T*)malloc(_nb_of_elem*sizeof(T));
    T *w=ret;
    for(int i=0;i<nbOfComp;i++)
      for(std::size_t j=0;j<nbOfTuples;j++,w++)
        *w=pt[j*nbOfComp+i];
    return ret;
  }

  template<class T>
  void MemArray<T>::sort(bool asc)
  {
    T *pt=_pointer.getPointer();
    if(asc)
      std::sort(pt,pt+_nb_of_elem);
    else
      {
        typename std::reverse_iterator<T *> it1(pt+_nb_of_elem);
        typename std::reverse_iterator<T *> it2(pt);
        std::sort(it1,it2);
      }
  }

  template<class T>
  void MemArray<T>::reverse(int nbOfComp)
  {
    if(nbOfComp<1)
      throw INTERP_KERNEL::Exception("MemArray<T>::reverse : only supported with 'this' array with ONE or more than ONE component !");
    T *pt=_pointer.getPointer();
    if(nbOfComp==1)
      {
        std::reverse(pt,pt+_nb_of_elem);
        return ;
      }
    else
      {
        T *pt2=pt+_nb_of_elem-nbOfComp;
        std::size_t nbOfTuples=_nb_of_elem/nbOfComp;
        for(std::size_t i=0;i<nbOfTuples/2;i++,pt+=nbOfComp,pt2-=nbOfComp)
          {
            for(int j=0;j<nbOfComp;j++)
              std::swap(pt[j],pt2[j]);
          }
      }
  }

  template<class T>
  void MemArray<T>::alloc(std::size_t nbOfElements)
  {
    destroy();
    _nb_of_elem=nbOfElements;
    _nb_of_elem_alloc=nbOfElements;
    _pointer.setInternal((T*)malloc(_nb_of_elem_alloc*sizeof(T)));
    _ownership=true;
    _dealloc=CDeallocator;
  }

  /*!
   * This method performs systematically an allocation of \a newNbOfElements elements in \a this.
   * \a _nb_of_elem and \a _nb_of_elem_alloc will \b NOT be systematically equal (contrary to MemArray<T>::reAlloc method.
   * So after the call of this method \a _nb_of_elem will be equal tostd::min<std::size_t>(_nb_of_elem,newNbOfElements) and \a _nb_of_elem_alloc equal to 
   * \a newNbOfElements. This method is typically used to perform a pushBack to avoid systematic allocations-copy-deallocation.
   * So after the call of this method the accessible content is perfectly set.
   * 
   * So this method should not be confused with MemArray<T>::reserve that is close to MemArray<T>::reAlloc but not same.
   */
  template<class T>
  void MemArray<T>::reserve(std::size_t newNbOfElements)
  {
    if(_nb_of_elem_alloc==newNbOfElements)
      return ;
    T *pointer=(T*)malloc(newNbOfElements*sizeof(T));
    std::copy(_pointer.getConstPointer(),_pointer.getConstPointer()+std::min<std::size_t>(_nb_of_elem,newNbOfElements),pointer);
    if(_ownership)
      DestroyPointer(const_cast<T *>(_pointer.getConstPointer()),_dealloc,_param_for_deallocator);//Do not use getPointer because in case of _external
    _pointer.setInternal(pointer);
    _nb_of_elem=std::min<std::size_t>(_nb_of_elem,newNbOfElements);
    _nb_of_elem_alloc=newNbOfElements;
    _ownership=true;
    _dealloc=CDeallocator;
    _param_for_deallocator=0;
  }

  /*!
   * This method performs systematically an allocation of \a newNbOfElements elements in \a this.
   * \a _nb_of_elem and \a _nb_of_elem_alloc will be equal even if only std::min<std::size_t>(_nb_of_elem,newNbOfElements) come from the .
   * The remaing part of the new allocated chunk are available but not set previouly !
   * 
   * So this method should not be confused with MemArray<T>::reserve that is close to MemArray<T>::reAlloc but not same.
   */
  template<class T>
  void MemArray<T>::reAlloc(std::size_t newNbOfElements)
  {
    if(_nb_of_elem==newNbOfElements)
      return ;
    T *pointer=(T*)malloc(newNbOfElements*sizeof(T));
    std::copy(_pointer.getConstPointer(),_pointer.getConstPointer()+std::min<std::size_t>(_nb_of_elem,newNbOfElements),pointer);
    if(_ownership)
      DestroyPointer(const_cast<T *>(_pointer.getConstPointer()),_dealloc,_param_for_deallocator);//Do not use getPointer because in case of _external
    _pointer.setInternal(pointer);
    _nb_of_elem=newNbOfElements;
    _nb_of_elem_alloc=newNbOfElements;
    _ownership=true;
    _dealloc=CDeallocator;
    _param_for_deallocator=0;
  }

  template<class T>
  void MemArray<T>::CPPDeallocator(void *pt, void *param)
  {
    delete [] reinterpret_cast<T*>(pt);
  }

  template<class T>
  void MemArray<T>::CDeallocator(void *pt, void *param)
  {
    free(pt);
  }

  template<class T>
  typename MemArray<T>::Deallocator MemArray<T>::BuildFromType(DeallocType type)
  {
    switch(type)
    {
      case CPP_DEALLOC:
        return CPPDeallocator;
      case C_DEALLOC:
        return CDeallocator;
      default:
        throw INTERP_KERNEL::Exception("Invalid deallocation requested ! Unrecognized enum DeallocType !");
    }
  }

  template<class T>
  void MemArray<T>::DestroyPointer(T *pt, typename MemArray<T>::Deallocator dealloc, void *param)
  {
    if(dealloc)
      dealloc(pt,param);
  }

  template<class T>
  void MemArray<T>::destroy()
  {
    if(_ownership)
      DestroyPointer(const_cast<T *>(_pointer.getConstPointer()),_dealloc,_param_for_deallocator);//Do not use getPointer because in case of _external
    _pointer.null();
    _ownership=false;
    _dealloc=NULL;
    _param_for_deallocator=NULL;
    _nb_of_elem=0;
    _nb_of_elem_alloc=0;
  }

  template<class T>
  MemArray<T> &MemArray<T>::operator=(const MemArray<T>& other)
  {
    alloc(other._nb_of_elem);
    std::copy(other._pointer.getConstPointer(),other._pointer.getConstPointer()+_nb_of_elem,_pointer.getPointer());
    return *this;
  }

  //////////////////////////////////
  
  template<class T>
  std::size_t DataArrayTemplate<T>::getHeapMemorySizeWithoutChildren() const
  {
    std::size_t sz(_mem.getNbOfElemAllocated());
    sz*=sizeof(T);
    return DataArray::getHeapMemorySizeWithoutChildren()+sz;
  }
  
  /*!
   * Allocates the raw data in memory. If the memory was already allocated, then it is
   * freed and re-allocated. See an example of this method use
   * \ref MEDCouplingArraySteps1WC "here".
   *  \param [in] nbOfTuple - number of tuples of data to allocate.
   *  \param [in] nbOfCompo - number of components of data to allocate.
   *  \throw If \a nbOfTuple < 0 or \a nbOfCompo < 0.
   */
  template<class T>
  void DataArrayTemplate<T>::alloc(int nbOfTuple, int nbOfCompo)
  {
    if(nbOfTuple<0 || nbOfCompo<0)
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::alloc : request for negative length of data !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    _info_on_compo.resize(nbOfCompo);
    _mem.alloc(nbOfCompo*(std::size_t)nbOfTuple);
    declareAsNew();
  }

  /*!
   * Sets a C array to be used as raw data of \a this. The previously set info
   *  of components is retained and re-sized. 
   * For more info see \ref MEDCouplingArraySteps1.
   *  \param [in] array - the C array to be used as raw data of \a this.
   *  \param [in] ownership - if \a true, \a array will be deallocated at destruction of \a this.
   *  \param [in] type - specifies how to deallocate \a array. If \a type == MEDCoupling::CPP_DEALLOC,
   *                     \c delete [] \c array; will be called. If \a type == MEDCoupling::C_DEALLOC,
   *                     \c free(\c array ) will be called.
   *  \param [in] nbOfTuple - new number of tuples in \a this.
   *  \param [in] nbOfCompo - new number of components in \a this.
   */
  template<class T>
  void DataArrayTemplate<T>::useArray(const T *array, bool ownership, DeallocType type, int nbOfTuple, int nbOfCompo)
  {
    _info_on_compo.resize(nbOfCompo);
    _mem.useArray(array,ownership,type,(std::size_t)nbOfTuple*nbOfCompo);
    declareAsNew();
  }
  
  template<class T>
  void DataArrayTemplate<T>::useExternalArrayWithRWAccess(const T *array, int nbOfTuple, int nbOfCompo)
  {
    _info_on_compo.resize(nbOfCompo);
    _mem.useExternalArrayWithRWAccess(array,(std::size_t)nbOfTuple*nbOfCompo);
    declareAsNew();
  }

  /*!
   * Returns a value located at specified tuple and component.
   * This method is equivalent to DataArrayTemplate<T>::getIJ() except that validity of
   * parameters is checked. So this method is safe but expensive if used to go through
   * all values of \a this.
   *  \param [in] tupleId - index of tuple of interest.
   *  \param [in] compoId - index of component of interest.
   *  \return double - value located by \a tupleId and \a compoId.
   *  \throw If \a this is not allocated.
   *  \throw If condition <em>( 0 <= tupleId < this->getNumberOfTuples() )</em> is violated.
   *  \throw If condition <em>( 0 <= compoId < this->getNumberOfComponents() )</em> is violated.
   */
  template<class T>
  T DataArrayTemplate<T>::getIJSafe(int tupleId, int compoId) const
  {
    checkAllocated();
    if(tupleId<0 || tupleId>=getNumberOfTuples())
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::getIJSafe : request for tupleId " << tupleId << " should be in [0," << getNumberOfTuples() << ") !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    if(compoId<0 || compoId>=getNumberOfComponents())
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::getIJSafe : request for compoId " << compoId << " should be in [0," << getNumberOfComponents() << ") !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    return _mem[tupleId*_info_on_compo.size()+compoId];
  }

  /*!
   * This method \b do \b not modify content of \a this. It only modify its memory footprint if the allocated memory is to high regarding real data to store.
   *
   * \sa DataArray::getHeapMemorySizeWithoutChildren, DataArrayTemplate<T>::reserve
   */
  template<class T>
  void DataArrayTemplate<T>::pack() const
  {
    _mem.pack();
  }

  /*!
   * Checks if raw data is allocated. Read more on the raw data
   * in \ref MEDCouplingArrayBasicsTuplesAndCompo "DataArrays infos" for more information.
   *  \return bool - \a true if the raw data is allocated, \a false else.
   */
  template<class T>
  bool DataArrayTemplate<T>::isAllocated() const
  {
    return getConstPointer()!=0;
  }
  
  /*!
   * Checks if raw data is allocated and throws an exception if it is not the case.
   *  \throw If the raw data is not allocated.
   */
  template<class T>
  void DataArrayTemplate<T>::checkAllocated() const
  {
    if(!isAllocated())
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::checkAllocated : Array is defined but not allocated ! Call alloc or setValues method first !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  }
  
  /*!
   * This method desallocated \a this without modification of informations relative to the components.
   * After call of this method, DataArrayDouble::isAllocated will return false.
   * If \a this is already not allocated, \a this is let unchanged.
   */
  template<class T>
  void DataArrayTemplate<T>::desallocate()
  {
    _mem.destroy();
  }

  /*!
   * This method reserve nbOfElems elements in memory ( nbOfElems*8 bytes ) \b without impacting the number of tuples in \a this.
   * If \a this has already been allocated, this method checks that \a this has only one component. If not an INTERP_KERNEL::Exception will be thrown.
   * If \a this has not already been allocated, number of components is set to one.
   * This method allows to reduce number of reallocations on invokation of DataArrayDouble::pushBackSilent and DataArrayDouble::pushBackValsSilent on \a this.
   * 
   * \sa DataArrayDouble::pack, DataArrayDouble::pushBackSilent, DataArrayDouble::pushBackValsSilent
   */
  template<class T>
  void DataArrayTemplate<T>::reserve(std::size_t nbOfElems)
  {
    int nbCompo(getNumberOfComponents());
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
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::reserve : not available for DataArrayDouble with number of components different than 1 !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  }
  
  /*!
   * This method adds at the end of \a this the single value \a val. This method do \b not update its time label to avoid useless incrementation
   * of counter. So the caller is expected to call TimeLabel::declareAsNew on \a this at the end of the push session.
   *
   * \param [in] val the value to be added in \a this
   * \throw If \a this has already been allocated with number of components different from one.
   * \sa DataArrayDouble::pushBackValsSilent
   */
  template<class T>
  void DataArrayTemplate<T>::pushBackSilent(T val)
  {
    int nbCompo(getNumberOfComponents());
    if(nbCompo==1)
      _mem.pushBack(val);
    else if(nbCompo==0)
      {
        _info_on_compo.resize(1);
        _mem.pushBack(val);
      }
    else
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::pushBackSilent : not available for DataArrayDouble with number of components different than 1 !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  }
  
  /*!
   * This method adds at the end of \a this a serie of values [\c valsBg,\c valsEnd). This method do \b not update its time label to avoid useless incrementation
   * of counter. So the caller is expected to call TimeLabel::declareAsNew on \a this at the end of the push session.
   *
   *  \param [in] valsBg - an array of values to push at the end of \c this.
   *  \param [in] valsEnd - specifies the end of the array \a valsBg, so that
   *              the last value of \a valsBg is \a valsEnd[ -1 ].
   * \throw If \a this has already been allocated with number of components different from one.
   * \sa DataArrayDouble::pushBackSilent
   */
  template<class T>
  void DataArrayTemplate<T>::pushBackValsSilent(const T *valsBg, const T *valsEnd)
  {
    int nbCompo(getNumberOfComponents());
    if(nbCompo==1)
      _mem.insertAtTheEnd(valsBg,valsEnd);
    else if(nbCompo==0)
      {
        _info_on_compo.resize(1);
        _mem.insertAtTheEnd(valsBg,valsEnd);
      }
    else
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::pushBackValsSilent : not available for DataArrayDouble with number of components different than 1 !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  }
  
  /*!
   * This method returns silently ( without updating time label in \a this ) the last value, if any and suppress it.
   * \throw If \a this is already empty.
   * \throw If \a this has number of components different from one.
   */
  template<class T>
  T DataArrayTemplate<T>::popBackSilent()
  {
    if(getNumberOfComponents()==1)
      return _mem.popBack();
    else
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::popBackSilent : not available for DataArrayDouble with number of components different than 1 !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
  }
  
  /*!
   * Allocates the raw data in memory. If exactly same memory as needed already
   * allocated, it is not re-allocated.
   *  \param [in] nbOfTuple - number of tuples of data to allocate.
   *  \param [in] nbOfCompo - number of components of data to allocate.
   *  \throw If \a nbOfTuple < 0 or \a nbOfCompo < 0.
   */
  template<class T>
  void DataArrayTemplate<T>::allocIfNecessary(int nbOfTuple, int nbOfCompo)
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
   * Checks the number of tuples.
   *  \return bool - \a true if getNumberOfTuples() == 0, \a false else.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  bool DataArrayTemplate<T>::empty() const
  {
    checkAllocated();
    return getNumberOfTuples()==0;
  }

  /*!
   * Copies all the data from another DataArrayDouble. For more info see
   * \ref MEDCouplingArrayBasicsCopyDeepAssign.
   *  \param [in] other - another instance of DataArrayDouble to copy data from.
   *  \throw If the \a other is not allocated.
   */
  template<class T>
  void DataArrayTemplate<T>::deepCopyFrom(const DataArrayTemplate<T>& other)
  {
    other.checkAllocated();
    int nbOfTuples(other.getNumberOfTuples()),nbOfComp(other.getNumberOfComponents());
    allocIfNecessary(nbOfTuples,nbOfComp);
    std::size_t nbOfElems((std::size_t)nbOfTuples*nbOfComp);
    T *pt(getPointer());
    const T *ptI(other.begin());
    for(std::size_t i=0;i<nbOfElems;i++)
      pt[i]=ptI[i];
    copyStringInfoFrom(other);
  }

  /*!
   * Reverse the array values.
   *  \throw If \a this->getNumberOfComponents() < 1.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  void DataArrayTemplate<T>::reverse()
  {
    checkAllocated();
    _mem.reverse(getNumberOfComponents());
    declareAsNew();
  }

  /*!
   * Assign \a val to all values in \a this array. To know more on filling arrays see
   * \ref MEDCouplingArrayFill.
   *  \param [in] val - the value to fill with.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  void DataArrayTemplate<T>::fillWithValue(T val)
  {
    checkAllocated();
    _mem.fillWithValue(val);
    declareAsNew();
  }

  /*!
   * Changes number of tuples in the array. If the new number of tuples is smaller
   * than the current number the array is truncated, otherwise the array is extended.
   *  \param [in] nbOfTuples - new number of tuples. 
   *  \throw If \a this is not allocated.
   *  \throw If \a nbOfTuples is negative.
   */
  template<class T>
  void DataArrayTemplate<T>::reAlloc(int nbOfTuples)
  {
    if(nbOfTuples<0)
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::reAlloc : input new number of tuples should be >=0 !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    checkAllocated();
    _mem.reAlloc(getNumberOfComponents()*(std::size_t)nbOfTuples);
    declareAsNew();
  }

  /*!
   * Permutes values of \a this array as required by \a old2New array. The values are
   * permuted so that \c new[ \a old2New[ i ]] = \c old[ i ]. Number of tuples remains
   * the same as in \c this one.
   * If a permutation reduction is needed, subArray() or selectByTupleId() should be used.
   * For more info on renumbering see \ref numbering.
   *  \param [in] old2New - C array of length equal to \a this->getNumberOfTuples()
   *     giving a new position for i-th old value.
   */
  template<class T>
  void DataArrayTemplate<T>::renumberInPlace(const int *old2New)
  {
    checkAllocated();
    int nbTuples(getNumberOfTuples()),nbOfCompo(getNumberOfComponents());
    T *tmp(new T[nbTuples*nbOfCompo]);
    const T *iptr(begin());
    for(int i=0;i<nbTuples;i++)
      {
        int v=old2New[i];
        if(v>=0 && v<nbTuples)
          std::copy(iptr+nbOfCompo*i,iptr+nbOfCompo*(i+1),tmp+nbOfCompo*v);
        else
          {
            std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::renumberInPlace : At place #" << i << " value is " << v << " ! Should be in [0," << nbTuples << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    std::copy(tmp,tmp+nbTuples*nbOfCompo,getPointer());
    delete [] tmp;
    declareAsNew();
  }


  /*!
   * Permutes values of \a this array as required by \a new2Old array. The values are
   * permuted so that \c new[ i ] = \c old[ \a new2Old[ i ]]. Number of tuples remains
   * the same as in \c this one.
   * For more info on renumbering see \ref numbering.
   *  \param [in] new2Old - C array of length equal to \a this->getNumberOfTuples()
   *     giving a previous position of i-th new value.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   */
  template<class T>
  void DataArrayTemplate<T>::renumberInPlaceR(const int *new2Old)
  {
    checkAllocated();
    int nbTuples(getNumberOfTuples()),nbOfCompo(getNumberOfComponents());
    T *tmp(new T[nbTuples*nbOfCompo]);
    const T *iptr(begin());
    for(int i=0;i<nbTuples;i++)
      {
        int v=new2Old[i];
        if(v>=0 && v<nbTuples)
          std::copy(iptr+nbOfCompo*v,iptr+nbOfCompo*(v+1),tmp+nbOfCompo*i);
        else
          {
            std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::renumberInPlaceR : At place #" << i << " value is " << v << " ! Should be in [0," << nbTuples << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    std::copy(tmp,tmp+nbTuples*nbOfCompo,getPointer());
    delete [] tmp;
    declareAsNew();
  }

  /*!
   * Sorts values of the array.
   *  \param [in] asc - \a true means ascending order, \a false, descending.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   */
  template<class T>
  void DataArrayTemplate<T>::sort(bool asc)
  {
    checkAllocated();
    if(getNumberOfComponents()!=1)
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::sort : only supported with 'this' array with ONE component !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    _mem.sort(asc);
    declareAsNew();
  }

  /*!
   * Returns a copy of \a this array with values permuted as required by \a old2New array.
   * The values are permuted so that  \c new[ \a old2New[ i ]] = \c old[ i ].
   * Number of tuples in the result array remains the same as in \c this one.
   * If a permutation reduction is needed, renumberAndReduce() should be used.
   * For more info on renumbering see \ref numbering.
   *  \param [in] old2New - C array of length equal to \a this->getNumberOfTuples()
   *          giving a new position for i-th old value.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::renumber(const int *old2New) const
  {
    checkAllocated();
    int nbTuples(getNumberOfTuples()),nbOfCompo(getNumberOfComponents());
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    ret->alloc(nbTuples,nbOfCompo);
    ret->copyStringInfoFrom(*this);
    const T *iptr(begin());
    T *optr(ret->getPointer());
    for(int i=0;i<nbTuples;i++)
      std::copy(iptr+nbOfCompo*i,iptr+nbOfCompo*(i+1),optr+nbOfCompo*old2New[i]);
    ret->copyStringInfoFrom(*this);
    return ret.retn();
  }

  /*!
   * Returns a copy of \a this array with values permuted as required by \a new2Old array.
   * The values are permuted so that  \c new[ i ] = \c old[ \a new2Old[ i ]]. Number of
   * tuples in the result array remains the same as in \c this one.
   * If a permutation reduction is needed, subArray() or selectByTupleId() should be used.
   * For more info on renumbering see \ref numbering.
   *  \param [in] new2Old - C array of length equal to \a this->getNumberOfTuples()
   *     giving a previous position of i-th new value.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::renumberR(const int *new2Old) const
  {
    checkAllocated();
    int nbTuples(getNumberOfTuples()),nbOfCompo(getNumberOfComponents());
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    ret->alloc(nbTuples,nbOfCompo);
    ret->copyStringInfoFrom(*this);
    const T *iptr(getConstPointer());
    T *optr(ret->getPointer());
    for(int i=0;i<nbTuples;i++)
      std::copy(iptr+nbOfCompo*new2Old[i],iptr+nbOfCompo*(new2Old[i]+1),optr+i*nbOfCompo);
    ret->copyStringInfoFrom(*this);
    return ret.retn();
  }

  /*!
   * Returns a shorten and permuted copy of \a this array. The new DataArrayDouble is
   * of size \a newNbOfTuple and it's values are permuted as required by \a old2New array.
   * The values are permuted so that  \c new[ \a old2New[ i ]] = \c old[ i ] for all
   * \a old2New[ i ] >= 0. In other words every i-th tuple in \a this array, for which 
   * \a old2New[ i ] is negative, is missing from the result array.
   * For more info on renumbering see \ref numbering.
   *  \param [in] old2New - C array of length equal to \a this->getNumberOfTuples()
   *     giving a new position for i-th old tuple and giving negative position for
   *     for i-th old tuple that should be omitted.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::renumberAndReduce(const int *old2New, int newNbOfTuple) const
  {
    checkAllocated();
    int nbTuples(getNumberOfTuples()),nbOfCompo(getNumberOfComponents());
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    ret->alloc(newNbOfTuple,nbOfCompo);
    const T *iptr=getConstPointer();
    T *optr=ret->getPointer();
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
   * Returns a shorten and permuted copy of \a this array. The new DataArrayDouble is
   * of size \a new2OldEnd - \a new2OldBg and it's values are permuted as required by
   * \a new2OldBg array.
   * The values are permuted so that  \c new[ i ] = \c old[ \a new2OldBg[ i ]].
   * This method is equivalent to renumberAndReduce() except that convention in input is
   * \c new2old and \b not \c old2new.
   * For more info on renumbering see \ref numbering.
   *  \param [in] new2OldBg - pointer to the beginning of a permutation array that gives a
   *              tuple index in \a this array to fill the i-th tuple in the new array.
   *  \param [in] new2OldEnd - specifies the end of the permutation array that starts at
   *              \a new2OldBg, so that pointer to a tuple index (\a pi) varies as this:
   *              \a new2OldBg <= \a pi < \a new2OldEnd.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::mySelectByTupleId(const int *new2OldBg, const int *new2OldEnd) const
  {
    checkAllocated();
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    int nbComp(getNumberOfComponents());
    ret->alloc((int)std::distance(new2OldBg,new2OldEnd),nbComp);
    ret->copyStringInfoFrom(*this);
    T *pt(ret->getPointer());
    const T *srcPt(getConstPointer());
    int i(0);
    for(const int *w=new2OldBg;w!=new2OldEnd;w++,i++)
      std::copy(srcPt+(*w)*nbComp,srcPt+((*w)+1)*nbComp,pt+i*nbComp);
    ret->copyStringInfoFrom(*this);
    return ret.retn();
  }

  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::mySelectByTupleId(const DataArrayInt& di) const
  {
    return DataArrayTemplate<T>::mySelectByTupleId(di.begin(),di.end());
  }
  
  /*!
   * Returns a shorten and permuted copy of \a this array. The new DataArrayDouble is
   * of size \a new2OldEnd - \a new2OldBg and it's values are permuted as required by
   * \a new2OldBg array.
   * The values are permuted so that  \c new[ i ] = \c old[ \a new2OldBg[ i ]].
   * This method is equivalent to renumberAndReduce() except that convention in input is
   * \c new2old and \b not \c old2new.
   * This method is equivalent to selectByTupleId() except that it prevents coping data
   * from behind the end of \a this array.
   * For more info on renumbering see \ref numbering.
   *  \param [in] new2OldBg - pointer to the beginning of a permutation array that gives a
   *              tuple index in \a this array to fill the i-th tuple in the new array.
   *  \param [in] new2OldEnd - specifies the end of the permutation array that starts at
   *              \a new2OldBg, so that pointer to a tuple index (\a pi) varies as this:
   *              \a new2OldBg <= \a pi < \a new2OldEnd.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   *  \throw If \a new2OldEnd - \a new2OldBg > \a this->getNumberOfTuples().
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::mySelectByTupleIdSafe(const int *new2OldBg, const int *new2OldEnd) const
  {
    checkAllocated();
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    int nbComp(getNumberOfComponents()),oldNbOfTuples(getNumberOfTuples());
    ret->alloc((int)std::distance(new2OldBg,new2OldEnd),nbComp);
    ret->copyStringInfoFrom(*this);
    T *pt(ret->getPointer());
    const T *srcPt(getConstPointer());
    int i(0);
    for(const int *w=new2OldBg;w!=new2OldEnd;w++,i++)
      if(*w>=0 && *w<oldNbOfTuples)
        std::copy(srcPt+(*w)*nbComp,srcPt+((*w)+1)*nbComp,pt+i*nbComp);
      else
        {
          std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::selectByTupleIdSafe : some ids has been detected to be out of [0,this->getNumberOfTuples) !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    ret->copyStringInfoFrom(*this);
    return ret.retn();
  }

  /*!
   * Changes the number of components within \a this array so that its raw data **does
   * not** change, instead splitting this data into tuples changes.
   *  \warning This method erases all (name and unit) component info set before!
   *  \param [in] newNbOfComp - number of components for \a this array to have.
   *  \throw If \a this is not allocated
   *  \throw If getNbOfElems() % \a newNbOfCompo != 0.
   *  \throw If \a newNbOfCompo is lower than 1.
   *  \throw If the rearrange method would lead to a number of tuples higher than 2147483647 (maximal capacity of int32 !).
   *  \warning This method erases all (name and unit) component info set before!
   */
  template<class T>
  void DataArrayTemplate<T>::rearrange(int newNbOfCompo)
  {
    checkAllocated();
    if(newNbOfCompo<1)
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::rearrange : input newNbOfCompo must be > 0 !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    std::size_t nbOfElems=getNbOfElems();
    if(nbOfElems%newNbOfCompo!=0)
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::rearrange : nbOfElems%newNbOfCompo!=0 !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    if(nbOfElems/newNbOfCompo>(std::size_t)std::numeric_limits<int>::max())
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::rearrange : the rearrangement leads to too high number of tuples (> 2147483647) !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    _info_on_compo.clear();
    _info_on_compo.resize(newNbOfCompo);
    declareAsNew();
  }

  /*!
   * Changes the number of components within \a this array to be equal to its number
   * of tuples, and inversely its number of tuples to become equal to its number of 
   * components. So that its raw data **does not** change, instead splitting this
   * data into tuples changes.
   *  \warning This method erases all (name and unit) component info set before!
   *  \warning Do not confuse this method with fromNoInterlace() and toNoInterlace()!
   *  \throw If \a this is not allocated.
   *  \sa rearrange()
   */
  template<class T>
  void DataArrayTemplate<T>::transpose()
  {
    checkAllocated();
    int nbOfTuples(getNumberOfTuples());
    rearrange(nbOfTuples);
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
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::changeNbOfComponents(int newNbOfComp, T dftValue) const
  {
    checkAllocated();
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    ret->alloc(getNumberOfTuples(),newNbOfComp);
    const T *oldc(getConstPointer());
    T *nc(ret->getPointer());
    int nbOfTuples(getNumberOfTuples()),oldNbOfComp(getNumberOfComponents());
    int dim(std::min(oldNbOfComp,newNbOfComp));
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
   * The new DataArrayDouble has the same number of tuples but includes components
   * specified by \a compoIds parameter. So that getNbOfElems() of the result array
   * can be either less, same or more than \a this->getNbOfElems().
   *  \param [in] compoIds - sequence of zero based indices of components to include
   *              into the new array.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   *  \throw If \a this is not allocated.
   *  \throw If a component index (\a i) is not valid: 
   *         \a i < 0 || \a i >= \a this->getNumberOfComponents().
   *
   *  \if ENABLE_EXAMPLES
   *  \ref py_mcdataarraydouble_KeepSelectedComponents "Here is a Python example".
   *  \endif
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::myKeepSelectedComponents(const std::vector<int>& compoIds) const
  {
    checkAllocated();
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    std::size_t newNbOfCompo(compoIds.size());
    int oldNbOfCompo(getNumberOfComponents());
    for(std::vector<int>::const_iterator it=compoIds.begin();it!=compoIds.end();it++)
      if((*it)<0 || (*it)>=oldNbOfCompo)
        {
          std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::keepSelectedComponents : invalid requested component : " << *it << " whereas it should be in [0," << oldNbOfCompo << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    int nbOfTuples(getNumberOfTuples());
    ret->alloc(nbOfTuples,(int)newNbOfCompo);
    ret->copyPartOfStringInfoFrom(*this,compoIds);
    const T *oldc(getConstPointer());
    T *nc(ret->getPointer());
    for(int i=0;i<nbOfTuples;i++)
      for(std::size_t j=0;j<newNbOfCompo;j++,nc++)
        *nc=oldc[i*oldNbOfCompo+compoIds[j]];
    return ret.retn();
  }

  /*!
   * Returns a shorten copy of \a this array. The new DataArrayDouble contains all
   * tuples starting from the \a tupleIdBg-th tuple and including all tuples located before
   * the \a tupleIdEnd-th one. This methods has a similar behavior as std::string::substr().
   * This method is a specialization of selectByTupleIdSafeSlice().
   *  \param [in] tupleIdBg - index of the first tuple to copy from \a this array.
   *  \param [in] tupleIdEnd - index of the tuple before which the tuples to copy are located.
   *          If \a tupleIdEnd == -1, all the tuples till the end of \a this array are copied.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   *  \throw If \a tupleIdBg < 0.
   *  \throw If \a tupleIdBg > \a this->getNumberOfTuples().
   *  \throw If \a tupleIdEnd != -1 && \a tupleIdEnd < \a this->getNumberOfTuples().
   *  \sa DataArrayDouble::selectByTupleIdSafeSlice
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::subArray(int tupleIdBg, int tupleIdEnd) const
  {
    checkAllocated();
    int nbt(getNumberOfTuples());
    if(tupleIdBg<0)
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::subArray : The tupleIdBg parameter must be greater than 0 !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    if(tupleIdBg>nbt)
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << ":subArray : The tupleIdBg parameter is greater than number of tuples !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    int trueEnd=tupleIdEnd;
    if(tupleIdEnd!=-1)
      {
        if(tupleIdEnd>nbt)
          {
            std::ostringstream oss; oss << Traits<T>::ArrayTypeName << ":subArray : The tupleIdBg parameter is greater than number of tuples !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    else
      trueEnd=nbt;
    int nbComp(getNumberOfComponents());
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    ret->alloc(trueEnd-tupleIdBg,nbComp);
    ret->copyStringInfoFrom(*this);
    std::copy(getConstPointer()+tupleIdBg*nbComp,getConstPointer()+trueEnd*nbComp,ret->getPointer());
    return ret.retn();
  }

  /*!
   * Returns a shorten copy of \a this array. The new DataArrayDouble contains every
   * (\a bg + \c i * \a step)-th tuple of \a this array located before the \a end2-th
   * tuple. Indices of the selected tuples are the same as ones returned by the Python
   * command \c range( \a bg, \a end2, \a step ).
   * This method is equivalent to selectByTupleIdSafe() except that the input array is
   * not constructed explicitly.
   * For more info on renumbering see \ref numbering.
   *  \param [in] bg - index of the first tuple to copy from \a this array.
   *  \param [in] end2 - index of the tuple before which the tuples to copy are located.
   *  \param [in] step - index increment to get index of the next tuple to copy.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   *  \sa DataArrayDouble::subArray.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::mySelectByTupleIdSafeSlice(int bg, int end2, int step) const
  {
    checkAllocated();
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    int nbComp(getNumberOfComponents());
    std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::selectByTupleIdSafeSlice : ";
    int newNbOfTuples(GetNumberOfItemGivenBESRelative(bg,end2,step,oss.str()));
    ret->alloc(newNbOfTuples,nbComp);
    T *pt(ret->getPointer());
    const T *srcPt(getConstPointer()+bg*nbComp);
    for(int i=0;i<newNbOfTuples;i++,srcPt+=step*nbComp)
      std::copy(srcPt,srcPt+nbComp,pt+i*nbComp);
    ret->copyStringInfoFrom(*this);
    return ret.retn();
  }
  
  /*!
   * Copy all values from another DataArrayDouble into specified tuples and components
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
   *  \if ENABLE_EXAMPLES
   *  \ref py_mcdataarraydouble_setpartofvalues1 "Here is a Python example".
   *  \endif
   */
  template<class T>
  void DataArrayTemplate<T>::setPartOfValues1(const typename Traits<T>::ArrayType *a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare)
  {
    if(!a)
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::setPartOfValues1 : input DataArrayDouble is NULL !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    const char msg[]="DataArrayTemplate::setPartOfValues1";
    checkAllocated();
    a->checkAllocated();
    int newNbOfTuples(DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg));
    int newNbOfComp(DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg));
    int nbComp(getNumberOfComponents()),nbOfTuples(getNumberOfTuples());
    DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
    DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
    bool assignTech(true);
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
    const T *srcPt(a->getConstPointer());
    T *pt(getPointer()+bgTuples*nbComp+bgComp);
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
            const T*srcPt2=srcPt;
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
 *            for \c this array.
 *
 *  \if ENABLE_EXAMPLES
 *  \ref py_mcdataarraydouble_setpartofvaluessimple1 "Here is a Python example".
 *  \endif
 */
  template<class T>
  void DataArrayTemplate<T>::setPartOfValuesSimple1(T a, int bgTuples, int endTuples, int stepTuples, int bgComp, int endComp, int stepComp)
  {
    const char msg[]="DataArrayTemplate::setPartOfValuesSimple1";
    checkAllocated();
    int newNbOfTuples(DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg));
    int newNbOfComp(DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg));
    int nbComp(getNumberOfComponents()),nbOfTuples(getNumberOfTuples());
    DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
    DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
    T *pt=getPointer()+bgTuples*nbComp+bgComp;
    for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
      for(int j=0;j<newNbOfComp;j++)
        pt[j*stepComp]=a;
  }
  
  /*!
   * Copy all values from another DataArrayDouble (\a a) into specified tuples and 
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
   *  \if ENABLE_EXAMPLES
   *  \ref py_mcdataarraydouble_setpartofvalues2 "Here is a Python example".
   *  \endif
   */
  template<class T>
  void DataArrayTemplate<T>::setPartOfValues2(const typename Traits<T>::ArrayType *a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp, bool strictCompoCompare)
  {
    if(!a)
      throw INTERP_KERNEL::Exception("DataArrayDouble::setPartOfValues2 : input DataArrayDouble is NULL !");
    const char msg[]="DataArrayTemplate::setPartOfValues2";
    checkAllocated();
    a->checkAllocated();
    int nbComp(getNumberOfComponents()),nbOfTuples(getNumberOfTuples());
    for(const int *z=bgComp;z!=endComp;z++)
      DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
    int newNbOfTuples((int)std::distance(bgTuples,endTuples));
    int newNbOfComp((int)std::distance(bgComp,endComp));
    bool assignTech(true);
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
    T *pt(getPointer());
    const T *srcPt(a->getConstPointer());
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
            const T *srcPt2=srcPt;
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
   *  \if ENABLE_EXAMPLES
   *  \ref py_mcdataarraydouble_setpartofvaluessimple2 "Here is a Python example".
   *  \endif
   */
  template<class T>
  void DataArrayTemplate<T>::setPartOfValuesSimple2(T a, const int *bgTuples, const int *endTuples, const int *bgComp, const int *endComp)
  {
    checkAllocated();
    int nbComp(getNumberOfComponents()),nbOfTuples(getNumberOfTuples());
    for(const int *z=bgComp;z!=endComp;z++)
      DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
    T *pt(getPointer());
    for(const int *w=bgTuples;w!=endTuples;w++)
      for(const int *z=bgComp;z!=endComp;z++)
        {
          DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
          pt[(std::size_t)(*w)*nbComp+(*z)]=a;
        }
  }
  
  /*!
   * Copy all values from another DataArrayDouble (\a a) into specified tuples and 
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
   *            for \c this array.
   *
   *  \if ENABLE_EXAMPLES
   *  \ref py_mcdataarraydouble_setpartofvalues3 "Here is a Python example".
   *  \endif
   */
  template<class T>
  void DataArrayTemplate<T>::setPartOfValues3(const typename Traits<T>::ArrayType *a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp, bool strictCompoCompare)
  {
    if(!a)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setPartOfValues3 : input DataArrayDouble is NULL !");
    const char msg[]="DataArrayTemplate::setPartOfValues3";
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
    T *pt(getPointer()+bgComp);
    const T *srcPt(a->getConstPointer());
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
            const T *srcPt2=srcPt;
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
   *            for \c this array.
   *
   *  \if ENABLE_EXAMPLES
   *  \ref py_mcdataarraydouble_setpartofvaluessimple3 "Here is a Python example".
   *  \endif
   */
  template<class T>
  void DataArrayTemplate<T>::setPartOfValuesSimple3(T a, const int *bgTuples, const int *endTuples, int bgComp, int endComp, int stepComp)
  {
    const char msg[]="DataArrayTemplate::setPartOfValuesSimple3";
    checkAllocated();
    int newNbOfComp(DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg));
    int nbComp(getNumberOfComponents()),nbOfTuples(getNumberOfTuples());
    DataArray::CheckValueInRangeEx(nbComp,bgComp,endComp,"invalid component value");
    T *pt(getPointer()+bgComp);
    for(const int *w=bgTuples;w!=endTuples;w++)
      for(int j=0;j<newNbOfComp;j++)
        {
          DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
          pt[(std::size_t)(*w)*nbComp+j*stepComp]=a;
        }
  }

  /*!
   * Copy all values from another DataArrayDouble into specified tuples and components
   * of \a this array. Textual data is not copied.
   * The tree parameters defining set of indices of tuples and components are similar to
   * the tree parameters of the Python function \c range(\c start,\c stop,\c step).
   *  \param [in] a - the array to copy values from.
   *  \param [in] bgTuples - index of the first tuple of \a this array to assign values to.
   *  \param [in] endTuples - index of the tuple before which the tuples to assign to
   *              are located.
   *  \param [in] stepTuples - index increment to get index of the next tuple to assign to.
   *  \param [in] bgComp - pointer to an array of component indices of \a this array to
   *              assign \a a to.
   *  \param [in] endComp - specifies the end of the array \a bgTuples, so that
   *              pointer to a component index (\a pi) varies as this: 
   *              \a bgComp <= \a pi < \a endComp.
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
   */
  template<class T>
  void DataArrayTemplate<T>::setPartOfValues4(const typename Traits<T>::ArrayType *a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp, bool strictCompoCompare)
  {if(!a)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setPartOfValues4 : input DataArrayTemplate is NULL !");
    const char msg[]="DataArrayTemplate::setPartOfValues4";
    checkAllocated();
    a->checkAllocated();
    int newNbOfTuples(DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg));
    int newNbOfComp((int)std::distance(bgComp,endComp));
    int nbComp(getNumberOfComponents());
    for(const int *z=bgComp;z!=endComp;z++)
      DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
    int nbOfTuples(getNumberOfTuples());
    DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
    bool assignTech(true);
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
    const T *srcPt(a->getConstPointer());
    T *pt(getPointer()+bgTuples*nbComp);
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
          const T *srcPt2(srcPt);
          for(const int *z=bgComp;z!=endComp;z++,srcPt2++)
            pt[*z]=*srcPt2;
        }
      }
  }

  template<class T>
  void DataArrayTemplate<T>::setPartOfValuesSimple4(T a, int bgTuples, int endTuples, int stepTuples, const int *bgComp, const int *endComp)
  {
    const char msg[]="DataArrayTemplate::setPartOfValuesSimple4";
    checkAllocated();
    int newNbOfTuples(DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg));
    int nbComp(getNumberOfComponents());
    for(const int *z=bgComp;z!=endComp;z++)
      DataArray::CheckValueInRange(nbComp,*z,"invalid component id");
    int nbOfTuples(getNumberOfTuples());
    DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
    T *pt=getPointer()+bgTuples*nbComp;
    for(int i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
      for(const int *z=bgComp;z!=endComp;z++)
        pt[*z]=a;
  }
  
  /*!
   * Copy some tuples from another DataArrayDouble into specified tuples
   * of \a this array. Textual data is not copied. Both arrays must have equal number of
   * components.
   * Both the tuples to assign and the tuples to assign to are defined by a DataArrayInt.
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
  template<class T>
  void DataArrayTemplate<T>::setPartOfValuesAdv(const typename Traits<T>::ArrayType *a, const DataArrayInt *tuplesSelec)
  {
    if(!a || !tuplesSelec)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setPartOfValuesAdv : input DataArrayTemplate is NULL !");
    checkAllocated();
    a->checkAllocated();
    tuplesSelec->checkAllocated();
    int nbOfComp=getNumberOfComponents();
    if(nbOfComp!=a->getNumberOfComponents())
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setPartOfValuesAdv : This and a do not have the same number of components !");
    if(tuplesSelec->getNumberOfComponents()!=2)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setPartOfValuesAdv : Expecting to have a tuple selector DataArrayInt instance with exactly 2 components !");
    int thisNt(getNumberOfTuples());
    int aNt(a->getNumberOfTuples());
    T *valsToSet(getPointer());
    const T *valsSrc(a->getConstPointer());
    for(const int *tuple=tuplesSelec->begin();tuple!=tuplesSelec->end();tuple+=2)
    {
      if(tuple[1]>=0 && tuple[1]<aNt)
        {
          if(tuple[0]>=0 && tuple[0]<thisNt)
            std::copy(valsSrc+nbOfComp*tuple[1],valsSrc+nbOfComp*(tuple[1]+1),valsToSet+nbOfComp*tuple[0]);
          else
            {
              std::ostringstream oss; oss << "DataArrayTemplate::setPartOfValuesAdv : Tuple #" << std::distance(tuplesSelec->begin(),tuple)/2;
              oss << " of 'tuplesSelec' request of tuple id #" << tuple[0] << " in 'this' ! It should be in [0," << thisNt << ") !";
              throw INTERP_KERNEL::Exception(oss.str().c_str());
            }
        }
      else
        {
          std::ostringstream oss; oss << "DataArrayTemplate::setPartOfValuesAdv : Tuple #" << std::distance(tuplesSelec->begin(),tuple)/2;
          oss << " of 'tuplesSelec' request of tuple id #" << tuple[1] << " in 'a' ! It should be in [0," << aNt << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    }
  }
  
  /*!
   * Copy some tuples from another DataArrayDouble (\a aBase) into contiguous tuples
   * of \a this array. Textual data is not copied. Both arrays must have equal number of
   * components.
   * The tuples to assign to are defined by index of the first tuple, and
   * their number is defined by \a tuplesSelec->getNumberOfTuples().
   * The tuples to copy are defined by values of a DataArrayInt.
   * All components of selected tuples are copied.
   *  \param [in] tupleIdStart - index of the first tuple of \a this array to assign
   *              values to.
   *  \param [in] aBase - the array to copy values from.
   *  \param [in] tuplesSelec - the array specifying tuples of \a a to copy.
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
  template<class T>
  void DataArrayTemplate<T>::setContigPartOfSelectedValues(int tupleIdStart, const DataArray *aBase, const DataArrayInt *tuplesSelec)
  {
    if(!aBase || !tuplesSelec)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValues : input DataArray is NULL !");
    const typename Traits<T>::ArrayType *a(dynamic_cast<const typename Traits<T>::ArrayType *>(aBase));
    if(!a)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValues : input DataArray aBase is not a DataArrayDouble !");
    checkAllocated();
    a->checkAllocated();
    tuplesSelec->checkAllocated();
    int nbOfComp(getNumberOfComponents());
    if(nbOfComp!=a->getNumberOfComponents())
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValues : This and a do not have the same number of components !");
    if(tuplesSelec->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValues : Expecting to have a tuple selector DataArrayInt instance with exactly 1 component !");
    int thisNt(getNumberOfTuples());
    int aNt(a->getNumberOfTuples());
    int nbOfTupleToWrite(tuplesSelec->getNumberOfTuples());
    T *valsToSet(getPointer()+tupleIdStart*nbOfComp);
    if(tupleIdStart+nbOfTupleToWrite>thisNt)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValues : invalid number range of values to write !");
    const T *valsSrc=a->getConstPointer();
    for(const int *tuple=tuplesSelec->begin();tuple!=tuplesSelec->end();tuple++,valsToSet+=nbOfComp)
      {
        if(*tuple>=0 && *tuple<aNt)
          {
            std::copy(valsSrc+nbOfComp*(*tuple),valsSrc+nbOfComp*(*tuple+1),valsToSet);
          }
        else
          {
            std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::setContigPartOfSelectedValues : Tuple #" << std::distance(tuplesSelec->begin(),tuple);
            oss << " of 'tuplesSelec' request of tuple id #" << *tuple << " in 'a' ! It should be in [0," << aNt << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
  }
  
  /*!
   * Copy some tuples from another DataArrayDouble (\a aBase) into contiguous tuples
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
  template<class T>
  void DataArrayTemplate<T>::setContigPartOfSelectedValuesSlice(int tupleIdStart, const DataArray *aBase, int bg, int end2, int step)
  {
    if(!aBase)
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::setContigPartOfSelectedValuesSlice : input DataArray is NULL !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    const typename Traits<T>::ArrayType *a(dynamic_cast<const typename Traits<T>::ArrayType *>(aBase));
    if(!a)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValuesSlice : input DataArray aBase is not a DataArrayDouble !");
    checkAllocated();
    a->checkAllocated();
    int nbOfComp(getNumberOfComponents());
    const char msg[]="DataArrayDouble::setContigPartOfSelectedValuesSlice";
    int nbOfTupleToWrite(DataArray::GetNumberOfItemGivenBES(bg,end2,step,msg));
    if(nbOfComp!=a->getNumberOfComponents())
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValuesSlice : This and a do not have the same number of components !");
    int thisNt(getNumberOfTuples()),aNt(a->getNumberOfTuples());
    T *valsToSet(getPointer()+tupleIdStart*nbOfComp);
    if(tupleIdStart+nbOfTupleToWrite>thisNt)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValuesSlice : invalid number range of values to write !");
    if(end2>aNt)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValuesSlice : invalid range of values to read !");
    const T *valsSrc(a->getConstPointer()+bg*nbOfComp);
    for(int i=0;i<nbOfTupleToWrite;i++,valsToSet+=nbOfComp,valsSrc+=step*nbOfComp)
      {
        std::copy(valsSrc,valsSrc+nbOfComp,valsToSet);
      }
  }

  /*!
   * Returns a shorten copy of \a this array. The new DataArrayDouble contains ranges
   * of tuples specified by \a ranges parameter.
   * For more info on renumbering see \ref numbering.
   *  \param [in] ranges - std::vector of std::pair's each of which defines a range
   *              of tuples in [\c begin,\c end) format.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   *  \throw If \a end < \a begin.
   *  \throw If \a end > \a this->getNumberOfTuples().
   *  \throw If \a this is not allocated.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::mySelectByTupleRanges(const std::vector<std::pair<int,int> >& ranges) const
  {
    checkAllocated();
    int nbOfComp(getNumberOfComponents()),nbOfTuplesThis(getNumberOfTuples());
    if(ranges.empty())
      {
        MCAuto<DataArray> ret0(buildNewEmptyInstance());
        MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
        ret->alloc(0,nbOfComp);
        ret->copyStringInfoFrom(*this);
        return ret.retn();
      }
    int ref(ranges.front().first),nbOfTuples(0);
    bool isIncreasing(true);
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
                std::ostringstream oss; oss << "DataArrayTemplate::selectByTupleRanges : on range #" << std::distance(ranges.begin(),it);
                oss << " (" << (*it).first << "," << (*it).second << ") is greater than number of tuples of this :" << nbOfTuples << " !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayTemplate::selectByTupleRanges : on range #" << std::distance(ranges.begin(),it);
            oss << " (" << (*it).first << "," << (*it).second << ") end is before begin !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    if(isIncreasing && nbOfTuplesThis==nbOfTuples)
      return static_cast<typename Traits<T>::ArrayType *>(deepCopy());
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    ret->alloc(nbOfTuples,nbOfComp);
    ret->copyStringInfoFrom(*this);
    const T *src(getConstPointer());
    T *work(ret->getPointer());
    for(std::vector<std::pair<int,int> >::const_iterator it=ranges.begin();it!=ranges.end();it++)
      work=std::copy(src+(*it).first*nbOfComp,src+(*it).second*nbOfComp,work);
    return ret.retn();
  }

  /*!
   * Returns the first value of \a this. 
   *  \return double - the last value of \a this array.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If \a this->getNumberOfTuples() < 1.
   */
  template<class T>
  T DataArrayTemplate<T>::front() const
  {
    checkAllocated();
    if(getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::front : number of components not equal to one !");
    int nbOfTuples(getNumberOfTuples());
    if(nbOfTuples<1)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::front : number of tuples must be >= 1 !");
    return *(getConstPointer());
  }
  
  /*!
   * Returns the last value of \a this. 
   *  \return double - the last value of \a this array.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If \a this->getNumberOfTuples() < 1.
   */
  template<class T>
  T DataArrayTemplate<T>::back() const
  {
    checkAllocated();
    if(getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::back : number of components not equal to one !");
    int nbOfTuples(getNumberOfTuples());
    if(nbOfTuples<1)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::back : number of tuples must be >= 1 !");
    return *(getConstPointer()+nbOfTuples-1);
  }
  
  /*!
   * Returns the maximal value and its location within \a this one-dimensional array.
   *  \param [out] tupleId - index of the tuple holding the maximal value.
   *  \return double - the maximal value among all values of \a this array.
   *  \throw If \a this->getNumberOfComponents() != 1
   *  \throw If \a this->getNumberOfTuples() < 1
   */
  template<class T>
  T DataArrayTemplate<T>::getMaxValue(int& tupleId) const
  {
    checkAllocated();
    if(getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxValue : must be applied on DataArrayDouble with only one component, you can call 'rearrange' method before or call 'getMaxValueInArray' method !");
    int nbOfTuples(getNumberOfTuples());
    if(nbOfTuples<=0)
      throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxValue : array exists but number of tuples must be > 0 !");
    const T *vals(getConstPointer());
    const T *loc(std::max_element(vals,vals+nbOfTuples));
    tupleId=(int)std::distance(vals,loc);
    return *loc;
  }
  
  /*!
   * Returns the maximal value within \a this array that is allowed to have more than
   *  one component.
   *  \return double - the maximal value among all values of \a this array.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  T DataArrayTemplate<T>::getMaxValueInArray() const
  {
    checkAllocated();
    const T *loc(std::max_element(begin(),end()));
    return *loc;
  }
  
  /*!
   * Returns the minimal value and its location within \a this one-dimensional array.
   *  \param [out] tupleId - index of the tuple holding the minimal value.
   *  \return double - the minimal value among all values of \a this array.
   *  \throw If \a this->getNumberOfComponents() != 1
   *  \throw If \a this->getNumberOfTuples() < 1
   */
  template<class T>
  T DataArrayTemplate<T>::getMinValue(int& tupleId) const
  {
    checkAllocated();
    if(getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayDouble::getMinValue : must be applied on DataArrayDouble with only one component, you can call 'rearrange' method before call 'getMinValueInArray' method !");
    int nbOfTuples(getNumberOfTuples());
    if(nbOfTuples<=0)
      throw INTERP_KERNEL::Exception("DataArrayDouble::getMinValue : array exists but number of tuples must be > 0 !");
    const T *vals(getConstPointer());
    const T *loc(std::min_element(vals,vals+nbOfTuples));
    tupleId=(int)std::distance(vals,loc);
    return *loc;
  }
  
  /*!
   * Returns the minimal value within \a this array that is allowed to have more than
   *  one component.
   *  \return double - the minimal value among all values of \a this array.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  T DataArrayTemplate<T>::getMinValueInArray() const
  {
    checkAllocated();
    const T *loc=std::min_element(begin(),end());
    return *loc;
  }
  
}

#endif
