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
}

#endif
