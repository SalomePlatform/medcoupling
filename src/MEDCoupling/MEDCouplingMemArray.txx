// Copyright (C) 2007-2025  CEA, EDF
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
// Author : Anthony Geay (EDF R&D)

#ifndef __PARAMEDMEM_MEDCOUPLINGMEMARRAY_TXX__
#define __PARAMEDMEM_MEDCOUPLINGMEMARRAY_TXX__

#include "MEDCouplingMemArray.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"
#include "InterpolationUtils.hxx"
#include "MEDCouplingPartDefinition.hxx"
#include "InterpKernelAutoPtr.hxx"
#include "MCAuto.hxx"
#include "MEDCouplingMap.txx"
#include "BBTreeDiscrete.txx"

#include <set>
#include <sstream>
#include <cstdlib>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <fstream>

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
        useArray(pointer,true,DeallocType::C_DEALLOC,other._nb_of_elem);
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
  bool MemArray<T>::reprHeader(mcIdType sl, std::ostream& stream) const
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
  void MemArray<T>::repr(mcIdType sl, std::ostream& stream) const
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
  void MemArray<T>::reprZip(mcIdType sl, std::ostream& stream) const
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
  void MemArray<T>::reprNotTooLong(mcIdType sl, std::ostream& stream) const
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
  T *MemArray<T>::fromNoInterlace(std::size_t nbOfComp) const
  {
    if(nbOfComp<1)
      throw INTERP_KERNEL::Exception("MemArray<T>::fromNoInterlace : number of components must be > 0 !");
    const T *pt=_pointer.getConstPointer();
    std::size_t nbOfTuples=_nb_of_elem/nbOfComp;
    T *ret=(T*)malloc(_nb_of_elem*sizeof(T));
    T *w=ret;
    for(std::size_t i=0;i<nbOfTuples;i++)
      for(std::size_t j=0;j<nbOfComp;j++,w++)
        *w=pt[j*nbOfTuples+i];
    return ret;
  }

  template<class T>
  T *MemArray<T>::toNoInterlace(std::size_t nbOfComp) const
  {
    if(nbOfComp<1)
      throw INTERP_KERNEL::Exception("MemArray<T>::toNoInterlace : number of components must be > 0 !");
    const T *pt=_pointer.getConstPointer();
    std::size_t nbOfTuples=_nb_of_elem/nbOfComp;
    T *ret=(T*)malloc(_nb_of_elem*sizeof(T));
    T *w=ret;
    for(std::size_t i=0;i<nbOfComp;i++)
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
  void MemArray<T>::reverse(std::size_t nbOfComp)
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
            for(std::size_t j=0;j<nbOfComp;j++)
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
   * The remaining part of the new allocated chunk are available but not set previously !
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
  void MemArray<T>::CPPDeallocator(void *pt, void * /*param*/)
  {
    delete [] reinterpret_cast<T*>(pt);
  }

  template<class T>
  void MemArray<T>::CDeallocator(void *pt, void * /*param*/)
  {
    free(pt);
  }

  template<class T>
  void MemArray<T>::COffsetDeallocator(void *pt, void *param)
  {
    int64_t *offset(reinterpret_cast<int64_t *>(param));
    char *ptcast(reinterpret_cast<char *>(pt));
    free(ptcast+*offset);
  }

  template<class T>
  typename MemArray<T>::Deallocator MemArray<T>::BuildFromType(DeallocType type)
  {
    switch(type)
    {
      case DeallocType::CPP_DEALLOC:
        return CPPDeallocator;
      case DeallocType::C_DEALLOC:
        return CDeallocator;
      case DeallocType::C_DEALLOC_WITH_OFFSET:
        return COffsetDeallocator;
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
  DataArrayIterator<T>::DataArrayIterator(typename Traits<T>::ArrayType *da):_da(da),_tuple_id(0),_nb_comp(0),_nb_tuple(0)
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

  template<class T>
  DataArrayIterator<T>::~DataArrayIterator()
  {
    if(_da)
      _da->decrRef();
  }

  template<class T>
  typename Traits<T>::ArrayTuple *DataArrayIterator<T>::nextt()
  {
    if(_tuple_id<_nb_tuple)
      {
        _tuple_id++;
        typename Traits<T>::ArrayTuple *ret=new typename Traits<T>::ArrayTuple(_pt,_nb_comp);
        _pt+=_nb_comp;
        return ret;
      }
    else
      return 0;
  }

  //////////////////////////////////

  template<class T>
  DataArrayTuple<T>::DataArrayTuple(T *pt, std::size_t nbOfComp):_pt(pt),_nb_of_compo(nbOfComp)
  {
  }

  template<class T>
  T DataArrayTuple<T>::zeValue() const
  {
    if(_nb_of_compo==1)
      return *_pt;
    throw INTERP_KERNEL::Exception("DataArrayTuple<T>::zeValue : DataArrayTuple instance has not exactly 1 component -> Not possible to convert it into a single value !");
  }

  template<class T>
  typename Traits<T>::ArrayType *DataArrayTuple<T>::buildDA(std::size_t nbOfTuples, std::size_t nbOfCompo) const
  {
    if((_nb_of_compo==nbOfCompo && nbOfTuples==1) || (_nb_of_compo==nbOfTuples && nbOfCompo==1))
    {
      typename Traits<T>::ArrayType *ret=Traits<T>::ArrayType::New();
      ret->useExternalArrayWithRWAccess(_pt,nbOfTuples,nbOfCompo);
      return ret;
    }
  else
    {
      std::ostringstream oss; oss << "DataArrayTuple<T>::buildDA : unable to build a requested DataArrayDouble instance with nbofTuple=" << nbOfTuples << " and nbOfCompo=" << nbOfCompo;
      oss << ".\nBecause the number of elements in this is " << _nb_of_compo << " !";
      throw INTERP_KERNEL::Exception(oss.str().c_str());
    }
  }

  //////////////////////////////////

  /*!
   * This method is useful to slice work among a pool of threads or processes. \a begin, \a end \a step is the input whole slice of work to perform,
   * typically it is a whole slice of tuples of DataArray or cells, nodes of a mesh...
   *
   * The input \a sliceId should be an id in [0, \a nbOfSlices) that specifies the slice of work.
   *
   * \param [in] start - the start of the input slice of the whole work to perform split into slices.
   * \param [in] stop - the stop of the input slice of the whole work to perform split into slices.
   * \param [in] step - the step (that can be <0) of the input slice of the whole work to perform split into slices.
   * \param [in] sliceId - the slice id considered
   * \param [in] nbOfSlices - the number of slices (typically the number of cores on which the work is expected to be sliced)
   * \param [out] startSlice - the start of the slice considered
   * \param [out] stopSlice - the stop of the slice consided
   *
   * \throw If \a step == 0
   * \throw If \a nbOfSlices not > 0
   * \throw If \a sliceId not in [0,nbOfSlices)
   */
  template <class T>
  void DataArrayTools<T>::GetSlice(T start, T stop, T step, mcIdType sliceId, mcIdType nbOfSlices, T& startSlice, T& stopSlice)
  {
    if(nbOfSlices<=0)
      {
        std::ostringstream oss; oss << "DataArray::GetSlice : nbOfSlices (" << nbOfSlices << ") must be > 0 !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    if(sliceId<0 || sliceId>=nbOfSlices)
      {
        std::ostringstream oss; oss << "DataArray::GetSlice : sliceId (" << nbOfSlices << ") must be in [0 , nbOfSlices (" << nbOfSlices << ") ) !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    mcIdType nbElems=DataArrayTools<T>::GetNumberOfItemGivenBESRelative(start,stop,step,"DataArray::GetSlice");
    mcIdType minNbOfElemsPerSlice=nbElems/nbOfSlices;
    startSlice=start+FromIdType<T>(minNbOfElemsPerSlice)*step*FromIdType<T>(sliceId);
    if(sliceId<nbOfSlices-1)
      stopSlice=start+FromIdType<T>(minNbOfElemsPerSlice)*step*(FromIdType<T>(sliceId+1));
    else
      stopSlice=stop;
  }

  template <class T>
  mcIdType DataArrayTools<T>::GetNumberOfItemGivenBES(T begin, T end, T step, const std::string& msg)
  {
    if(end<begin)
      {
        std::ostringstream oss; oss << msg << " : end before begin !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    if(end==begin)
      return 0;
    if(step<=0)
      {
        std::ostringstream oss; oss << msg << " : invalid step should be > 0 !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    return ToIdType((end-1-begin)/step+1);
  }

  template <class T>
  mcIdType DataArrayTools<T>::GetNumberOfItemGivenBESRelative(T begin, T end, T step, const std::string& msg)
  {
    if(step==0)
      throw INTERP_KERNEL::Exception("DataArray::GetNumberOfItemGivenBES : step=0 is not allowed !");
    if(end<begin && step>0)
      {
        std::ostringstream oss; oss << msg << " : end before begin whereas step is positive !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    if(begin<end && step<0)
      {
        std::ostringstream oss; oss << msg << " : invalid step should be > 0 !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    if(begin!=end)
      return ToIdType((std::max(begin,end)-1-std::min(begin,end))/std::abs(step)+1);
    else
      return 0;
  }

  template <class T>
  mcIdType DataArrayTools<T>::GetPosOfItemGivenBESRelativeNoThrow(T value, T begin, T end, T step)
  {
    if (step == 0)
      return -1;

    if((step>0 && begin<=value && value<end) ||
       (step<0 && begin>=value && value>end))
    {
      mcIdType id = ToIdType((value-begin)/step);
      if (begin + step * id == value)
        return id;
      else
        return -1;
    }
    else
      return -1;
  }

  //////////////////////////////////

  template<class T>
  MCAuto< typename Traits<T>::ArrayTypeCh > DataArrayTemplate<T>::NewFromStdVector(const typename std::vector<T>& v)
  {
    std::size_t sz(v.size());
    MCAuto< typename Traits<T>::ArrayTypeCh > ret(Traits<T>::ArrayTypeCh::New());
    ret->alloc(sz,1);
    T *pt(ret->getPointer());
    std::copy(v.begin(),v.end(),pt);
    return ret;
  }

  /*!
   * Returns a newly created array containing a copy of the input array defined by [ \a arrBegin, \a arrEnd )
   */
  template<class T>
  MCAuto< typename Traits<T>::ArrayTypeCh > DataArrayTemplate<T>::NewFromArray(const T *arrBegin, const T *arrEnd)
  {
    using DataArrayT = typename Traits<T>::ArrayTypeCh;
    MCAuto< DataArrayT > ret(DataArrayT::New());
    std::size_t nbElts(std::distance(arrBegin,arrEnd));
    ret->alloc(nbElts,1);
    std::copy(arrBegin,arrEnd,ret->getPointer());
    return ret;
  }

  template<class T>
  std::vector< MCAuto< typename Traits<T>::ArrayTypeCh > > DataArrayTemplate<T>::explodeComponents() const
  {
    checkAllocated();
    std::size_t sz(getNumberOfComponents());
    mcIdType nbTuples(getNumberOfTuples());
    std::string name(getName());
    std::vector<std::string> compNames(getInfoOnComponents());
    std::vector< MCAuto< typename Traits<T>::ArrayTypeCh > > ret(sz);
    const T *thisPt(begin());
    for(std::size_t i=0;i<sz;i++)
      {
        MCAuto< typename Traits<T>::ArrayTypeCh > part(Traits<T>::ArrayTypeCh::New());
        part->alloc(nbTuples,1);
        part->setName(name);
        part->setInfoOnComponent(0,compNames[i]);
        T *otherPt(part->getPointer());
        for(mcIdType j=0;j<nbTuples;j++)
          otherPt[j]=thisPt[sz*j+i];
        ret[i]=part;
      }
    return ret;
  }

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
  void DataArrayTemplate<T>::alloc(std::size_t nbOfTuple, std::size_t nbOfCompo)
  {
    _info_on_compo.resize(nbOfCompo);
    _mem.alloc(nbOfCompo*nbOfTuple);
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
  void DataArrayTemplate<T>::useArray(const T *array, bool ownership, DeallocType type, std::size_t nbOfTuple, std::size_t nbOfCompo)
  {
    _info_on_compo.resize(nbOfCompo);
    _mem.useArray(array,ownership,type,nbOfTuple*nbOfCompo);
    declareAsNew();
  }

  template<class T>
  void DataArrayTemplate<T>::useExternalArrayWithRWAccess(const T *array, std::size_t nbOfTuple, std::size_t nbOfCompo)
  {
    _info_on_compo.resize(nbOfCompo);
    _mem.useExternalArrayWithRWAccess(array,nbOfTuple*nbOfCompo);
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
  T DataArrayTemplate<T>::getIJSafe(std::size_t tupleId, std::size_t compoId) const
  {
    checkAllocated();
    if(ToIdType(tupleId)>=getNumberOfTuples())
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::getIJSafe : request for tupleId " << tupleId << " should be in [0," << getNumberOfTuples() << ") !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    if(compoId>=getNumberOfComponents())
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
   * This method deallocated \a this without modification of information relative to the components.
   * After call of this method, DataArrayDouble::isAllocated will return false.
   * If \a this is already not allocated, \a this is let unchanged.
   */
  template<class T>
  void DataArrayTemplate<T>::desallocate()
  {
    _mem.destroy();
  }

  /*!
   * This method reserve nbOfElems elements in memory ( nbOfElems*sizeof(T) ) \b without impacting the number of tuples in \a this.
   * If \a this has already been allocated, this method checks that \a this has only one component. If not an INTERP_KERNEL::Exception will be thrown.
   * If \a this has not already been allocated, number of components is set to one.
   * This method allows to reduce number of reallocations on invocation of DataArrayDouble::pushBackSilent and DataArrayDouble::pushBackValsSilent on \a this.
   *
   * \sa DataArrayDouble::pack, DataArrayDouble::pushBackSilent, DataArrayDouble::pushBackValsSilent
   */
  template<class T>
  void DataArrayTemplate<T>::reserve(std::size_t nbOfElems)
  {
    std::size_t nbCompo(getNumberOfComponents());
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
    std::size_t nbCompo(getNumberOfComponents());
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
  void DataArrayTemplate<T>::allocIfNecessary(std::size_t nbOfTuple, std::size_t nbOfCompo)
  {
    if(isAllocated())
      {
        if(ToIdType(nbOfTuple)!=getNumberOfTuples() || nbOfCompo!=getNumberOfComponents())
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
    mcIdType nbOfTuples(other.getNumberOfTuples());
    std::size_t nbOfComp(other.getNumberOfComponents());
    allocIfNecessary(nbOfTuples,nbOfComp);
    std::size_t nbOfElems(nbOfTuples*nbOfComp);
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
  void DataArrayTemplate<T>::reAlloc(std::size_t nbOfTuples)
  {
    checkAllocated();
    _mem.reAlloc(getNumberOfComponents()*nbOfTuples);
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
  void DataArrayTemplate<T>::renumberInPlace(const mcIdType *old2New)
  {
    checkAllocated();
    mcIdType nbTuples(getNumberOfTuples());
    std::size_t nbOfCompo(getNumberOfComponents());
    T *tmp(new T[nbTuples*nbOfCompo]);
    const T *iptr(begin());
    for(mcIdType i=0;i<nbTuples;i++)
      {
        mcIdType v=old2New[i];
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
   */
  template<class T>
  void DataArrayTemplate<T>::renumberInPlaceR(const mcIdType *new2Old)
  {
    checkAllocated();
    mcIdType nbTuples(getNumberOfTuples());
    std::size_t nbOfCompo(getNumberOfComponents());
    T *tmp(new T[nbTuples*nbOfCompo]);
    const T *iptr(begin());
    for(mcIdType i=0;i<nbTuples;i++)
      {
        mcIdType v=new2Old[i];
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
   * Sorts values of the array. \b Warning, this method is not const, it alterates \a this content.
   *
   *  \param [in] asc - \a true means ascending order, \a false, descending.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \sa copySorted
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
    * Sorts value within every tuple of \a this array.
    *  \param [in] asc - if \a true, the values are sorted in ascending order, else,
    *              in descending order.
    *  \throw If \a this is not allocated.
    */
    template<class T>
    void DataArrayTemplate<T>::sortPerTuple(bool asc)
    {
      this->checkAllocated();
      T *pt( this->getPointer() );
      mcIdType nbOfTuple(this->getNumberOfTuples());
      std::size_t nbOfComp(this->getNumberOfComponents());
      if(asc)
        for(mcIdType i=0;i<nbOfTuple;i++,pt+=nbOfComp)
          std::sort(pt,pt+nbOfComp);
      else
        for(mcIdType i=0;i<nbOfTuple;i++,pt+=nbOfComp)
          std::sort(pt,pt+nbOfComp,std::greater<double>());
      this->declareAsNew();
    }

  /*!
   * Sorts values of the array and put the result in a newly allocated returned array.
   * This method does not alterate \a this content.
   *
   *  \param [in] asc - \a true means ascending order, \a false, descending.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \sa sort
   */
  template<class T>
  typename Traits<T>::ArrayTypeCh *DataArrayTemplate<T>::copySortedImpl(bool asc) const
  {
    MCAuto<typename Traits<T>::ArrayTypeCh> ret(static_cast<typename Traits<T>::ArrayTypeCh *>(this->deepCopy()));
    ret->sort(asc);
    return ret.retn();
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
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::renumber(const mcIdType *old2New) const
  {
    checkAllocated();
    mcIdType nbTuples(getNumberOfTuples());
    std::size_t nbOfCompo(getNumberOfComponents());
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    ret->alloc(nbTuples,nbOfCompo);
    ret->copyStringInfoFrom(*this);
    const T *iptr(begin());
    T *optr(ret->getPointer());
    for(mcIdType i=0;i<nbTuples;i++)
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
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::renumberR(const mcIdType *new2Old) const
  {
    checkAllocated();
    mcIdType nbTuples(getNumberOfTuples());
    std::size_t nbOfCompo(getNumberOfComponents());
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    ret->alloc(nbTuples,nbOfCompo);
    ret->copyStringInfoFrom(*this);
    const T *iptr(getConstPointer());
    T *optr(ret->getPointer());
    for(mcIdType i=0;i<nbTuples;i++)
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
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::renumberAndReduce(const mcIdType *old2New, mcIdType newNbOfTuple) const
  {
    checkAllocated();
    mcIdType nbTuples(getNumberOfTuples());
    std::size_t nbOfCompo(getNumberOfComponents());
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    ret->alloc(newNbOfTuple,nbOfCompo);
    const T *iptr=getConstPointer();
    T *optr=ret->getPointer();
    for(mcIdType i=0;i<nbTuples;i++)
      {
        mcIdType w=old2New[i];
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
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::mySelectByTupleId(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const
  {
    checkAllocated();
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    std::size_t nbComp(getNumberOfComponents());
    ret->alloc(std::distance(new2OldBg,new2OldEnd),nbComp);
    ret->copyStringInfoFrom(*this);
    T *pt(ret->getPointer());
    const T *srcPt(getConstPointer());
    std::size_t i(0);
    for(const mcIdType *w=new2OldBg;w!=new2OldEnd;w++,i++)
      std::copy(srcPt+(*w)*nbComp,srcPt+((*w)+1)*nbComp,pt+i*nbComp);
    ret->copyStringInfoFrom(*this);
    return ret.retn();
  }

  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::mySelectByTupleId(const DataArrayIdType& di) const
  {
    return this->mySelectByTupleId(di.begin(),di.end());
  }

  template<class T>
  MCAuto<typename Traits<T>::ArrayTypeCh> DataArrayTemplate<T>::selectPartDef(const PartDefinition *pd) const
  {
    if(!pd)
      throw INTERP_KERNEL::Exception("DataArrayTemplate<T>::selectPartDef : null input pointer !");
    MCAuto<typename Traits<T>::ArrayTypeCh> ret(Traits<T>::ArrayTypeCh::New());
    const SlicePartDefinition *spd(dynamic_cast<const SlicePartDefinition *>(pd));
    if(spd)
      {
        mcIdType a,b,c;
        spd->getSlice(a,b,c);
        if(a==0 && b==getNumberOfTuples() && c==1)
          {
            DataArrayTemplate<T> *directRet(const_cast<DataArrayTemplate<T> *>(this));
            directRet->incrRef();
            MCAuto<DataArrayTemplate<T> > ret2(directRet);
            return DynamicCastSafe<DataArrayTemplate<T>,typename Traits<T>::ArrayTypeCh>(ret2);
          }
        else
          {
            MCAuto<DataArray> ret2(selectByTupleIdSafeSlice(a,b,c));
            return DynamicCastSafe<DataArray,typename Traits<T>::ArrayTypeCh>(ret2);
          }
      }
    const DataArrayPartDefinition *dpd(dynamic_cast<const DataArrayPartDefinition *>(pd));
    if(dpd)
      {
        MCAuto<DataArrayIdType> arr(dpd->toDAI());
        MCAuto<DataArray> ret2(selectByTupleIdSafe(arr->begin(),arr->end()));
        return DynamicCastSafe<DataArray,typename Traits<T>::ArrayTypeCh>(ret2);

      }
    throw INTERP_KERNEL::Exception("DataArrayTemplate<T>::selectPartDef : unrecognized part def !");
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
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::mySelectByTupleIdSafe(const mcIdType *new2OldBg, const mcIdType *new2OldEnd) const
  {
    checkAllocated();
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    std::size_t nbComp(getNumberOfComponents());
    mcIdType oldNbOfTuples(getNumberOfTuples());
    ret->alloc(std::distance(new2OldBg,new2OldEnd),nbComp);
    ret->copyStringInfoFrom(*this);
    T *pt(ret->getPointer());
    const T *srcPt(getConstPointer());
    mcIdType i(0);
    for(const mcIdType *w=new2OldBg;w!=new2OldEnd;w++,i++)
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
   *  \param [in] newNbOfCompo - number of components for \a this array to have.
   *  \throw If \a this is not allocated
   *  \throw If getNbOfElems() % \a newNbOfCompo != 0.
   *  \throw If \a newNbOfCompo is lower than 1.
   *  \throw If the rearrange method would lead to a number of tuples higher than 2147483647 (maximal capacity of int32 !).
   *  \warning This method erases all (name and unit) component info set before!
   */
  template<class T>
  void DataArrayTemplate<T>::rearrange(std::size_t newNbOfCompo)
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
    if(nbOfElems/newNbOfCompo>(std::size_t)std::numeric_limits<mcIdType>::max())
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
    rearrange(getNumberOfTuples());
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
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::changeNbOfComponents(std::size_t newNbOfComp, T dftValue) const
  {
    checkAllocated();
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    ret->alloc(getNumberOfTuples(),newNbOfComp);
    const T *oldc(getConstPointer());
    T *nc(ret->getPointer());
    mcIdType nbOfTuples=getNumberOfTuples();
    std::size_t oldNbOfComp=getNumberOfComponents();
    std::size_t dim(std::min(oldNbOfComp,newNbOfComp));
    for(mcIdType i=0;i<nbOfTuples;i++)
      {
        std::size_t j=0;
        for(;j<dim;j++)
          nc[newNbOfComp*i+j]=oldc[i*oldNbOfComp+j];
        for(;j<newNbOfComp;j++)
          nc[newNbOfComp*i+j]=dftValue;
      }
    ret->setName(getName());
    for(std::size_t i=0;i<dim;i++)
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
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::myKeepSelectedComponents(const std::vector<std::size_t>& compoIds) const
  {
    checkAllocated();
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    std::size_t newNbOfCompo=compoIds.size();
    std::size_t oldNbOfCompo=getNumberOfComponents();
    for(std::vector<std::size_t>::const_iterator it=compoIds.begin();it!=compoIds.end();it++)
      if((*it)>=oldNbOfCompo)  // (*it) >= 0 (it is a size_t)
        {
          std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::keepSelectedComponents : invalid requested component : " << *it << " whereas it should be in [0," << oldNbOfCompo << ") !";
          throw INTERP_KERNEL::Exception(oss.str().c_str());
        }
    mcIdType nbOfTuples(getNumberOfTuples());
    ret->alloc(nbOfTuples,newNbOfCompo);
    ret->copyPartOfStringInfoFrom(*this,compoIds);
    const T *oldc(getConstPointer());
    T *nc(ret->getPointer());
    for(mcIdType i=0;i<nbOfTuples;i++)
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
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::subArray(mcIdType tupleIdBg, mcIdType tupleIdEnd) const
  {
    checkAllocated();
    mcIdType nbt=getNumberOfTuples();
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
    mcIdType trueEnd=tupleIdEnd;
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
    std::size_t nbComp=getNumberOfComponents();
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
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::mySelectByTupleIdSafeSlice(mcIdType bg, mcIdType end2, mcIdType step) const
  {
    checkAllocated();
    MCAuto<DataArray> ret0(buildNewEmptyInstance());
    MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
    std::size_t nbComp(getNumberOfComponents());
    std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::selectByTupleIdSafeSlice : ";
    mcIdType newNbOfTuples(GetNumberOfItemGivenBESRelative(bg,end2,step,oss.str()));
    ret->alloc(newNbOfTuples,nbComp);
    T *pt(ret->getPointer());
    const T *srcPt(getConstPointer()+bg*nbComp);
    for(mcIdType i=0;i<newNbOfTuples;i++,srcPt+=step*nbComp)
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
  void DataArrayTemplate<T>::setPartOfValues1(const typename Traits<T>::ArrayType *a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp, bool strictCompoCompare)
  {
    if(!a)
      {
        std::ostringstream oss; oss << Traits<T>::ArrayTypeName << "::setPartOfValues1 : input DataArrayDouble is NULL !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    const char msg[]="DataArrayTemplate::setPartOfValues1";
    checkAllocated();
    a->checkAllocated();
    mcIdType newNbOfTuples(DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg));
    mcIdType newNbOfComp(DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg));
    std::size_t nbComp(getNumberOfComponents());
    mcIdType nbOfTuples(getNumberOfTuples());
    DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
    DataArray::CheckValueInRangeEx(ToIdType(nbComp),bgComp,endComp,"invalid component value");
    bool assignTech(true);
    if(a->getNbOfElems()==newNbOfTuples*newNbOfComp)
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
        for(mcIdType i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
          for(mcIdType j=0;j<newNbOfComp;j++,srcPt++)
            pt[j*stepComp]=*srcPt;
      }
    else
      {
        for(mcIdType i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
          {
            const T*srcPt2=srcPt;
            for(mcIdType j=0;j<newNbOfComp;j++,srcPt2++)
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
  void DataArrayTemplate<T>::setPartOfValuesSimple1(T a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp)
  {
    const char msg[]="DataArrayTemplate::setPartOfValuesSimple1";
    checkAllocated();
    mcIdType newNbOfTuples(DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg));
    mcIdType newNbOfComp(DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg));
    std::size_t nbComp(getNumberOfComponents());
    mcIdType nbOfTuples(getNumberOfTuples());
    DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
    DataArray::CheckValueInRangeEx(ToIdType(nbComp),bgComp,endComp,"invalid component value");
    T *pt=getPointer()+bgTuples*nbComp+bgComp;
    for(mcIdType i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
      for(mcIdType j=0;j<newNbOfComp;j++)
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
  void DataArrayTemplate<T>::setPartOfValues2(const typename Traits<T>::ArrayType *a, const mcIdType *bgTuples, const mcIdType *endTuples, const mcIdType *bgComp, const mcIdType *endComp, bool strictCompoCompare)
  {
    if(!a)
      throw INTERP_KERNEL::Exception("DataArrayDouble::setPartOfValues2 : input DataArrayDouble is NULL !");
    const char msg[]="DataArrayTemplate::setPartOfValues2";
    checkAllocated();
    a->checkAllocated();
    std::size_t nbComp(getNumberOfComponents());
    mcIdType nbOfTuples(getNumberOfTuples());
    for(const mcIdType *z=bgComp;z!=endComp;z++)
      DataArray::CheckValueInRange(ToIdType(nbComp),*z,"invalid component id");
    mcIdType newNbOfTuples(ToIdType(std::distance(bgTuples,endTuples)));
    mcIdType newNbOfComp(ToIdType(std::distance(bgComp,endComp)));
    bool assignTech(true);
    if(a->getNbOfElems()==newNbOfTuples*newNbOfComp)
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
        for(const mcIdType *w=bgTuples;w!=endTuples;w++)
          {
            DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
            for(const mcIdType *z=bgComp;z!=endComp;z++,srcPt++)
              {
                pt[(std::size_t)(*w)*nbComp+(*z)]=*srcPt;
              }
          }
      }
    else
      {
        for(const mcIdType *w=bgTuples;w!=endTuples;w++)
          {
            const T *srcPt2=srcPt;
            DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
            for(const mcIdType *z=bgComp;z!=endComp;z++,srcPt2++)
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
  void DataArrayTemplate<T>::setPartOfValuesSimple2(T a, const mcIdType *bgTuples, const mcIdType *endTuples, const mcIdType *bgComp, const mcIdType *endComp)
  {
    checkAllocated();
    std::size_t nbComp=getNumberOfComponents();
    mcIdType nbOfTuples=getNumberOfTuples();
    for(const mcIdType *z=bgComp;z!=endComp;z++)
      DataArray::CheckValueInRange(ToIdType(nbComp),*z,"invalid component id");
    T *pt(getPointer());
    for(const mcIdType *w=bgTuples;w!=endTuples;w++)
      for(const mcIdType *z=bgComp;z!=endComp;z++)
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
  void DataArrayTemplate<T>::setPartOfValues3(const typename Traits<T>::ArrayType *a, const mcIdType *bgTuples, const mcIdType *endTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp, bool strictCompoCompare)
  {
    if(!a)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setPartOfValues3 : input DataArrayDouble is NULL !");
    const char msg[]="DataArrayTemplate::setPartOfValues3";
    checkAllocated();
    a->checkAllocated();
    mcIdType newNbOfComp=DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg);
    std::size_t nbComp(getNumberOfComponents());
    mcIdType nbOfTuples(getNumberOfTuples());
    DataArray::CheckValueInRangeEx(ToIdType(nbComp),bgComp,endComp,"invalid component value");
    mcIdType newNbOfTuples=ToIdType(std::distance(bgTuples,endTuples));
    bool assignTech=true;
    if(a->getNbOfElems()==newNbOfTuples*newNbOfComp)
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
        for(const mcIdType *w=bgTuples;w!=endTuples;w++)
          for(mcIdType j=0;j<newNbOfComp;j++,srcPt++)
            {
              DataArray::CheckValueInRange(nbOfTuples,*w,"invalid tuple id");
              pt[(std::size_t)(*w)*nbComp+j*stepComp]=*srcPt;
            }
      }
    else
      {
        for(const mcIdType *w=bgTuples;w!=endTuples;w++)
          {
            const T *srcPt2=srcPt;
            for(mcIdType j=0;j<newNbOfComp;j++,srcPt2++)
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
  void DataArrayTemplate<T>::setPartOfValuesSimple3(T a, const mcIdType *bgTuples, const mcIdType *endTuples, mcIdType bgComp, mcIdType endComp, mcIdType stepComp)
  {
    const char msg[]="DataArrayTemplate::setPartOfValuesSimple3";
    checkAllocated();
    std::size_t newNbOfComp(DataArray::GetNumberOfItemGivenBES(bgComp,endComp,stepComp,msg));
    std::size_t nbComp(getNumberOfComponents());
    mcIdType nbOfTuples(getNumberOfTuples());
    DataArray::CheckValueInRangeEx(ToIdType(nbComp),bgComp,endComp,"invalid component value");
    T *pt(getPointer()+bgComp);
    for(const mcIdType *w=bgTuples;w!=endTuples;w++)
      for(std::size_t j=0;j<newNbOfComp;j++)
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
  void DataArrayTemplate<T>::setPartOfValues4(const typename Traits<T>::ArrayType *a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, const mcIdType *bgComp, const mcIdType *endComp, bool strictCompoCompare)
  {if(!a)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setPartOfValues4 : input DataArrayTemplate is NULL !");
    const char msg[]="DataArrayTemplate::setPartOfValues4";
    checkAllocated();
    a->checkAllocated();
    mcIdType newNbOfTuples(DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg));
    std::size_t newNbOfComp(std::distance(bgComp,endComp));
    std::size_t nbComp(getNumberOfComponents());
    for(const mcIdType *z=bgComp;z!=endComp;z++)
      DataArray::CheckValueInRange(ToIdType(nbComp),*z,"invalid component id");
    mcIdType nbOfTuples(getNumberOfTuples());
    DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
    bool assignTech(true);
    if(a->getNbOfElems()==ToIdType(newNbOfTuples*newNbOfComp))
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
        for(mcIdType i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
          for(const mcIdType *z=bgComp;z!=endComp;z++,srcPt++)
            pt[*z]=*srcPt;
      }
    else
      {
      for(mcIdType i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
        {
          const T *srcPt2(srcPt);
          for(const mcIdType *z=bgComp;z!=endComp;z++,srcPt2++)
            pt[*z]=*srcPt2;
        }
      }
  }

  template<class T>
  void DataArrayTemplate<T>::setPartOfValuesSimple4(T a, mcIdType bgTuples, mcIdType endTuples, mcIdType stepTuples, const mcIdType *bgComp, const mcIdType *endComp)
  {
    const char msg[]="DataArrayTemplate::setPartOfValuesSimple4";
    checkAllocated();
    mcIdType newNbOfTuples(DataArray::GetNumberOfItemGivenBES(bgTuples,endTuples,stepTuples,msg));
    std::size_t nbComp(getNumberOfComponents());
    for(const mcIdType *z=bgComp;z!=endComp;z++)
      DataArray::CheckValueInRange(ToIdType(nbComp),*z,"invalid component id");
    mcIdType nbOfTuples(getNumberOfTuples());
    DataArray::CheckValueInRangeEx(nbOfTuples,bgTuples,endTuples,"invalid tuple value");
    T *pt=getPointer()+bgTuples*nbComp;
    for(mcIdType i=0;i<newNbOfTuples;i++,pt+=stepTuples*nbComp)
      for(const mcIdType *z=bgComp;z!=endComp;z++)
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
  void DataArrayTemplate<T>::setPartOfValuesAdv(const typename Traits<T>::ArrayType *a, const DataArrayIdType *tuplesSelec)
  {
    if(!a || !tuplesSelec)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setPartOfValuesAdv : input DataArrayTemplate is NULL !");
    checkAllocated();
    a->checkAllocated();
    tuplesSelec->checkAllocated();
    std::size_t nbOfComp(getNumberOfComponents());
    if(nbOfComp!=a->getNumberOfComponents())
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setPartOfValuesAdv : This and a do not have the same number of components !");
    if(tuplesSelec->getNumberOfComponents()!=2)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setPartOfValuesAdv : Expecting to have a tuple selector DataArrayInt instance with exactly 2 components !");
    mcIdType thisNt(getNumberOfTuples());
    mcIdType aNt(a->getNumberOfTuples());
    T *valsToSet(getPointer());
    const T *valsSrc(a->getConstPointer());
    for(const mcIdType *tuple=tuplesSelec->begin();tuple!=tuplesSelec->end();tuple+=2)
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
  void DataArrayTemplate<T>::setContigPartOfSelectedValues(mcIdType tupleIdStart, const DataArray *aBase, const DataArrayIdType *tuplesSelec)
  {
    if(!aBase || !tuplesSelec)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValues : input DataArray is NULL !");
    const typename Traits<T>::ArrayType *a(dynamic_cast<const typename Traits<T>::ArrayType *>(aBase));
    if(!a)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValues : input DataArray aBase is not a DataArrayDouble !");
    checkAllocated();
    a->checkAllocated();
    tuplesSelec->checkAllocated();
    std::size_t nbOfComp(getNumberOfComponents());
    if(nbOfComp!=a->getNumberOfComponents())
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValues : This and a do not have the same number of components !");
    if(tuplesSelec->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValues : Expecting to have a tuple selector DataArrayInt instance with exactly 1 component !");
    mcIdType thisNt(getNumberOfTuples());
    mcIdType aNt(a->getNumberOfTuples());
    mcIdType nbOfTupleToWrite(tuplesSelec->getNumberOfTuples());
    T *valsToSet(getPointer()+tupleIdStart*nbOfComp);
    if(tupleIdStart+nbOfTupleToWrite>thisNt)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValues : invalid number range of values to write !");
    const T *valsSrc=a->getConstPointer();
    for(const mcIdType *tuple=tuplesSelec->begin();tuple!=tuplesSelec->end();tuple++,valsToSet+=nbOfComp)
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
  void DataArrayTemplate<T>::setContigPartOfSelectedValuesSlice(mcIdType tupleIdStart, const DataArray *aBase, mcIdType bg, mcIdType end2, mcIdType step)
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
    std::size_t nbOfComp(getNumberOfComponents());
    const char msg[]="DataArrayDouble::setContigPartOfSelectedValuesSlice";
    mcIdType nbOfTupleToWrite(DataArray::GetNumberOfItemGivenBES(bg,end2,step,msg));
    if(nbOfComp!=a->getNumberOfComponents())
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValuesSlice : This and a do not have the same number of components !");
    mcIdType thisNt(getNumberOfTuples());
    mcIdType aNt(a->getNumberOfTuples());
    T *valsToSet(getPointer()+tupleIdStart*nbOfComp);
    if(tupleIdStart+nbOfTupleToWrite>thisNt)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValuesSlice : invalid number range of values to write !");
    if(end2>aNt)
      throw INTERP_KERNEL::Exception("DataArrayTemplate::setContigPartOfSelectedValuesSlice : invalid range of values to read !");
    const T *valsSrc(a->getConstPointer()+bg*nbOfComp);
    for(mcIdType i=0;i<nbOfTupleToWrite;i++,valsToSet+=nbOfComp,valsSrc+=step*nbOfComp)
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
  typename Traits<T>::ArrayType *DataArrayTemplate<T>::mySelectByTupleRanges(const std::vector<std::pair<mcIdType,mcIdType> >& ranges) const
  {
    checkAllocated();
    std::size_t nbOfComp(getNumberOfComponents());
    mcIdType nbOfTuplesThis(getNumberOfTuples());
    if(ranges.empty())
      {
        MCAuto<DataArray> ret0(buildNewEmptyInstance());
        MCAuto< typename Traits<T>::ArrayType > ret(DynamicCastSafe<DataArray,typename Traits<T>::ArrayType>(ret0));
        ret->alloc(0,nbOfComp);
        ret->copyStringInfoFrom(*this);
        return ret.retn();
      }
    mcIdType ref(ranges.front().first),nbOfTuples(0);
    bool isIncreasing(true);
    for(std::vector<std::pair<mcIdType,mcIdType> >::const_iterator it=ranges.begin();it!=ranges.end();it++)
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
    for(std::vector<std::pair<mcIdType,mcIdType> >::const_iterator it=ranges.begin();it!=ranges.end();it++)
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
    mcIdType nbOfTuples=getNumberOfTuples();
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
    mcIdType nbOfTuples=getNumberOfTuples();
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
   *  \sa getMaxAbsValue, getMinValue
   */
  template<class T>
  T DataArrayTemplate<T>::getMaxValue(mcIdType& tupleId) const
  {
    checkAllocated();
    if(getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxValue : must be applied on DataArrayDouble with only one component, you can call 'rearrange' method before or call 'getMaxValueInArray' method !");
    mcIdType nbOfTuples=getNumberOfTuples();
    if(nbOfTuples<=0)
      throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxValue : array exists but number of tuples must be > 0 !");
    const T *vals(getConstPointer());
    const T *loc(std::max_element(vals,vals+nbOfTuples));
    tupleId=ToIdType(std::distance(vals,loc));
    return *loc;
  }

  /*!
   * Returns the maximal value within \a this array that is allowed to have more than
   *  one component.
   *  \return double - the maximal value among all values of \a this array.
   *  \throw If \a this is not allocated.
   *         If \a this is empty.
   *  \sa getMaxAbsValueInArray, getMinValueInArray
   */
  template<class T>
  T DataArrayTemplate<T>::getMaxValueInArray() const
  {
    checkAllocated();
    if( empty() )
      THROW_IK_EXCEPTION("getMaxValueInArray : this is empty !");
    const T *loc(std::max_element(begin(),end()));
    return *loc;
  }

  /*!
   * Returns the maximal absolute value in \a this and the first occurrence location associated to it.
   * \return the element in this (positive or negative) having the max abs value in \a this.
   *  \throw If \a this is not allocated.
   *  \throw If \a this is non one component array.
   *  \throw If \a this is empty.
   */
  template<class T>
  T DataArrayTemplate<T>::getMaxAbsValue(std::size_t& tupleId) const
  {
    checkAllocated();
    if(getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayDouble::getMaxAbsValue : must be applied on DataArrayDouble with only one component, you can call 'rearrange' method before or call 'getMaxValueInArray' method !");
    mcIdType nbTuples(this->getNumberOfTuples());
    if(nbTuples==0)
      throw INTERP_KERNEL::Exception("DataArrayTemplate<T>::getMaxAbsValue : empty array !");
    T ret((T)-1);
    tupleId=0;
    const T *pt(begin());
    for(mcIdType i=0;i<nbTuples;i++,pt++)
      {
        T cand((T)std::abs(*pt));
        if(cand>ret)
          {
            ret=cand;
            tupleId=i;
          }
      }
    return this->getIJ(ToIdType(tupleId),0);
  }

  /*!
   * Returns the maximal absolute value in \a this.
   *  \throw If \a this is not allocated.
   *  \throw If \a this is non one component array.
   *  \throw If \a this is empty.
   */
  template<class T>
  T DataArrayTemplate<T>::getMaxAbsValueInArray() const
  {
    std::size_t dummy;
    return getMaxAbsValue(dummy);
  }

  /*!
   * Returns the minimal value and its location within \a this one-dimensional array.
   *  \param [out] tupleId - index of the tuple holding the minimal value.
   *  \return double - the minimal value among all values of \a this array.
   *  \throw If \a this->getNumberOfComponents() != 1
   *  \throw If \a this->getNumberOfTuples() < 1
   */
  template<class T>
  T DataArrayTemplate<T>::getMinValue(mcIdType& tupleId) const
  {
    checkAllocated();
    if(getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayDouble::getMinValue : must be applied on DataArrayDouble with only one component, you can call 'rearrange' method before call 'getMinValueInArray' method !");
    mcIdType nbOfTuples=getNumberOfTuples();
    if(nbOfTuples<=0)
      throw INTERP_KERNEL::Exception("DataArrayDouble::getMinValue : array exists but number of tuples must be > 0 !");
    const T *vals(getConstPointer());
    const T *loc(std::min_element(vals,vals+nbOfTuples));
    tupleId=ToIdType(std::distance(vals,loc));
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

  template<class T>
  void DataArrayTemplate<T>::circularPermutation(mcIdType nbOfShift)
  {
    checkAllocated();
    std::size_t nbOfCompo(getNumberOfComponents());
    mcIdType nbTuples(getNumberOfTuples());
    mcIdType effNbSh(EffectiveCircPerm(nbOfShift,nbTuples));
    if(effNbSh==0)
      return ;
    T *work(getPointer());
    if(effNbSh<nbTuples-effNbSh)
      {
        typename INTERP_KERNEL::AutoPtr<T> buf(new T[effNbSh*nbOfCompo]);
        std::copy(work,work+effNbSh*nbOfCompo,(T *)buf);
        std::copy(work+effNbSh*nbOfCompo,work+nbTuples*nbOfCompo,work);// ze big shift
        std::copy((T *)buf,(T *)buf+effNbSh*nbOfCompo,work+(nbTuples-effNbSh)*nbOfCompo);
      }
    else
      {
        typename INTERP_KERNEL::AutoPtr<T> buf(new T[(nbTuples-effNbSh)*nbOfCompo]);
        std::copy(work+effNbSh*nbOfCompo,work+nbTuples*nbOfCompo,(T *)buf);
        std::copy(work,work+effNbSh*nbOfCompo,work+(nbTuples-effNbSh)*nbOfCompo);// ze big shift
        std::copy((T*)buf,(T *)buf+(nbTuples-effNbSh)*nbOfCompo,work);
      }
  }

  template<class T>
  void DataArrayTemplate<T>::circularPermutationPerTuple(mcIdType nbOfShift)
  {
    checkAllocated();
    std::size_t nbOfCompo(getNumberOfComponents());
    mcIdType nbTuples(getNumberOfTuples());
    mcIdType effNbSh(EffectiveCircPerm(nbOfShift,ToIdType(nbOfCompo)));
    if(effNbSh==0)
      return ;
    T *work(getPointer());
    if(effNbSh<ToIdType(nbOfCompo)-effNbSh)
      {
        typename INTERP_KERNEL::AutoPtr<T> buf(new T[effNbSh]);
        for(mcIdType i=0;i<nbTuples;i++,work+=nbOfCompo)
          {
            std::copy(work,work+effNbSh,(T *)buf);
            std::copy(work+effNbSh,work+nbOfCompo,work);// ze big shift
            std::copy((T *)buf,(T *)buf+effNbSh,work+(nbOfCompo-effNbSh));
          }
      }
    else
      {
        typename INTERP_KERNEL::AutoPtr<T> buf(new T[nbOfCompo-effNbSh]);
        for(mcIdType i=0;i<nbTuples;i++,work+=nbOfCompo)
          {
            std::copy(work+effNbSh,work+nbOfCompo,(T *)buf);
            std::copy(work,work+effNbSh,work+(nbOfCompo-effNbSh));// ze big shift
            std::copy((T*)buf,(T *)buf+(nbOfCompo-effNbSh),work);
          }
      }
    std::vector<std::string> sts(nbOfCompo);
    for(std::size_t i=0;i<nbOfCompo;i++)
      sts[i]=_info_on_compo[(i+effNbSh)%nbOfCompo];
    setInfoOnComponents(sts);
  }

  template<class T>
  void DataArrayTemplate<T>::reversePerTuple()
  {
    checkAllocated();
    std::size_t nbOfCompo(getNumberOfComponents());
    mcIdType nbTuples(getNumberOfTuples());
    if(nbOfCompo<=1)
      return ;
    T *work(getPointer());
    for(mcIdType i=0;i<nbTuples;i++,work+=nbOfCompo)
      std::reverse(work,work+nbOfCompo);
    std::reverse(_info_on_compo.begin(),_info_on_compo.end());
  }

  /*!
   * Assign pointer to one array to a pointer to another appay. Reference counter of
   * \a arrayToSet is incremented / decremented.
   *  \param [in] newArray - the pointer to array to assign to \a arrayToSet.
   *  \param [in,out] arrayToSet - the pointer to array to assign to.
   */
  template<class T>
  void DataArrayTemplate<T>::SetArrayIn(typename Traits<T>::ArrayType *newArray, typename Traits<T>::ArrayType* &arrayToSet)
  {
    if(newArray!=arrayToSet)
      {
        if(arrayToSet)
          arrayToSet->decrRef();
        arrayToSet=newArray;
        if(arrayToSet)
          arrayToSet->incrRef();
      }
  }

    /*!
   * Assign zero to all values in \a this array. To know more on filling arrays see
   * \ref MEDCouplingArrayFill.
   * \throw If \a this is not allocated.
   */
  template<class T>
  void DataArrayTemplate<T>::fillWithZero()
  {
    fillWithValue((T)0);
  }

  /*!
  * Write an array in a file for debugging purpose. Do not use it for long term storage. Use it with care.
  */
  template<class T>
  void DataArrayTemplate<T>::writeForDbg(const std::string& fileName) const
  {
    checkAllocated();
    mcIdType nbTuples( getNumberOfTuples() );
    std::size_t nbCompo( getNumberOfComponents() );
    std::ofstream outFile(fileName, std::ios::binary);
    char szOfData( sizeof(Type) );
    outFile.exceptions(std::ios::failbit | std::ios::badbit);
    {
      outFile.write(&szOfData,1);
      outFile.write(reinterpret_cast<const char*>(&nbTuples),sizeof(mcIdType));
      outFile.write(reinterpret_cast<const char*>(&nbCompo),sizeof(std::size_t));
      outFile.write(reinterpret_cast<const char*>(begin()),nbCompo*nbTuples*sizeof(Type));
    }
  }

  template<class T>
  MCAuto<typename Traits<T>::ArrayTypeCh> DataArrayTemplate<T>::LoadForDbg(const std::string& fileName)
  {
    MCAuto< typename Traits<T>::ArrayTypeCh > ret(Traits<T>::ArrayTypeCh::New());
    std::ifstream inFile(fileName, std::ios::binary);
    inFile.exceptions(std::ios::failbit | std::ios::badbit);
    inFile.seekg(0, std::ios::end);
    mcIdType fileSize( ToIdType( inFile.tellg() ) );
    if( fileSize < 1 + ToIdType( sizeof(std::size_t) ) + ToIdType( sizeof(mcIdType) ) )
      THROW_IK_EXCEPTION( "Input file \"" << fileName << "\" is invalid !" );
    inFile.seekg(0, std::ios::beg);
    char szOfData(0);
    inFile.read(&szOfData,1);
    if( szOfData != sizeof(Type) )
      THROW_IK_EXCEPTION( "Input file \"" << fileName << "\" contain atomic data of size " << (int) szOfData << " this atomic data size is " << sizeof(Type) );
    mcIdType nbTuples(0);
    inFile.read(reinterpret_cast<char *>(&nbTuples),sizeof(mcIdType));
    std::size_t nbCompo(0);
    inFile.read(reinterpret_cast<char *>(&nbCompo),sizeof(std::size_t));
    ret->alloc(nbTuples,nbCompo);
    mcIdType sizeOfDataToRead( ToIdType(nbCompo)*nbTuples*szOfData );
    mcIdType sizeExpected( 1 + ToIdType( sizeof(std::size_t) ) + ToIdType( sizeof(mcIdType) ) + sizeOfDataToRead );
    if( fileSize != sizeExpected )
      THROW_IK_EXCEPTION( "Input file \"" << fileName << "\" length is invalid : Size expected " << sizeExpected << " actual : " << fileSize );
    inFile.read(reinterpret_cast<char *>(ret->getPointer()),sizeOfDataToRead);
    return ret;
  }

  //////////////////////////////

  namespace
  {
    // local static function to copy arrays without warnings
    template <class TIn, class TOut>
    static void copyCast (const TIn *begin, const TIn *end, TOut* dest)
    {
      for (const TIn *src = begin; src != end; ++src, ++dest)
        *dest=static_cast<TOut>(*src);
    }
  }

  template<class T>
  template<class U>
  MCAuto< typename Traits<U>::ArrayType > DataArrayTemplateClassic<T>::convertToOtherTypeOfArr() const
  {
    this->checkAllocated();
    MCAuto<typename Traits<U>::ArrayType> ret(Traits<U>::ArrayType::New());
    ret->alloc(this->getNumberOfTuples(),this->getNumberOfComponents());
    std::size_t nbOfVals(this->getNbOfElems());
    const T *src(this->begin());
    U *dest(ret->getPointer());
    // to make Visual C++ happy : instead of std::size_t nbOfVals=getNbOfElems(); std::copy(src,src+nbOfVals,dest);
    copyCast(src, src+nbOfVals, dest);
    //std::copy(src,src+nbOfVals,dest);
    ret->copyStringInfoFrom(*this);
    return ret;
  }

  /*!
   * Creates a new DataArrayDouble and assigns all (textual and numerical) data of \a this
   * array to the new one.
   *  \return DataArrayDouble * - the new instance of DataArrayInt.
   */
  template<class T>
  MCAuto<DataArrayDouble> DataArrayTemplateClassic<T>::convertToDblArr() const
  {
    return convertToOtherTypeOfArr<double>();
  }

  /*!
   * Creates a new DataArrayInt and assigns all (textual and numerical) data of \a this
   * array to the new one.
   *  \return DataArrayInt * - the new instance of DataArrayInt.
   */
  template<class T>
  MCAuto<DataArrayInt> DataArrayTemplateClassic<T>::convertToIntArr() const
  {
    return convertToOtherTypeOfArr<int>();
  }

  /*!
   * Creates a new DataArrayInt64 and assigns all (textual and numerical) data of \a this
   * array to the new one.
   *  \return DataArrayInt * - the new instance of DataArrayInt64.
   */
  template<class T>
  MCAuto<DataArrayInt64> DataArrayTemplateClassic<T>::convertToInt64Arr() const
  {
    return convertToOtherTypeOfArr<Int64>();
  }
  /*!
   * Creates a new DataArrayFloat and assigns all (textual and numerical) data of \a this
   * array to the new one.
   *  \return DataArrayFloat * - the new instance of DataArrayInt.
   */
  template<class T>
  MCAuto<DataArrayFloat> DataArrayTemplateClassic<T>::convertToFloatArr() const
  {
    return convertToOtherTypeOfArr<float>();
  }

  /*!
   * Apply a linear function to a given component of \a this array, so that
   * an array element <em>(x)</em> becomes \f$ a * x + b \f$.
   *  \param [in] a - the first coefficient of the function.
   *  \param [in] b - the second coefficient of the function.
   *  \param [in] compoId - the index of component to modify.
   *  \throw If \a this is not allocated, or \a compoId is not in [0,\c this->getNumberOfComponents() ).
   */
  template<class T>
  void DataArrayTemplateClassic<T>::applyLin(T a, T b, std::size_t compoId)
  {
    this->checkAllocated();
    std::size_t nbOfComp=this->getNumberOfComponents();
    if(compoId>=nbOfComp)
      {
        std::ostringstream oss; oss << "DataArrayDouble::applyLin : The compoId requested (" << compoId << ") is not valid ! Must be in [0," << nbOfComp << ") !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    T *ptr(this->getPointer()+compoId);
    mcIdType nbOfTuple=this->getNumberOfTuples();
    for(mcIdType i=0;i<nbOfTuple;i++,ptr+=nbOfComp)
      *ptr=a*(*ptr)+b;
    this->declareAsNew();
  }

  /*!
   * Apply a linear function to all elements of \a this array, so that
   * an element _x_ becomes \f$ a * x + b \f$.
   *  \param [in] a - the first coefficient of the function.
   *  \param [in] b - the second coefficient of the function.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  void DataArrayTemplateClassic<T>::applyLin(T a, T b)
  {
    this->checkAllocated();
    T *ptr(this->getPointer());
    std::size_t nbOfElems(this->getNbOfElems());
    for(std::size_t i=0;i<nbOfElems;i++,ptr++)
      *ptr=a*(*ptr)+b;
    this->declareAsNew();
  }

  /*!
   * Returns a full copy of \a this array except that sign of all elements is reversed.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble containing the
   *          same number of tuples and component as \a this array.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::negate() const
  {
    this->checkAllocated();
    MCAuto<typename Traits<T>::ArrayType> newArr(Traits<T>::ArrayType::New());
    mcIdType nbOfTuples(this->getNumberOfTuples());
    std::size_t nbOfComp(this->getNumberOfComponents());
    newArr->alloc(nbOfTuples,nbOfComp);
    const T *cptr(this->begin());
    std::transform(cptr,cptr+nbOfTuples*nbOfComp,newArr->getPointer(),std::negate<T>());
    newArr->copyStringInfoFrom(*this);
    return newArr.retn();
  }

  template<class T>
  template<class FCT>
  void DataArrayTemplateClassic<T>::somethingEqual(const typename Traits<T>::ArrayType *other)
  {
    if(!other)
      throw INTERP_KERNEL::Exception("DataArray<T>::SomethingEqual : input DataArray<T> instance is NULL !");
    const char *msg="Nb of tuples mismatch for DataArrayDouble::multiplyEqual !";
    this->checkAllocated();
    other->checkAllocated();
    mcIdType nbOfTuple(this->getNumberOfTuples());
    mcIdType nbOfTuple2(other->getNumberOfTuples());
    std::size_t nbOfComp(this->getNumberOfComponents());
    std::size_t nbOfComp2(other->getNumberOfComponents());
    if(nbOfTuple==nbOfTuple2)
      {
        if(nbOfComp==nbOfComp2)
          {
            std::transform(this->begin(),this->end(),other->begin(),this->getPointer(),FCT());
          }
        else if(nbOfComp2==1)
          {
            T *ptr(this->getPointer());
            const T *ptrc(other->begin());
            for(mcIdType i=0;i<nbOfTuple;i++)
              std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptr+i*nbOfComp,std::bind(FCT(),std::placeholders::_1,*ptrc++));
          }
        else
          throw INTERP_KERNEL::Exception(msg);
      }
    else if(nbOfTuple2==1)
      {
        if(nbOfComp2==nbOfComp)
          {
            T *ptr(this->getPointer());
            const T *ptrc(other->begin());
            for(mcIdType i=0;i<nbOfTuple;i++)
              std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptrc,ptr+i*nbOfComp,FCT());
          }
        else
          throw INTERP_KERNEL::Exception(msg);
      }
    else
      throw INTERP_KERNEL::Exception(msg);
    this->declareAsNew();
  }

  /*!
   * Adds values of another DataArrayDouble to values of \a this one. There are 3
   * valid cases.
   * 1.  The arrays have same number of tuples and components. Then each value of
   *   \a other array is added to the corresponding value of \a this array, i.e.:
   *   _a_ [ i, j ] += _other_ [ i, j ].
   * 2.  The arrays have same number of tuples and \a other array has one component. Then
   *   _a_ [ i, j ] += _other_ [ i, 0 ].
   * 3.  The arrays have same number of components and \a other array has one tuple. Then
   *   _a_ [ i, j ] += _a2_ [ 0, j ].
   *
   *  \param [in] other - an array to add to \a this one.
   *  \throw If \a other is NULL.
   *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
   *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
   *         \a other has number of both tuples and components not equal to 1.
   */
  template<class T>
  void DataArrayTemplateClassic<T>::addEqual(const typename Traits<T>::ArrayType *other)
  {
    this->somethingEqual< std::plus<T> >(other);
  }

  /*!
   * Subtract values of another DataArrayDouble from values of \a this one. There are 3
   * valid cases.
   * 1.  The arrays have same number of tuples and components. Then each value of
   *   \a other array is subtracted from the corresponding value of \a this array, i.e.:
   *   _a_ [ i, j ] -= _other_ [ i, j ].
   * 2.  The arrays have same number of tuples and \a other array has one component. Then
   *   _a_ [ i, j ] -= _other_ [ i, 0 ].
   * 3.  The arrays have same number of components and \a other array has one tuple. Then
   *   _a_ [ i, j ] -= _a2_ [ 0, j ].
   *
   *  \param [in] other - an array to subtract from \a this one.
   *  \throw If \a other is NULL.
   *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
   *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
   *         \a other has number of both tuples and components not equal to 1.
   */
  template<class T>
  void DataArrayTemplateClassic<T>::substractEqual(const typename Traits<T>::ArrayType *other)
  {
    this->somethingEqual< std::minus<T> >(other);
  }

  /*!
   * Multiply values of another DataArrayDouble to values of \a this one. There are 3
   * valid cases.
   * 1.  The arrays have same number of tuples and components. Then each value of
   *   \a other array is multiplied to the corresponding value of \a this array, i.e.
   *   _this_ [ i, j ] *= _other_ [ i, j ].
   * 2.  The arrays have same number of tuples and \a other array has one component. Then
   *   _this_ [ i, j ] *= _other_ [ i, 0 ].
   * 3.  The arrays have same number of components and \a other array has one tuple. Then
   *   _this_ [ i, j ] *= _a2_ [ 0, j ].
   *
   *  \param [in] other - an array to multiply to \a this one.
   *  \throw If \a other is NULL.
   *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
   *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
   *         \a other has number of both tuples and components not equal to 1.
   */
  template<class T>
  void DataArrayTemplateClassic<T>::multiplyEqual(const typename Traits<T>::ArrayType *other)
  {
    this->somethingEqual< std::multiplies<T> >(other);
  }

  /*!
   * Divide values of \a this array by values of another DataArrayDouble. There are 3
   * valid cases.
   * 1.  The arrays have same number of tuples and components. Then each value of
   *    \a this array is divided by the corresponding value of \a other one, i.e.:
   *   _a_ [ i, j ] /= _other_ [ i, j ].
   * 2.  The arrays have same number of tuples and \a other array has one component. Then
   *   _a_ [ i, j ] /= _other_ [ i, 0 ].
   * 3.  The arrays have same number of components and \a other array has one tuple. Then
   *   _a_ [ i, j ] /= _a2_ [ 0, j ].
   *
   *  \warning No check of division by zero is performed!
   *  \param [in] other - an array to divide \a this one by.
   *  \throw If \a other is NULL.
   *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
   *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
   *         \a other has number of both tuples and components not equal to 1.
   */
  template<class T>
  void DataArrayTemplateClassic<T>::divideEqual(const typename Traits<T>::ArrayType *other)
  {
    this->somethingEqual< std::divides<T> >(other);
  }

  template<class T, class FCT>
  typename Traits<T>::ArrayType *DivSub(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2)
  {
    if(!a1 || !a2)
      throw INTERP_KERNEL::Exception("DivSub : input DataArrayDouble instance is NULL !");
    mcIdType nbOfTuple1(a1->getNumberOfTuples());
    mcIdType nbOfTuple2(a2->getNumberOfTuples());
    std::size_t nbOfComp1(a1->getNumberOfComponents());
    std::size_t nbOfComp2(a2->getNumberOfComponents());
    if(nbOfTuple2==nbOfTuple1)
      {
        if(nbOfComp1==nbOfComp2)
          {
            MCAuto<typename Traits<T>::ArrayType> ret(Traits<T>::ArrayType::New());
            ret->alloc(nbOfTuple2,nbOfComp1);
            std::transform(a1->begin(),a1->end(),a2->begin(),ret->getPointer(),FCT());
            ret->copyStringInfoFrom(*a1);
            return ret.retn();
          }
        else if(nbOfComp2==1)
          {
            MCAuto<typename Traits<T>::ArrayType> ret(Traits<T>::ArrayType::New());
            ret->alloc(nbOfTuple1,nbOfComp1);
            const T *a2Ptr(a2->begin()),*a1Ptr(a1->begin());
            T *res(ret->getPointer());
            for(mcIdType i=0;i<nbOfTuple1;i++)
              res=std::transform(a1Ptr+i*nbOfComp1,a1Ptr+(i+1)*nbOfComp1,res,std::bind(FCT(),std::placeholders::_1,a2Ptr[i]));
            ret->copyStringInfoFrom(*a1);
            return ret.retn();
          }
        else
          {
            a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Divide !");
            return 0;
          }
      }
    else if(nbOfTuple2==1)
      {
        a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Divide !");
        MCAuto<typename Traits<T>::ArrayType> ret(Traits<T>::ArrayType::New());
        ret->alloc(nbOfTuple1,nbOfComp1);
        const T *a1ptr=a1->begin(),*a2ptr(a2->begin());
        T *pt(ret->getPointer());
        for(mcIdType i=0;i<nbOfTuple1;i++)
          pt=std::transform(a1ptr+i*nbOfComp1,a1ptr+(i+1)*nbOfComp1,a2ptr,pt,FCT());
        ret->copyStringInfoFrom(*a1);
        return ret.retn();
      }
    else
      {
        a1->checkNbOfTuples(nbOfTuple2,"Nb of tuples mismatch for array Divide !");//will always throw an exception
        return 0;
      }
  }

  /*!
   * Returns a new DataArrayDouble that is a subtraction of two given arrays. There are 3
   * valid cases.
   * 1.  The arrays have same number of tuples and components. Then each value of
   *   the result array (_a_) is a subtraction of the corresponding values of \a a1 and
   *   \a a2, i.e.: _a_ [ i, j ] = _a1_ [ i, j ] - _a2_ [ i, j ].
   * 2.  The arrays have same number of tuples and one array, say _a2_, has one
   *   component. Then
   *   _a_ [ i, j ] = _a1_ [ i, j ] - _a2_ [ i, 0 ].
   * 3.  The arrays have same number of components and one array, say _a2_, has one
   *   tuple. Then
   *   _a_ [ i, j ] = _a1_ [ i, j ] - _a2_ [ 0, j ].
   *
   * Info on components is copied either from the first array (in the first case) or from
   * the array with maximal number of elements (getNbOfElems()).
   *  \param [in] a1 - an array to subtract from.
   *  \param [in] a2 - an array to subtract.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If either \a a1 or \a a2 is NULL.
   *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
   *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
   *         none of them has number of tuples or components equal to 1.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::Substract(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2)
  {
    return DivSub< T,std::minus<T> >(a1,a2);
  }

  /*!
   * Returns a new DataArrayDouble that is a division of two given arrays. There are 3
   * valid cases.
   * 1.  The arrays have same number of tuples and components. Then each value of
   *   the result array (_a_) is a division of the corresponding values of \a a1 and
   *   \a a2, i.e.: _a_ [ i, j ] = _a1_ [ i, j ] / _a2_ [ i, j ].
   * 2.  The arrays have same number of tuples and one array, say _a2_, has one
   *   component. Then
   *   _a_ [ i, j ] = _a1_ [ i, j ] / _a2_ [ i, 0 ].
   * 3.  The arrays have same number of components and one array, say _a2_, has one
   *   tuple. Then
   *   _a_ [ i, j ] = _a1_ [ i, j ] / _a2_ [ 0, j ].
   *
   * Info on components is copied either from the first array (in the first case) or from
   * the array with maximal number of elements (getNbOfElems()).
   *  \warning No check of division by zero is performed!
   *  \param [in] a1 - a numerator array.
   *  \param [in] a2 - a denominator array.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If either \a a1 or \a a2 is NULL.
   *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
   *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
   *         none of them has number of tuples or components equal to 1.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::Divide(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2)
  {
    return DivSub< T,std::divides<T> >(a1,a2);
  }

  template<class T, class FCT>
  typename Traits<T>::ArrayType *MulAdd(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2)
  {
    if(!a1 || !a2)
      throw INTERP_KERNEL::Exception("DataArrayDouble::MulAdd : input DataArrayDouble instance is NULL !");
    mcIdType nbOfTuple(a1->getNumberOfTuples());
    mcIdType nbOfTuple2(a2->getNumberOfTuples());
    std::size_t nbOfComp(a1->getNumberOfComponents());
    std::size_t nbOfComp2(a2->getNumberOfComponents());
    MCAuto<typename Traits<T>::ArrayType> ret=0;
    if(nbOfTuple==nbOfTuple2)
      {
        if(nbOfComp==nbOfComp2)
          {
            ret=Traits<T>::ArrayType::New();
            ret->alloc(nbOfTuple,nbOfComp);
            std::transform(a1->begin(),a1->end(),a2->begin(),ret->getPointer(),FCT());
            ret->copyStringInfoFrom(*a1);
          }
        else
          {
            std::size_t nbOfCompMin,nbOfCompMax;
            const typename Traits<T>::ArrayType *aMin, *aMax;
            if(nbOfComp>nbOfComp2)
              {
                nbOfCompMin=nbOfComp2; nbOfCompMax=nbOfComp;
                aMin=a2; aMax=a1;
              }
            else
              {
                nbOfCompMin=nbOfComp; nbOfCompMax=nbOfComp2;
                aMin=a1; aMax=a2;
              }
            if(nbOfCompMin==1)
              {
                ret=Traits<T>::ArrayType::New();
                ret->alloc(nbOfTuple,nbOfCompMax);
                const T *aMinPtr(aMin->begin());
                const T *aMaxPtr(aMax->begin());
                T *res=ret->getPointer();
                for(mcIdType i=0;i<nbOfTuple;i++)
                  res=std::transform(aMaxPtr+i*nbOfCompMax,aMaxPtr+(i+1)*nbOfCompMax,res,std::bind(FCT(),std::placeholders::_1,aMinPtr[i]));
                ret->copyStringInfoFrom(*aMax);
              }
            else
              throw INTERP_KERNEL::Exception("Nb of components mismatch for array MulAdd !");
          }
      }
    else if((nbOfTuple==1 && nbOfTuple2>1) || (nbOfTuple>1 && nbOfTuple2==1))
      {
        if(nbOfComp==nbOfComp2)
          {
            mcIdType nbOfTupleMax=std::max(nbOfTuple,nbOfTuple2);
            const typename Traits<T>::ArrayType *aMin(nbOfTuple>nbOfTuple2?a2:a1);
            const typename Traits<T>::ArrayType *aMax(nbOfTuple>nbOfTuple2?a1:a2);
            const T *aMinPtr(aMin->begin()),*aMaxPtr(aMax->begin());
            ret=Traits<T>::ArrayType::New();
            ret->alloc(nbOfTupleMax,nbOfComp);
            T *res(ret->getPointer());
            for(mcIdType i=0;i<nbOfTupleMax;i++)
              res=std::transform(aMaxPtr+i*nbOfComp,aMaxPtr+(i+1)*nbOfComp,aMinPtr,res,FCT());
            ret->copyStringInfoFrom(*aMax);
          }
        else
          throw INTERP_KERNEL::Exception("Nb of components mismatch for array MulAdd !");
      }
    else
      throw INTERP_KERNEL::Exception("Nb of tuples mismatch for array MulAdd !");
    return ret.retn();
  }

  /*!
   * Returns a new DataArrayDouble that is a product of two given arrays. There are 3
   * valid cases.
   * 1.  The arrays have same number of tuples and components. Then each value of
   *   the result array (_a_) is a product of the corresponding values of \a a1 and
   *   \a a2, i.e. _a_ [ i, j ] = _a1_ [ i, j ] * _a2_ [ i, j ].
   * 2.  The arrays have same number of tuples and one array, say _a2_, has one
   *   component. Then
   *   _a_ [ i, j ] = _a1_ [ i, j ] * _a2_ [ i, 0 ].
   * 3.  The arrays have same number of components and one array, say _a2_, has one
   *   tuple. Then
   *   _a_ [ i, j ] = _a1_ [ i, j ] * _a2_ [ 0, j ].
   *
   * Info on components is copied either from the first array (in the first case) or from
   * the array with maximal number of elements (getNbOfElems()).
   *  \param [in] a1 - a factor array.
   *  \param [in] a2 - another factor array.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If either \a a1 or \a a2 is NULL.
   *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
   *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
   *         none of them has number of tuples or components equal to 1.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::Multiply(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2)
  {
    return MulAdd< T , std::multiplies<T> >(a1,a2);
  }

  /*!
   * Returns a new DataArrayDouble that is a sum of two given arrays. There are 3
   * valid cases.
   * 1.  The arrays have same number of tuples and components. Then each value of
   *   the result array (_a_) is a sum of the corresponding values of \a a1 and \a a2,
   *   i.e.: _a_ [ i, j ] = _a1_ [ i, j ] + _a2_ [ i, j ].
   * 2.  The arrays have same number of tuples and one array, say _a2_, has one
   *   component. Then
   *   _a_ [ i, j ] = _a1_ [ i, j ] + _a2_ [ i, 0 ].
   * 3.  The arrays have same number of components and one array, say _a2_, has one
   *   tuple. Then
   *   _a_ [ i, j ] = _a1_ [ i, j ] + _a2_ [ 0, j ].
   *
   * Info on components is copied either from the first array (in the first case) or from
   * the array with maximal number of elements (getNbOfElems()).
   *  \param [in] a1 - an array to sum up.
   *  \param [in] a2 - another array to sum up.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If either \a a1 or \a a2 is NULL.
   *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
   *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
   *         none of them has number of tuples or components equal to 1.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::Add(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2)
  {
    return MulAdd< T , std::plus<T> >(a1,a2);
  }

  /*!
   * Returns either a \a deep or \a shallow copy of this array. For more info see
   * \ref MEDCouplingArrayBasicsCopyDeep and \ref MEDCouplingArrayBasicsCopyShallow.
   *  \param [in] dCpy - if \a true, a deep copy is returned, else, a shallow one.
   *  \return DataArrayDouble * - either a new instance of DataArrayDouble (if \a dCpy
   *          == \a true) or \a this instance (if \a dCpy == \a false).
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::PerformCopyOrIncrRef(bool dCpy, const typename Traits<T>::ArrayType& self)
  {
    if(dCpy)
      return self.deepCopy();
    else
      {
        self.incrRef();
        return const_cast<typename Traits<T>::ArrayType *>(&self);
      }
  }

  template<class T>
  struct GreatEqual
  {
    GreatEqual(T v):_v(v) { }
    bool operator()(T v) const { return v>=_v; }
    T _v;
  };

  template<class T>
  struct GreaterThan
  {
    GreaterThan(T v):_v(v) { }
    bool operator()(T v) const { return v>_v; }
    T _v;
  };

  template<class T>
  struct LowerEqual
  {
    LowerEqual(T v):_v(v) { }
    bool operator()(T v) const { return v<=_v; }
    T _v;
  };

  template<class T>
  struct LowerThan
  {
    LowerThan(T v):_v(v) { }
    bool operator()(T v) const { return v<_v; }
    T _v;
  };

  template<class T>
  struct InRange
  {
    InRange(T a, T b):_a(a),_b(b) { }
    bool operator()(T v) const { return v>=_a && v<_b; }
    T _a,_b;
  };

template<class T>
struct NotInRange
{
  NotInRange(T a, T b):_a(a),_b(b) { }
  bool operator()(T v) const { return v<_a || v>=_b; }
  T _a,_b;
};

  /*!
   * This method works only on data array with one component. This method returns a newly allocated array storing stored ascendantly of tuple ids in \a this so that this[id]<0.
   *
   * \return a newly allocated data array that the caller should deal with.
   * \sa DataArrayInt::findIdsInRange
   */
  template<class T>
  DataArrayIdType *DataArrayTemplateClassic<T>::findIdsStrictlyNegative() const
  {
    LowerThan<T> lt((T)0);
    MCAuto<DataArrayIdType> ret(findIdsAdv(lt));
    return ret.retn();
  }

  template<class T>
  MCAuto<DataArrayIdType> DataArrayTemplateClassic<T>::findIdsGreaterOrEqualTo(T val) const
  {
    GreatEqual<T> ge(val);
    return findIdsAdv(ge);
  }

  template<class T>
  MCAuto<DataArrayIdType> DataArrayTemplateClassic<T>::findIdsGreaterThan(T val) const
  {
    GreaterThan<T> gt(val);
    return findIdsAdv(gt);
  }

  template<class T>
  MCAuto<DataArrayIdType> DataArrayTemplateClassic<T>::findIdsLowerOrEqualTo(T val) const
  {
    LowerEqual<T> le(val);
    return findIdsAdv(le);
  }

  template<class T>
  MCAuto<DataArrayIdType> DataArrayTemplateClassic<T>::findIdsLowerThan(T val) const
  {
    LowerThan<T> lt(val);
    return findIdsAdv(lt);
  }

  /*!
   * Returns a new DataArrayDouble by aggregating two given arrays, so that (1) the number
   * of components in the result array is a sum of the number of components of given arrays
   * and (2) the number of tuples in the result array is same as that of each of given
   * arrays. In other words the i-th tuple of result array includes all components of
   * i-th tuples of all given arrays.
   * Number of tuples in the given arrays must be  the same.
   *  \param [in] a1 - an array to include in the result array.
   *  \param [in] a2 - another array to include in the result array.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If both \a a1 and \a a2 are NULL.
   *  \throw If any given array is not allocated.
   *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::Meld(const typename Traits<T>::ArrayType *a1, const typename Traits<T>::ArrayType *a2)
  {
    std::vector<const typename Traits<T>::ArrayType *> arr(2);
    arr[0]=a1; arr[1]=a2;
    return Meld(arr);
  }

  /*!
   * Returns a new DataArrayDouble by aggregating all given arrays, so that (1) the number
   * of components in the result array is a sum of the number of components of given arrays
   * and (2) the number of tuples in the result array is same as that of each of given
   * arrays. In other words the i-th tuple of result array includes all components of
   * i-th tuples of all given arrays.
   * Number of tuples in the given arrays must be  the same.
   *  \param [in] arr - a sequence of arrays to include in the result array.
   *  \return DataArrayDouble * - the new instance of DataArrayDouble.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If all arrays within \a arr are NULL.
   *  \throw If any given array is not allocated.
   *  \throw If getNumberOfTuples() of arrays within \a arr is different.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::Meld(const std::vector<const typename Traits<T>::ArrayType *>& arr)
  {
    std::vector<const typename Traits<T>::ArrayType *> a;
    for(typename std::vector<const typename Traits<T>::ArrayType *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
      if(*it4)
        a.push_back(*it4);
    if(a.empty())
      throw INTERP_KERNEL::Exception("DataArrayDouble::Meld : input list must contain at least one NON EMPTY DataArrayDouble !");
    typename std::vector<const typename Traits<T>::ArrayType *>::const_iterator it;
    for(it=a.begin();it!=a.end();it++)
      (*it)->checkAllocated();
    it=a.begin();
    mcIdType nbOfTuples((*it)->getNumberOfTuples());
    std::vector<std::size_t> nbc(a.size());
    std::vector<const T *> pts(a.size());
    nbc[0]=(*it)->getNumberOfComponents();
    pts[0]=(*it++)->getConstPointer();
    for(mcIdType i=1;it!=a.end();it++,i++)
      {
        if(nbOfTuples!=(*it)->getNumberOfTuples())
          throw INTERP_KERNEL::Exception("DataArrayDouble::Meld : mismatch of number of tuples !");
        nbc[i]=(*it)->getNumberOfComponents();
        pts[i]=(*it)->getConstPointer();
      }
    std::size_t totalNbOfComp=std::accumulate(nbc.begin(),nbc.end(),(std::size_t)0);
    typename Traits<T>::ArrayType *ret(Traits<T>::ArrayType::New());
    ret->alloc(nbOfTuples,totalNbOfComp);
    T *retPtr(ret->getPointer());
    for(mcIdType i=0;i<nbOfTuples;i++)
      for(std::size_t j=0;j<a.size();j++)
        {
          retPtr=std::copy(pts[j],pts[j]+nbc[j],retPtr);
          pts[j]+=nbc[j];
        }
    std::size_t k=0;
    for(std::size_t i=0;i<a.size();i++)
      for(std::size_t j=0;j<nbc[i];j++,k++)
        ret->setInfoOnComponent(k,a[i]->getInfoOnComponent(j));
    return ret;
  }

  /*!
   * Returns a new DataArrayDouble holding the same values as \a this array but differently
   * arranged in memory. If \a this array holds 2 components of 3 values:
   * \f$ x_0,x_1,x_2,y_0,y_1,y_2 \f$, then the result array holds these values arranged
   * as follows: \f$ x_0,y_0,x_1,y_1,x_2,y_2 \f$.
   *  \warning Do not confuse this method with transpose()!
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::fromNoInterlace() const
  {
    if(this->_mem.isNull())
      throw INTERP_KERNEL::Exception("DataArrayDouble::fromNoInterlace : Not defined array !");
    T *tab(this->_mem.fromNoInterlace(this->getNumberOfComponents()));
    MCAuto<typename Traits<T>::ArrayType> ret(Traits<T>::ArrayType::New());
    ret->useArray(tab,true,DeallocType::C_DEALLOC,this->getNumberOfTuples(),this->getNumberOfComponents());
    return ret.retn();
  }

  /*!
   * Returns a new DataArrayDouble holding the same values as \a this array but differently
   * arranged in memory. If \a this array holds 2 components of 3 values:
   * \f$ x_0,y_0,x_1,y_1,x_2,y_2 \f$, then the result array holds these values arranged
   * as follows: \f$ x_0,x_1,x_2,y_0,y_1,y_2 \f$.
   *  \warning Do not confuse this method with transpose()!
   *  \return DataArrayDouble * - the new instance of DataArrayDouble that the caller
   *          is to delete using decrRef() as it is no more needed.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::toNoInterlace() const
  {
    if(this->_mem.isNull())
      throw INTERP_KERNEL::Exception("DataArrayDouble::toNoInterlace : Not defined array !");
    T *tab(this->_mem.toNoInterlace(this->getNumberOfComponents()));
    MCAuto<typename Traits<T>::ArrayType> ret(Traits<T>::ArrayType::New());
    ret->useArray(tab,true,DeallocType::C_DEALLOC,this->getNumberOfTuples(),this->getNumberOfComponents());
    return ret.retn();
  }

  /*!
   * Appends components of another array to components of \a this one, tuple by tuple.
   * So that the number of tuples of \a this array remains the same and the number of
   * components increases.
   *  \param [in] other - the DataArrayDouble to append to \a this one.
   *  \throw If \a this is not allocated.
   *  \throw If \a this and \a other arrays have different number of tuples.
   *
   *  \if ENABLE_EXAMPLES
   *  \ref cpp_mcdataarraydouble_meldwith "Here is a C++ example".
   *
   *  \ref py_mcdataarraydouble_meldwith "Here is a Python example".
   *  \endif
   */
  template<class T>
  void DataArrayTemplateClassic<T>::meldWith(const typename Traits<T>::ArrayType *other)
  {
    this->checkAllocated();
    other->checkAllocated();
    mcIdType nbOfTuples(this->getNumberOfTuples());
    if(nbOfTuples!=other->getNumberOfTuples())
      throw INTERP_KERNEL::Exception("DataArrayDouble::meldWith : mismatch of number of tuples !");
    std::size_t nbOfComp1=this->getNumberOfComponents();
    std::size_t nbOfComp2=other->getNumberOfComponents();
    T *newArr=(T *)malloc((nbOfTuples*(nbOfComp1+nbOfComp2))*sizeof(T));
    T *w=newArr;
    const T *inp1(this->begin()),*inp2(other->begin());
    for(mcIdType i=0;i<nbOfTuples;i++,inp1+=nbOfComp1,inp2+=nbOfComp2)
      {
        w=std::copy(inp1,inp1+nbOfComp1,w);
        w=std::copy(inp2,inp2+nbOfComp2,w);
      }
    this->useArray(newArr,true,DeallocType::C_DEALLOC,nbOfTuples,nbOfComp1+nbOfComp2);
    std::vector<std::size_t> compIds(nbOfComp2);
    for(std::size_t i=0;i<nbOfComp2;i++)
      compIds[i]=nbOfComp1+i;
    this->copyPartOfStringInfoFrom2(compIds,*other);
  }

  /*!
   *
   * \param [in] nbTimes specifies the nb of times each tuples in \a this will be duplicated contiguouly in returned DataArrayDouble instance.
   *             \a nbTimes  should be at least equal to 1.
   * \return a newly allocated DataArrayDouble having one component and number of tuples equal to \a nbTimes * \c this->getNumberOfTuples.
   * \throw if \a this is not allocated or if \a this has not number of components set to one or if \a nbTimes is lower than 1.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::duplicateEachTupleNTimes(mcIdType nbTimes) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayDouble::duplicateEachTupleNTimes : this should have only one component !");
    if(nbTimes<1)
      throw INTERP_KERNEL::Exception("DataArrayDouble::duplicateEachTupleNTimes : nb times should be >= 1 !");
    mcIdType nbTuples=this->getNumberOfTuples();
    const T *inPtr(this->begin());
    MCAuto<typename Traits<T>::ArrayType> ret(Traits<T>::ArrayType::New()); ret->alloc(nbTimes*nbTuples,1);
    T *retPtr(ret->getPointer());
    for(mcIdType i=0;i<nbTuples;i++,inPtr++)
      {
        T val(*inPtr);
        for(mcIdType j=0;j<nbTimes;j++,retPtr++)
          *retPtr=val;
      }
    ret->copyStringInfoFrom(*this);
    return ret.retn();
  }

  template<class T>
  MCAuto<typename Traits<T>::ArrayType> DataArrayTemplateClassic<T>::duplicateNTimes(mcIdType nbTimes) const
  {
    this->checkAllocated();
    std::size_t nbCompo(this->getNumberOfComponents());
    mcIdType nbElems(this->getNbOfElems()),nbTuples(this->getNumberOfTuples());
    MCAuto<typename Traits<T>::ArrayType> ret(Traits<T>::ArrayType::New()); ret->alloc(nbTimes*nbTuples,nbCompo);
    ret->copyStringInfoFrom(*this);
    T *retPtr(ret->getPointer());
    const T *inPtr(this->begin());
    for(mcIdType i=0;i<nbTimes;i++)
    {
      retPtr = std::copy(inPtr,inPtr+nbElems,retPtr);
    }
    return ret;
  }

  template<class T>
  void DataArrayTemplateClassic<T>::aggregate(const typename Traits<T>::ArrayType *other)
  {
    if(!other)
      throw INTERP_KERNEL::Exception("DataArrayDouble::aggregate : null pointer !");
    if(this->getNumberOfComponents()!=other->getNumberOfComponents())
      throw INTERP_KERNEL::Exception("DataArrayDouble::aggregate : mismatch number of components !");
    this->_mem.insertAtTheEnd(other->begin(),other->end());
  }

  /*!
   * Converts every value of \a this array to its absolute value.
   * \b WARNING this method is non const. If a new DataArrayDouble instance should be built containing the result of abs DataArrayDouble::computeAbs
   * should be called instead.
   *
   * \throw If \a this is not allocated.
   * \sa DataArrayDouble::computeAbs
   */
  template<class T>
  void DataArrayTemplateClassic<T>::abs()
  {
    this->checkAllocated();
    T *ptr(this->getPointer());
    std::size_t nbOfElems(this->getNbOfElems());
    std::transform(ptr,ptr+nbOfElems,ptr,[](T c){return std::abs(c);});
    this->declareAsNew();
  }

  /*!
   * This method builds a new instance of \a this object containing the result of std::abs applied of all elements in \a this.
   * This method is a const method (that do not change any values in \a this) contrary to  DataArrayDouble::abs method.
   *
   * \return DataArrayDouble * - the new instance of DataArrayDouble containing the
   *         same number of tuples and component as \a this array.
   *         The caller is to delete this result array using decrRef() as it is no more
   *         needed.
   * \throw If \a this is not allocated.
   * \sa DataArrayDouble::abs
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::computeAbs() const
  {
    this->checkAllocated();
    MCAuto<typename Traits<T>::ArrayType> newArr(Traits<T>::ArrayType::New());
    mcIdType nbOfTuples(this->getNumberOfTuples());
    std::size_t nbOfComp(this->getNumberOfComponents());
    newArr->alloc(nbOfTuples,nbOfComp);
    std::transform(this->begin(),this->end(),newArr->getPointer(),[](T c){return std::abs(c);});
    newArr->copyStringInfoFrom(*this);
    return newArr.retn();
  }

  /*!
   * Returns either a \a deep or \a shallow copy of this array. For more info see
   * \ref MEDCouplingArrayBasicsCopyDeep and \ref MEDCouplingArrayBasicsCopyShallow.
   *  \param [in] dCpy - if \a true, a deep copy is returned, else, a shallow one.
   *  \return DataArrayDouble * - either a new instance of DataArrayDouble (if \a dCpy
   *          == \a true) or \a this instance (if \a dCpy == \a false).
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::performCopyOrIncrRef(bool dCpy) const
  {
    const typename Traits<T>::ArrayType *thisC(static_cast<const typename Traits<T>::ArrayType *>(this));
    return DataArrayTemplateClassic<T>::PerformCopyOrIncrRef(dCpy,*thisC);
  }

  /*!
   * Computes for each tuple the sum of number of components values in the tuple and return it.
   *
   * \return DataArrayDouble * - the new instance of DataArrayDouble containing the
   *          same number of tuples as \a this array and one component.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayTemplateClassic<T>::sumPerTuple() const
  {
    this->checkAllocated();
    std::size_t nbOfComp(this->getNumberOfComponents());
    mcIdType nbOfTuple(this->getNumberOfTuples());
    MCAuto<typename Traits<T>::ArrayType> ret(Traits<T>::ArrayType::New());
    ret->alloc(nbOfTuple,1);
    const T *src(this->begin());
    T *dest(ret->getPointer());
    for(mcIdType i=0;i<nbOfTuple;i++,dest++,src+=nbOfComp)
      *dest=std::accumulate(src,src+nbOfComp,(T)0);
    return ret.retn();
  }

  /*!
   * Set all values in \a this array so that the i-th element equals to \a init + i
   * (i starts from zero). To know more on filling arrays see \ref MEDCouplingArrayFill.
   *  \param [in] init - value to assign to the first element of array.
   *  \throw If \a this->getNumberOfComponents() != 1
   *  \throw If \a this is not allocated.
   */
  template<class T>
  void DataArrayTemplateClassic<T>::iota(T init)
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayDouble::iota : works only for arrays with only one component, you can call 'rearrange' method before !");
    T *ptr(this->getPointer());
    mcIdType ntuples(this->getNumberOfTuples());
    for(mcIdType i=0;i<ntuples;i++)
      ptr[i]=init+(T)i;
    this->declareAsNew();
  }

  template<class T>
  struct ImplReprTraits { static void SetPrecision(std::ostream& /*oss*/) { } };

  template<>
  struct ImplReprTraits<double> {  static void SetPrecision(std::ostream& oss) { oss.precision(17); } };

  template<>
  struct ImplReprTraits<float> {  static void SetPrecision(std::ostream& oss) { oss.precision(7); } };

  template<class T>
  void DataArrayTemplateClassic<T>::reprStream(std::ostream& stream) const
  {
    stream << "Name of " << Traits<T>::ReprStr << " array : \"" << this->_name << "\"\n";
    reprWithoutNameStream(stream);
  }

  template<class T>
  void DataArrayTemplateClassic<T>::reprZipStream(std::ostream& stream) const
  {
    stream << "Name of " << Traits<T>::ReprStr << " array : \"" << this->_name << "\"\n";
    reprZipWithoutNameStream(stream);
  }

  template<class T>
  void DataArrayTemplateClassic<T>::reprNotTooLongStream(std::ostream& stream) const
  {
    stream << "Name of "<< Traits<T>::ReprStr << " array : \"" << this->_name << "\"\n";
    reprNotTooLongWithoutNameStream(stream);
  }

  template<class T>
  void DataArrayTemplateClassic<T>::reprWithoutNameStream(std::ostream& stream) const
  {
    DataArray::reprWithoutNameStream(stream);
    ImplReprTraits<T>::SetPrecision(stream);
    this->_mem.repr(ToIdType(this->getNumberOfComponents()),stream);
  }

  template<class T>
  void DataArrayTemplateClassic<T>::reprZipWithoutNameStream(std::ostream& stream) const
  {
    DataArray::reprWithoutNameStream(stream);
    ImplReprTraits<T>::SetPrecision(stream);
    this->_mem.reprZip(ToIdType(this->getNumberOfComponents()),stream);
  }

  template<class T>
  void DataArrayTemplateClassic<T>::reprNotTooLongWithoutNameStream(std::ostream& stream) const
  {
    DataArray::reprWithoutNameStream(stream);
    ImplReprTraits<T>::SetPrecision(stream);
    this->_mem.reprNotTooLong(ToIdType(this->getNumberOfComponents()),stream);
  }

  /*!
   * This method is close to repr method except that when \a this has more than 1000 tuples, all tuples are not
   * printed out to avoid to consume too much space in interpretor.
   * \sa repr
   */
  template<class T>
  std::string DataArrayTemplateClassic<T>::reprNotTooLong() const
  {
    std::ostringstream ret;
    reprNotTooLongStream(ret);
    return ret.str();
  }

  /*!
   * Returns a textual and human readable representation of \a this instance of
   * DataArrayInt. This text is shown when a DataArrayInt is printed in Python.
   * \return std::string - text describing \a this DataArrayInt.
   *
   * \sa reprNotTooLong, reprZip
   */
  template<class T>
  std::string DataArrayTemplateClassic<T>::repr() const
  {
    std::ostringstream ret;
    DataArrayTemplateClassic<T>::reprStream(ret);
    return ret.str();
  }

  template<class T>
  std::string DataArrayTemplateClassic<T>::reprZip() const
  {
    std::ostringstream ret;
    DataArrayTemplateClassic<T>::reprZipStream(ret);
    return ret.str();
  }

  /////////////////////////////////

  /*!
   * Checks if all values in \a this array are equal to \a val at precision \a eps.
   *  \param [in] val - value to check equality of array values to.
   *  \param [in] eps - precision to check the equality.
   *  \return bool - \a true if all values are in range (_val_ - _eps_; _val_ + _eps_),
   *                 \a false else.
   *  \throw If \a this->getNumberOfComponents() != 1
   *  \throw If \a this is not allocated.
   */
  template<class T>
  bool DataArrayTemplateFP<T>::isUniform(T val, T eps) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayDouble::isUniform : must be applied on DataArrayDouble with only one component, you can call 'rearrange' method before !");
    const T *w(this->begin()),*end2(this->end());
    const T vmin(val-eps),vmax(val+eps);
    for(;w!=end2;w++)
      if(*w<vmin || *w>vmax)
        return false;
    return true;
  }

  /////////////////////////////////

  /*!
   * Returns the only one value in \a this, if and only if number of elements
   * (nb of tuples * nb of components) is equal to 1, and that \a this is allocated.
   *  \return double - the sole value stored in \a this array.
   *  \throw If at least one of conditions stated above is not fulfilled.
   */
  template <class T>
  T DataArrayDiscrete<T>::intValue() const
  {
    if(this->isAllocated())
      {
        if(this->getNbOfElems()==1)
          {
            return *this->getConstPointer();
          }
        else
          throw INTERP_KERNEL::Exception("DataArrayInt::intValue : DataArrayInt instance is allocated but number of elements is not equal to 1 !");
      }
    else
      throw INTERP_KERNEL::Exception("DataArrayInt::intValue : DataArrayInt instance is not allocated !");
  }

  /*!
   * Equivalent to DataArrayInt::isEqual except that if false the reason of
   * mismatch is given.
   *
   * \param [in] other the instance to be compared with \a this
   * \param [out] reason In case of inequality returns the reason.
   * \sa DataArrayInt::isEqual
   */
  template<class T>
  bool DataArrayDiscrete<T>::isEqualIfNotWhy(const DataArrayDiscrete<T>& other, std::string& reason) const
  {
    if(!this->areInfoEqualsIfNotWhy(other,reason))
      return false;
    return this->_mem.isEqual(other._mem,0,reason);
  }

  /*!
   * Checks if \a this and another DataArrayInt are fully equal. For more info see
   * \ref MEDCouplingArrayBasicsCompare.
   *  \param [in] other - an instance of DataArrayInt to compare with \a this one.
   *  \return bool - \a true if the two arrays are equal, \a false else.
   */
  template<class T>
  bool DataArrayDiscrete<T>::isEqual(const DataArrayDiscrete<T>& other) const
  {
    std::string tmp;
    return isEqualIfNotWhy(other,tmp);
  }

  /*!
   * Returns a new instance of DataArrayInt. The caller is to delete this array
   * using decrRef() as it is no more needed.
   */
  template<class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::New()
  {
    return new typename Traits<T>::ArrayType;
  }

  /*!
   * Checks if values of \a this and another DataArrayInt are equal. For more info see
   * \ref MEDCouplingArrayBasicsCompare.
   *  \param [in] other - an instance of DataArrayInt to compare with \a this one.
   *  \return bool - \a true if the values of two arrays are equal, \a false else.
   */
  template<class T>
  bool DataArrayDiscrete<T>::isEqualWithoutConsideringStr(const DataArrayDiscrete<T>& other) const
  {
    std::string tmp;
    return this->_mem.isEqual(other._mem,0,tmp);
  }

  /*!
   * Checks if values of \a this and another DataArrayInt are equal. Comparison is
   * performed on sorted value sequences.
   * For more info see\ref MEDCouplingArrayBasicsCompare.
   *  \param [in] other - an instance of DataArrayInt to compare with \a this one.
   *  \return bool - \a true if the sorted values of two arrays are equal, \a false else.
   */
  template<class T>
  bool DataArrayDiscrete<T>::isEqualWithoutConsideringStrAndOrder(const typename Traits<T>::ArrayType& other) const
  {
    MCAuto<typename Traits<T>::ArrayType> a((static_cast<const typename Traits<T>::ArrayType *>(this))->deepCopy());
    MCAuto<typename Traits<T>::ArrayType> b((static_cast<const typename Traits<T>::ArrayType *>(&other))->deepCopy());
    a->sort();
    b->sort();
    return a->isEqualWithoutConsideringStr(*b);
  }

  template<class T>
  template<class ALG>
  void DataArrayDiscrete<T>::switchOnTupleAlg(T val, std::vector<bool>& vec, ALG algo) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::switchOnTupleEqualTo : number of components of this should be equal to one !");
    mcIdType nbOfTuples(this->getNumberOfTuples());
    if(nbOfTuples!=ToIdType(vec.size()))
      throw INTERP_KERNEL::Exception("DataArrayInt::switchOnTupleEqualTo : number of tuples of this should be equal to size of input vector of bool !");
    const T *pt(this->begin());
    for(mcIdType i=0;i<nbOfTuples;i++)
      if(algo(pt[i],val))
        vec[i]=true;
  }

  /*!
   * This method assumes that \a this has one component and is allocated. This method scans all tuples in \a this and for all tuple equal to \a val
   * put True to the corresponding entry in \a vec.
   * \a vec is expected to be with the same size than the number of tuples of \a this.
   *
   *  \sa DataArrayInt::switchOnTupleNotEqualTo.
   */
  template<class T>
  void DataArrayDiscrete<T>::switchOnTupleEqualTo(T val, std::vector<bool>& vec) const
  {
    switchOnTupleAlg(val,vec,std::equal_to<T>());
  }

  /*!
   * This method assumes that \a this has one component and is allocated. This method scans all tuples in \a this and for all tuple different from \a val
   * put True to the corresponding entry in \a vec.
   * \a vec is expected to be with the same size than the number of tuples of \a this.
   *
   *  \sa DataArrayInt::switchOnTupleEqualTo.
   */
  template<class T>
  void DataArrayDiscrete<T>::switchOnTupleNotEqualTo(T val, std::vector<bool>& vec) const
  {
    switchOnTupleAlg(val,vec,std::not_equal_to<T>());
  }

  /*!
   * Compute for each element in \a this the occurence rank of that element. This method is typically useful of one-component array having a same element
   * appearing several times. If each element in \a this appears once an 1 component array containing only 0 will be returned.
   *
   * \b Example:
   * - \a this : [5, 3, 2, 1, 4, 5, 2, 1, 0, 11, 5, 4]
   * - \a return is : [0, 0, 0, 0, 0, 1, 1, 1, 0, 0, 2, 1] because at pos #0 of \a this (ie value 5) is the first occurrence ->0. At pos #10 of \a this (ie value 5 also) is the third occurrence of 5 -> 2.
   *
   * \return DataArrayInt * - a new instance of DataArrayInt with same number of tuples than \a this. The caller is to delete this
   *          array using decrRef() as it is no more needed.
   * \throw If either this not allocated or not with one component.
   *
   * \sa DataArrayInt::FindPermutationFromFirstToSecond
   */
  template<class T>
  DataArrayIdType *DataArrayDiscrete<T>::occurenceRankInThis() const
  {
    constexpr char MSG0[] = "occurenceRankInThis :";
    this->checkAllocated();
    this->checkNbOfComps(1,MSG0);
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New());
    ret->alloc(this->getNumberOfTuples(),1);
    mcIdType *retPtr(ret->getPointer());
    std::map<T,mcIdType> m;
    for(const T *pt = this->begin() ; pt != this->end() ; ++pt, ++retPtr )
    {
      auto it = m.find(*pt);
      if( it == m.end() )
      {
        *retPtr = 0;
        m[*pt] = 1;
      }
      else
      {
        *retPtr = (*it).second++;
      }
    }
    return ret.retn();
  }

  /*!
   * Creates a new one-dimensional DataArrayInt of the same size as \a this and a given
   * one-dimensional arrays that must be of the same length. The result array describes
   * correspondence between \a this and \a other arrays, so that
   * <em> other.getIJ(i,0) == this->getIJ(ret->getIJ(i),0)</em>. If such a permutation is
   * not possible because some element in \a other is not in \a this, an exception is thrown.
   *  \param [in] other - an array to compute permutation to.
   *  \return DataArrayInt * - a new instance of DataArrayInt, which is a permutation array
   * from \a this to \a other. The caller is to delete this array using decrRef() as it is
   * no more needed.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If \a other->getNumberOfComponents() != 1.
   *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples().
   *  \throw If \a other includes a value which is not in \a this array.
   *
   *  \if ENABLE_EXAMPLES
   *  \ref cpp_mcdataarrayint_buildpermutationarr "Here is a C++ example".
   *
   *  \ref py_mcdataarrayint_buildpermutationarr "Here is a Python example".
   *  \endif
   */
  template<class T>
  DataArrayIdType *DataArrayDiscrete<T>::buildPermutationArr(const DataArrayDiscrete<T>& other) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1 || other.getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildPermutationArr : 'this' and 'other' have to have exactly ONE component !");
    mcIdType nbTuple(this->getNumberOfTuples());
    other.checkAllocated();
    if(nbTuple!=other.getNumberOfTuples())
      throw INTERP_KERNEL::Exception("DataArrayInt::buildPermutationArr : 'this' and 'other' must have the same number of tuple !");
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New());
    ret->alloc(nbTuple,1);
    ret->fillWithValue(-1);
    const T *pt(this->begin());
    std::map<mcIdType,mcIdType> mm;
    for(mcIdType i=0;i<nbTuple;i++)
      mm[ToIdType(pt[i])]=i;
    pt=other.begin();
    mcIdType *retToFill(ret->getPointer());
    for(mcIdType i=0;i<nbTuple;i++)
      {
        std::map<mcIdType,mcIdType>::const_iterator it=mm.find(ToIdType(pt[i]));
        if(it==mm.end())
          {
            std::ostringstream oss; oss << "DataArrayInt::buildPermutationArr : Arrays mismatch : element (" << pt[i] << ") in 'other' not findable in 'this' !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
        retToFill[i]=(*it).second;
      }
    return ret.retn();
  }

  /*!
   * Elements of \a partOfThis are expected to be included in \a this.
   * The returned array \a ret is so that this[ret]==partOfThis
   *
   * For example, if \a this array contents are [9,10,0,6,4,11,3,8] and if \a partOfThis contains [6,0,11,8]
   * the return array will contain [3,2,5,7].
   *
   * \a this is expected to be a 1 compo allocated array.
   * \param [in] partOfThis - A 1 compo allocated array
   * \return - A newly allocated array to be dealed by caller having the same number of tuples than \a partOfThis.
   * \throw if two same element is present twice in \a this
   * \throw if an element in \a partOfThis is \b NOT in \a this.
   */
  template<class T>
  DataArrayIdType *DataArrayDiscrete<T>::indicesOfSubPart(const DataArrayDiscrete<T>& partOfThis) const
  {
    if(this->getNumberOfComponents()!=1 || partOfThis.getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::indicesOfSubPart : this and input array must be one component array !");
    this->checkAllocated(); partOfThis.checkAllocated();
    mcIdType thisNbTuples(this->getNumberOfTuples()),nbTuples(partOfThis.getNumberOfTuples());
    const T *thisPt(this->begin()),*pt(partOfThis.begin());
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New());
    ret->alloc(nbTuples,1);
    mcIdType *retPt(ret->getPointer());
    std::map<mcIdType,mcIdType> m;
    for(mcIdType i=0;i<thisNbTuples;i++,thisPt++)
      m[ToIdType(*thisPt)]=i;
    if(ToIdType(m.size())!=thisNbTuples)
      throw INTERP_KERNEL::Exception("DataArrayInt::indicesOfSubPart : some elements appears more than once !");
    for(mcIdType i=0;i<nbTuples;i++,retPt++,pt++)
      {
        std::map<mcIdType,mcIdType>::const_iterator it(m.find(ToIdType(*pt)));
        if(it!=m.end())
          *retPt=(*it).second;
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::indicesOfSubPart : At pos #" << i << " of input array value is " << *pt << " not in this !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
    return ret.retn();
  }

  /*!
   * Checks that \a this array is consistently **increasing** or **decreasing** in value.
   * If not an exception is thrown.
   *  \param [in] increasing - if \a true, the array values should be increasing.
   *  \throw If sequence of values is not strictly monotonic in agreement with \a
   *         increasing arg.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  void DataArrayDiscrete<T>::checkMonotonic(bool increasing) const
  {
    if(!isMonotonic(increasing))
      {
        if (increasing)
          throw INTERP_KERNEL::Exception("DataArrayInt::checkMonotonic : 'this' is not INCREASING monotonic !");
        else
          throw INTERP_KERNEL::Exception("DataArrayInt::checkMonotonic : 'this' is not DECREASING monotonic !");
      }
  }

  /*!
   * Checks that \a this array is consistently **increasing** or **decreasing** in value.
   *  \param [in] increasing - if \a true, array values should be increasing.
   *  \return bool - \a true if values change in accordance with \a increasing arg.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  bool DataArrayDiscrete<T>::isMonotonic(bool increasing) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::isMonotonic : only supported with 'this' array with ONE component !");
    std::size_t nbOfElements(this->getNumberOfTuples());
    const T *ptr(this->begin());
    if(nbOfElements==0)
      return true;
    T ref(ptr[0]);
    if(increasing)
      {
        for(std::size_t i=1;i<nbOfElements;i++)
          {
            if(ptr[i]>=ref)
              ref=ptr[i];
            else
              return false;
          }
      }
    else
      {
        for(std::size_t i=1;i<nbOfElements;i++)
          {
            if(ptr[i]<=ref)
              ref=ptr[i];
            else
              return false;
          }
      }
    return true;
  }

  /*!
   * This method check that array consistently INCREASING or DECREASING in value.
   */
  template<class T>
  bool DataArrayDiscrete<T>::isStrictlyMonotonic(bool increasing) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::isStrictlyMonotonic : only supported with 'this' array with ONE component !");
    std::size_t nbOfElements(this->getNumberOfTuples());
    const T *ptr(this->begin());
    if(nbOfElements==0)
      return true;
    T ref(ptr[0]);
    if(increasing)
      {
        for(std::size_t i=1;i<nbOfElements;i++)
          {
            if(ptr[i]>ref)
              ref=ptr[i];
            else
              return false;
          }
      }
    else
      {
        for(std::size_t i=1;i<nbOfElements;i++)
          {
            if(ptr[i]<ref)
              ref=ptr[i];
            else
              return false;
          }
      }
    return true;
  }

  /*!
   * This method check that array consistently INCREASING or DECREASING in value.
   */
  template<class T>
  void DataArrayDiscrete<T>::checkStrictlyMonotonic(bool increasing) const
  {
    if(!isStrictlyMonotonic(increasing))
      {
        if (increasing)
          throw INTERP_KERNEL::Exception("DataArrayInt::checkStrictlyMonotonic : 'this' is not strictly INCREASING monotonic !");
        else
          throw INTERP_KERNEL::Exception("DataArrayInt::checkStrictlyMonotonic : 'this' is not strictly DECREASING monotonic !");
      }
  }

  /*!
   * Returns an integer value characterizing \a this array, which is useful for a quick
   * comparison of many instances of DataArrayInt.
   *  \return mcIdType - the hash value.
   *  \throw If \a this is not allocated.
   */
  template<class T>
  mcIdType DataArrayDiscrete<T>::getHashCode() const
  {
    this->checkAllocated();
    mcIdType nbOfElems=ToIdType(this->getNbOfElems());
    mcIdType ret=nbOfElems*65536;
    mcIdType delta=3;
    if(nbOfElems>48)
      delta=nbOfElems/8;
    T ret0(0);
    const T *pt(this->begin());
    for(mcIdType i=0;i<nbOfElems;i+=delta)
      ret0+=pt[i] & 0x1FFF;
    return ToIdType(ret+ret0);
  }

  template<class T>
  void DataArrayDiscrete<T>::reprCppStream(const std::string& varName, std::ostream& stream) const
  {
    mcIdType nbTuples(this->getNumberOfTuples());
    std::size_t nbComp(this->getNumberOfComponents());
    const T *data(this->getConstPointer());
    stream << Traits<T>::ArrayTypeName << " *" << varName << "=" << Traits<T>::ArrayTypeName << "::New();" << std::endl;
    if(nbTuples*nbComp>=1)
      {
        stream << "const mcIdType " << varName << "Data[" << nbTuples*nbComp << "]={";
        std::copy(data,data+nbTuples*nbComp-1,std::ostream_iterator<T>(stream,","));
        stream << data[nbTuples*nbComp-1] << "};" << std::endl;
        stream << varName << "->useArray(" << varName << "Data,false,CPP_DEALLOC," << nbTuples << "," << nbComp << ");" << std::endl;
      }
    else
      stream << varName << "->alloc(" << nbTuples << "," << nbComp << ");" << std::endl;
    stream << varName << "->setName(\"" << this->getName() << "\");" << std::endl;
  }

  /*!
   * Method that gives a quick overvien of \a this for python.
   */
  template<class T>
  void DataArrayDiscrete<T>::reprQuickOverview(std::ostream& stream) const
  {
    static const std::size_t MAX_NB_OF_BYTE_IN_REPR=300;
    stream << Traits<T>::ArrayTypeName << " C++ instance at " << this << ". ";
    if(this->isAllocated())
      {
        std::size_t nbOfCompo(this->getNumberOfComponents());
        if(nbOfCompo>=1)
          {
            mcIdType nbOfTuples(this->getNumberOfTuples());
            stream << "Number of tuples : " << nbOfTuples << ". Number of components : " << nbOfCompo << "." << std::endl;
            reprQuickOverviewData(stream,MAX_NB_OF_BYTE_IN_REPR);
          }
        else
          stream << "Number of components : 0.";
      }
    else
      stream << "*** No data allocated ****";
  }

  template<class T>
  void DataArrayDiscrete<T>::reprQuickOverviewData(std::ostream& stream, std::size_t maxNbOfByteInRepr) const
  {
    const T *data(this->begin());
    mcIdType nbOfTuples(this->getNumberOfTuples());
    std::size_t nbOfCompo(this->getNumberOfComponents());
    std::ostringstream oss2; oss2 << "[";
    std::string oss2Str(oss2.str());
    bool isFinished=true;
    for(mcIdType i=0;i<nbOfTuples && isFinished;i++)
      {
        if(nbOfCompo>1)
          {
            oss2 << "(";
            for(std::size_t j=0;j<nbOfCompo;j++,data++)
              {
                oss2 << *data;
                if(j!=nbOfCompo-1) oss2 << ", ";
              }
            oss2 << ")";
          }
        else
          oss2 << *data++;
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

  template<class T>
  void DataArrayDiscrete<T>::writeVTK(std::ostream& ofs, mcIdType indent, const std::string& type, const std::string& nameInFile, DataArrayByte *byteArr) const
  {
    static const char SPACE[4]={' ',' ',' ',' '};
    this->checkAllocated();
    std::string idt(indent,' ');
    ofs << idt << "<DataArray type=\"" << type << "\" Name=\"" << nameInFile << "\" NumberOfComponents=\"" << this->getNumberOfComponents() << "\"";
    if(byteArr)
      {
        ofs << " format=\"appended\" offset=\"" << byteArr->getNumberOfTuples() << "\">";
        if(std::string(type)==Traits<T>::VTKReprStr)
          {
            const char *data(reinterpret_cast<const char *>(this->begin()));
            std::size_t sz(this->getNbOfElems()*sizeof(T));
            byteArr->insertAtTheEnd(data,data+sz);
            byteArr->insertAtTheEnd(SPACE,SPACE+4);
          }
        else if(std::string(type)=="Int8")
          {
            INTERP_KERNEL::AutoPtr<char> tmp(new char[this->getNbOfElems()]);
            copyCast(this->begin(),this->end(),(char *)tmp);
            byteArr->insertAtTheEnd((char *)tmp,(char *)tmp+this->getNbOfElems());
            byteArr->insertAtTheEnd(SPACE,SPACE+4);
          }
        else if(std::string(type)=="UInt8")
          {
            INTERP_KERNEL::AutoPtr<unsigned char> tmp(new unsigned char[this->getNbOfElems()]);
            copyCast(this->begin(),this->end(),(unsigned char *)tmp);
            byteArr->insertAtTheEnd((unsigned char *)tmp,(unsigned char *)tmp+this->getNbOfElems());
            byteArr->insertAtTheEnd(SPACE,SPACE+4);
          }
        else
          {
            std::ostringstream oss;
            oss << Traits<T>::ArrayTypeName << "::writeVTK : Only " << Traits<T>::VTKReprStr << ", Int8 and UInt8 supported !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
    else
      {
        ofs << " RangeMin=\"" << this->getMinValueInArray() << "\" RangeMax=\"" << this->getMaxValueInArray() << "\" format=\"ascii\">\n" << idt;
        std::copy(this->begin(),this->end(),std::ostream_iterator<T>(ofs," "));
      }
      ofs << std::endl << idt << "</DataArray>\n";
  }

  /*!
   * Modifies in place \a this one-dimensional array so that each value \a v = \a indArrBg[ \a v ],
   * i.e. a current value is used as in index to get a new value from \a indArrBg.
   *  \param [in] indArrBg - pointer to the first element of array of new values to assign
   *         to \a this array.
   *  \param [in] indArrEnd - specifies the end of the array \a indArrBg, so that
   *              the last value of \a indArrBg is \a indArrEnd[ -1 ].
   *  \throw If \a this->getNumberOfComponents() != 1
   *  \throw If any value of \a this can't be used as a valid index for
   *         [\a indArrBg, \a indArrEnd).
   *
   *  \sa changeValue, findIdForEach
   */
  template<class T>
  void DataArrayDiscrete<T>::transformWithIndArr(const T *indArrBg, const T *indArrEnd)
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("Call transformWithIndArr method on DataArrayInt with only one component, you can call 'rearrange' method before !");
    mcIdType nbElemsIn=ToIdType(std::distance(indArrBg,indArrEnd));
    mcIdType nbOfTuples(this->getNumberOfTuples());
    T *pt(this->getPointer());
    for(mcIdType i=0;i<nbOfTuples;i++,pt++)
      {
        if(*pt>=0 && *pt<nbElemsIn)
          *pt=indArrBg[*pt];
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::transformWithIndArr : error on tuple #" << i << " of this value is " << *pt << ", should be in [0," << nbElemsIn << ") !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
    this->declareAsNew();
  }

  template<class T>
  void DataArrayDiscrete<T>::transformWithIndArr(const MapKeyVal<T, T>& m)
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("Call transformWithIndArr method on DataArrayInt with only one component, you can call 'rearrange' method before !");
    const typename std::map<T,T>& dat(m.data());
    mcIdType nbOfTuples(this->getNumberOfTuples());
    T *pt(this->getPointer());
    for(mcIdType i=0;i<nbOfTuples;i++,pt++)
      {
        typename std::map<T,T>::const_iterator it(dat.find(*pt));
        if(it!=dat.end())
          *pt=(*it).second;
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::transformWithIndArr : error on tuple #" << i << " of this value is " << *pt << " not in map !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
    this->declareAsNew();
  }

  template<class T>
  template<int SPACEDIM>
  void DataArrayDiscrete<T>::findCommonTuplesAlg(const T *bbox, mcIdType nbNodes, mcIdType limitNodeId, DataArrayIdType *c, DataArrayIdType *cI) const
  {
    const T *coordsPtr(this->begin());
    BBTreeDiscrete<SPACEDIM,T,mcIdType> myTree(bbox,nullptr,0,nbNodes);
    std::vector<bool> isDone(nbNodes);
    for(mcIdType i=0;i<nbNodes;i++)
    {
      if(!isDone[i])
        {
          std::vector<mcIdType> intersectingElems;
          myTree.getIntersectingElems(coordsPtr+i*SPACEDIM,intersectingElems);
          if(intersectingElems.size()>1)
            {
              std::vector<mcIdType> commonNodes;
              for(std::vector<mcIdType>::const_iterator it=intersectingElems.begin();it!=intersectingElems.end();it++)
                if(*it!=i)
                  if(*it>=limitNodeId)
                    {
                      commonNodes.push_back(*it);
                      isDone[*it]=true;
                    }
              if(!commonNodes.empty())
                {
                  cI->pushBackSilent(cI->back()+ToIdType(commonNodes.size())+1);
                  c->pushBackSilent(i);
                  c->insertAtTheEnd(commonNodes.begin(),commonNodes.end());
                }
            }
        }
    }
  }

  template<class T>
  void DataArrayDiscrete<T>::findCommonTuples(mcIdType limitTupleId, MCAuto<DataArrayIdType> &comm, MCAuto<DataArrayIdType>& commIndex) const
  {
    this->checkAllocated();
    std::size_t nbOfCompo( this->getNumberOfComponents() );
    if ((nbOfCompo<1) || (nbOfCompo>4)) //test before work
      throw INTERP_KERNEL::Exception("DataArrayDiscrete::findCommonTuples : Unexpected spacedim of coords. Must be 1, 2, 3 or 4.");

    mcIdType nbOfTuples( this->getNumberOfTuples() );
    //
    comm = DataArrayIdType::New(); commIndex = DataArrayIdType::New(); comm->alloc(0,1); commIndex->pushBackSilent(0);
    switch(nbOfCompo)
    {
      case 4:
        findCommonTuplesAlg<4>(this->begin(),nbOfTuples,limitTupleId,comm,commIndex);
        break;
      case 3:
        findCommonTuplesAlg<3>(this->begin(),nbOfTuples,limitTupleId,comm,commIndex);
        break;
      case 2:
        findCommonTuplesAlg<2>(this->begin(),nbOfTuples,limitTupleId,comm,commIndex);
        break;
      case 1:
        findCommonTuplesAlg<1>(this->begin(),nbOfTuples,limitTupleId,comm,commIndex);
        break;
      default:
        throw INTERP_KERNEL::Exception("DataArrayDiscrete::findCommonTuples : nb of components managed are 1,2,3 and 4 ! not implemented for other number of components !");
    }
  }

  /*!
   * Creates a new DataArrayInt containing IDs (indices) of tuples holding value equal to a
   * given one. The ids are sorted in the ascending order.
   *  \param [in] val - the value to find within \a this.
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *          array using decrRef() as it is no more needed.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \sa DataArrayInt::findIdsEqualTuple
   */
  template<class T>
  DataArrayIdType *DataArrayDiscrete<T>::findIdsEqual(T val) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdsEqual : the array must have only one component, you can call 'rearrange' method before !");
    const T *cptr(this->getConstPointer());
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(0,1);
    mcIdType nbOfTuples(this->getNumberOfTuples());
    for(mcIdType i=0;i<nbOfTuples;i++,cptr++)
      if(*cptr==val)
        ret->pushBackSilent(ToIdType(i));
    return ret.retn();
  }

  /*!
   * Creates a one-dimensional DataArrayInt (\a res) whose contents are computed from
   * values of \a this (\a a) and the given (\a indArr) arrays as follows:
   * \a res[ \a indArr[ \a a[ i ]]] = i. I.e. for each value in place i \a v = \a a[ i ],
   * new value in place \a indArr[ \a v ] is i.
   *  \param [in] indArrBg - the array holding indices within the result array to assign
   *         indices of values of \a this array pointing to values of \a indArrBg.
   *  \param [in] indArrEnd - specifies the end of the array \a indArrBg, so that
   *              the last value of \a indArrBg is \a indArrEnd[ -1 ].
   *  \return DataArrayInt * - the new instance of DataArrayInt.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If any value of \a this array is not a valid index for \a indArrBg array.
   *  \throw If any value of \a indArrBg is not a valid index for \a this array.
   */
  template<class T>
  DataArrayIdType *DataArrayDiscrete<T>::transformWithIndArrR(const T *indArrBg, const T *indArrEnd) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("Call transformWithIndArrR method on DataArrayInt with only one component, you can call 'rearrange' method before !");
    mcIdType nbElemsIn=ToIdType(std::distance(indArrBg,indArrEnd));
    mcIdType nbOfTuples(this->getNumberOfTuples());
    const T *pt=this->getConstPointer();
    MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
    ret->alloc(nbOfTuples,1);
    ret->fillWithValue(-1);
    mcIdType *tmp=ret->getPointer();
    for(mcIdType i=0;i<nbOfTuples;i++,pt++)
      {
        if(*pt>=0 && *pt<nbElemsIn)
          {
            T pos=indArrBg[*pt];
            if(pos>=0 && pos<nbOfTuples)
              tmp[ToIdType(pos)]=i;
            else
              {
                std::ostringstream oss; oss << "DataArrayInt::transformWithIndArrR : error on tuple #" << i << " value of new pos is " << pos << " ( indArrBg[" << *pt << "]) ! Should be in [0," << nbOfTuples << ") !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::transformWithIndArrR : error on tuple #" << i << " value is " << *pt << " and indirectionnal array as a size equal to " << nbElemsIn << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    return ret.retn();
  }

  /*!
   * Computes distribution of values of \a this one-dimensional array between given value
   * ranges (casts). This method is typically useful for entity number splitting by types,
   * for example.
   *  \warning The values contained in \a arrBg should be sorted ascendently. No
   *           check of this is be done. If not, the result is not warranted.
   *  \param [in] arrBg - the array of ascending values defining the value ranges. The i-th
   *         value of \a arrBg (\a arrBg[ i ]) gives the lowest value of the i-th range,
   *         and the greatest value of the i-th range equals to \a arrBg[ i+1 ] - 1. \a
   *         arrBg containing \a n values defines \a n-1 ranges. The last value of \a arrBg
   *         should be more than every value in \a this array.
   *  \param [in] arrEnd - specifies the end of the array \a arrBg, so that
   *              the last value of \a arrBg is \a arrEnd[ -1 ].
   *  \param [out] castArr - a new instance of DataArrayInt, of same size as \a this array
   *         (same number of tuples and components), the caller is to delete
   *         using decrRef() as it is no more needed.
   *         This array contains indices of ranges for every value of \a this array. I.e.
   *         the i-th value of \a castArr gives the index of range the i-th value of \a this
   *         belongs to. Or, in other words, this parameter contains for each tuple in \a
   *         this in which cast it holds.
   *  \param [out] rankInsideCast - a new instance of DataArrayInt, of same size as \a this
   *         array, the caller is to delete using decrRef() as it is no more needed.
   *         This array contains ranks of values of \a this array within ranges
   *         they belongs to. I.e. the i-th value of \a rankInsideCast gives the rank of
   *         the i-th value of \a this array within the \a castArr[ i ]-th range, to which
   *         the i-th value of \a this belongs to. Or, in other words, this param contains
   *         for each tuple its rank inside its cast. The rank is computed as difference
   *         between the value and the lowest value of range.
   *  \param [out] castsPresent - a new instance of DataArrayInt, containing indices of
   *         ranges (casts) to which at least one value of \a this array belongs.
   *         Or, in other words, this param contains the casts that \a this contains.
   *         The caller is to delete this array using decrRef() as it is no more needed.
   *
   * \b Example: If \a this contains [6,5,0,3,2,7,8,1,4] and \a arrBg contains [0,4,9] then
   *            the output of this method will be :
   * - \a castArr       : [1,1,0,0,0,1,1,0,1]
   * - \a rankInsideCast: [2,1,0,3,2,3,4,1,0]
   * - \a castsPresent  : [0,1]
   *
   * I.e. values of \a this array belong to 2 ranges: #0 and #1. Value 6 belongs to the
   * range #1 and its rank within this range is 2; etc.
   *
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If \a arrEnd - arrBg < 2.
   *  \throw If any value of \a this is not less than \a arrEnd[-1].
   */
  template<class T>
  void DataArrayDiscrete<T>::splitByValueRange(const T *arrBg, const T *arrEnd,
                                               DataArrayType *& castArr, DataArrayType *& rankInsideCast, DataArrayType *& castsPresent) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("Call splitByValueRange  method on DataArrayInt with only one component, you can call 'rearrange' method before !");
    mcIdType nbOfTuples=this->getNumberOfTuples();
    std::size_t nbOfCast=std::distance(arrBg,arrEnd);
    if(nbOfCast<2)
      throw INTERP_KERNEL::Exception("DataArrayInt::splitByValueRange : The input array giving the cast range values should be of size >=2 !");
    nbOfCast--;
    const T *work=this->getConstPointer();
    typedef std::reverse_iterator<const T *> rintstart;
    rintstart bg(arrEnd);//OK no problem because size of 'arr' is greater or equal 2
    rintstart end2(arrBg);
    MCAuto<DataArrayType> ret1=DataArrayType::New();
    MCAuto<DataArrayType> ret2=DataArrayType::New();
    MCAuto<DataArrayType> ret3=DataArrayType::New();
    ret1->alloc(nbOfTuples,1);
    ret2->alloc(nbOfTuples,1);
    T *ret1Ptr=ret1->getPointer();
    T *ret2Ptr=ret2->getPointer();
    std::set<T> castsDetected;
    for(mcIdType i=0;i<nbOfTuples;i++)
      {
        rintstart res=std::find_if(bg,end2,std::bind(std::less_equal<T>(),std::placeholders::_1,work[i]));
        std::size_t pos=std::distance(bg,res);
        std::size_t pos2=nbOfCast-pos;
        if(pos2<nbOfCast)
          {
            ret1Ptr[i]=static_cast<T>(pos2);
            ret2Ptr[i]=work[i]-arrBg[pos2];
            castsDetected.insert(ret1Ptr[i]);
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::splitByValueRange : At rank #" << i << " the value is " << work[i] << " should be in [0," << *bg << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    ret3->alloc(castsDetected.size(),1);
    std::copy(castsDetected.begin(),castsDetected.end(),ret3->getPointer());
    castArr=ret1.retn();
    rankInsideCast=ret2.retn();
    castsPresent=ret3.retn();
  }

  /*!
   * This method look at \a this if it can be considered as a range defined by the 3-tuple ( \a strt , \a sttoopp , \a stteepp ).
   * If false is returned the tuple must be ignored. If true is returned \a this can be considered by a range( \a strt , \a sttoopp , \a stteepp ).
   * This method works only if \a this is allocated and single component. If not an exception will be thrown.
   *
   * \param [out] strt - the start of the range (included) if true is returned.
   * \param [out] sttoopp - the end of the range (not included) if true is returned.
   * \param [out] stteepp - the step of the range if true is returned.
   * \return the verdict of the check.
   *
   * \sa DataArray::GetNumberOfItemGivenBES
   */
  template<class T>
  bool DataArrayDiscrete<T>::isRange(T& strt, T& sttoopp, T& stteepp) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::isRange : this must be single component array !");
    mcIdType nbTuples(this->getNumberOfTuples());
    if(nbTuples==0)
      { strt=0; sttoopp=0; stteepp=1; return true; }
    const T *pt(this->begin());
    strt=*pt;
    if(nbTuples==1)
      { sttoopp=strt+1; stteepp=1; return true; }
    strt=*pt; sttoopp=pt[nbTuples-1];
    if(strt==sttoopp)
      return false;
    if(sttoopp>strt)
      {
        sttoopp++;
        T a(sttoopp-1-strt),tmp(strt);
        if(a%(nbTuples-1)!=0)
          return false;
        stteepp=a/(FromIdType<T>(nbTuples)-1);
        for(mcIdType i=0;i<nbTuples;i++,tmp+=stteepp)
          if(pt[i]!=tmp)
            return false;
        return true;
      }
    else
      {
        sttoopp--;
        T a(strt-sttoopp-1),tmp(strt);
        if(a%(nbTuples-1)!=0)
          return false;
        stteepp=-(a/(FromIdType<T>(nbTuples)-1));
        for(mcIdType i=0;i<nbTuples;i++,tmp+=stteepp)
          if(pt[i]!=tmp)
            return false;
        return true;
      }
  }

  /*!
   * Creates a one-dimensional DataArrayInt of given length, whose contents are computed
   * from values of \a this array, which is supposed to contain a renumbering map in
   * "Old to New" mode. The result array contains a renumbering map in "New to Old" mode.
   * To know how to use the renumbering maps see \ref numbering.
   *  \param [in] newNbOfElem - the number of tuples in the result array.
   *  \return DataArrayInt * - the new instance of DataArrayInt.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *
   *  \if ENABLE_EXAMPLES
   *  \ref cpp_mcdataarrayint_invertarrayo2n2n2o "Here is a C++ example".<br>
   *  \ref py_mcdataarrayint_invertarrayo2n2n2o  "Here is a Python example".
   *  \endif
   */
  template<class T>
  DataArrayIdType * DataArrayDiscrete<T>::invertArrayO2N2N2O(mcIdType newNbOfElem) const
  {
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New());
    ret->alloc(newNbOfElem,1);
    mcIdType nbOfOldNodes(this->getNumberOfTuples());
    const T *old2New(this->begin());
    mcIdType *pt(ret->getPointer());
    for(mcIdType i=0;i!=nbOfOldNodes;i++)
      {
        T newp(old2New[i]);
        if(newp!=-1)
          {
            if(newp>=0 && newp<newNbOfElem)
              pt[newp]=i;
            else
              {
                std::ostringstream oss; oss << "DataArrayInt::invertArrayO2N2N2O : At place #" << i << " the newplace is " << newp << " must be in [0," << newNbOfElem << ") !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
          }
      }
    return ret.retn();
  }

  /*!
   * Creates a one-dimensional DataArrayInt of given length, whose contents are computed
   * from values of \a this array, which is supposed to contain a renumbering map in
   * "New to Old" mode. The result array contains a renumbering map in "Old to New" mode.
   * To know how to use the renumbering maps see \ref numbering.
   *  \param [in] oldNbOfElem - the number of tuples in the result array.
   *  \return DataArrayInt * - the new instance of DataArrayInt.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *
   *  \if ENABLE_EXAMPLES
   *  \ref cpp_mcdataarrayint_invertarrayn2o2o2n "Here is a C++ example".
   *
   *  \ref py_mcdataarrayint_invertarrayn2o2o2n "Here is a Python example".
   *  \sa invertArrayN2O2O2NOptimized
   *  \endif
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::invertArrayN2O2O2N(mcIdType oldNbOfElem) const
  {
    this->checkAllocated();
    MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
    ret->alloc(oldNbOfElem,1);
    const T *new2Old=this->getConstPointer();
    mcIdType *pt=ret->getPointer();
    std::fill(pt,pt+oldNbOfElem,-1);
    mcIdType nbOfNewElems(this->getNumberOfTuples());
    for(mcIdType i=0;i<nbOfNewElems;i++)
      {
        T v(new2Old[i]);
        if(v>=0 && v<oldNbOfElem)
          pt[v]=i;
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::invertArrayN2O2O2N : in new id #" << i << " old value is " << v << " expected to be in [0," << oldNbOfElem << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    return ret.retn();
  }

  /*!
   * This method is similar to DataArrayInt::invertArrayO2N2N2O except that
   * Example : If \a this contains [0,1,2,0,3,4,5,4,6,4] this method will return [0,1,2,4,5,6,8] whereas DataArrayInt::invertArrayO2N2N2O returns [3,1,2,4,9,6,8]
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::invertArrayO2N2N2OBis(mcIdType newNbOfElem) const
  {
    MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
    ret->alloc(newNbOfElem,1);
    mcIdType nbOfOldNodes(this->getNumberOfTuples());
    const T *old2New=this->getConstPointer();
    mcIdType *pt=ret->getPointer();
    for(mcIdType i=nbOfOldNodes-1;i>=0;i--)
      {
        T newp(old2New[i]);
        if(newp!=-1)
          {
            if(newp>=0 && newp<newNbOfElem)
              pt[newp]=i;
            else
              {
                std::ostringstream oss; oss << "DataArrayInt::invertArrayO2N2N2OBis : At place #" << i << " the newplace is " << newp << " must be in [0," << newNbOfElem << ") !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
          }
      }
    return ret.retn();
  }

  /*!
   * Creates a map, whose contents are computed
   * from values of \a this array, which is supposed to contain a renumbering map in
   * "New to Old" mode. The result array contains a renumbering map in "Old to New" mode.
   * To know how to use the renumbering maps see \ref numbering.
   *  \return MapII  - the new instance of Map.
   *
   *  \if ENABLE_EXAMPLES
   *  \ref cpp_mcdataarrayint_invertarrayn2o2o2n "Here is a C++ example".
   *
   *  \ref py_mcdataarrayint_invertarrayn2o2o2n "Here is a Python example".
   *  \sa invertArrayN2O2O2N, giveN2OOptimized, MEDCouplingPointSet::renumberNodesInConn
   *  \endif
   */
  template <class T>
  MCAuto< MapKeyVal<T, mcIdType> > DataArrayDiscrete<T>::invertArrayN2O2O2NOptimized() const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::invertArrayN2O2O2NOptimized : single component expected !");
    MCAuto< MapKeyVal<T, mcIdType> > ret(MapKeyVal<T, mcIdType>::New());
    std::map<T, mcIdType>& m(ret->data());
    const T *new2Old(this->begin());
    mcIdType nbOfNewElems(this->getNumberOfTuples());
    for(mcIdType i=0;i<nbOfNewElems;i++)
      {
        T v(new2Old[i]);
        m[v]=i;
      }
    return ret;
  }

  /*!
   * Creates a map, whose contents are computed
   * from values of \a this array, which is supposed to contain a renumbering map in
   * "New to Old" mode. The result array contains a renumbering map in "New to Old" mode as C++ map for performance reasons.
   *
   * \sa invertArrayN2O2O2NOptimized, MEDCouplingPointSet::renumberNodesInConn
   */
  template <class T>
  MCAuto< MapKeyVal<mcIdType, T> > DataArrayDiscrete<T>::giveN2OOptimized() const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::giveN2OOptimized : single component expected !");
    MCAuto< MapKeyVal<mcIdType, T> > ret(MapKeyVal<mcIdType, T>::New());
    std::map<mcIdType,T>& m(ret->data());
    const T *new2Old(this->begin());
    mcIdType nbOfNewElems(this->getNumberOfTuples());
    for(mcIdType i=0;i<nbOfNewElems;i++)
      {
        T v(new2Old[i]);
        m[i]=v;
      }
    return ret;
  }

  /*!
   * This method finds for each element \a ELT in [valsBg,valsEnd) elements in \a this equal to it. Associated to ELT
   * this method will return the tuple id of last element found. If there is no element in \a this equal to ELT
   * an exception will be thrown.
   *
   * In case of success this[ret]==vals. Samely ret->transformWithIndArr(this->begin(),this->end())==vals.
   * Where \a vals is the [valsBg,valsEnd) array and \a ret the array returned by this method.
   * This method can be seen as an extension of FindPermutationFromFirstToSecond.
   * <br>
   * \b Example: <br>
   * - \a this: [17,27,2,10,-4,3,12,27,16]
   * - \a val : [3,16,-4,27,17]
   * - result: [5,8,4,7,0]
   *
   * \return - An array of size std::distance(valsBg,valsEnd)
   *
   * \sa DataArrayInt::FindPermutationFromFirstToSecond , DataArrayInt::FindPermutationFromFirstToSecondDuplicate
   */
  template <class T>
  MCAuto<DataArrayIdType> DataArrayDiscrete<T>::findIdForEach(const T *valsBg, const T *valsEnd) const
  {
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New());
    std::size_t nbOfTuplesOut(std::distance(valsBg,valsEnd));
    ret->alloc(nbOfTuplesOut,1);
    MCAuto< MapKeyVal<T, mcIdType> > zeMap(this->invertArrayN2O2O2NOptimized());
    const std::map<T, mcIdType>& dat(zeMap->data());
    mcIdType *ptToFeed(ret->getPointer());
    for(const T *pt=valsBg;pt!=valsEnd;pt++)
      {
        typename std::map<T,mcIdType>::const_iterator it(dat.find(*pt));
        if(it!=dat.end())
          *ptToFeed++=(*it).second;
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::findIdForEach : error for element at place " << std::distance(valsBg,pt);
            oss << " of input array value is " << *pt << " which is not in this !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
    return ret;
  }

  /*!
   * Returns a new DataArrayInt containing a renumbering map in "Old to New" mode.
   * This map, if applied to \a this array, would make it sorted. For example, if
   * \a this array contents are [9,10,0,6,4,11,3,7] then the contents of the result array
   * are [5,6,0,3,2,7,1,4]; if this result array (\a res) is used as an argument in call
   * \a this->renumber(\a res) then the returned array contains [0,3,4,6,7,9,10,11].
   * This method is useful for renumbering (in MED file for example). For more info
   * on renumbering see \ref numbering.
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *          array using decrRef() as it is no more needed.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If there are equal values in \a this array.
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::checkAndPreparePermutation() const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::checkAndPreparePermutation : number of components must == 1 !");
    mcIdType nbTuples(this->getNumberOfTuples());
    const T *pt=this->getConstPointer();
    mcIdType *pt2=this->CheckAndPreparePermutation(pt,pt+nbTuples);
    DataArrayIdType *ret=DataArrayIdType::New();
    ret->useArray(pt2,true,DeallocType::C_DEALLOC,nbTuples,1);
    return ret;
  }

  /*!
   * Returns two arrays describing a surjective mapping from \a this set of values (\a A)
   * onto a set of values of size \a targetNb (\a B). The surjective function is
   * \a B[ \a A[ i ]] = i. That is to say that for each \a id in [0,\a targetNb), where \a
   * targetNb < \a this->getNumberOfTuples(), there exists at least one tupleId (\a tid) so
   * that <em> this->getIJ( tid, 0 ) == id</em>. <br>
   * The first of out arrays returns indices of elements of \a this array, grouped by their
   * place in the set \a B. The second out array is the index of the first one; it shows how
   * many elements of \a A are mapped into each element of \a B. <br>
   * For more info on
   * mapping and its usage in renumbering see \ref numbering. <br>
   * \b Example:
   * - \a this: [0,3,2,3,2,2,1,2]
   * - \a targetNb: 4
   * - \a arr:  [0,  6,  2,4,5,7,  1,3]
   * - \a arrI: [0,1,2,6,8]
   *
   * This result means: <br>
   * the element of \a B 0 encounters within \a A once (\a arrI[ 0+1 ] - \a arrI[ 0 ]) and
   * its index within \a A is 0 ( \a arr[ 0:1 ] == \a arr[ \a arrI[ 0 ] : \a arrI[ 0+1 ]]);<br>
   * the element of \a B 2 encounters within \a A 4 times (\a arrI[ 2+1 ] - \a arrI[ 2 ]) and
   * its indices within \a A are [2,4,5,7] ( \a arr[ 2:6 ] == \a arr[ \a arrI[ 2 ] :
   * \a arrI[ 2+1 ]]); <br> etc.
   *  \param [in] targetNb - the size of the set \a B. \a targetNb must be equal or more
   *         than the maximal value of \a A.
   *  \param [out] arr - a new instance of DataArrayInt returning indices of
   *         elements of \a this, grouped by their place in the set \a B. The caller is to delete
   *         this array using decrRef() as it is no more needed.
   *  \param [out] arrI - a new instance of DataArrayInt returning size of groups of equal
   *         elements of \a this. The caller is to delete this array using decrRef() as it
   *         is no more needed.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If any value in \a this is more or equal to \a targetNb.
   */
  template <class T>
  void DataArrayDiscrete<T>::changeSurjectiveFormat(T targetNb, DataArrayIdType *&arr, DataArrayIdType *&arrI) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::changeSurjectiveFormat : number of components must == 1 !");
    mcIdType nbOfTuples(this->getNumberOfTuples());
    const T *input=this->getConstPointer();
    std::vector< std::vector<mcIdType> > tmp(targetNb);
    for(mcIdType i=0;i<nbOfTuples;i++)
      {
        T tmp2=input[i];
        if(tmp2>=0 && tmp2<targetNb)
          tmp[tmp2].push_back(i);
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::changeSurjectiveFormat : At pos " << i << " presence of element " << tmp2 << " ! should be in [0," << targetNb << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }

    MCAuto<DataArrayIdType> retI(DataArrayIdType::New());
    retI->alloc(targetNb+1,1);
    mcIdType *retIPtr=retI->getPointer();
    *retIPtr=0;
    for(std::vector< std::vector<mcIdType> >::const_iterator it1=tmp.begin();it1!=tmp.end();it1++,retIPtr++)
      retIPtr[1]=retIPtr[0]+ToIdType((*it1).size());
    if(nbOfTuples!=retI->getIJ(ToIdType(targetNb),0))
      throw INTERP_KERNEL::Exception("DataArrayInt::changeSurjectiveFormat : big problem should never happen !");
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New());
    ret->alloc(nbOfTuples,1);
    mcIdType *retPtr=ret->getPointer();
    for(std::vector< std::vector<mcIdType> >::const_iterator it1=tmp.begin();it1!=tmp.end();it1++)
      retPtr=std::copy((*it1).begin(),(*it1).end(),retPtr);
    arr=ret.retn();
    arrI=retI.retn();
  }

  /*!
   * Returns a new DataArrayInt containing a renumbering map in "New to Old" mode,
   * which if applied to \a this array would make it sorted ascendingly.
   * For more info on renumbering see \ref numbering. <br>
   * \b Example: <br>
   * - \a this: [2,0,1,1,0,1,2,0,1,1,0,0]
   * - result: [10,0,5,6,1,7,11,2,8,9,3,4]
   * - after applying result to \a this: [0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 2, 2]
   *
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *          array using decrRef() as it is no more needed.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::buildPermArrPerLevel() const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildPermArrPerLevel : number of components must == 1 !");
    mcIdType nbOfTuples=this->getNumberOfTuples();
    const T *pt=this->getConstPointer();
    std::map<T,mcIdType> m;
    MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
    ret->alloc(nbOfTuples,1);
    mcIdType *opt=ret->getPointer();
    for(mcIdType i=0;i<nbOfTuples;i++,pt++,opt++)
      {
        T val=*pt;
        typename std::map<T,mcIdType>::iterator it=m.find(val);
        if(it!=m.end())
          {
            *opt=(*it).second;
            (*it).second++;
          }
        else
          {
            *opt=0;
            m.insert(std::pair<T,mcIdType>(val,1));
          }
      }
    mcIdType sum=0;
    for(typename std::map<T,mcIdType>::iterator it=m.begin();it!=m.end();it++)
      {
        mcIdType vt=(*it).second;
        (*it).second=sum;
        sum+=vt;
      }
    pt=this->getConstPointer();
    opt=ret->getPointer();
    for(mcIdType i=0;i<nbOfTuples;i++,pt++,opt++)
      *opt+=m[*pt];
    //
    return ret.retn();
  }

  /*!
   * Checks if \a this array has the given size, and if its contents is equal to an array filled with
   * iota(). This method is particularly useful for DataArrayInt instances that represent
   * a renumbering array, to check if there is a real need in renumbering.
   * This method checks than \a this can be considered as an identity mapping
   * of a set having \a sizeExpected elements into itself.
   *
   *  \param [in] sizeExpected - The number of elements expected.
   *  \return bool - \a true if \a this array contents == \a range( \a this->getNumberOfTuples())
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   */
  template <class T>
  bool DataArrayDiscrete<T>::isIota(mcIdType sizeExpected) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      return false;
    mcIdType nbOfTuples(this->getNumberOfTuples());
    if(nbOfTuples!=sizeExpected)
      return false;
    const T *pt=this->getConstPointer();
    for(mcIdType i=0;i<nbOfTuples;i++,pt++)
      if(*pt!=i)
        return false;
    return true;
  }

  /*!
   * Checks if all values in \a this array are equal to \a val.
   *  \param [in] val - value to check equality of array values to.
   *  \return bool - \a true if all values are \a val.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1
   *  \sa DataArrayInt::checkUniformAndGuess
   */
  template <class T>
  bool DataArrayDiscrete<T>::isUniform(T val) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::isUniform : must be applied on DataArrayInt with only one component, you can call 'rearrange' method before !");
    const T *w(this->begin()),*end2(this->end());
    for(;w!=end2;w++)
      if(*w!=val)
        return false;
    return true;
  }

  /*!
   * This method checks that \a this is uniform. If not and exception will be thrown.
   * In case of uniformity the corresponding value is returned.
   *
   * \return mcIdType - the unique value contained in this
   * \throw If \a this is not allocated.
   * \throw If \a this->getNumberOfComponents() != 1
   * \throw If \a this is not uniform.
   * \sa DataArrayInt::isUniform
   */
  template <class T>
  T DataArrayDiscrete<T>::checkUniformAndGuess() const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::checkUniformAndGuess : must be applied on DataArrayInt with only one component, you can call 'rearrange' method before !");
    if(this->empty())
      throw INTERP_KERNEL::Exception("DataArrayInt::checkUniformAndGuess : this is empty !");
    const T *w(this->begin()),*end2(this->end());
    T ret(*w);
    for(;w!=end2;w++)
      if(*w!=ret)
        throw INTERP_KERNEL::Exception("DataArrayInt::checkUniformAndGuess : this is not uniform !");
    return ret;
  }

  /*!
   * Checks if all values in \a this array are unique.
   *  \return bool - \a true if condition above is true
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1
   */
  template <class T>
  bool DataArrayDiscrete<T>::hasUniqueValues() const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::hasOnlyUniqueValues: must be applied on DataArrayInt with only one component, you can call 'rearrange' method before !");
    std::size_t nbOfElements(this->getNumberOfTuples());
    std::set<T> s(this->begin(),this->end());  // in C++11, should use unordered_set (O(1) complexity)
    if (s.size() != nbOfElements)
      return false;
    return true;
  }

  /*!
   * Copy all components in a specified order from another DataArrayInt.
   * The specified components become the first ones in \a this array.
   * Both numerical and textual data is copied. The number of tuples in \a this and
   * the other array can be different.
   *  \param [in] a - the array to copy data from.
   *  \param [in] compoIds - sequence of zero based indices of components, data of which is
   *              to be copied.
   *  \throw If \a a is NULL.
   *  \throw If \a compoIds.size() != \a a->getNumberOfComponents().
   *  \throw If \a compoIds[i] < 0 or \a compoIds[i] > \a this->getNumberOfComponents().
   *
   *  \if ENABLE_EXAMPLES
   *  \ref py_mcdataarrayint_setselectedcomponents "Here is a Python example".
   *  \endif
   */
  template <class T>
  void DataArrayDiscrete<T>::setSelectedComponents(const DataArrayType *a, const std::vector<std::size_t>& compoIds)
  {
    if(!a)
      throw INTERP_KERNEL::Exception("DataArrayInt::setSelectedComponents : input DataArrayInt is NULL !");
    this->checkAllocated();
    a->checkAllocated();
    this->copyPartOfStringInfoFrom2(compoIds,*a);
    std::size_t partOfCompoSz=compoIds.size();
    std::size_t nbOfCompo = this->getNumberOfComponents();
    mcIdType nbOfTuples=std::min(this->getNumberOfTuples(),a->getNumberOfTuples());
    const T *ac=a->getConstPointer();
    T *nc=this->getPointer();
    for(mcIdType i=0;i<nbOfTuples;i++)
      for(std::size_t j=0;j<partOfCompoSz;j++,ac++)
        nc[nbOfCompo*i+compoIds[j]]=*ac;
  }

  /*!
   * This method searches each value in \a valToSearchIntoTuples among values in \a this given the corresponding tuple to find into.
   * If the value at the corresponding tuple is not found in the tuple an exception will be thrown.
   * If the value is found the corresponding component id is returned.
   *
   *  \param [in] valToSearchIntoTuples - a one component array containing the values to find in \a this
   *  \param [in] tupleIdHint - a one component array having same size than \a valToSearchIntoTuples giving for each value the tuple to find into
   *  \return DataArrayInt * - A newly allocated array having same size than \a valToSearchIntoTuples with one component
   *
   *   \b Example: <br>
   *          - \a this: [(0, 1, 2), (3, 4, 5), (6, 2, 3), (7, 8, 9), (9, 0, 10), (11, 12, 13), (14, 5, 11), (15, 16, 17)]
   *          - \a valToSearchIntoTuples: [1, 4, 6, 8, 10, 12, 14, 16, 17]
   *          - \a tupleIdHint: [0, 1, 2, 3, 4, 5, 6, 7, 7]
   *          - result array: [1, 1, 0, 1, 2, 1, 0, 1, 2] == <br>
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::locateComponentId(const DataArrayType *valToSearchIntoTuples, const DataArrayIdType *tupleIdHint) const
  {
    if(!valToSearchIntoTuples || !tupleIdHint)
      THROW_IK_EXCEPTION("DataArrayInt::locateComponentId : valToSearchIntoTuples and tupleIdHint must be not nullptr !");
    valToSearchIntoTuples->checkAllocated(); tupleIdHint->checkAllocated();
    this->checkAllocated();
    constexpr char MSG1[] = "DataArrayInt::locateComponentId : single component array expected";
    valToSearchIntoTuples->checkNbOfComps(1,MSG1); tupleIdHint->checkNbOfComps(1,MSG1);
    std::size_t nbOfCompo( this->getNumberOfComponents() );
    mcIdType thisNbTuples( this->getNumberOfTuples() );
    mcIdType nbOfTuples( valToSearchIntoTuples->getNumberOfTuples() );
    tupleIdHint->checkNbOfTuples(nbOfTuples,"Number of tuples of input arrays must be the same.");
    const T *cPtr(this->begin()),*valSearchPt(valToSearchIntoTuples->begin());
    const mcIdType *tHintPtr(tupleIdHint->begin());
    MCAuto<DataArrayIdType> ret = DataArrayIdType::New();
    ret->alloc(nbOfTuples,1);
    mcIdType *retPtr(ret->getPointer());
    for(auto i = 0 ; i < nbOfTuples ; ++i)
    {
      if( tHintPtr[i] >=0 && tHintPtr[i] < thisNbTuples )
      {
        auto strtSearch(cPtr+nbOfCompo*tHintPtr[i]),endSearch(cPtr+nbOfCompo*(tHintPtr[i]+1));
        auto pos = std::find(strtSearch,endSearch,valSearchPt[i]);
        if(pos != endSearch)
          *retPtr++ = ToIdType( std::distance(strtSearch,pos) );
        else
          THROW_IK_EXCEPTION("At pos " << i << " value " << valSearchPt[i] << " is not present at tuple " << tHintPtr[i]);
      }
      else
        THROW_IK_EXCEPTION("At pos " << i << " hint tuple is " << tHintPtr[i] << " not in [0," << thisNbTuples << ")");
    }
    return ret.retn();
  }

  /*!
   * Creates a new DataArrayInt containing IDs (indices) of tuples holding value \b not
   * equal to a given one.
   *  \param [in] val - the value to ignore within \a this.
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *          array using decrRef() as it is no more needed.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::findIdsNotEqual(T val) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdsNotEqual : the array must have only one component, you can call 'rearrange' method before !");
    const T *cptr(this->getConstPointer());
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New());
    ret->alloc(0,1);
    mcIdType nbOfTuples(this->getNumberOfTuples());
    for(mcIdType i=0;i<nbOfTuples;i++,cptr++)
      if(*cptr!=val)
        ret->pushBackSilent(i);
    return ret.retn();
  }

  /*!
   * Creates a new DataArrayInt containing IDs (indices) of tuples holding tuple equal to those defined by [ \a tupleBg , \a tupleEnd )
   * This method is an extension of  DataArrayInt::findIdsEqual method.
   *
   *  \param [in] tupleBg - the begin (included) of the input tuple to find within \a this.
   *  \param [in] tupleEnd - the end (excluded) of the input tuple to find within \a this.
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *          array using decrRef() as it is no more needed.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != std::distance(tupleBg,tupleEnd).
   * \throw If \a this->getNumberOfComponents() is equal to 0.
   * \sa DataArrayInt::findIdsEqual
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::findIdsEqualTuple(const T *tupleBg, const T *tupleEnd) const
  {
    std::size_t nbOfCompoExp=std::distance(tupleBg,tupleEnd);
    this->checkAllocated();
    if(this->getNumberOfComponents()!=nbOfCompoExp)
      {
        std::ostringstream oss; oss << "DataArrayInt::findIdsEqualTuple : mismatch of number of components. Input tuple has " << nbOfCompoExp << " whereas this array has " << this->getNumberOfComponents() << " components !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    if(nbOfCompoExp==0)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdsEqualTuple : number of components should be > 0 !");
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New());
    ret->alloc(0,1);
    const T *bg(this->begin()),*end2(this->end()),*work(this->begin());
    while(work!=end2)
      {
        work=std::search(work,end2,tupleBg,tupleEnd);
        if(work!=end2)
          {
            std::ptrdiff_t pos=std::distance(bg,work);
            if(pos%nbOfCompoExp==0)
              ret->pushBackSilent(ToIdType(pos/nbOfCompoExp));
            work++;
          }
      }
    return ret.retn();
  }

  /*!
   * Creates a new DataArrayInt containing IDs (indices) of tuples holding value equal to
   * one of given values.
   *  \param [in] valsBg - an array of values to find within \a this array.
   *  \param [in] valsEnd - specifies the end of the array \a valsBg, so that
   *              the last value of \a valsBg is \a valsEnd[ -1 ].
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *          array using decrRef() as it is no more needed.
   *  \throw If \a this->getNumberOfComponents() != 1.
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::findIdsEqualList(const T *valsBg, const T *valsEnd) const
  {
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdsEqualList : the array must have only one component, you can call 'rearrange' method before !");
    std::set<T> vals2(valsBg,valsEnd);
    const T *cptr(this->getConstPointer());
    mcIdType nbOfTuples(this->getNumberOfTuples());
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(0,1);
    for(mcIdType i=0;i<nbOfTuples;i++,cptr++)
      if(vals2.find(*cptr)!=vals2.end())
        ret->pushBackSilent(i);
    return ret.retn();
  }

  /*!
   * Creates a new DataArrayInt containing IDs (indices) of tuples holding values \b not
   * equal to any of given values.
   *  \param [in] valsBg - an array of values to ignore within \a this array.
   *  \param [in] valsEnd - specifies the end of the array \a valsBg, so that
   *              the last value of \a valsBg is \a valsEnd[ -1 ].
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *          array using decrRef() as it is no more needed.
   *  \throw If \a this->getNumberOfComponents() != 1.
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::findIdsNotEqualList(const T *valsBg, const T *valsEnd) const
  {
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdsNotEqualList : the array must have only one component, you can call 'rearrange' method before !");
    std::set<T> vals2(valsBg,valsEnd);
    const T *cptr=this->getConstPointer();
    mcIdType nbOfTuples(this->getNumberOfTuples());
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(0,1);
    for(mcIdType i=0;i<nbOfTuples;i++,cptr++)
      if(vals2.find(*cptr)==vals2.end())
        ret->pushBackSilent(i);
    return ret.retn();
  }

  /*!
   * This method expects to be called when number of components of this is equal to one.
   * This method returns the tuple id, if it exists, of the first tuple equal to \b value.
   * If not any tuple contains \b value -1 is returned.
   * \sa DataArrayInt::presenceOfValue
   */
  template <class T>
  mcIdType DataArrayDiscrete<T>::findIdFirstEqual(T value) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::presenceOfValue : the array must have only one component, you can call 'rearrange' method before !");
    const T *cptr=this->getConstPointer();
    mcIdType nbOfTuples(this->getNumberOfTuples());
    const T *ret=std::find(cptr,cptr+nbOfTuples,value);
    if(ret!=cptr+nbOfTuples)
      return ToIdType(std::distance(cptr,ret));
    return -1;
  }

  /*!
   * This method expects to be called when number of components of this is equal to one.
   * This method returns the tuple id, if it exists, of the first tuple so that the value is contained in \b vals.
   * If not any tuple contains one of the values contained in 'vals' -1 is returned.
   * \sa DataArrayInt::presenceOfValue
   */
  template <class T>
  mcIdType DataArrayDiscrete<T>::findIdFirstEqual(const std::vector<T>& vals) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::presenceOfValue : the array must have only one component, you can call 'rearrange' method before !");
    std::set<T> vals2(vals.begin(),vals.end());
    const T *cptr=this->getConstPointer();
    mcIdType nbOfTuples(this->getNumberOfTuples());
    for(const T *w=cptr;w!=cptr+nbOfTuples;w++)
      if(vals2.find(*w)!=vals2.end())
        return ToIdType(std::distance(cptr,w));
    return -1;
  }

  /*!
   * This method is an extension of DataArrayInt::findIdFirstEqual method because this method works for DataArrayInt with
   * any number of components excepted 0 (an INTERP_KERNEL::Exception is thrown in this case).
   * This method searches in \b this is there is a tuple that matched the input parameter \b tupl.
   * If any the tuple id is returned. If not -1 is returned.
   *
   * This method throws an INTERP_KERNEL::Exception if the number of components in \b this mismatches with the size of
   * the input vector. An INTERP_KERNEL::Exception is thrown too if \b this is not allocated.
   *
   * \return tuple id where \b tupl is. -1 if no such tuple exists in \b this.
   * \sa DataArrayInt::findIdSequence, DataArrayInt::presenceOfTuple.
   */
  template <class T>
  mcIdType DataArrayDiscrete<T>::findIdFirstEqualTuple(const std::vector<T>& tupl) const
  {
    this->checkAllocated();
    std::size_t nbOfCompo(this->getNumberOfComponents());
    if(nbOfCompo==0)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdFirstEqualTuple : 0 components in 'this' !");
    if(nbOfCompo!=tupl.size())
      {
        std::ostringstream oss; oss << "DataArrayInt::findIdFirstEqualTuple : 'this' contains " << nbOfCompo << " components and searching for a tuple of length " << tupl.size() << " !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    const T *cptr=this->getConstPointer();
    std::size_t nbOfVals=this->getNbOfElems();
    for(const T *work=cptr;work!=cptr+nbOfVals;)
      {
        work=std::search(work,cptr+nbOfVals,tupl.begin(),tupl.end());
        if(work!=cptr+nbOfVals)
          {
            if(std::distance(cptr,work)%nbOfCompo!=0)
              work++;
            else
              return ToIdType (std::distance(cptr,work)/nbOfCompo);
          }
      }
    return -1;
  }

  /*!
   * This method searches the sequence specified in input parameter \b vals in \b this.
   * This works only for DataArrayInt having number of components equal to one (if not an INTERP_KERNEL::Exception will be thrown).
   * This method differs from DataArrayInt::findIdFirstEqualTuple in that the position is internal raw data is not considered here contrary to DataArrayInt::findIdFirstEqualTuple.
   * \sa DataArrayInt::findIdFirstEqualTuple
   */
  template <class T>
  mcIdType DataArrayDiscrete<T>::findIdSequence(const std::vector<T>& vals) const
  {
    this->checkAllocated();
    std::size_t nbOfCompo=this->getNumberOfComponents();
    if(nbOfCompo!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdSequence : works only for DataArrayInt instance with one component !");
    const T *cptr=this->getConstPointer();
    std::size_t nbOfVals=this->getNbOfElems();
    const T *loc=std::search(cptr,cptr+nbOfVals,vals.begin(),vals.end());
    if(loc!=cptr+nbOfVals)
      return ToIdType(std::distance(cptr,loc));
    return -1;
  }

  /*!
   * Assigns \a newValue to all elements holding \a oldValue within \a this
   * one-dimensional array.
   *  \param [in] oldValue - the value to replace.
   *  \param [in] newValue - the value to assign.
   *  \return mcIdType - number of replacements performed.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   */
  template <class T>
  mcIdType DataArrayDiscrete<T>::changeValue(T oldValue, T newValue)
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::changeValue : the array must have only one component, you can call 'rearrange' method before !");
    if(oldValue==newValue)
      return 0;
    T *start(this->getPointer()),*end2(start+this->getNbOfElems());
    mcIdType ret(0);
    for(T *val=start;val!=end2;val++)
      {
        if(*val==oldValue)
          {
            *val=newValue;
            ret++;
          }
      }
    if(ret>0)
      this->declareAsNew();
    return ret;
  }

  /*!
   * This method returns the number of values in \a this that are equals to input parameter \a value.
   * This method only works for single component array.
   *
   * \return a value in [ 0, \c this->getNumberOfTuples() )
   *
   * \throw If \a this is not allocated
   *
   */
  template <class T>
  mcIdType DataArrayDiscrete<T>::count(T value) const
  {
    mcIdType ret=0;
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::count : must be applied on DataArrayInt with only one component, you can call 'rearrange' method before !");
    const T *vals=this->begin();
    std::size_t nbOfElements=this->getNumberOfTuples();
    for(std::size_t i=0;i<nbOfElements;i++,vals++)
      if(*vals==value)
        ret++;
    return ret;
  }

  /*!
   * This method is an extension of DataArrayInt::presenceOfValue method because this method works for DataArrayInt with
   * any number of components excepted 0 (an INTERP_KERNEL::Exception is thrown in this case).
   * This method searches in \b this is there is a tuple that matched the input parameter \b tupl.
   * This method throws an INTERP_KERNEL::Exception if the number of components in \b this mismatches with the size of
   * the input vector. An INTERP_KERNEL::Exception is thrown too if \b this is not allocated.
   * \sa DataArrayInt::findIdFirstEqualTuple
   */
  template <class T>
  bool DataArrayDiscrete<T>::presenceOfTuple(const std::vector<T>& tupl) const
  {
    return this->findIdFirstEqualTuple(tupl)!=-1;
  }


  /*!
   * Returns \a true if a given value is present within \a this one-dimensional array.
   *  \param [in] value - the value to find within \a this array.
   *  \return bool - \a true in case if \a value is present within \a this array.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \sa findIdFirstEqual()
   */
  template <class T>
  bool DataArrayDiscrete<T>::presenceOfValue(T value) const
  {
    return this->findIdFirstEqual(value)!=-1;
  }

  /*!
   * This method expects to be called when number of components of this is equal to one.
   * This method returns true if it exists a tuple so that the value is contained in \b vals.
   * If not any tuple contains one of the values contained in 'vals' false is returned.
   * \sa DataArrayInt::findIdFirstEqual
   */
  template <class T>
  bool DataArrayDiscrete<T>::presenceOfValue(const std::vector<T>& vals) const
  {
    return this->findIdFirstEqual(vals)!=-1;
  }

  /*!
   * Accumulates values of each component of \a this array.
   *  \param [out] res - an array of length \a this->getNumberOfComponents(), allocated
   *         by the caller, that is filled by this method with sum value for each
   *         component.
   *  \throw If \a this is not allocated.
   */
  template <class T>
  void DataArrayDiscrete<T>::accumulate(T *res) const
  {
    this->checkAllocated();
    const T *ptr=this->getConstPointer();
    mcIdType nbTuple(this->getNumberOfTuples());
    std::size_t nbComps(this->getNumberOfComponents());
    std::fill(res,res+nbComps,0);
    for(mcIdType i=0;i<nbTuple;i++)
      std::transform(ptr+i*nbComps,ptr+(i+1)*nbComps,res,res,std::plus<T>());
  }

  template <class T>
  T DataArrayDiscrete<T>::accumulate(std::size_t compId) const
  {
    this->checkAllocated();
    const T *ptr=this->getConstPointer();
    mcIdType nbTuple(this->getNumberOfTuples());
    std::size_t nbComps(this->getNumberOfComponents());
    if(compId>=nbComps) // compId >= 0 (it is a size_t)
      throw INTERP_KERNEL::Exception("DataArrayInt::accumulate : Invalid compId specified : No such nb of components !");
    T ret=0;
    for(mcIdType i=0;i<nbTuple;i++)
      ret+=ptr[i*nbComps+compId];
    return ret;
  }

  /*!
   * This method accumulate using addition tuples in \a this using input index array [ \a bgOfIndex, \a endOfIndex ).
   * The returned array will have same number of components than \a this and number of tuples equal to
   * \c std::distance(bgOfIndex,endOfIndex) \b minus \b one.
   *
   * The input index array is expected to be ascendingly sorted in which the all referenced ids should be in [0, \c this->getNumberOfTuples).
   *
   * \param [in] bgOfIndex - begin (included) of the input index array.
   * \param [in] endOfIndex - end (excluded) of the input index array.
   * \return DataArrayInt * - the new instance having the same number of components than \a this.
   *
   * \throw If bgOfIndex or end is NULL.
   * \throw If input index array is not ascendingly sorted.
   * \throw If there is an id in [ \a bgOfIndex, \a endOfIndex ) not in [0, \c this->getNumberOfTuples).
   * \throw If std::distance(bgOfIndex,endOfIndex)==0.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::accumulatePerChunck(const mcIdType *bgOfIndex, const mcIdType *endOfIndex) const
  {
    if(!bgOfIndex || !endOfIndex)
      throw INTERP_KERNEL::Exception("DataArrayInt::accumulatePerChunck : input pointer NULL !");
    this->checkAllocated();
    std::size_t nbCompo(this->getNumberOfComponents());
    mcIdType nbOfTuples(this->getNumberOfTuples());
    mcIdType sz=ToIdType(std::distance(bgOfIndex,endOfIndex));
    if(sz<1)
      throw INTERP_KERNEL::Exception("DataArrayInt::accumulatePerChunck : invalid size of input index array !");
    sz--;
    MCAuto<DataArrayType> ret=DataArrayType::New(); ret->alloc(sz,nbCompo);
    const mcIdType *w=bgOfIndex;
    if(*w<0 || *w>=nbOfTuples)
      throw INTERP_KERNEL::Exception("DataArrayInt::accumulatePerChunck : The first element of the input index not in [0,nbOfTuples) !");
    const T *srcPt=this->begin()+(*w)*nbCompo;
    T *tmp=ret->getPointer();
    for(mcIdType i=0;i<sz;i++,tmp+=nbCompo,w++)
      {
        std::fill(tmp,tmp+nbCompo,0);
        if(w[1]>=w[0])
          {
            for(mcIdType j=w[0];j<w[1];j++,srcPt+=nbCompo)
              {
                if(j>=0 && j<nbOfTuples)
                  std::transform(srcPt,srcPt+nbCompo,tmp,tmp,std::plus<T>());
                else
                  {
                    std::ostringstream oss; oss << "DataArrayInt::accumulatePerChunck : At rank #" << i << " the input index array points to id " << j << " should be in [0," << nbOfTuples << ") !";
                    throw INTERP_KERNEL::Exception(oss.str().c_str());
                  }
              }
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::accumulatePerChunck : At rank #" << i << " the input index array is not in ascendingly sorted.";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    ret->copyStringInfoFrom(*this);
    return ret.retn();
  }

  /*!
   * Returns in a single walk in \a this the min value and the max value in \a this.
   * \a this is expected to be single component array.
   *
   * \param [out] minValue - the min value in \a this.
   * \param [out] maxValue - the max value in \a this.
   *
   * \sa getMinValueInArray, getMinValue, getMaxValueInArray, getMaxValue
   */
  template <class T>
  void DataArrayDiscrete<T>::getMinMaxValues(T& minValue, T& maxValue) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::getMinMaxValues : must be applied on DataArrayInt with only one component !");
    std::size_t nbElements(this->getNumberOfTuples());
    const T *pt(this->begin());
    minValue=std::numeric_limits<T>::max(); maxValue=-std::numeric_limits<T>::max();
    for(std::size_t i=0;i<nbElements;i++,pt++)
      {
        if(*pt<minValue)
          minValue=*pt;
        if(*pt>maxValue)
          maxValue=*pt;
      }
  }

  /*!
   * Modify all elements of \a this array, so that
   * an element _x_ becomes \f$ numerator / x \f$.
   *  \warning If an exception is thrown because of presence of 0 element in \a this
   *           array, all elements processed before detection of the zero element remain
   *           modified.
   *  \param [in] numerator - the numerator used to modify array elements.
   *  \throw If \a this is not allocated.
   *  \throw If there is an element equal to 0 in \a this array.
   */
  template <class T>
  void DataArrayDiscrete<T>::applyInv(T numerator)
  {
    this->checkAllocated();
    T *ptr=this->getPointer();
    std::size_t nbOfElems=this->getNbOfElems();
    for(std::size_t i=0;i<nbOfElems;i++,ptr++)
      {
        if(*ptr!=0)
          {
            *ptr=numerator/(*ptr);
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::applyInv : presence of null value in tuple #" << i/(this->getNumberOfComponents()) << " component #" << i%(this->getNumberOfComponents());
            oss << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    this->declareAsNew();
  }

  /*!
   * Modify all elements of \a this array, so that
   * an element _x_ becomes \f$ x / val \f$.
   *  \param [in] val - the denominator used to modify array elements.
   *  \throw If \a this is not allocated.
   *  \throw If \a val == 0.
   */
  template <class T>
  void DataArrayDiscrete<T>::applyDivideBy(T val)
  {
    if(val==0)
      throw INTERP_KERNEL::Exception("DataArrayInt::applyDivideBy : Trying to divide by 0 !");
    this->checkAllocated();
    T *ptr=this->getPointer();
    std::size_t nbOfElems=this->getNbOfElems();
    std::transform(ptr,ptr+nbOfElems,ptr,std::bind(std::divides<T>(),std::placeholders::_1,val));
    this->declareAsNew();
  }

  /*!
   * Modify all elements of \a this array, so that
   * an element _x_ becomes  <em> x % val </em>.
   *  \param [in] val - the divisor used to modify array elements.
   *  \throw If \a this is not allocated.
   *  \throw If \a val <= 0.
   */
  template <class T>
  void DataArrayDiscrete<T>::applyModulus(T val)
  {
    if(val<=0)
      throw INTERP_KERNEL::Exception("DataArrayInt::applyDivideBy : Trying to operate modulus on value <= 0 !");
    this->checkAllocated();
    T *ptr=this->getPointer();
    std::size_t nbOfElems=this->getNbOfElems();
    std::transform(ptr,ptr+nbOfElems,ptr,std::bind(std::modulus<T>(),std::placeholders::_1,val));
    this->declareAsNew();
  }

  /*!
   * Modify all elements of \a this array, so that
   * an element _x_ becomes <em> val % x </em>.
   *  \warning If an exception is thrown because of presence of an element <= 0 in \a this
   *           array, all elements processed before detection of the zero element remain
   *           modified.
   *  \param [in] val - the divident used to modify array elements.
   *  \throw If \a this is not allocated.
   *  \throw If there is an element equal to or less than 0 in \a this array.
   */
  template <class T>
  void DataArrayDiscrete<T>::applyRModulus(T val)
  {
    this->checkAllocated();
    T *ptr=this->getPointer();
    std::size_t nbOfElems=this->getNbOfElems();
    for(std::size_t i=0;i<nbOfElems;i++,ptr++)
      {
        if(*ptr>0)
          {
            *ptr=val%(*ptr);
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::applyRModulus : presence of value <=0 in tuple #" << i/(this->getNumberOfComponents()) << " component #" << i%(this->getNumberOfComponents());
            oss << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    this->declareAsNew();
  }

  /*!
   * Modify all elements of \a this array, so that
   * an element _x_ becomes <em> val ^ x </em>.
   *  \param [in] val - the value used to apply pow on all array elements.
   *  \throw If \a this is not allocated.
   *  \throw If \a val < 0.
   */
  template <class T>
  void DataArrayDiscrete<T>::applyPow(T val)
  {
    this->checkAllocated();
    if(val<0)
      throw INTERP_KERNEL::Exception("DataArrayInt::applyPow : input pow in < 0 !");
    T *ptr=this->getPointer();
    std::size_t nbOfElems=this->getNbOfElems();
    if(val==0)
      {
        std::fill(ptr,ptr+nbOfElems,1);
        return ;
      }
    for(std::size_t i=0;i<nbOfElems;i++,ptr++)
      {
        T tmp=1;
        for(T j=0;j<val;j++)
          tmp*=*ptr;
        *ptr=tmp;
      }
    this->declareAsNew();
  }

  /*!
   * Modify all elements of \a this array, so that
   * an element _x_ becomes \f$ val ^ x \f$.
   *  \param [in] val - the value used to apply pow on all array elements.
   *  \throw If \a this is not allocated.
   *  \throw If there is an element < 0 in \a this array.
   *  \warning If an exception is thrown because of presence of 0 element in \a this
   *           array, all elements processed before detection of the zero element remain
   *           modified.
   */
  template <class T>
  void DataArrayDiscrete<T>::applyRPow(T val)
  {
    this->checkAllocated();
    T *ptr=this->getPointer();
    std::size_t nbOfElems=this->getNbOfElems();
    for(std::size_t i=0;i<nbOfElems;i++,ptr++)
      {
        if(*ptr>=0)
          {
            T tmp=1;
            for(T j=0;j<*ptr;j++)
              tmp*=val;
            *ptr=tmp;
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::applyRPow : presence of negative value in tuple #" << i/(this->getNumberOfComponents()) << " component #" << i%(this->getNumberOfComponents());
            oss << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    this->declareAsNew();
  }

  /*!
   * This method works only on data array with one component.
   * This method returns a newly allocated array storing stored ascendantly tuple ids in \b this so that
   * this[*id] in [\b vmin,\b vmax)
   *
   * \param [in] vmin begin of range. This value is included in range (included).
   * \param [in] vmax end of range. This value is \b not included in range (excluded).
   * \return a newly allocated data array that the caller should deal with.
   *
   * \sa DataArrayInt::findIdsNotInRange , DataArrayInt::findIdsStricltyNegative
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::findIdsInRange(T vmin, T vmax) const
  {
    InRange<T> ir(vmin,vmax);
    MCAuto<DataArrayIdType> ret(this->findIdsAdv(ir));
    return ret.retn();
  }

  /*!
   * This method works only on data array with one component.
   * This method returns a newly allocated array storing stored ascendantly tuple ids in \b this so that
   * this[*id] \b not in [\b vmin,\b vmax)
   *
   * \param [in] vmin begin of range. This value is \b not included in range (excluded).
   * \param [in] vmax end of range. This value is included in range (included).
   * \return a newly allocated data array that the caller should deal with.
   *
   * \sa DataArrayInt::findIdsInRange , DataArrayInt::findIdsStricltyNegative
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::findIdsNotInRange(T vmin, T vmax) const
  {
    NotInRange<T> nir(vmin,vmax);
    MCAuto<DataArrayIdType> ret(this->findIdsAdv(nir));
    return ret.retn();
  }

  /*!
   * This method works only on data array with one component.
   * This method checks that all ids in \b this are in [ \b vmin, \b vmax ). If there is at least one element in \a this not in [ \b vmin, \b vmax ) an exception will be thrown.
   *
   * \param [in] vmin begin of range. This value is included in range (included).
   * \param [in] vmax end of range. This value is \b not included in range (excluded).
   * \return if all ids in \a this are so that (*this)[i]==i for all i in [ 0, \c this->getNumberOfTuples() ). */
  template <class T>
  bool DataArrayDiscrete<T>::checkAllIdsInRange(T vmin, T vmax) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::checkAllIdsInRange : this must have exactly one component !");
    mcIdType nbOfTuples(this->getNumberOfTuples());
    bool ret=true;
    const T *cptr=this->getConstPointer();
    for(mcIdType i=0;i<nbOfTuples;i++,cptr++)
      {
        if(*cptr>=vmin && *cptr<vmax)
          { ret=ret && *cptr==i; }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::checkAllIdsInRange : tuple #" << i << " has value " << *cptr << " should be in [" << vmin << "," << vmax << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    return ret;
  }

  /*!
   * Returns a new DataArrayInt which contains a complement of elements of \a this
   * one-dimensional array. I.e. the result array contains all elements from the range [0,
   * \a nbOfElement) not present in \a this array.
   *  \param [in] nbOfElement - maximal size of the result array.
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *         array using decrRef() as it is no more needed.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If any element \a x of \a this array violates condition ( 0 <= \a x < \a
   *         nbOfElement ).
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::buildComplement(mcIdType nbOfElement) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildComplement : only single component allowed !");
    std::vector<bool> tmp(nbOfElement);
    const T *pt=this->getConstPointer();
    std::size_t nbOfElements=this->getNumberOfTuples();
    for(const T *w=pt;w!=pt+nbOfElements;w++)
      if(*w>=0 && *w<nbOfElement)
        tmp[*w]=true;
      else
        throw INTERP_KERNEL::Exception("DataArrayInt::buildComplement : an element is not in valid range : [0,nbOfElement) !");
    std::size_t nbOfRetVal=std::count(tmp.begin(),tmp.end(),false);
    DataArrayIdType *ret=DataArrayIdType::New();
    ret->alloc(nbOfRetVal,1);
    mcIdType j=0;
    mcIdType *retPtr=ret->getPointer();
    for(mcIdType i=0;i<nbOfElement;i++)
      if(!tmp[i])
        retPtr[j++]=i;
    return ret;
  }

  /*!
   * Returns a new DataArrayInt containing elements of \a this one-dimensional missing
   * from an \a other one-dimensional array.
   *  \param [in] other - a DataArrayInt containing elements not to include in the result array.
   *  \return DataArrayInt * - a new instance of DataArrayInt with one component. The
   *         caller is to delete this array using decrRef() as it is no more needed.
   *  \throw If \a other is NULL.
   *  \throw If \a other is not allocated.
   *  \throw If \a other->getNumberOfComponents() != 1.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \sa DataArrayInt::buildSubstractionOptimized()
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::buildSubstraction(const DataArrayType *other) const
  {
    if(!other)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildSubstraction : DataArrayInt pointer in input is NULL !");
    this->checkAllocated();
    other->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildSubstraction : only single component allowed !");
    if(other->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildSubstraction : only single component allowed for other type !");
    const T *pt=this->getConstPointer();
    std::size_t nbOfElements=this->getNumberOfTuples();
    std::set<T> s1(pt,pt+nbOfElements);
    pt=other->getConstPointer();
    nbOfElements=other->getNumberOfTuples();
    std::set<T> s2(pt,pt+nbOfElements);
    std::vector<T> r;
    std::set_difference(s1.begin(),s1.end(),s2.begin(),s2.end(),std::back_insert_iterator< std::vector<T> >(r));
    DataArrayType *ret=DataArrayType::New();
    ret->alloc(r.size(),1);
    std::copy(r.begin(),r.end(),ret->getPointer());
    return ret;
  }

  /*!
   * \a this is expected to have one component and to be sorted ascendingly (as for \a other).
   * \a other is expected to be a part of \a this. If not DataArrayInt::buildSubstraction should be called instead.
   *
   * \param [in] other an array with one component and expected to be sorted ascendingly.
   * \ret list of ids in \a this but not in \a other.
   * \sa DataArrayInt::buildSubstraction
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::buildSubstractionOptimized(const DataArrayType *other) const
  {
    static const char *MSG="DataArrayInt::buildSubstractionOptimized : only single component allowed !";
    if(!other) throw INTERP_KERNEL::Exception("DataArrayInt::buildSubstractionOptimized : NULL input array !");
    this->checkAllocated(); other->checkAllocated();
    if(this->getNumberOfComponents()!=1) throw INTERP_KERNEL::Exception(MSG);
    if(other->getNumberOfComponents()!=1) throw INTERP_KERNEL::Exception(MSG);
    const T *pt1Bg(this->begin()),*pt1End(this->end()),*pt2Bg(other->begin()),*pt2End(other->end());
    const T *work1(pt1Bg),*work2(pt2Bg);
    MCAuto<DataArrayType> ret(DataArrayType::New()); ret->alloc(0,1);
    for(;work1!=pt1End;work1++)
      {
        if(work2!=pt2End && *work1==*work2)
          work2++;
        else
          ret->pushBackSilent(*work1);
      }
    return ret.retn();
  }

  /*!
   * Returns a new DataArrayInt which contains all elements of \a this and a given
   * one-dimensional arrays. The result array does not contain any duplicates
   * and its values are sorted in ascending order.
   *  \param [in] other - an array to unite with \a this one.
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *         array using decrRef() as it is no more needed.
   *  \throw If \a this or \a other is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If \a other->getNumberOfComponents() != 1.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::buildUnion(const DataArrayType *other) const
  {
    std::vector<const DataArrayType *>arrs(2);
    arrs[0]=dynamic_cast<const DataArrayType *>(this); arrs[1]=other;
    return DataArrayDiscrete<T>::BuildUnion(arrs);
  }

  /*!
   * Returns a new DataArrayInt which contains elements present in both \a this and a given
   * one-dimensional arrays. The result array does not contain any duplicates
   * and its values are sorted in ascending order.
   *  \param [in] other - an array to intersect with \a this one.
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *         array using decrRef() as it is no more needed.
   *  \throw If \a this or \a other is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If \a other->getNumberOfComponents() != 1.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::buildIntersection(const DataArrayType *other) const
  {
    std::vector<const DataArrayType *>arrs(2);
    arrs[0]=dynamic_cast<const DataArrayType *>(this); arrs[1]=other;
    return DataArrayDiscrete<T>::BuildIntersection(arrs);
  }
  /*!
   * This method can be applied on allocated with one component DataArrayInt instance.
   * Locate groups of all consecutive same values in \a this and return them into an indexed array of positions pointing to \a this starting with 0.
   * Number of tuples of returned array is equal to size of \a this->buildUnique() + 1.
   * Last value of returned array is equal to \a this->getNumberOfTuples()
   *
   * \b Example:
   * - \a this : [0, 1, 1, 2, 2, 3, 4, 4, 5, 5, 5, 11]
   * - \a return is : [0, 1, 3, 5, 6, 8, 11, 12]
   *
   * \return a newly allocated array containing the indexed array format of groups by same consecutive value.
   * \throw if \a this is not allocated or if \a this has not exactly one component.
   * \sa DataArrayInt::buildUnique, MEDCouplingSkyLineArray::groupPacks
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::indexOfSameConsecutiveValueGroups() const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::indexOfSameConsecutiveValueGroups : only single component allowed !");
    const T *pt(this->begin());
    const T *const ptEnd(this->end()) , * const ptBg(this->begin());
    // first find nb of different values in this
    std::size_t nbOfTuplesOut(0);
    while( pt != ptEnd )
    {
      T val(*pt);
      const T *endOfPack(std::find_if(pt+1,ptEnd,[val](T elt){ return val!=elt; }));
      pt = endOfPack;
      ++nbOfTuplesOut;
    }
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(nbOfTuplesOut+1,1);
    mcIdType *retPtr(ret->getPointer()); *retPtr++ = 0;
    pt = this->begin();
    while( pt != ptEnd )
    {
      T val(*pt);
      const T *endOfPack(std::find_if(pt+1,ptEnd,[val](T elt){ return val!=elt; }));
      *retPtr++ = ToIdType( std::distance(ptBg,endOfPack) );
      pt = endOfPack;
      ++nbOfTuplesOut;
    }
    return ret.retn();
  }

  /*!
   * This method can be applied on allocated with one component DataArrayInt instance.
   * This method is typically relevant for sorted arrays. All consecutive duplicated items in \a this will appear only once in returned DataArrayInt instance.
   * Example : if \a this contains [1,2,2,3,3,3,3,4,5,5,7,7,7,19] the returned array will contain [1,2,3,4,5,7,19]
   *
   * \return a newly allocated array that contain the result of the unique operation applied on \a this.
   * \throw if \a this is not allocated or if \a this has not exactly one component.
   * \sa DataArrayInt::buildUniqueNotSorted, DataArrayInt::indexOfSameConsecutiveValueGroups
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::buildUnique() const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildUnique : only single component allowed !");
    std::size_t nbOfElements(this->getNumberOfTuples());
    MCAuto<DataArrayType> tmp(DataArrayType::New());
    tmp->deepCopyFrom(*this);
    T *data(tmp->getPointer());
    T *last(std::unique(data,data+nbOfElements));
    MCAuto<DataArrayType> ret(DataArrayType::New());
    ret->alloc(std::distance(data,last),1);
    std::copy(data,last,ret->getPointer());
    return ret.retn();
  }

  /*!
   * This method can be applied on allocated with one component DataArrayInt instance.
   * This method keep elements only once by keeping the same order in \a this that is not expected to be sorted.
   *
   * \return a newly allocated array that contain the result of the unique operation applied on \a this.
   *
   * \throw if \a this is not allocated or if \a this has not exactly one component.
   *
   * \sa DataArrayInt::buildUnique
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::buildUniqueNotSorted() const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildUniqueNotSorted : only single component allowed !");
    T minVal,maxVal;
    this->getMinMaxValues(minVal,maxVal);
    std::vector<bool> b(maxVal-minVal+1,false);
    const T *ptBg(this->begin()),*endBg(this->end());
    MCAuto<DataArrayType> ret(DataArrayType::New()); ret->alloc(0,1);
    for(const T *pt=ptBg;pt!=endBg;pt++)
      {
        if(!b[*pt-minVal])
          {
            ret->pushBackSilent(*pt);
            b[*pt-minVal]=true;
          }
      }
    ret->copyStringInfoFrom(*this);
    return ret.retn();
  }

  /*!
   * Returns a new DataArrayInt which contains size of every of groups described by \a this
   * "index" array. Such "index" array is returned for example by
   * \ref MEDCoupling::MEDCouplingUMesh::buildDescendingConnectivity
   * "MEDCouplingUMesh::buildDescendingConnectivity" and
   * \ref MEDCoupling::MEDCouplingUMesh::getNodalConnectivityIndex
   * "MEDCouplingUMesh::getNodalConnectivityIndex" etc.
   * This method performs the reverse operation of DataArrayInt::computeOffsetsFull.
   *  \return DataArrayInt * - a new instance of DataArrayInt, whose number of tuples
   *          equals to \a this->getNumberOfComponents() - 1, and number of components is 1.
   *          The caller is to delete this array using decrRef() as it is no more needed.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If \a this->getNumberOfTuples() < 1.
   *
   *  \b Example: <br>
   *         - this contains [1,3,6,7,7,9,15]
   *         - result array contains [2,3,1,0,2,6],
   *          where 2 = 3 - 1, 3 = 6 - 3, 1 = 7 - 6 etc.
   *
   * \sa DataArrayInt::computeOffsetsFull
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::deltaShiftIndex() const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::deltaShiftIndex : only single component allowed !");
    std::size_t nbOfElements=this->getNumberOfTuples();
    if(nbOfElements<1)
      throw INTERP_KERNEL::Exception("DataArrayInt::deltaShiftIndex : 2 tuples at least must be present in 'this' !");
    const T *ptr=this->getConstPointer();
    DataArrayType *ret=DataArrayType::New();
    ret->alloc(nbOfElements-1,1);
    T *out=ret->getPointer();
    std::transform(ptr+1,ptr+nbOfElements,ptr,out,std::minus<T>());
    return ret;
  }

  /*!
   * Modifies \a this one-dimensional array so that value of each element \a x
   * of \a this array (\a a) is computed as \f$ x_i = \sum_{j=0}^{i-1} a[ j ] \f$.
   * Or: for each i>0 new[i]=new[i-1]+old[i-1] for i==0 new[i]=0. Number of tuples
   * and components remains the same.<br>
   * This method is useful for allToAllV in MPI with contiguous policy. This method
   * differs from computeOffsetsFull() in that the number of tuples is \b not changed by
   * this one.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *
   *  \b Example: <br>
   *          - Before \a this contains [3,5,1,2,0,8]
   *          - After \a this contains  [0,3,8,9,11,11]<br>
   *          Note that the last element 19 = 11 + 8 is missing because size of \a this
   *          array is retained and thus there is no space to store the last element.
   */
  template <class T>
  void DataArrayDiscrete<T>::computeOffsets()
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::computeOffsets : only single component allowed !");
    std::size_t nbOfElements=this->getNumberOfTuples();
    if(nbOfElements==0)
      return ;
    T *work=this->getPointer();
    T tmp=work[0];
    work[0]=0;
    for(std::size_t i=1;i<nbOfElements;i++)
      {
        T tmp2=work[i];
        work[i]=work[i-1]+tmp;
        tmp=tmp2;
      }
    this->declareAsNew();
  }

  /*!
   * Modifies \a this one-dimensional array so that value of each element \a x
   * of \a this array (\a a) is computed as \f$ x_i = \sum_{j=0}^{i-1} a[ j ] \f$.
   * Or: for each i>0 new[i]=new[i-1]+old[i-1] for i==0 new[i]=0. Number
   * components remains the same and number of tuples is inceamented by one.<br>
   * This method is useful for allToAllV in MPI with contiguous policy. This method
   * differs from computeOffsets() in that the number of tuples is changed by this one.
   * This method performs the reverse operation of DataArrayInt::deltaShiftIndex.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *
   *  \b Example: <br>
   *          - Before \a this contains [3,5,1,2,0,8]
   *          - After \a this contains  [0,3,8,9,11,11,19]<br>
   * \sa DataArrayInt::deltaShiftIndex
   */
  template <class T>
  void DataArrayDiscrete<T>::computeOffsetsFull()
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::computeOffsetsFull : only single component allowed !");
    std::size_t nbOfElements=this->getNumberOfTuples();
    T *ret=(T *)malloc((nbOfElements+1)*sizeof(T));
    const T *work=this->getConstPointer();
    ret[0]=0;
    for(std::size_t i=0;i<nbOfElements;i++)
      ret[i+1]=work[i]+ret[i];
    this->useArray(ret,true,DeallocType::C_DEALLOC,nbOfElements+1,1);
    this->declareAsNew();
  }

  /*!
   * Returns two new DataArrayInt instances whose contents is computed from that of \a this and \a listOfIds arrays as follows.
   * \a this is expected to be an offset format ( as returned by DataArrayInt::computeOffsetsFull ) that is to say with one component
   * and ** sorted strictly increasingly **. \a listOfIds is expected to be sorted ascendingly (not strictly needed for \a listOfIds).
   * This methods searches in \a this, considered as a set of contiguous \c this->getNumberOfComponents() ranges, all ids in \a listOfIds
   * filling completely one of the ranges in \a this.
   *
   * \param [in] listOfIds a list of ids that has to be sorted ascendingly.
   * \param [out] rangeIdsFetched the range ids fetched
   * \param [out] idsInInputListThatFetch contains the list of ids in \a listOfIds that are \b fully included in a range in \a this. So
   *              \a idsInInputListThatFetch is a part of input \a listOfIds.
   *
   * \sa DataArrayInt::computeOffsetsFull
   *
   *  \b Example: <br>
   *          - \a this : [0,3,7,9,15,18]
   *          - \a listOfIds contains  [0,1,2,3,7,8,15,16,17]
   *          - \a rangeIdsFetched result array: [0,2,4]
   *          - \a idsInInputListThatFetch result array: [0,1,2,7,8,15,16,17]
   * In this example id 3 in input \a listOfIds is alone so it do not appear in output \a idsInInputListThatFetch.
   * <br>
   */
  template <class T>
  void DataArrayDiscrete<T>::findIdsRangesInListOfIds(const DataArrayType *listOfIds, DataArrayIdType *& rangeIdsFetched, DataArrayType *& idsInInputListThatFetch) const
  {
    if(!listOfIds)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdsRangesInListOfIds : input list of ids is null !");
    listOfIds->checkAllocated(); this->checkAllocated();
    if(listOfIds->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdsRangesInListOfIds : input list of ids must have exactly one component !");
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdsRangesInListOfIds : this must have exactly one component !");
    MCAuto<DataArrayIdType> ret0=DataArrayIdType::New(); ret0->alloc(0,1);
    MCAuto<DataArrayType> ret1=DataArrayType::New(); ret1->alloc(0,1);
    const T *tupPtr(listOfIds->begin()), *tupEnd(listOfIds->end());
    const T *offBg(this->begin()),*offEnd(this->end()-1);
    const T *offPtr(offBg);
    while(tupPtr!=tupEnd && offPtr!=offEnd)
      {
        if(*tupPtr==*offPtr)
          {
            T i=offPtr[0];
            while(i<offPtr[1] && *tupPtr==i && tupPtr!=tupEnd) { i++; tupPtr++; }
            if(i==offPtr[1])
              {
                ret0->pushBackSilent(ToIdType(std::distance(offBg,offPtr)));
                ret1->pushBackValsSilent(tupPtr-(offPtr[1]-offPtr[0]),tupPtr);
                offPtr++;
              }
          }
        else
          { if(*tupPtr<*offPtr) tupPtr++; else offPtr++; }
      }
    rangeIdsFetched=ret0.retn();
    idsInInputListThatFetch=ret1.retn();
  }

  /*!
   * Returns a new DataArrayInt whose contents is computed from that of \a this and \a
   * offsets arrays as follows. \a offsets is a one-dimensional array considered as an
   * "index" array of a "iota" array, thus, whose each element gives an index of a group
   * beginning within the "iota" array. And \a this is a one-dimensional array
   * considered as a selector of groups described by \a offsets to include into the result array.
   *  \throw If \a offsets is NULL.
   *  \throw If \a offsets is not allocated.
   *  \throw If \a offsets->getNumberOfComponents() != 1.
   *  \throw If \a offsets is not monotonically increasing.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If any element of \a this is not a valid index for \a offsets array.
   *
   *  \b Example: <br>
   *          - \a this: [0,2,3]
   *          - \a offsets: [0,3,6,10,14,20]
   *          - result array: [0,1,2,6,7,8,9,10,11,12,13] == <br>
   *            \c range(0,3) + \c range(6,10) + \c range(10,14) ==<br>
   *            \c range( \a offsets[ \a this[0] ], offsets[ \a this[0]+1 ]) +
   *            \c range( \a offsets[ \a this[1] ], offsets[ \a this[1]+1 ]) +
   *            \c range( \a offsets[ \a this[2] ], offsets[ \a this[2]+1 ])
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::buildExplicitArrByRanges(const DataArrayType *offsets) const
  {
    if(!offsets)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildExplicitArrByRanges : DataArrayInt pointer in input is NULL !");
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildExplicitArrByRanges : only single component allowed !");
    offsets->checkAllocated();
    if(offsets->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildExplicitArrByRanges : input array should have only single component !");
    mcIdType othNbTuples=offsets->getNumberOfTuples()-1;
    mcIdType nbOfTuples=this->getNumberOfTuples();
    T retNbOftuples=0;
    const T *work=this->getConstPointer();
    const T *offPtr=offsets->getConstPointer();
    for(mcIdType i=0;i<nbOfTuples;i++)
      {
        T val=work[i];
        if(val>=0 && val<othNbTuples)
          {
            T delta=offPtr[val+1]-offPtr[val];
            if(delta>=0)
              retNbOftuples+=delta;
            else
              {
                std::ostringstream oss; oss << "DataArrayInt::buildExplicitArrByRanges : Tuple #" << val << " of offset array has a delta < 0 !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::buildExplicitArrByRanges : Tuple #" << i << " in this contains " << val;
            oss << " whereas offsets array is of size " << othNbTuples+1 << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    MCAuto<DataArrayType> ret=DataArrayType::New();
    ret->alloc(retNbOftuples,1);
    T *retPtr=ret->getPointer();
    for(mcIdType i=0;i<nbOfTuples;i++)
      {
        T val=work[i];
        T start=offPtr[val];
        T off=offPtr[val+1]-start;
        for(T j=0;j<off;j++,retPtr++)
          *retPtr=start+j;
      }
    return ret.retn();
  }

  /*!
   * Returns a new DataArrayInt whose contents is computed using \a this that must be a
   * scaled array (monotonically increasing).
  from that of \a this and \a
   * offsets arrays as follows. \a offsets is a one-dimensional array considered as an
   * "index" array of a "iota" array, thus, whose each element gives an index of a group
   * beginning within the "iota" array. And \a this is a one-dimensional array
   * considered as a selector of groups described by \a offsets to include into the result array.
   *  \throw If \a  is NULL.
   *  \throw If \a this is not allocated.
   *  \throw If \a this->getNumberOfComponents() != 1.
   *  \throw If \a this->getNumberOfTuples() == 0.
   *  \throw If \a this is not monotonically increasing.
   *  \throw If any element of ids in ( \a bg \a stop \a step ) points outside the scale in \a this.
   *
   *  \b Example: <br>
   *          - \a bg , \a stop and \a step : (0,5,2)
   *          - \a this: [0,3,6,10,14,20]
   *          - result array: [0,0,0, 2,2,2,2, 4,4,4,4,4,4] == <br>
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::buildExplicitArrOfSliceOnScaledArr(T bg, T stop, T step) const
  {
    if(!this->isAllocated())
      throw INTERP_KERNEL::Exception("DataArrayInt::buildExplicitArrOfSliceOnScaledArr : not allocated array !");
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildExplicitArrOfSliceOnScaledArr : number of components is expected to be equal to one !");
    mcIdType nbOfTuples(this->getNumberOfTuples());
    if(nbOfTuples==0)
      throw INTERP_KERNEL::Exception("DataArrayInt::buildExplicitArrOfSliceOnScaledArr : number of tuples must be != 0 !");
    const T *ids(this->begin());
    mcIdType nbOfEltsInSlc=DataArrayTools<T>::GetNumberOfItemGivenBESRelative(bg,stop,step,"DataArrayInt::buildExplicitArrOfSliceOnScaledArr");
    T sz(0),pos(bg);
    for(mcIdType i=0;i<nbOfEltsInSlc;i++,pos+=step)
      {
        if(pos>=0 && pos<nbOfTuples-1)
          {
            T delta(ids[pos+1]-ids[pos]);
            sz+=delta;
            if(delta<0)
              {
                std::ostringstream oss; oss << "DataArrayInt::buildExplicitArrOfSliceOnScaledArr : At pos #" << i << " of input slice, value is " << pos << " and at this pos this is not monotonically increasing !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::buildExplicitArrOfSliceOnScaledArr : At pos #" << i << " of input slice, value is " << pos << " should be in [0," << nbOfTuples-1 << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    MCAuto<DataArrayType> ret(DataArrayType::New()); ret->alloc(sz,1);
    T *retPtr(ret->getPointer());
    pos=bg;
    for(mcIdType i=0;i<nbOfEltsInSlc;i++,pos+=step)
      {
        T delta(ids[pos+1]-ids[pos]);
        for(T j=0;j<delta;j++,retPtr++)
          *retPtr=pos;
      }
    return ret.retn();
  }

  /*!
   * Given in input ranges \a ranges, it returns a newly allocated DataArrayInt instance having one component and the same number of tuples than \a this.
   * For each tuple at place **i** in \a this it tells which is the first range in \a ranges that contains value \c this->getIJ(i,0) and put the result
   * in tuple **i** of returned DataArrayInt.
   * If ranges overlapped (in theory it should not) this method do not detect it and always returns the first range.
   *
   * For example if \a this contains : [1,24,7,8,10,17] and \a ranges contains [(0,3),(3,8),(8,15),(15,22),(22,30)]
   * The return DataArrayInt will contain : **[0,4,1,2,2,3]**
   *
   * \param [in] ranges typically come from output of MEDCouplingUMesh::ComputeRangesFromTypeDistribution. Each range is specified like this : 1st component is
   *             for lower value included and 2nd component is the upper value of corresponding range **excluded**.
   * \throw If offsets is a null pointer or does not have 2 components or if \a this is not allocated or \a this do not have exactly one component. To finish an exception
   *        is thrown if no ranges in \a ranges contains value in \a this.
   *
   * \sa DataArrayInt::findIdInRangeForEachTuple
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::findRangeIdForEachTuple(const DataArrayType *ranges) const
  {
    if(!ranges)
      throw INTERP_KERNEL::Exception("DataArrayInt::findRangeIdForEachTuple : null input pointer !");
    if(ranges->getNumberOfComponents()!=2)
      throw INTERP_KERNEL::Exception("DataArrayInt::findRangeIdForEachTuple : input DataArrayInt instance should have 2 components !");
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::findRangeIdForEachTuple : this should have only one component !");
    mcIdType nbTuples(this->getNumberOfTuples());
    MCAuto<DataArrayIdType> ret=DataArrayIdType::New(); ret->alloc(nbTuples,1);
    mcIdType nbOfRanges(ranges->getNumberOfTuples());
    const T *rangesPtr=ranges->getConstPointer();
    mcIdType *retPtr=ret->getPointer();
    const T *inPtr=this->getConstPointer();
    for(mcIdType i=0;i<nbTuples;i++,retPtr++)
      {
        T val=inPtr[i];
        bool found=false;
        for(mcIdType j=0;j<nbOfRanges && !found;j++)
          if(val>=rangesPtr[2*j] && val<rangesPtr[2*j+1])
            { *retPtr=j; found=true; }
        if(found)
          continue;
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::findRangeIdForEachTuple : tuple #" << i << " not found by any ranges !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    return ret.retn();
  }

  /*!
   * Given in input ranges \a ranges, it returns a newly allocated DataArrayInt instance having one component and the same number of tuples than \a this.
   * For each tuple at place **i** in \a this it tells which is the sub position of the first range in \a ranges that contains value \c this->getIJ(i,0) and put the result
   * in tuple **i** of returned DataArrayInt.
   * If ranges overlapped (in theory it should not) this method do not detect it and always returns the sub position of the first range.
   *
   * For example if \a this contains : [1,24,7,8,10,17] and \a ranges contains [(0,3),(3,8),(8,15),(15,22),(22,30)]
   * The return DataArrayInt will contain : **[1,2,4,0,2,2]**
   * This method is often called in pair with DataArrayInt::findRangeIdForEachTuple method.
   *
   * \param [in] ranges typically come from output of MEDCouplingUMesh::ComputeRangesFromTypeDistribution. Each range is specified like this : 1st component is
   *             for lower value included and 2nd component is the upper value of corresponding range **excluded**.
   * \throw If offsets is a null pointer or does not have 2 components or if \a this is not allocated or \a this do not have exactly one component. To finish an exception
   *        is thrown if no ranges in \a ranges contains value in \a this.
   * \sa DataArrayInt::findRangeIdForEachTuple
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::findIdInRangeForEachTuple(const DataArrayType *ranges) const
  {
    if(!ranges)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdInRangeForEachTuple : null input pointer !");
    if(ranges->getNumberOfComponents()!=2)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdInRangeForEachTuple : input DataArrayInt instance should have 2 components !");
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::findIdInRangeForEachTuple : this should have only one component !");
    mcIdType nbTuples=this->getNumberOfTuples();
    MCAuto<DataArrayType> ret=DataArrayType::New(); ret->alloc(nbTuples,1);
    mcIdType nbOfRanges=ranges->getNumberOfTuples();
    const T *rangesPtr=ranges->getConstPointer();
    T *retPtr=ret->getPointer();
    const T *inPtr=this->getConstPointer();
    for(mcIdType i=0;i<nbTuples;i++,retPtr++)
      {
        T val=inPtr[i];
        bool found=false;
        for(mcIdType j=0;j<nbOfRanges && !found;j++)
          if(val>=rangesPtr[2*j] && val<rangesPtr[2*j+1])
            { *retPtr=val-rangesPtr[2*j]; found=true; }
        if(found)
          continue;
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::findIdInRangeForEachTuple : tuple #" << i << " not found by any ranges !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    return ret.retn();
  }

  /*!
   * \b WARNING this method is a \b non \a const \b method. This method works tuple by tuple. Each tuple is expected to be pairs (number of components must be equal to 2).
   * This method rearrange each pair in \a this so that, tuple with id \b tid will be after the call \c this->getIJ(tid,0)==this->getIJ(tid-1,1) and \c this->getIJ(tid,1)==this->getIJ(tid+1,0).
   * If it is impossible to reach such condition an exception will be thrown ! \b WARNING In case of throw \a this can be partially modified !
   * If this method has correctly worked, \a this will be able to be considered as a linked list.
   * This method does nothing if number of tuples is lower of equal to 1.
   *
   * This method is useful for users having an unstructured mesh having only SEG2 to rearrange internally the connectivity without any coordinates consideration.
   *
   * \sa MEDCouplingUMesh::orderConsecutiveCells1D, DataArrayInt::sortToHaveConsecutivePairs(), DataArrayInt::fromLinkedListOfPairToList
   */
  template <class T>
  void DataArrayDiscrete<T>::sortEachPairToMakeALinkedList()
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=2)
      throw INTERP_KERNEL::Exception("DataArrayInt::sortEachPairToMakeALinkedList : Only works on DataArrayInt instance with nb of components equal to 2 !");
    mcIdType nbOfTuples(this->getNumberOfTuples());
    if(nbOfTuples<=1)
      return ;
    T *conn(this->getPointer());
    for(mcIdType i=1;i<nbOfTuples;i++,conn+=2)
      {
        if(i>1)
          {
            if(conn[2]==conn[3])
              {
                std::ostringstream oss; oss << "DataArrayInt::sortEachPairToMakeALinkedList : In the tuple #" << i << " presence of a pair filled with same ids !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
            if(conn[2]!=conn[1] && conn[3]!=conn[1])
              {
                std::ostringstream oss; oss << "DataArrayInt::sortEachPairToMakeALinkedList : There is no common noeud between the tuple #" << i << " and the tuple #" << i + 1 << ". Need call DataArrayInt::sortToHaveConsecutivePairs() ";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
            if(conn[2]!=conn[1] && conn[3]==conn[1] && conn[2]!=conn[0])
              std::swap(conn[2],conn[3]);
            //not(conn[2]==conn[1] && conn[3]!=conn[1] && conn[3]!=conn[0])
            if(conn[2]!=conn[1] || conn[3]==conn[1] || conn[3]==conn[0])
              {
                std::ostringstream oss; oss << "DataArrayInt::sortEachPairToMakeALinkedList : In the tuple #" << i << " something is invalid !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
          }
        else
          {
            if(conn[0]==conn[1] || conn[2]==conn[3])
              throw INTERP_KERNEL::Exception("DataArrayInt::sortEachPairToMakeALinkedList : In the 2 first tuples presence of a pair filled with same ids !");
            T tmp[4];
            std::set<T> s;
            s.insert(conn,conn+4);
            if(s.size()!=3)
              throw INTERP_KERNEL::Exception("DataArrayInt::sortEachPairToMakeALinkedList : This can't be considered as a linked list regarding 2 first tuples !");
            if(std::count(conn,conn+4,conn[0])==2)
              {
                tmp[0]=conn[1];
                tmp[1]=conn[0];
                tmp[2]=conn[0];
                if(conn[2]==conn[0])
                  { tmp[3]=conn[3]; }
                else
                  { tmp[3]=conn[2];}
                std::copy(tmp,tmp+4,conn);
              }
            else
              {//here we are sure to have (std::count(conn,conn+4,conn[1])==2)
                if(conn[1]==conn[3])
                  std::swap(conn[2],conn[3]);
              }
          }
      }
  }

  /*!
   * This method is the improvement from the method sortEachPairToMakeALinkedList().
   *
   * \sa MEDCouplingUMesh::orderConsecutiveCells1D, DataArrayInt::sortEachPairToMakeALinkedList(), DataArrayInt::fromLinkedListOfPairToList
   */
  template <class T>
  void DataArrayDiscrete<T>::sortToHaveConsecutivePairs()
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=2)
      throw INTERP_KERNEL::Exception("DataArrayInt::sortToHaveConsecutivePairs : Only works on DataArrayInt instance with nb of components equal to 2 !");
    mcIdType nbOfTuples(this->getNumberOfTuples());
    T *thisPtr(this->getPointer());
    mcIdType idOfLastTuple = 0;
    std::pair<T*,T*> tmp;
    if(thisPtr[0]==thisPtr[1])
      {
        THROW_IK_EXCEPTION("DataArrayInt::sortToHaveConsecutivePairs : In the first tuple presence of a pair filled with same ids !");
      }
    while(idOfLastTuple < nbOfTuples-1)
    {
      mcIdType i = idOfLastTuple+1;
      tmp.first = &thisPtr[2*i];
      tmp.second = &thisPtr[2*i+1];
      if(std::get<0>(tmp)[0] == std::get<1>(tmp)[0])
        {
          THROW_IK_EXCEPTION("DataArrayInt::sortToHaveConsecutivePairs : In the tuple #" << i << " presence of a pair filled with same ids !");
        }
      while((this->getIJ(idOfLastTuple,1) != std::get<0>(tmp)[0] &&
            this->getIJ(idOfLastTuple,1) != std::get<1>(tmp)[0]) &&
            i < nbOfTuples-1)
        {
          std::swap(std::get<0>(tmp)[0],thisPtr[2*(i+1)]);
          std::swap(std::get<1>(tmp)[0],thisPtr[2*(i+1)+1]);
          i++;
        }
      if(i < nbOfTuples-1 ||
         this->getIJ(idOfLastTuple,1) == std::get<0>(tmp)[0] ||
         this->getIJ(idOfLastTuple,1) == std::get<1>(tmp)[0])
      {
        if(this->getIJ(idOfLastTuple,1) != std::get<0>(tmp)[0])
          std::swap(std::get<0>(tmp)[0],std::get<1>(tmp)[0]);
        idOfLastTuple++;
      }
      else
      {
        THROW_IK_EXCEPTION("DataArrayInt::sortToHaveConsecutivePairs : not found the tuple which have the common noeud = " << this->getIJ(idOfLastTuple,1));
      }
    }
  }

  /*!
   * \a this is expected to be a correctly linked list of pairs.
   *
   * \sa DataArrayInt::sortEachPairToMakeALinkedList
   */
  template <class T>
  MCAuto<typename Traits<T>::ArrayType> DataArrayDiscrete<T>::fromLinkedListOfPairToList() const
  {
    this->checkAllocated();
    this->checkNbOfComps(2,"DataArrayInt::fromLinkedListOfPairToList : this is expected to have 2 components");
    mcIdType nbTuples(this->getNumberOfTuples());
    if(nbTuples<1)
      throw INTERP_KERNEL::Exception("DataArrayInt::fromLinkedListOfPairToList : no tuples in this ! Not a linked list !");
    MCAuto<DataArrayType> ret(DataArrayType::New()); ret->alloc(nbTuples+1,1);
    const T *thisPtr(this->begin());
    T *retPtr(ret->getPointer());
    retPtr[0]=thisPtr[0];
    for(mcIdType i=0;i<nbTuples;i++)
      {
        retPtr[i+1]=thisPtr[2*i+1];
        if(i<nbTuples-1)
          if(thisPtr[2*i+1]!=thisPtr[2*(i+1)+0])
            {
              std::ostringstream oss; oss << "DataArrayInt::fromLinkedListOfPairToList : this is not a proper linked list of pair. The link is broken between tuple #" << i << " and tuple #" << i+1 << " ! Call sortEachPairToMakeALinkedList ?";
              throw INTERP_KERNEL::Exception(oss.str());
            }
      }
    return ret;
  }

  /*!
   * This method returns all different values found in \a this. This method throws if \a this has not been allocated.
   * But the number of components can be different from one.
   * \return a newly allocated array (that should be dealt by the caller) containing different values in \a this.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::getDifferentValues() const
  {
    this->checkAllocated();
    std::set<T> ret;
    ret.insert(this->begin(),this->end());
    MCAuto<DataArrayType> ret2=DataArrayType::New();
    ret2->alloc(ret.size(),1);
    std::copy(ret.begin(),ret.end(),ret2->getPointer());
    return ret2.retn();
  }

  template<class T>
  class PartitionCfg : public std::set<T>
  {
    public:
      PartitionCfg() = default;
      void setID(T id) { _id = id; }
      T getID() const { return _id; }
      bool operator<(const PartitionCfg<T>& other) { return std::set<T>::operator<( other._data ); }
      bool operator>(const PartitionCfg<T>& other) { return std::set<T>::operator>( other._data ); }
      bool operator==(const PartitionCfg<T>& other) { return std::set<T>::operator==( other._data ); }
    private:
      T _id = 0;
  };

  /*!
   * \a this is considered as an array defining a partition. It means that values in \a this define partition-ids. All tuples in \a this with same value can be considered as partition.
   * Typically MED stores families uses this compact storage method.
   *
   * This method computes and returns the partition stored in \a this on reduced space using a reduction operation specified by pairs \a commonEntities and \a commonEntitiesIndex parameters.
   * The reduction operation is the consequence of a fusion of enties. It explains the storage format defining reduction.
   *
   * An another way to consider this method : This method can be seen as an interpolation on integer field (\a this). For fused entities it returns in compact way all combinations of partition configuration.
   *
   * \param [out] partitionsToBeModified - For all partitions impacted by the reduction a n-uplet is returned (n is deduced thanks to \a partitionsToBeModifiedIndex)
   *                                       Each n-uplet represents a list of partitions impacted by reduction. The first element of n-uplet represents the ID to locate entities (using for example findIdsEqual in returned array)
   *                                       Remaining IDs in n-uplet are Partition IDs impacted.
   * \param [out] partitionsToBeModifiedIndex - Index attached to \a partitionsToBeModified to interprete.
   *
   * \return - The partition array that is the output of the reduction of \a this.
   *
   * \sa MEDCouplingUMesh::findCommonCells, DataArrayDouble::findCommonTuples
   */
  template <class T>
  MCAuto< typename DataArrayDiscrete<T>::DataArrayType > DataArrayDiscrete<T>::forThisAsPartitionBuildReduction(const MCAuto<DataArrayIdType>& commonEntities, const MCAuto<DataArrayIdType>& commonEntitiesIndex, MCAuto<DataArrayType>& partitionsToBeModified, MCAuto<DataArrayIdType>& partitionsToBeModifiedIndex) const
  {
    constexpr char MSG[] = "this should have only one component";
    this->checkAllocated(); this->checkNbOfComps(1,MSG);
    commonEntities->checkAllocated(); commonEntities->checkNbOfComps(1,MSG);
    commonEntitiesIndex->checkAllocated(); commonEntitiesIndex->checkNbOfComps(1,MSG);
    std::size_t initialSpaceSz( this->getNumberOfTuples() );
    mcIdType returnedSpaceSz( 0 );// store size of reducted returned size
    //
    MCAuto<DataArrayIdType> sizeOfPacks( commonEntitiesIndex->deltaShiftIndex() );
    //
    MCAuto<DataArrayIdType> o2n( DataArrayIdType::ConvertIndexArrayToO2N(initialSpaceSz,commonEntities->begin(),commonEntitiesIndex->begin(),commonEntitiesIndex->end(),returnedSpaceSz) );
    MCAuto< DataArrayType > ret( DataArrayDiscrete<T>::New() );
    ret->alloc( returnedSpaceSz, 1 );
    ret->fillWithValue( std::numeric_limits<T>::max() );
    // First deal with entities not fused.
    MCAuto<DataArrayIdType> eltsNotFused( commonEntities->copySorted() );
    eltsNotFused = eltsNotFused->buildComplement( initialSpaceSz );
    MCAuto<DataArrayType> partionIdsOfNotFused = this->mySelectByTupleIdSafe(eltsNotFused->begin(),eltsNotFused->end());
    MCAuto<DataArrayIdType> tupleIdsToSelect = o2n->selectByTupleIdSafe(eltsNotFused->begin(),eltsNotFused->end());
    ret->setPartOfValues3(partionIdsOfNotFused,tupleIdsToSelect->begin(),tupleIdsToSelect->end(),0,1,1);
    T *retPt( ret->getPointer() );
    const mcIdType *o2nPt( o2n->begin() );
    //
    partitionsToBeModified = DataArrayType::New(); partitionsToBeModified->alloc(0,1);
    partitionsToBeModifiedIndex = DataArrayIdType::New(); partitionsToBeModifiedIndex->alloc(1,1); partitionsToBeModifiedIndex->setIJSilent(0,0,0);
    if( !sizeOfPacks->empty() )// if empty -> no fusion -> ret is already ready at this point -> nothing to do.
    {// ready to work
      mcIdType maxSizeOfPacks = sizeOfPacks->getMaxValueInArray();// store the max size of common entities
      T newValueInThis = this->getMaxValueInArray() + 1;
      const T *ptOfThisData( this->begin() );
      const mcIdType *ceBg( commonEntities->begin() ), *ceiBg( commonEntitiesIndex->begin() );
      for(mcIdType szOfPack = 2 ; szOfPack <= maxSizeOfPacks ; ++szOfPack)
      {
        MCAuto<DataArrayIdType> idsInThisWithSamePackSz = sizeOfPacks->findIdsEqual( FromIdType<T>( szOfPack ) );
        std::set< PartitionCfg<T> > partitionCfgHolder;
        for( const mcIdType *idsInThisWithSamePackSzIt = idsInThisWithSamePackSz->begin() ; idsInThisWithSamePackSzIt != idsInThisWithSamePackSz->end() ; ++idsInThisWithSamePackSzIt )
        {
          PartitionCfg<T> partitionCfg;
          std::transform(ceBg + ceiBg[*idsInThisWithSamePackSzIt],ceBg + ceiBg[*idsInThisWithSamePackSzIt + 1], std::inserter(partitionCfg,partitionCfg.end()),[ptOfThisData](mcIdType elt) { return ptOfThisData[elt]; });
          auto existCfg = partitionCfgHolder.find( partitionCfg );
          if( existCfg != partitionCfgHolder.end() )
          {//partition already exist by a previous pack -> reuse it !
            T newPartitionID = existCfg->getID();
            retPt[ o2nPt[ ceBg [ ceiBg[*idsInThisWithSamePackSzIt] ] ] ] = newPartitionID;// hypothesis that o2n is so that all o2n[ceBg + ceiBg[*idsInThisWithSamePackSzIt],ceBg + ceiBg[*idsInThisWithSamePackSzIt + 1]) points to the same point in new renumbering
          }
          else
          {//partition does not exist yet -> create it !
            partitionCfg.setID( newValueInThis++ );
            partitionCfgHolder.insert( partitionCfg );
            retPt[ o2nPt[ ceBg [ ceiBg[*idsInThisWithSamePackSzIt] ] ] ] = partitionCfg.getID();
            partitionsToBeModified->pushBackSilent( partitionCfg.getID() );
            partitionsToBeModified->pushBackValsSilent(partitionCfg.begin(),partitionCfg.end());
            partitionsToBeModifiedIndex->pushBackSilent( partitionsToBeModifiedIndex->back() + partitionCfg.size() + 1 );
          }
        }
      }
    }
    return ret;
  }

  /*!
   * This method constructs for a list of pairs of ids in this to a pair arr,arrIndex (typically followed by a ConvertIndexArrayToO2N call). This method is useful to 
   * convert in context of // computing to move to global ids approach.
   * 
   * \a this is expected to be a 2 components array.
   * \b Warning, the pairs are expected to be sorted. (3,5), (0,7), (5,7), (0,3) wil. If not the output can
   * 
   * 
   *\b Example: <br>
   * - \a this          : [0,3, 5,4, 0,7]
   * - \a arrOut        : [ 0,3,7, 4,5 ]
   * - \a arrIndexOut       : [0,3,5]
   * 
   * \sa ConvertIndexArrayToO2N
   */
  template <class T>
  void DataArrayDiscrete<T>::fromListOfPairsToIndexArray(MCAuto<DataArrayType>& arrOut, MCAuto<DataArrayIdType>& arrIndexOut) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=2)
      throw INTERP_KERNEL::Exception("DataArrayInt::fromListOfPairsToIndexArray : this should have exactly two components !");
    std::map< T,std::set<T> > workStruct;
    mcIdType nbTuples( this->getNumberOfTuples() );
    const T *curTuple( this->begin() );
    for( mcIdType i = 0 ; i < nbTuples ; ++i )
    {
      T zeMin ( std::min(curTuple[0],curTuple[1]) ), zeMax( std::max( curTuple[0],curTuple[1] ) );
      auto it = workStruct.find( zeMin );
      if( it != workStruct.end() )
      {
        (*it).second.insert( zeMax );
      }
      else
      {
        auto it2 = workStruct.find( zeMax );
        if( it2 == workStruct.end() )
        {
          workStruct[zeMin] = {zeMin, zeMax};
        }
        else
        {
          std::set<T> tmp( std::move( (*it2).second ) );
          tmp.insert( zeMin );
          workStruct.erase( it2 );
          workStruct[ zeMin ] = tmp;
        }
      }
      curTuple += 2;
    }
    // put result in arrays
    mcIdType nbOfGroups( ToIdType( workStruct.size() ) );
    arrIndexOut = DataArrayIdType::New(); arrIndexOut->alloc(nbOfGroups + 1,1);
    mcIdType *arrIndexOutPt = arrIndexOut->getPointer(); *arrIndexOutPt = 0;
    for( const auto& it : workStruct )
    {
      mcIdType nbElemsInCurGrp( ToIdType( it.second.size() ) );
      arrIndexOutPt[1] = arrIndexOutPt[0] + nbElemsInCurGrp;
      arrIndexOutPt++;
    }
    arrOut = DataArrayType::New(); arrOut->alloc(*arrIndexOutPt,1);
    T *arrOutPt( arrOut->getPointer() );
    for( const auto& it : workStruct )
    {
      arrOutPt = std::copy(it.second.begin(),it.second.end(),arrOutPt);
    }
  }

  /*!
   * This method is a refinement of DataArrayInt::getDifferentValues because it returns not only different values in \a this but also, for each of
   * them it tells which tuple id have this id.
   * This method works only on arrays with one component (if it is not the case call DataArrayInt::rearrange(1) ).
   * This method returns two arrays having same size.
   * The instances of DataArrayInt in the returned vector have be specially allocated and computed by this method. Each of them should be dealt by the caller of this method.
   * Example : if this is equal to [1,0,1,2,0,2,2,-3,2] -> differentIds=[-3,0,1,2] and returned array will be equal to [[7],[1,4],[0,2],[3,5,6,8]]
   */
  template <class T>
  std::vector<DataArrayIdType *> DataArrayDiscrete<T>::partitionByDifferentValues(std::vector<T>& differentIds) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::partitionByDifferentValues : this should have only one component !");
    mcIdType id=0;
    std::map<T,mcIdType> m,m2,m3;
    for(const T *w=this->begin();w!=this->end();w++)
      m[*w]++;
    differentIds.resize(m.size());
    std::vector<DataArrayIdType *> ret(m.size());
    std::vector<mcIdType *> retPtr(m.size());
    for(typename std::map<T,mcIdType>::const_iterator it=m.begin();it!=m.end();it++,id++)
      {
        m2[(*it).first]=id;
        ret[id]=DataArrayIdType::New();
        ret[id]->alloc((*it).second,1);
        retPtr[id]=ret[id]->getPointer();
        differentIds[id]=(*it).first;
      }
    id=0;
    for(const T *w=this->begin();w!=this->end();w++,id++)
      {
        retPtr[m2[*w]][m3[*w]++]=id;
      }
    return ret;
  }

  /*!
   * This method split ids in [0, \c this->getNumberOfTuples() ) using \a this array as a field of weight (>=0 each).
   * The aim of this method is to return a set of \a nbOfSlices chunk of contiguous ids as balanced as possible.
   *
   * \param [in] nbOfSlices - number of slices expected.
   * \return - a vector having a size equal to \a nbOfSlices giving the start (included) and the stop (excluded) of each chunks.
   *
   * \sa DataArray::GetSlice
   * \throw If \a this is not allocated or not with exactly one component.
   * \throw If an element in \a this if < 0.
   */
  template <class T>
  std::vector< std::pair<mcIdType,mcIdType> > DataArrayDiscrete<T>::splitInBalancedSlices(mcIdType nbOfSlices) const
  {
    if(!this->isAllocated() || this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::splitInBalancedSlices : this array should have number of components equal to one and must be allocated !");
    if(nbOfSlices<=0)
      throw INTERP_KERNEL::Exception("DataArrayInt::splitInBalancedSlices : number of slices must be >= 1 !");
    T sum(this->accumulate((std::size_t)0));
    mcIdType nbOfTuples(this->getNumberOfTuples());
    T sumPerSlc(sum/FromIdType<T>(nbOfSlices));
    mcIdType pos(0);
    const T *w(this->begin());
    std::vector< std::pair<mcIdType,mcIdType> > ret(nbOfSlices);
    for(mcIdType i=0;i<nbOfSlices;i++)
      {
        std::pair<mcIdType, mcIdType> p(pos,-1);
        T locSum(0);
        while(locSum<sumPerSlc && pos<nbOfTuples) { pos++; locSum+=*w++; }
        if(i!=nbOfSlices-1)
          p.second=pos;
        else
          p.second=nbOfTuples;
        ret[i]=p;
      }
    return ret;
  }

  /*!
   * Modify \a this array so that each value becomes a modulus of division of this value by
   * a value of another DataArrayInt. There are 3 valid cases.
   * 1.  The arrays have same number of tuples and components. Then each value of
   *    \a this array is divided by the corresponding value of \a other one, i.e.:
   *   _a_ [ i, j ] %= _other_ [ i, j ].
   * 2.  The arrays have same number of tuples and \a other array has one component. Then
   *   _a_ [ i, j ] %= _other_ [ i, 0 ].
   * 3.  The arrays have same number of components and \a other array has one tuple. Then
   *   _a_ [ i, j ] %= _a2_ [ 0, j ].
   *
   *  \warning No check of division by zero is performed!
   *  \param [in] other - a divisor array.
   *  \throw If \a other is NULL.
   *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples() and
   *         \a this->getNumberOfComponents() != \a other->getNumberOfComponents() and
   *         \a other has number of both tuples and components not equal to 1.
   */
  template <class T>
  void DataArrayDiscrete<T>::modulusEqual(const DataArrayType *other)
  {
    if(!other)
      throw INTERP_KERNEL::Exception("DataArrayInt::modulusEqual : input DataArrayInt instance is NULL !");
    const char *msg="Nb of tuples mismatch for DataArrayInt::modulusEqual !";
    this->checkAllocated(); other->checkAllocated();
    mcIdType nbOfTuple(this->getNumberOfTuples());
    mcIdType nbOfTuple2(other->getNumberOfTuples());
    std::size_t nbOfComp(this->getNumberOfComponents());
    std::size_t nbOfComp2(other->getNumberOfComponents());
    if(nbOfTuple==nbOfTuple2)
      {
        if(nbOfComp==nbOfComp2)
          {
            std::transform(this->begin(),this->end(),other->begin(),this->getPointer(),std::modulus<T>());
          }
        else if(nbOfComp2==1)
          {
            if(nbOfComp2==nbOfComp)
              {
                T *ptr=this->getPointer();
                const T *ptrc=other->getConstPointer();
                for(mcIdType i=0;i<nbOfTuple;i++)
                  std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptr+i*nbOfComp,std::bind(std::modulus<T>(),std::placeholders::_1,*ptrc++));
              }
            else
              throw INTERP_KERNEL::Exception(msg);
          }
        else
          throw INTERP_KERNEL::Exception(msg);
      }
    else if(nbOfTuple2==1)
      {
        T *ptr=this->getPointer();
        const T *ptrc=other->getConstPointer();
        for(mcIdType i=0;i<nbOfTuple;i++)
          std::transform(ptr+i*nbOfComp,ptr+(i+1)*nbOfComp,ptrc,ptr+i*nbOfComp,std::modulus<T>());
      }
    else
      throw INTERP_KERNEL::Exception(msg);
    this->declareAsNew();
  }

  /*!
   * Apply pow on values of another DataArrayInt to values of \a this one.
   *
   *  \param [in] other - an array to pow to \a this one.
   *  \throw If \a other is NULL.
   *  \throw If \a this->getNumberOfTuples() != \a other->getNumberOfTuples()
   *  \throw If \a this->getNumberOfComponents() != 1 or \a other->getNumberOfComponents() != 1
   *  \throw If there is a negative value in \a other.
   */
  template <class T>
  void DataArrayDiscrete<T>::powEqual(const DataArrayType *other)
  {
    if(!other)
      throw INTERP_KERNEL::Exception("DataArrayInt::powEqual : input instance is null !");
    mcIdType nbOfTuple=this->getNumberOfTuples();
    mcIdType nbOfTuple2=other->getNumberOfTuples();
    std::size_t nbOfComp=this->getNumberOfComponents();
    std::size_t nbOfComp2=other->getNumberOfComponents();
    if(nbOfTuple!=nbOfTuple2)
      throw INTERP_KERNEL::Exception("DataArrayInt::powEqual : number of tuples mismatches !");
    if(nbOfComp!=1 || nbOfComp2!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::powEqual : number of components of both arrays must be equal to 1 !");
    T *ptr=this->getPointer();
    const T *ptrc=other->begin();
    for(mcIdType i=0;i<nbOfTuple;i++,ptrc++,ptr++)
      {
        if(*ptrc>=0)
          {
            T tmp=1;
            for(T j=0;j<*ptrc;j++)
              tmp*=*ptr;
            *ptr=tmp;
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::powEqual : on tuple #" << i << " of other value is < 0 (" << *ptrc << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    this->declareAsNew();
  }

  ////////////////////////////////////
  /*!
   * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
   * Server side.
   */
  template <class T>
  void DataArrayDiscrete<T>::getTinySerializationIntInformation(std::vector<mcIdType>& tinyInfo) const
  {
    tinyInfo.resize(2);
    if(this->isAllocated())
      {
        tinyInfo[0]=this->getNumberOfTuples();
        tinyInfo[1]=ToIdType(this->getNumberOfComponents());
      }
    else
      {
        tinyInfo[0]=-1;
        tinyInfo[1]=-1;
      }
  }

  /*!
   * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
   * Server side.
   */
  template <class T>
  void DataArrayDiscrete<T>::getTinySerializationStrInformation(std::vector<std::string>& tinyInfo) const
  {
    if(this->isAllocated())
      {
        std::size_t nbOfCompo(this->getNumberOfComponents());
        tinyInfo.resize(nbOfCompo+1);
        tinyInfo[0]=this->getName();
        for(std::size_t i=0;i<nbOfCompo;i++)
          tinyInfo[i+1]=this->getInfoOnComponent(i);
      }
    else
      {
        tinyInfo.resize(1);
        tinyInfo[0]=this->getName();
      }
  }

  /*!
   * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
   * This method returns if a feeding is needed.
   */
  template <class T>
  bool DataArrayDiscrete<T>::resizeForUnserialization(const std::vector<mcIdType>& tinyInfoI)
  {
    mcIdType nbOfTuple=tinyInfoI[0];
    mcIdType nbOfComp=tinyInfoI[1];
    if(nbOfTuple!=-1 || nbOfComp!=-1)
      {
        this->alloc(nbOfTuple,nbOfComp);
        return true;
      }
    return false;
  }

  /*!
   * Useless method for end user. Only for MPI/Corba/File serialsation for multi arrays class.
   * This method returns if a feeding is needed.
   */
  template <class T>
  void DataArrayDiscrete<T>::finishUnserialization(const std::vector<mcIdType>& tinyInfoI, const std::vector<std::string>& tinyInfoS)
  {
    this->setName(tinyInfoS[0]);
    if(this->isAllocated())
      {
        mcIdType nbOfCompo=tinyInfoI[1];
        for(mcIdType i=0;i<nbOfCompo;i++)
          this->setInfoOnComponent(i,tinyInfoS[i+1]);
      }
  }

  ////////////////////////////////////

  /*!
   * Returns a new DataArrayInt that is the result of pow of two given arrays. There are 3
   * valid cases.
   *
   *  \param [in] a1 - an array to pow up.
   *  \param [in] a2 - another array to sum up.
   *  \return DataArrayInt * - the new instance of DataArrayInt.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If either \a a1 or \a a2 is NULL.
   *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples()
   *  \throw If \a a1->getNumberOfComponents() != 1 or \a a2->getNumberOfComponents() != 1.
   *  \throw If there is a negative value in \a a2.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::Pow(const DataArrayType *a1, const DataArrayType *a2)
  {
    if(!a1 || !a2)
      throw INTERP_KERNEL::Exception("DataArrayInt::Pow : at least one of input instances is null !");
    mcIdType nbOfTuple=a1->getNumberOfTuples();
    mcIdType nbOfTuple2=a2->getNumberOfTuples();
    std::size_t nbOfComp=a1->getNumberOfComponents();
    std::size_t nbOfComp2=a2->getNumberOfComponents();
    if(nbOfTuple!=nbOfTuple2)
      throw INTERP_KERNEL::Exception("DataArrayInt::Pow : number of tuples mismatches !");
    if(nbOfComp!=1 || nbOfComp2!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::Pow : number of components of both arrays must be equal to 1 !");
    MCAuto<DataArrayType> ret=DataArrayType::New(); ret->alloc(nbOfTuple,1);
    const T *ptr1(a1->begin()),*ptr2(a2->begin());
    T *ptr=ret->getPointer();
    for(mcIdType i=0;i<nbOfTuple;i++,ptr1++,ptr2++,ptr++)
      {
        if(*ptr2>=0)
          {
            T tmp=1;
            for(T j=0;j<*ptr2;j++)
              tmp*=*ptr1;
            *ptr=tmp;
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::Pow : on tuple #" << i << " of a2 value is < 0 (" << *ptr2 << ") !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    return ret.retn();
  }

  /*!
   * Returns a new DataArrayInt that is a modulus of two given arrays. There are 3
   * valid cases.
   * 1.  The arrays have same number of tuples and components. Then each value of
   *   the result array (_a_) is a division of the corresponding values of \a a1 and
   *   \a a2, i.e.: _a_ [ i, j ] = _a1_ [ i, j ] % _a2_ [ i, j ].
   * 2.  The arrays have same number of tuples and one array, say _a2_, has one
   *   component. Then
   *   _a_ [ i, j ] = _a1_ [ i, j ] % _a2_ [ i, 0 ].
   * 3.  The arrays have same number of components and one array, say _a2_, has one
   *   tuple. Then
   *   _a_ [ i, j ] = _a1_ [ i, j ] % _a2_ [ 0, j ].
   *
   * Info on components is copied either from the first array (in the first case) or from
   * the array with maximal number of elements (getNbOfElems()).
   *  \warning No check of division by zero is performed!
   *  \param [in] a1 - a dividend array.
   *  \param [in] a2 - a divisor array.
   *  \return DataArrayInt * - the new instance of DataArrayInt.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If either \a a1 or \a a2 is NULL.
   *  \throw If \a a1->getNumberOfTuples() != \a a2->getNumberOfTuples() and
   *         \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents() and
   *         none of them has number of tuples or components equal to 1.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::Modulus(const DataArrayType *a1, const DataArrayType *a2)
  {
    if(!a1 || !a2)
      throw INTERP_KERNEL::Exception("DataArrayInt::Modulus : input DataArrayInt instance is NULL !");
    mcIdType nbOfTuple1(a1->getNumberOfTuples());
    mcIdType nbOfTuple2(a2->getNumberOfTuples());
    std::size_t nbOfComp1(a1->getNumberOfComponents());
    std::size_t nbOfComp2(a2->getNumberOfComponents());
    if(nbOfTuple2==nbOfTuple1)
      {
        if(nbOfComp1==nbOfComp2)
          {
            MCAuto<DataArrayType> ret=DataArrayType::New();
            ret->alloc(nbOfTuple2,nbOfComp1);
            std::transform(a1->begin(),a1->end(),a2->begin(),ret->getPointer(),std::modulus<T>());
            ret->copyStringInfoFrom(*a1);
            return ret.retn();
          }
        else if(nbOfComp2==1)
          {
            MCAuto<DataArrayType> ret=DataArrayType::New();
            ret->alloc(nbOfTuple1,nbOfComp1);
            const T *a2Ptr=a2->getConstPointer();
            const T *a1Ptr=a1->getConstPointer();
            T *res=ret->getPointer();
            for(mcIdType i=0;i<nbOfTuple1;i++)
              res=std::transform(a1Ptr+i*nbOfComp1,a1Ptr+(i+1)*nbOfComp1,res,std::bind(std::modulus<T>(),std::placeholders::_1,a2Ptr[i]));
            ret->copyStringInfoFrom(*a1);
            return ret.retn();
          }
        else
          {
            a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Modulus !");
            return 0;
          }
      }
    else if(nbOfTuple2==1)
      {
        a1->checkNbOfComps(nbOfComp2,"Nb of components mismatch for array Modulus !");
        MCAuto<DataArrayType> ret=DataArrayType::New();
        ret->alloc(nbOfTuple1,nbOfComp1);
        const T *a1ptr=a1->getConstPointer(),*a2ptr=a2->getConstPointer();
        T *pt=ret->getPointer();
        for(mcIdType i=0;i<nbOfTuple1;i++)
          pt=std::transform(a1ptr+i*nbOfComp1,a1ptr+(i+1)*nbOfComp1,a2ptr,pt,std::modulus<T>());
        ret->copyStringInfoFrom(*a1);
        return ret.retn();
      }
    else
      {
        a1->checkNbOfTuples(nbOfTuple2,"Nb of tuples mismatch for array Modulus !");//will always throw an exception
        return 0;
      }
  }

  /*!
   * This method tries to find the permutation to apply to the first input \a ids1 to obtain the same array (without considering strings information) the second
   * input array \a ids2.
   * \a ids1 and \a ids2 are expected to be both a list of ids (both with number of components equal to one) not sorted and with values that can be negative.
   * This method will throw an exception is no such permutation array can be obtained. It is typically the case if there is some ids in \a ids1 not in \a ids2 or
   * inversely.
   * In case of success both assertion will be true (no throw) :
   * \c ids1->renumber(ret)->isEqual(ids2) where \a ret is the return of this method.
   * \c ret->transformWithIndArr(ids2)->isEqual(ids1)
   *
   * \b Example:
   * - \a ids1 : [3,1,103,4,6,10,-7,205]
   * - \a ids2 : [-7,1,205,10,6,3,103,4]
   * - \a return is : [5,1,6,7,4,3,0,2] because ids2[5]==ids1[0], ids2[1]==ids1[1], ids2[6]==ids1[2]...
   *
   * \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *          array using decrRef() as it is no more needed.
   * \throw If either ids1 or ids2 is null not allocated or not with one components.
   *
   * \sa DataArrayInt::findIdForEach, DataArrayInt::FindPermutationFromFirstToSecondDuplicate, DataArrayInt::rankOfElementInThis
   */
  template<class T>
  DataArrayIdType *DataArrayDiscrete<T>::FindPermutationFromFirstToSecond(const DataArrayType *ids1, const DataArrayType *ids2)
  {
    if(!ids1 || !ids2)
      throw INTERP_KERNEL::Exception("DataArrayInt::FindPermutationFromFirstToSecond : the two input arrays must be not null !");
    if(!ids1->isAllocated() || !ids2->isAllocated())
      throw INTERP_KERNEL::Exception("DataArrayInt::FindPermutationFromFirstToSecond : the two input arrays must be allocated !");
    if(ids1->getNumberOfComponents()!=1 || ids2->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::FindPermutationFromFirstToSecond : the two input arrays have exactly one component !");
    if(ids1->getNumberOfTuples()!=ids2->getNumberOfTuples())
      {
        std::ostringstream oss; oss << "DataArrayInt::FindPermutationFromFirstToSecond : first array has " << ids1->getNumberOfTuples() << " tuples and the second one " << ids2->getNumberOfTuples() << " tuples ! No chance to find a permutation between the 2 arrays !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    MCAuto<DataArrayType> c1(ids1->deepCopy());
    MCAuto<DataArrayType> c2(ids2->deepCopy());
    c1->sort(true); c2->sort(true);
    if(!c1->isEqualWithoutConsideringStr(*c2))
      throw INTERP_KERNEL::Exception("DataArrayInt::FindPermutationFromFirstToSecond : the two arrays are not lying on same ids ! Impossible to find a permutation between the 2 arrays !");
    MCAuto<DataArrayIdType> p1=ids1->checkAndPreparePermutation();
    MCAuto<DataArrayIdType> p2=ids2->checkAndPreparePermutation();
    p2=p2->invertArrayO2N2N2O(p2->getNumberOfTuples());
    p2=p2->selectByTupleIdSafe(p1->begin(),p1->end());
    return p2.retn();
  }

  /*!
   * This method tries to find the permutation to apply to the first input \a ids1 to obtain the same array (without considering strings information) the second
   * input array \a ids2.
   * \a ids1 and \a ids2 are expected to be both a list of ids (both with number of components equal to one) not sorted and with values that can be negative.
   * This method will throw an exception is no such permutation array can be obtained. It is typically the case if there is some ids in \a ids1 not in \a ids2 or
   * inversely.
   * The difference with DataArrayInt::FindPermutationFromFirstToSecond is that this method supports multiple same values in \a ids1 and \a ids2 whereas
   * DataArrayInt::FindPermutationFromFirstToSecond doesn't. It implies that this method my be slower than the DataArrayInt::FindPermutationFromFirstToSecond one.
   *
   * In case of success both assertion will be true (no throw) :
   * \c ids1->renumber(ret)->isEqual(ids2) where \a ret is the return of this method.
   * \c ret->transformWithIndArr(ids2)->isEqual(ids1)
   *
   * \b Example:
   * - \a ids1 : [5, 3, 2, 1, 4, 5, 2, 1, 0, 11, 5, 4]
   * - \a ids2 : [0, 1, 1, 2, 2, 3, 4, 4, 5, 5, 5, 11]
   * - \a return is : [8, 5, 3, 1, 6, 9, 4, 2, 0, 11, 10, 7] because ids2[8]==ids1[0], ids2[5]==ids1[1], ids2[3]==ids1[2], ids2[1]==ids1[3]...
   *
   * \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *          array using decrRef() as it is no more needed.
   * \throw If either ids1 or ids2 is null not allocated or not with one components.
   *
   * \sa DataArrayInt::findIdForEach, DataArrayInt::FindPermutationFromFirstToSecond, DataArrayInt::occurenceRankInThis
   */
  template<class T>
  DataArrayIdType *DataArrayDiscrete<T>::FindPermutationFromFirstToSecondDuplicate(const DataArrayType *ids1, const DataArrayType *ids2)
  {
    if(!ids1 || !ids2)
      throw INTERP_KERNEL::Exception("DataArrayInt::FindPermutationFromFirstToSecondDuplicate : the two input arrays must be not null !");
    constexpr char MSG0[] = "DataArrayInt::FindPermutationFromFirstToSecondDuplicate :";
    ids1->checkAllocated(); ids2->checkAllocated();
    ids1->checkNbOfComps(1,MSG0); ids2->checkNbOfComps(1,MSG0);
    mcIdType nbTuples(ids1->getNumberOfTuples());
    if(nbTuples != ids2->getNumberOfTuples())
      {
        std::ostringstream oss; oss << "DataArrayInt::FindPermutationFromFirstToSecondDuplicate : first array has " << ids1->getNumberOfTuples() << " tuples and the second one " << ids2->getNumberOfTuples() << " tuples ! No chance to find a permutation between the 2 arrays !";
        throw INTERP_KERNEL::Exception(oss.str().c_str());
      }
    MCAuto<DataArrayIdType> ret(DataArrayIdType::New()); ret->alloc(nbTuples,1);
    MCAuto<DataArrayIdType> oids2(ids2->occurenceRankInThis());
    std::map< std::pair<T,mcIdType>, mcIdType> m;
    mcIdType pos(0);
    const mcIdType *oids2Ptr(oids2->begin());
    for(const T * it2 = ids2->begin() ; it2 != ids2->end() ; ++it2, ++oids2Ptr, ++pos)
      m[{*it2,*oids2Ptr}] = pos;
    mcIdType *retPtr(ret->getPointer());
    //
    std::map<T,mcIdType> mOccurence1; // see DataArrayInt::occurenceRankInThis : avoid to compute additionnal temporary array
    //
    for(const T * it1 = ids1->begin() ; it1 != ids1->end() ; ++it1, ++retPtr)
    {
      auto it = mOccurence1.find(*it1);
      mcIdType occRk1;
      if( it == mOccurence1.end() )
      {
        occRk1 = 0;
        mOccurence1[*it1] = 1;
      }
      else
      {
        occRk1 = (*it).second++;
      }
      //
      auto it2 = m.find({*it1,occRk1});
      if(it2 != m.end())
      {
        *retPtr = (*it2).second;
      }
      else
      {
        std::ostringstream oss; oss << MSG0 << "At pos " << std::distance(ids1->begin(),it1) << " value is " << *it1 << " and occurence rank is " << occRk1 << ". No such item into second array !";
        throw INTERP_KERNEL::Exception(oss.str());
      }

    }
    return ret.retn();
  }

  /*!
   * Returns a C array which is a renumbering map in "Old to New" mode for the input array.
   * This map, if applied to \a start array, would make it sorted. For example, if
   * \a start array contents are [9,10,0,6,4,11,3,7] then the contents of the result array is
   * [5,6,0,3,2,7,1,4].
   *  \param [in] start - pointer to the first element of the array for which the
   *         permutation map is computed.
   *  \param [in] end - pointer specifying the end of the array \a start, so that
   *         the last value of \a start is \a end[ -1 ].
   *  \return mcIdType * - the result permutation array that the caller is to delete as it is no
   *         more needed.
   *  \throw If there are equal values in the input array.
   */
  template<class T>
  mcIdType *DataArrayDiscrete<T>::CheckAndPreparePermutation(const T *start, const T *end)
  {
    std::size_t sz=std::distance(start,end);
    mcIdType *ret=(mcIdType *)malloc(sz*sizeof(mcIdType));
    T *work=new T[sz];
    std::copy(start,end,work);
    std::sort(work,work+sz);
    if(std::unique(work,work+sz)!=work+sz)
      {
        delete [] work;
        free(ret);
        throw INTERP_KERNEL::Exception("Some elements are equals in the specified array !");
      }
    std::map<T,mcIdType> m;
    for(T *workPt=work;workPt!=work+sz;workPt++)
      m[*workPt]=ToIdType(std::distance(work,workPt));
    mcIdType *iter2=ret;
    for(const T *iter=start;iter!=end;iter++,iter2++)
      *iter2=m[*iter];
    delete [] work;
    return ret;
  }

  /*!
   * Returns a new DataArrayInt by concatenating two given arrays, so that (1) the number
   * of tuples in the result array is <em> a1->getNumberOfTuples() + a2->getNumberOfTuples() -
   * offsetA2</em> and (2)
   * the number of component in the result array is same as that of each of given arrays.
   * First \a offsetA2 tuples of \a a2 are skipped and thus are missing from the result array.
   * Info on components is copied from the first of the given arrays. Number of components
   * in the given arrays must be the same.
   *  \param [in] a1 - an array to include in the result array.
   *  \param [in] a2 - another array to include in the result array.
   *  \param [in] offsetA2 - number of tuples of \a a2 to skip.
   *  \return DataArrayInt * - the new instance of DataArrayInt.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If either \a a1 or \a a2 is NULL.
   *  \throw If \a a1->getNumberOfComponents() != \a a2->getNumberOfComponents().
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::Aggregate(const DataArrayType *a1, const DataArrayType *a2, T offsetA2)
  {
    if(!a1 || !a2)
      throw INTERP_KERNEL::Exception("DataArrayInt::Aggregate : input DataArrayInt instance is NULL !");
    std::size_t nbOfComp(a1->getNumberOfComponents());
    if(nbOfComp!=a2->getNumberOfComponents())
      throw INTERP_KERNEL::Exception("Nb of components mismatch for array Aggregation !");
    mcIdType nbOfTuple1(a1->getNumberOfTuples()),nbOfTuple2(a2->getNumberOfTuples());
    MCAuto<DataArrayType> ret(DataArrayType::New());
    ret->alloc(nbOfTuple1+nbOfTuple2-offsetA2,nbOfComp);
    T *pt=std::copy(a1->begin(),a1->end(),ret->getPointer());
    std::copy(a2->getConstPointer()+offsetA2*nbOfComp,a2->getConstPointer()+nbOfTuple2*nbOfComp,pt);
    ret->copyStringInfoFrom(*a1);
    return ret.retn();
  }

  /*!
   * Returns a new DataArrayInt by concatenating all given arrays, so that (1) the number
   * of tuples in the result array is a sum of the number of tuples of given arrays and (2)
   * the number of component in the result array is same as that of each of given arrays.
   * Info on components is copied from the first of the given arrays. Number of components
   * in the given arrays must be  the same.
   * If the number of non null of elements in \a arr is equal to one the returned object is a copy of it
   * not the object itself.
   *  \param [in] arr - a sequence of arrays to include in the result array.
   *  \return DataArrayInt * - the new instance of DataArrayInt.
   *          The caller is to delete this result array using decrRef() as it is no more
   *          needed.
   *  \throw If all arrays within \a arr are NULL.
   *  \throw If getNumberOfComponents() of arrays within \a arr.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::Aggregate(const std::vector<const DataArrayType *>& arr)
  {
    std::vector<const DataArrayType *> a;
    for(typename std::vector<const DataArrayType *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
      if(*it4)
        a.push_back(*it4);
    if(a.empty())
      throw INTERP_KERNEL::Exception("DataArrayInt::Aggregate : input list must be NON EMPTY !");
    typename std::vector<const DataArrayType *>::const_iterator it=a.begin();
    std::size_t nbOfComp((*it)->getNumberOfComponents());
    mcIdType nbt((*it++)->getNumberOfTuples());
    for(;it!=a.end();it++)
      {
        if((*it)->getNumberOfComponents()!=nbOfComp)
          throw INTERP_KERNEL::Exception("DataArrayInt::Aggregate : Nb of components mismatch for array aggregation !");
        nbt+=(*it)->getNumberOfTuples();
      }
    MCAuto<DataArrayType> ret=DataArrayType::New();
    ret->alloc(nbt,nbOfComp);
    T *pt=ret->getPointer();
    for(it=a.begin();it!=a.end();it++)
      pt=std::copy((*it)->getConstPointer(),(*it)->getConstPointer()+(*it)->getNbOfElems(),pt);
    ret->copyStringInfoFrom(*(a[0]));
    return ret.retn();
  }

  /*!
   * This method takes as input a list of DataArrayInt instances \a arrs that represent each a packed index arrays.
   * A packed index array is an allocated array with one component, and at least one tuple. The first element
   * of each array in \a arrs must be 0. Each array in \a arrs is expected to be increasingly monotonic.
   * This method is useful for users that want to aggregate a pair of DataArrayInt representing an indexed data (typically nodal connectivity index in unstructured meshes.
   *
   * \return DataArrayInt * - a new object to be managed by the caller.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::AggregateIndexes(const std::vector<const DataArrayType *>& arrs)
  {
    mcIdType retSz=1;
    for(typename std::vector<const DataArrayType *>::const_iterator it4=arrs.begin();it4!=arrs.end();it4++)
      {
        if(*it4)
          {
            (*it4)->checkAllocated();
            if((*it4)->getNumberOfComponents()!=1)
              {
                std::ostringstream oss; oss << "DataArrayInt::AggregateIndexes : presence of a DataArrayInt instance with nb of compo != 1 at pos " << std::distance(arrs.begin(),it4) << " !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
            mcIdType nbTupl((*it4)->getNumberOfTuples());
            if(nbTupl<1)
              {
                std::ostringstream oss; oss << "DataArrayInt::AggregateIndexes : presence of a DataArrayInt instance with nb of tuples < 1 at pos " << std::distance(arrs.begin(),it4) << " !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
            if((*it4)->front()!=0)
              {
                std::ostringstream oss; oss << "DataArrayInt::AggregateIndexes : presence of a DataArrayInt instance with front value != 0 at pos " << std::distance(arrs.begin(),it4) << " !";
                throw INTERP_KERNEL::Exception(oss.str().c_str());
              }
            retSz+=nbTupl-1;
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::AggregateIndexes : presence of a null instance at pos " << std::distance(arrs.begin(),it4) << " !";
            throw INTERP_KERNEL::Exception(oss.str().c_str());
          }
      }
    if(arrs.empty())
      throw INTERP_KERNEL::Exception("DataArrayInt::AggregateIndexes : input list must be NON EMPTY !");
    MCAuto<DataArrayType> ret=DataArrayType::New();
    ret->alloc(retSz,1);
    T *pt=ret->getPointer(); *pt++=0;
    for(typename std::vector<const DataArrayType *>::const_iterator it=arrs.begin();it!=arrs.end();it++)
      pt=std::transform((*it)->begin()+1,(*it)->end(),pt,std::bind(std::plus<T>(),std::placeholders::_1,pt[-1]));
    ret->copyStringInfoFrom(*(arrs[0]));
    return ret.retn();
  }

  /*!
   * Returns a new DataArrayInt which contains all elements of given one-dimensional
   * arrays. The result array does not contain any duplicates and its values
   * are sorted in ascending order.
   *  \param [in] arr - sequence of DataArrayInt's to unite.
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *         array using decrRef() as it is no more needed.
   *  \throw If any \a arr[i] is not allocated.
   *  \throw If \a arr[i]->getNumberOfComponents() != 1.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::BuildUnion(const std::vector<const DataArrayType *>& arr)
  {
    std::vector<const DataArrayType *> a;
    for(typename std::vector<const DataArrayType *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
      if(*it4)
        a.push_back(*it4);
    for(typename std::vector<const DataArrayType *>::const_iterator it=a.begin();it!=a.end();it++)
      {
        (*it)->checkAllocated();
        if((*it)->getNumberOfComponents()!=1)
          throw INTERP_KERNEL::Exception("DataArrayInt::BuildUnion : only single component allowed !");
      }
    //
    std::set<T> r;
    for(typename std::vector<const DataArrayType *>::const_iterator it=a.begin();it!=a.end();it++)
      {
        const T *pt=(*it)->getConstPointer();
        mcIdType nbOfTuples((*it)->getNumberOfTuples());
        r.insert(pt,pt+nbOfTuples);
      }
    DataArrayType *ret=DataArrayType::New();
    ret->alloc(r.size(),1);
    std::copy(r.begin(),r.end(),ret->getPointer());
    return ret;
  }

  /*!
   * Returns a new DataArrayInt which contains elements present in each of given one-dimensional
   * arrays. The result array does not contain any duplicates and its values
   * are sorted in ascending order.
   *  \param [in] arr - sequence of DataArrayInt's to intersect.
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *         array using decrRef() as it is no more needed.
   *  \throw If any \a arr[i] is not allocated.
   *  \throw If \a arr[i]->getNumberOfComponents() != 1.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::BuildIntersection(const std::vector<const DataArrayType *>& arr)
  {
    std::vector<const DataArrayType *> a;
    for(typename std::vector<const DataArrayType *>::const_iterator it4=arr.begin();it4!=arr.end();it4++)
      if(*it4)
        a.push_back(*it4);
    for(typename std::vector<const DataArrayType *>::const_iterator it=a.begin();it!=a.end();it++)
      {
        (*it)->checkAllocated();
        if((*it)->getNumberOfComponents()!=1)
          throw INTERP_KERNEL::Exception("DataArrayInt::BuildIntersection : only single component allowed !");
      }
    //
    if(a.size()==1)
      throw INTERP_KERNEL::Exception("DataArrayInt::BuildIntersection : only single not null element in array ! For safety reasons exception is raised !");
    //
    std::set<T> r;
    for(typename std::vector<const DataArrayType *>::const_iterator it=a.begin();it!=a.end();it++)
      {
        const T *pt=(*it)->getConstPointer();
        mcIdType nbOfTuples((*it)->getNumberOfTuples());
        std::set<T> s1(pt,pt+nbOfTuples);
        if(it!=a.begin())
          {
            std::set<T> r2;
            std::set_intersection(r.begin(),r.end(),s1.begin(),s1.end(),inserter(r2,r2.end()));
            r=r2;
          }
        else
          r=s1;
      }
    DataArrayType *ret(DataArrayType::New());
    ret->alloc(r.size(),1);
    std::copy(r.begin(),r.end(),ret->getPointer());
    return ret;
  }

  /*!
   * This method allows to put a vector of vector of integer into a more compact data structure (skyline).
   * This method is not available into python because no available optimized data structure available to map std::vector< std::vector<mcIdType> >.
   *
   * \param [in] v the input data structure to be translate into skyline format.
   * \param [out] data the first element of the skyline format. The user is expected to deal with newly allocated array.
   * \param [out] dataIndex the second element of the skyline format.
   */
  template <class T>
  void DataArrayDiscrete<T>::PutIntoToSkylineFrmt(const std::vector< std::vector<T> >& v, DataArrayType *& data, DataArrayIdType *& dataIndex)
  {
    std::size_t sz(v.size());
    MCAuto<DataArrayType> retDat(DataArrayType::New());
    MCAuto<DataArrayIdType> retIdx(DataArrayIdType::New());
    retIdx->alloc(sz+1,1);
    mcIdType *ptid(retIdx->getPointer()); *ptid=0;
    for(std::size_t i=0;i<sz;i++,ptid++)
      ptid[1]=ptid[0]+ToIdType(v[i].size());
    retDat->alloc(retIdx->back(),1);
    T *pt=retDat->getPointer();
    for(std::size_t i=0;i<sz;i++)
      pt=std::copy(v[i].begin(),v[i].end(),pt);
    data=retDat.retn(); dataIndex=retIdx.retn();
  }

  /*!
   * This method works on a pair input (\b arrIn, \b arrIndxIn) where \b arrIn indexes is in \b arrIndxIn
   * (\ref numbering-indirect).
   * This method returns the result of the extraction ( specified by a set of ids in [\b idsOfSelectBg , \b idsOfSelectEnd ) ).
   * The selection of extraction is done standardly in new2old format.
   * This method returns indexed arrays (\ref numbering-indirect) using 2 arrays (arrOut,arrIndexOut).
   *
   * \param [in] idsOfSelectBg begin of set of ids of the input extraction (included)
   * \param [in] idsOfSelectEnd end of set of ids of the input extraction (excluded)
   * \param [in] arrIn arr origin array from which the extraction will be done.
   * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
   * \param [out] arrOut the resulting array
   * \param [out] arrIndexOut the index array of the resulting array \b arrOut
   * \sa DataArrayInt::ExtractFromIndexedArraysSlice
   */
  template <class T>
  void DataArrayDiscrete<T>::ExtractFromIndexedArrays(const mcIdType *idsOfSelectBg, const mcIdType *idsOfSelectEnd,
                                                      const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                      DataArrayType* &arrOut, DataArrayIdType* &arrIndexOut)
  {
    if(!arrIn || !arrIndxIn)
      throw INTERP_KERNEL::Exception("DataArrayInt::ExtractFromIndexedArrays : input pointer is NULL !");
    arrIn->checkAllocated(); arrIndxIn->checkAllocated();
    if(arrIn->getNumberOfComponents()!=1 || arrIndxIn->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::ExtractFromIndexedArrays : input arrays must have exactly one component !");
    std::size_t sz=std::distance(idsOfSelectBg,idsOfSelectEnd);
    const T *arrInPtr=arrIn->begin();
    const mcIdType *arrIndxPtr=arrIndxIn->begin();
    mcIdType nbOfGrps=arrIndxIn->getNumberOfTuples()-1;
    if(nbOfGrps<0)
      throw INTERP_KERNEL::Exception("DataArrayInt::ExtractFromIndexedArrays : The format of \"arrIndxIn\" is invalid ! Its nb of tuples should be >=1 !");
    mcIdType maxSizeOfArr(arrIn->getNumberOfTuples());
    MCAuto<DataArrayType> arro=DataArrayType::New();
    MCAuto<DataArrayIdType> arrIo=DataArrayIdType::New();
    arrIo->alloc(sz+1,1);
    const mcIdType *idsIt=idsOfSelectBg;
    mcIdType *work=arrIo->getPointer();
    *work++=0;
    mcIdType lgth=0;
    for(std::size_t i=0;i<sz;i++,work++,idsIt++)
      {
        if(*idsIt>=0 && *idsIt<nbOfGrps)
          lgth+=arrIndxPtr[*idsIt+1]-arrIndxPtr[*idsIt];
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::ExtractFromIndexedArrays : id located on pos #" << i << " value is " << *idsIt << " ! Must be in [0," << nbOfGrps << ") !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
        if(lgth>=work[-1])
          *work=lgth;
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::ExtractFromIndexedArrays : id located on pos #" << i << " value is " << *idsIt << " and at this pos arrIndxIn[" << *idsIt;
            oss << "+1]-arrIndxIn[" << *idsIt << "] < 0 ! The input index array is bugged !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
    arro->alloc(lgth,1);
    T *data=arro->getPointer();
    idsIt=idsOfSelectBg;
    for(std::size_t i=0;i<sz;i++,idsIt++)
      {
        if(arrIndxPtr[*idsIt]>=0 && arrIndxPtr[*idsIt+1]<=maxSizeOfArr)
          data=std::copy(arrInPtr+arrIndxPtr[*idsIt],arrInPtr+arrIndxPtr[*idsIt+1],data);
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::ExtractFromIndexedArrays : id located on pos #" << i << " value is " << *idsIt << " arrIndx[" << *idsIt << "] must be >= 0 and arrIndx[";
            oss << *idsIt << "+1] <= " << maxSizeOfArr << " (the size of arrIn)!";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
    arrOut=arro.retn();
    arrIndexOut=arrIo.retn();
  }

  /*!
   * This method converts from VTK <= 9.3 polyhedra nodal connectivity to MED.
   *
   * \param [in] arrIn arr origin array from which the extraction will be done.
   * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
   * \param [out] arrOut the resulting array
   * \param [out] arrIndexOut the index array of the resulting array \b arrOut
   * \sa DataArrayDiscrete<T>::FromVTK93To94FacesInternaReprOfPolyedra
   */
  template <class T>
  void DataArrayDiscrete<T>::FromVTKInternalReprOfPolyedra(const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                           MCAuto<DataArrayType> &arrOut, MCAuto<DataArrayIdType> &arrIndexOut)
  {
    if(!arrIn || !arrIndxIn)
      throw INTERP_KERNEL::Exception("DataArrayInt::FromVTKInternalReprOfPolyedra : input pointer is NULL !");
    arrIn->checkAllocated(); arrIndxIn->checkAllocated();
    arrIn->checkNbOfComps(1,"1st array must have single component");
    arrIndxIn->checkNbOfComps(1,"2nd array must have single component");
    if(arrIndxIn->getNumberOfTuples()<1)
      THROW_IK_EXCEPTION("2nd input array must be of size >= 1");
    mcIdType nbCells(arrIndxIn->getNumberOfTuples()-1);
    const T *arrInPt(arrIn->begin());
    const mcIdType *arrIndxInPt(arrIndxIn->begin());
    arrIndexOut = DataArrayIdType::New(); arrIndexOut->alloc(arrIndxIn->getNumberOfTuples(),1);
    arrOut = DataArrayType::New(); arrOut->alloc(arrIn->getNumberOfTuples() - 2*nbCells,1);
    T *arrOutPt(arrOut->getPointer());
    mcIdType *arrOutIdxPt(arrIndexOut->getPointer()); *arrOutIdxPt = 0;
    for(auto i = 0 ; i < nbCells ; ++i)
    {
      T nbFaces = arrInPt[ arrIndxInPt[i] ];
      T *arrOutPtStart(arrOutPt);
      const T *facePtr = arrInPt + arrIndxInPt[i] + 1;
      for(T iFace = 0 ; iFace < nbFaces ; ++iFace)
      {
        T nbNodesInFace = *facePtr++;
        if(iFace>0)
        {
          *arrOutPt++ = -1;
        }
        arrOutPt = std::copy(facePtr,facePtr+nbNodesInFace,arrOutPt);
        facePtr += nbNodesInFace;
      }
      arrOutIdxPt[1] = arrOutIdxPt[0] + std::distance(arrOutPtStart,arrOutPt);
      ++arrOutIdxPt;
    }
  }

  /*!
   * This method converts from VTK >= 9.4 polyhedra nodal connectivity to MED.
   *
   * \param [in] arrIn arr origin array from which the extraction will be done.
   * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
   * \param [in] arrIndxIn2 is the input index array telling which is the faces id range in arrIndxIn of poly
   * \param [out] arrOut the resulting array
   * \param [out] arrIndexOut the index array of the resulting array \b arrOut
   * \sa DataArrayDiscrete<T>::FromVTK93To94FacesInternaReprOfPolyedra
   */
  template <class T>
  void DataArrayDiscrete<T>::FromVTK94InternalReprOfPolyedra(const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn, const DataArrayIdType *arrIndxIn2,
                                                             MCAuto<DataArrayType> &arrOut, MCAuto<DataArrayIdType> &arrIndexOut)
  {
    if(!arrIn || !arrIndxIn)
      throw INTERP_KERNEL::Exception("DataArrayInt::FromVTK94InternalReprOfPolyedra : input pointer is NULL !");
    arrIn->checkAllocated(); arrIndxIn->checkAllocated();
    arrIn->checkNbOfComps(1,"1st array must have single component");
    arrIndxIn->checkNbOfComps(1,"2nd array must have single component");
    arrIndxIn2->checkNbOfComps(1,"3nd array must have single component");
    if(arrIndxIn->getNumberOfTuples()<1 || arrIndxIn2->getNumberOfTuples()<1)
      THROW_IK_EXCEPTION("2nd and 3rd input array must be of size >= 1");
    mcIdType nbCells(arrIndxIn2->getNumberOfTuples()-1);
    const T *arrInPt(arrIn->begin());
    const mcIdType *arrIndxInPt(arrIndxIn->begin()), *arrIndx2InPt(arrIndxIn2->begin());
    arrIndexOut = DataArrayIdType::New(); arrIndexOut->alloc(arrIndxIn2->getNumberOfTuples(),1);
    arrOut = DataArrayType::New(); arrOut->alloc(arrIn->getNumberOfTuples() + ( arrIndxIn->getNumberOfTuples() - 1 ) - nbCells,1);
    T *arrOutPt(arrOut->getPointer());
    mcIdType *arrOutIdxPt(arrIndexOut->getPointer()); *arrOutIdxPt = 0;
    for(auto i = 0 ; i < nbCells ; ++i)
    {
      mcIdType nbFaces = arrIndx2InPt[i+1] - arrIndx2InPt[i];
      T *arrOutPtStart(arrOutPt);
      for(mcIdType iFace = 0 ; iFace < nbFaces ; ++iFace)
      {
        const T *facePtr = arrInPt + arrIndxInPt[ arrIndx2InPt[i] +  iFace ];
        mcIdType nbNodesInFace = arrIndxInPt[ arrIndx2InPt[i] +  iFace + 1 ] - arrIndxInPt[ arrIndx2InPt[i] +  iFace ];
        if(iFace>0)
        {
          *arrOutPt++ = -1;
        }
        arrOutPt = std::copy(facePtr,facePtr+nbNodesInFace,arrOutPt);
        facePtr += nbNodesInFace;
      }
      arrOutIdxPt[1] = arrOutIdxPt[0] + std::distance(arrOutPtStart,arrOutPt);
      ++arrOutIdxPt;
    }
  }

  /*!
   * This method converts from VTK 9.3 polyhedra nodal connectivity to faces connectivity compatible with MEDCoupling1DGTUMesh(NORM_POLYGON).
   *
   * \param [in] arrIn arr origin array from which the extraction will be done.
   * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
   * \param [out] arrOut the resulting array
   * \param [out] arrIndexOut the index array of the resulting array \b arrOut
   * \sa DataArrayDiscrete<T>::FromVTKInternalReprOfPolyedra
   */
  template <class T>
  void DataArrayDiscrete<T>::FromVTK93To94FacesInternaReprOfPolyedra(const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                                     MCAuto<DataArrayType> &arrOut, MCAuto<DataArrayIdType> &arrIndexOut)
  {
    if(!arrIn || !arrIndxIn)
      throw INTERP_KERNEL::Exception("DataArrayInt::FromVTK93FacesInternaReprOfPolyedra : input pointer is NULL !");
    arrIn->checkAllocated(); arrIndxIn->checkAllocated();
    arrIn->checkNbOfComps(1,"1st array must have single component");
    arrIndxIn->checkNbOfComps(1,"2nd array must have single component");
    if(arrIndxIn->getNumberOfTuples()<1)
      THROW_IK_EXCEPTION("2nd input array must be of size >= 1");
    mcIdType nbCells(arrIndxIn->getNumberOfTuples()-1);
    const T *arrInPt(arrIn->begin());
    const mcIdType *arrIndxInPt(arrIndxIn->begin());
    arrOut = DataArrayType::New(); arrOut->alloc(0,1);
    mcIdType curIndex(0);
    arrIndexOut = DataArrayIdType::New(); arrIndexOut->alloc(1,1); arrIndexOut->setIJSilent(0,0,curIndex);
    for( mcIdType i = 0 ; i <  nbCells ; ++i )
    {
      const T *ptOfCurPolyh = arrInPt + arrIndxInPt[i];
      const T *endOfCurPolyh = arrInPt + arrIndxInPt[ i+1 ];
      T nbOfFaces = *ptOfCurPolyh++;
      //
      for( T j = 0 ; j < nbOfFaces && ptOfCurPolyh <= endOfCurPolyh ; ++j )
      {
        T nbOfPtsInFace = *ptOfCurPolyh++;
        mcIdType newIndex = curIndex + ToIdType( nbOfPtsInFace );
        arrIndexOut->pushBackSilent( newIndex );
        curIndex = newIndex;
        const T *endOfFaceConn = ptOfCurPolyh + nbOfPtsInFace;
        arrOut->pushBackValsSilent( ptOfCurPolyh, endOfFaceConn);
        ptOfCurPolyh = endOfFaceConn;
      }
      if( ptOfCurPolyh != endOfCurPolyh )
        THROW_IK_EXCEPTION("For cell #" << i << " problem in input connectivity");
    }
    arrIndexOut->pack();
    arrOut->pack();
  }

  /*!
   * This method works on a pair input (\b arrIn, \b arrIndxIn) where \b arrIn indexes is in \b arrIndxIn
   * (\ref numbering-indirect).
   * This method returns the result of the extraction ( specified by a set of ids with a slice given by \a idsOfSelectStart, \a idsOfSelectStop and \a idsOfSelectStep ).
   * The selection of extraction is done standardly in new2old format.
   * This method returns indexed arrays (\ref numbering-indirect) using 2 arrays (arrOut,arrIndexOut).
   *
   * \param [in] idsOfSelectStart begin of set of ids of the input extraction (included)
   * \param [in] idsOfSelectStop end of set of ids of the input extraction (excluded)
   * \param [in] idsOfSelectStep step of set of ids of the input extraction
   * \param [in] arrIn arr origin array from which the extraction will be done.
   * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
   * \param [out] arrOut the resulting array
   * \param [out] arrIndexOut the index array of the resulting array \b arrOut
   * \sa DataArrayInt::ExtractFromIndexedArrays
   */
  template <class T>
  void DataArrayDiscrete<T>::ExtractFromIndexedArraysSlice(mcIdType idsOfSelectStart, mcIdType idsOfSelectStop, mcIdType idsOfSelectStep,
                                                           const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                           DataArrayType* &arrOut, DataArrayIdType* &arrIndexOut)
  {
    if(!arrIn || !arrIndxIn)
      throw INTERP_KERNEL::Exception("DataArrayInt::ExtractFromIndexedArraysSlice : input pointer is NULL !");
    arrIn->checkAllocated(); arrIndxIn->checkAllocated();
    if(arrIn->getNumberOfComponents()!=1 || arrIndxIn->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::ExtractFromIndexedArraysSlice : input arrays must have exactly one component !");
    mcIdType sz=DataArray::GetNumberOfItemGivenBESRelative(idsOfSelectStart,idsOfSelectStop,idsOfSelectStep,"MEDCouplingUMesh::ExtractFromIndexedArraysSlice : Input slice ");
    const T *arrInPtr=arrIn->begin();
    const mcIdType *arrIndxPtr=arrIndxIn->begin();
    mcIdType nbOfGrps=arrIndxIn->getNumberOfTuples()-1;
    if(nbOfGrps<0)
      throw INTERP_KERNEL::Exception("DataArrayInt::ExtractFromIndexedArraysSlice : The format of \"arrIndxIn\" is invalid ! Its nb of tuples should be >=1 !");
    mcIdType maxSizeOfArr(arrIn->getNumberOfTuples());
    MCAuto<DataArrayType> arro=DataArrayType::New();
    MCAuto<DataArrayIdType> arrIo=DataArrayIdType::New();
    arrIo->alloc(sz+1,1);
    mcIdType idsIt=idsOfSelectStart;
    mcIdType *work=arrIo->getPointer();
    *work++=0;
    mcIdType lgth=0;
    for(mcIdType i=0;i<sz;i++,work++,idsIt+=idsOfSelectStep)
      {
        if(idsIt>=0 && idsIt<nbOfGrps)
          lgth+=arrIndxPtr[idsIt+1]-arrIndxPtr[idsIt];
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::ExtractFromIndexedArraysSlice : id located on pos #" << i << " value is " << idsIt << " ! Must be in [0," << nbOfGrps << ") !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
        if(lgth>=work[-1])
          *work=lgth;
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::ExtractFromIndexedArraysSlice : id located on pos #" << i << " value is " << idsIt << " and at this pos arrIndxIn[" << idsIt;
            oss << "+1]-arrIndxIn[" << idsIt << "] < 0 ! The input index array is bugged !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
    arro->alloc(lgth,1);
    T *data=arro->getPointer();
    idsIt=idsOfSelectStart;
    for(mcIdType i=0;i<sz;i++,idsIt+=idsOfSelectStep)
      {
        if(arrIndxPtr[idsIt]>=0 && arrIndxPtr[idsIt+1]<=maxSizeOfArr)
          data=std::copy(arrInPtr+arrIndxPtr[idsIt],arrInPtr+arrIndxPtr[idsIt+1],data);
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::ExtractFromIndexedArraysSlice : id located on pos #" << i << " value is " << idsIt << " arrIndx[" << idsIt << "] must be >= 0 and arrIndx[";
            oss << idsIt << "+1] <= " << maxSizeOfArr << " (the size of arrIn)!";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
    arrOut=arro.retn();
    arrIndexOut=arrIo.retn();
  }

  /*!
   * This method works on an input pair (\b arrIn, \b arrIndxIn) where \b arrIn indexes is in \b arrIndxIn.
   * This method builds an output pair (\b arrOut,\b arrIndexOut) that is a copy from \b arrIn for all cell ids \b not \b in [ \b idsOfSelectBg , \b idsOfSelectEnd ) and for
   * cellIds \b in [ \b idsOfSelectBg , \b idsOfSelectEnd ) a copy coming from the corresponding values in input pair (\b srcArr, \b srcArrIndex).
   * This method is an generalization of MEDCouplingUMesh::SetPartOfIndexedArraysSameIdx that performs the same thing but by without building explicitly a result output arrays.
   *
   * \param [in] idsOfSelectBg begin of set of ids of the input extraction (included)
   * \param [in] idsOfSelectEnd end of set of ids of the input extraction (excluded)
   * \param [in] arrIn arr origin array from which the extraction will be done.
   * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
   * \param [in] srcArr input array that will be used as source of copy for ids in [ \b idsOfSelectBg, \b idsOfSelectEnd )
   * \param [in] srcArrIndex index array of \b srcArr
   * \param [out] arrOut the resulting array
   * \param [out] arrIndexOut the index array of the resulting array \b arrOut
   *
   * \sa DataArrayInt::SetPartOfIndexedArraysSameIdx
   */
  template <class T>
  void DataArrayDiscrete<T>::SetPartOfIndexedArrays(const mcIdType *idsOfSelectBg, const mcIdType *idsOfSelectEnd,
                                                    const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                    const DataArrayType *srcArr, const DataArrayIdType *srcArrIndex,
                                                    DataArrayType* &arrOut, DataArrayIdType* &arrIndexOut)
  {
    if(arrIn==0 || arrIndxIn==0 || srcArr==0 || srcArrIndex==0)
      throw INTERP_KERNEL::Exception("DataArrayInt::SetPartOfIndexedArrays : presence of null pointer in input parameter !");
    MCAuto<DataArrayType> arro=DataArrayType::New();
    MCAuto<DataArrayIdType> arrIo=DataArrayIdType::New();
    mcIdType nbOfTuples=arrIndxIn->getNumberOfTuples()-1;
    std::vector<bool> v(nbOfTuples,true);
    mcIdType offset=0;
    const mcIdType *arrIndxInPtr=arrIndxIn->begin();
    const mcIdType *srcArrIndexPtr=srcArrIndex->begin();
    for(const mcIdType *it=idsOfSelectBg;it!=idsOfSelectEnd;it++,srcArrIndexPtr++)
      {
        if(*it>=0 && *it<nbOfTuples)
          {
            v[*it]=false;
            offset+=(srcArrIndexPtr[1]-srcArrIndexPtr[0])-(arrIndxInPtr[*it+1]-arrIndxInPtr[*it]);
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::SetPartOfIndexedArrays : On pos #" << std::distance(idsOfSelectBg,it) << " value is " << *it << " not in [0," << nbOfTuples << ") !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
    srcArrIndexPtr=srcArrIndex->begin();
    arrIo->alloc(nbOfTuples+1,1);
    arro->alloc(arrIn->getNumberOfTuples()+offset,1);
    const T *arrInPtr=arrIn->begin();
    const T *srcArrPtr=srcArr->begin();
    mcIdType *arrIoPtr=arrIo->getPointer(); *arrIoPtr++=0;
    T *arroPtr=arro->getPointer();
    for(mcIdType ii=0;ii<nbOfTuples;ii++,arrIoPtr++)
      {
        if(v[ii])
          {
            arroPtr=std::copy(arrInPtr+arrIndxInPtr[ii],arrInPtr+arrIndxInPtr[ii+1],arroPtr);
            *arrIoPtr=arrIoPtr[-1]+(arrIndxInPtr[ii+1]-arrIndxInPtr[ii]);
          }
        else
          {
            std::size_t pos=std::distance(idsOfSelectBg,std::find(idsOfSelectBg,idsOfSelectEnd,ii));
            arroPtr=std::copy(srcArrPtr+srcArrIndexPtr[pos],srcArrPtr+srcArrIndexPtr[pos+1],arroPtr);
            *arrIoPtr=arrIoPtr[-1]+(srcArrIndexPtr[pos+1]-srcArrIndexPtr[pos]);
          }
      }
    arrOut=arro.retn();
    arrIndexOut=arrIo.retn();
  }

  /*!
   * This method works on an input pair (\b arrIn, \b arrIndxIn) where \b arrIn indexes is in \b arrIndxIn.
   * This method builds an output pair (\b arrOut,\b arrIndexOut) that is a copy from \b arrIn for all cell ids \b not \b in [ \b idsOfSelectBg , \b idsOfSelectEnd ) and for
   * cellIds \b in [\b idsOfSelectBg, \b idsOfSelectEnd) a copy coming from the corresponding values in input pair (\b srcArr, \b srcArrIndex).
   * This method is an generalization of DataArrayInt::SetPartOfIndexedArraysSameIdx that performs the same thing but by without building explicitly a result output arrays.
   *
   * \param [in] start begin of set of ids of the input extraction (included)
   * \param [in] end end of set of ids of the input extraction (excluded)
   * \param [in] step step of the set of ids in range mode.
   * \param [in] arrIn arr origin array from which the extraction will be done.
   * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
   * \param [in] srcArr input array that will be used as source of copy for ids in [\b idsOfSelectBg, \b idsOfSelectEnd)
   * \param [in] srcArrIndex index array of \b srcArr
   * \param [out] arrOut the resulting array
   * \param [out] arrIndexOut the index array of the resulting array \b arrOut
   *
   * \sa DataArrayInt::SetPartOfIndexedArraysSameIdx DataArrayInt::SetPartOfIndexedArrays
   */
  template <class T>
  void DataArrayDiscrete<T>::SetPartOfIndexedArraysSlice(mcIdType start, mcIdType end, mcIdType step,
                                                         const DataArrayType *arrIn, const DataArrayIdType *arrIndxIn,
                                                         const DataArrayType *srcArr, const DataArrayIdType *srcArrIndex,
                                                         DataArrayType* &arrOut, DataArrayIdType* &arrIndexOut)
  {
    if(arrIn==0 || arrIndxIn==0 || srcArr==0 || srcArrIndex==0)
      throw INTERP_KERNEL::Exception("DataArrayInt::SetPartOfIndexedArraysSlice : presence of null pointer in input parameter !");
    MCAuto<DataArrayType> arro=DataArrayType::New();
    MCAuto<DataArrayIdType> arrIo=DataArrayIdType::New();
    mcIdType nbOfTuples=arrIndxIn->getNumberOfTuples()-1;
    mcIdType offset=0;
    const mcIdType *arrIndxInPtr=arrIndxIn->begin();
    const mcIdType *srcArrIndexPtr=srcArrIndex->begin();
    mcIdType nbOfElemsToSet=DataArray::GetNumberOfItemGivenBESRelative(start,end,step,"DataArrayInt::SetPartOfIndexedArraysSlice : ");
    mcIdType it=start;
    for(mcIdType i=0;i<nbOfElemsToSet;i++,srcArrIndexPtr++,it+=step)
      {
        if(it>=0 && it<nbOfTuples)
          offset+=(srcArrIndexPtr[1]-srcArrIndexPtr[0])-(arrIndxInPtr[it+1]-arrIndxInPtr[it]);
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::SetPartOfIndexedArraysSlice : On pos #" << i << " value is " << it << " not in [0," << nbOfTuples << ") !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
    srcArrIndexPtr=srcArrIndex->begin();
    arrIo->alloc(nbOfTuples+1,1);
    arro->alloc(arrIn->getNumberOfTuples()+offset,1);
    const T *arrInPtr=arrIn->begin();
    const T *srcArrPtr=srcArr->begin();
    mcIdType *arrIoPtr=arrIo->getPointer(); *arrIoPtr++=0;
    T *arroPtr=arro->getPointer();
    for(mcIdType ii=0;ii<nbOfTuples;ii++,arrIoPtr++)
      {
        mcIdType pos=DataArray::GetPosOfItemGivenBESRelativeNoThrow(ii,start,end,step);
        if(pos<0)
          {
            arroPtr=std::copy(arrInPtr+arrIndxInPtr[ii],arrInPtr+arrIndxInPtr[ii+1],arroPtr);
            *arrIoPtr=arrIoPtr[-1]+(arrIndxInPtr[ii+1]-arrIndxInPtr[ii]);
          }
        else
          {
            arroPtr=std::copy(srcArrPtr+srcArrIndexPtr[pos],srcArrPtr+srcArrIndexPtr[pos+1],arroPtr);
            *arrIoPtr=arrIoPtr[-1]+(srcArrIndexPtr[pos+1]-srcArrIndexPtr[pos]);
          }
      }
    arrOut=arro.retn();
    arrIndexOut=arrIo.retn();
  }

  /*!
   * This method works on an input pair (\b arrIn, \b arrIndxIn) where \b arrIn indexes is in \b arrIndxIn.
   * This method is an specialization of MEDCouplingUMesh::SetPartOfIndexedArrays in the case of assignment do not modify the index in \b arrIndxIn.
   *
   * \param [in] idsOfSelectBg begin of set of ids of the input extraction (included)
   * \param [in] idsOfSelectEnd end of set of ids of the input extraction (excluded)
   * \param [in,out] arrInOut arr origin array from which the extraction will be done.
   * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
   * \param [in] srcArr input array that will be used as source of copy for ids in [ \b idsOfSelectBg , \b idsOfSelectEnd )
   * \param [in] srcArrIndex index array of \b srcArr
   *
   * \sa DataArrayInt::SetPartOfIndexedArrays
   */
  template <class T>
  void DataArrayDiscrete<T>::SetPartOfIndexedArraysSameIdx(const mcIdType *idsOfSelectBg, const mcIdType *idsOfSelectEnd,
                                                           DataArrayType *arrInOut, const DataArrayIdType *arrIndxIn,
                                                           const DataArrayType *srcArr, const DataArrayIdType *srcArrIndex)
  {
    if(arrInOut==0 || arrIndxIn==0 || srcArr==0 || srcArrIndex==0)
      throw INTERP_KERNEL::Exception("DataArrayInt::SetPartOfIndexedArraysSameIdx : presence of null pointer in input parameter !");
    mcIdType nbOfTuples=arrIndxIn->getNumberOfTuples()-1;
    const mcIdType *arrIndxInPtr=arrIndxIn->begin();
    const mcIdType *srcArrIndexPtr=srcArrIndex->begin();
    T *arrInOutPtr=arrInOut->getPointer();
    const T *srcArrPtr=srcArr->begin();
    for(const mcIdType *it=idsOfSelectBg;it!=idsOfSelectEnd;it++,srcArrIndexPtr++)
      {
        if(*it>=0 && *it<nbOfTuples)
          {
            if(srcArrIndexPtr[1]-srcArrIndexPtr[0]==arrIndxInPtr[*it+1]-arrIndxInPtr[*it])
              std::copy(srcArrPtr+srcArrIndexPtr[0],srcArrPtr+srcArrIndexPtr[1],arrInOutPtr+arrIndxInPtr[*it]);
            else
              {
                std::ostringstream oss; oss << "DataArrayInt::SetPartOfIndexedArraysSameIdx : On pos #" << std::distance(idsOfSelectBg,it) << " id (idsOfSelectBg[" << std::distance(idsOfSelectBg,it)<< "]) is " << *it << " arrIndxIn[id+1]-arrIndxIn[id]!=srcArrIndex[pos+1]-srcArrIndex[pos] !";
                throw INTERP_KERNEL::Exception(oss.str());
              }
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::SetPartOfIndexedArraysSameIdx : On pos #" << std::distance(idsOfSelectBg,it) << " value is " << *it << " not in [0," << nbOfTuples << ") !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
  }

  /*!
   * This method works on an input pair (\b arrIn, \b arrIndxIn) where \b arrIn indexes is in \b arrIndxIn.
   * This method is an specialization of MEDCouplingUMesh::SetPartOfIndexedArrays in the case of assignment do not modify the index in \b arrIndxIn.
   *
   * \param [in] start begin of set of ids of the input extraction (included)
   * \param [in] end end of set of ids of the input extraction (excluded)
   * \param [in] step step of the set of ids in range mode.
   * \param [in,out] arrInOut arr origin array from which the extraction will be done.
   * \param [in] arrIndxIn is the input index array allowing to walk into \b arrIn
   * \param [in] srcArr input array that will be used as source of copy for ids in [\b idsOfSelectBg, \b idsOfSelectEnd)
   * \param [in] srcArrIndex index array of \b srcArr
   *
   * \sa DataArrayInt::SetPartOfIndexedArraysSlice DataArrayInt::SetPartOfIndexedArraysSameIdx
   */
  template <class T>
  void DataArrayDiscrete<T>::SetPartOfIndexedArraysSameIdxSlice(mcIdType start, mcIdType end, mcIdType step,
                                                                DataArrayType *arrInOut, const DataArrayIdType *arrIndxIn,
                                                                const DataArrayType *srcArr, const DataArrayIdType *srcArrIndex)
  {
    if(arrInOut==0 || arrIndxIn==0 || srcArr==0 || srcArrIndex==0)
      throw INTERP_KERNEL::Exception("DataArrayInt::SetPartOfIndexedArraysSameIdxSlice : presence of null pointer in input parameter !");
    mcIdType nbOfTuples=arrIndxIn->getNumberOfTuples()-1;
    const mcIdType *arrIndxInPtr=arrIndxIn->begin();
    const mcIdType *srcArrIndexPtr=srcArrIndex->begin();
    T *arrInOutPtr=arrInOut->getPointer();
    const T *srcArrPtr=srcArr->begin();
    mcIdType nbOfElemsToSet=DataArray::GetNumberOfItemGivenBESRelative(start,end,step,"DataArrayInt::SetPartOfIndexedArraysSameIdxSlice : ");
    mcIdType it=start;
    for(mcIdType i=0;i<nbOfElemsToSet;i++,srcArrIndexPtr++,it+=step)
      {
        if(it>=0 && it<nbOfTuples)
          {
            if(srcArrIndexPtr[1]-srcArrIndexPtr[0]==arrIndxInPtr[it+1]-arrIndxInPtr[it])
              std::copy(srcArrPtr+srcArrIndexPtr[0],srcArrPtr+srcArrIndexPtr[1],arrInOutPtr+arrIndxInPtr[it]);
            else
              {
                std::ostringstream oss; oss << "DataArrayInt::SetPartOfIndexedArraysSameIdxSlice : On pos #" << i << " id (idsOfSelectBg[" << i << "]) is " << it << " arrIndxIn[id+1]-arrIndxIn[id]!=srcArrIndex[pos+1]-srcArrIndex[pos] !";
                throw INTERP_KERNEL::Exception(oss.str());
              }
          }
        else
          {
            std::ostringstream oss; oss << "DataArrayInt::SetPartOfIndexedArraysSameIdxSlice : On pos #" << i << " value is " << it << " not in [0," << nbOfTuples << ") !";
            throw INTERP_KERNEL::Exception(oss.str());
          }
      }
  }

  /*!
   * This method works on an input pair (\b arr, \b arrIndx) where \b arr indexes is in \b arrIndx.
   * This method will not impact the size of inout parameter \b arrIndx but the size of \b arr will be modified in case of suppression.
   *
   * \param [in] idsToRemoveBg begin of set of ids to remove in \b arr (included)
   * \param [in] idsToRemoveEnd end of set of ids to remove in \b arr (excluded)
   * \param [in,out] arr array in which the remove operation will be done.
   * \param [in,out] arrIndx array in the remove operation will modify
   * \param [in] offsetForRemoval (by default 0) offset so that for each i in [0,arrIndx->getNumberOfTuples()-1) removal process will be performed in the following range [arr+arrIndx[i]+offsetForRemoval,arr+arr[i+1])
   * \return true if \b arr and \b arrIndx have been modified, false if not.
   */
  template <class T>
  bool DataArrayDiscrete<T>::RemoveIdsFromIndexedArrays(const T *idsToRemoveBg, const T *idsToRemoveEnd,
                                                        DataArrayType *arr, DataArrayIdType *arrIndx, mcIdType offsetForRemoval)
  {
    if(!arrIndx || !arr)
      throw INTERP_KERNEL::Exception("DataArrayInt::RemoveIdsFromIndexedArrays : some input arrays are empty !");
    if(offsetForRemoval<0)
      throw INTERP_KERNEL::Exception("DataArrayInt::RemoveIdsFromIndexedArrays : offsetForRemoval should be >=0 !");
    std::set<T> s(idsToRemoveBg,idsToRemoveEnd);
    mcIdType nbOfGrps=arrIndx->getNumberOfTuples()-1;
    mcIdType *arrIPtr=arrIndx->getPointer();
    *arrIPtr++=0;
    mcIdType previousArrI=0;
    const T *arrPtr=arr->begin();
    std::vector<T> arrOut;//no utility to switch to DataArrayInt because copy always needed
    for(mcIdType i=0;i<nbOfGrps;i++,arrIPtr++)
      {
        if(*arrIPtr-previousArrI>offsetForRemoval)
          {
            for(const T *work=arrPtr+previousArrI+offsetForRemoval;work!=arrPtr+*arrIPtr;work++)
              {
                if(s.find(*work)==s.end())
                  arrOut.push_back(*work);
              }
          }
        previousArrI=*arrIPtr;
        *arrIPtr=ToIdType(arrOut.size());
      }
    if(arr->getNumberOfTuples()==ToIdType(arrOut.size()))
      return false;
    arr->alloc(arrOut.size(),1);
    std::copy(arrOut.begin(),arrOut.end(),arr->getPointer());
    return true;
  }

  /*!
   * Returns a new DataArrayInt containing an arithmetic progression
   * that is equal to the sequence returned by Python \c range(\a begin,\a  end,\a  step )
   * function.
   *  \param [in] begin - the start value of the result sequence.
   *  \param [in] end - limiting value, so that every value of the result array is less than
   *              \a end.
   *  \param [in] step - specifies the increment or decrement.
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *          array using decrRef() as it is no more needed.
   *  \throw If \a step == 0.
   *  \throw If \a end < \a begin && \a step > 0.
   *  \throw If \a end > \a begin && \a step < 0.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::Range(T begin, T end, T step)
  {
    mcIdType nbOfTuples=DataArrayTools<T>::GetNumberOfItemGivenBESRelative(begin,end,step,"DataArrayInt::Range");
    MCAuto<DataArrayType> ret=DataArrayType::New();
    ret->alloc(nbOfTuples,1);
    T *ptr=ret->getPointer();
    if(step>0)
      {
        for(T i=begin;i<end;i+=step,ptr++)
          *ptr=i;
      }
    else
      {
        for(T i=begin;i>end;i+=step,ptr++)
          *ptr=i;
      }
    return ret.retn();
  }

  /*!
   * Returns a new DataArrayInt containing a renumbering map in "Old to New" mode computed
   * from a zip representation of a surjective format (returned e.g. by
   * \ref MEDCoupling::DataArrayDouble::findCommonTuples() "DataArrayDouble::findCommonTuples()"
   * for example). The result array minimizes the permutation. <br>
   * For more info on renumbering see \ref numbering. <br>
   * \b Example: <br>
   * - \a nbOfOldTuples: 10
   * - \a arr          : [0,3, 5,7,9]
   * - \a arrIBg       : [0,2,5]
   * - \a newNbOfTuples: 7
   * - result array    : [0,1,2,0,3,4,5,4,6,4]
   *
   *  \param [in] nbOfOldTuples - number of tuples in the initial array \a arr.
   *  \param [in] arr - the array of tuple indices grouped by \a arrIBg array.
   *  \param [in] arrIBg - the array dividing all indices stored in \a arr into groups of
   *         (indices of) equal values. Its every element (except the last one) points to
   *         the first element of a group of equal values.
   *  \param [in] arrIEnd - specifies the end of \a arrIBg, so that the last element of \a
   *          arrIBg is \a arrIEnd[ -1 ].
   *  \param [out] newNbOfTuples - number of tuples after surjection application.
   *  \return DataArrayInt * - a new instance of DataArrayInt. The caller is to delete this
   *          array using decrRef() as it is no more needed.
   *  \throw If any value of \a arr breaks condition ( 0 <= \a arr[ i ] < \a nbOfOldTuples ).
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::ConvertIndexArrayToO2N(mcIdType nbOfOldTuples, const mcIdType *arr, const mcIdType *arrIBg, const mcIdType *arrIEnd, mcIdType &newNbOfTuples)
  {
    MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
    ret->alloc(nbOfOldTuples,1);
    mcIdType *pt=ret->getPointer();
    std::fill(pt,pt+nbOfOldTuples,-1);
    mcIdType nbOfGrps=ToIdType(std::distance(arrIBg,arrIEnd))-1;
    const mcIdType *cIPtr=arrIBg;
    for(mcIdType i=0;i<nbOfGrps;i++)
      pt[arr[cIPtr[i]]]=-(i+2);
    mcIdType newNb=0;
    for(mcIdType iNode=0;iNode<nbOfOldTuples;iNode++)
      {
        if(pt[iNode]<0)
          {
            if(pt[iNode]==-1)
              pt[iNode]=newNb++;
            else
              {
                mcIdType grpId=-(pt[iNode]+2);
                for(mcIdType j=cIPtr[grpId];j<cIPtr[grpId+1];j++)
                  {
                    if(arr[j]>=0 && arr[j]<nbOfOldTuples)
                      pt[arr[j]]=newNb;
                    else
                      {
                        std::ostringstream oss; oss << "DataArrayInt::ConvertIndexArrayToO2N : With element #" << j << " value is " << arr[j] << " should be in [0," << nbOfOldTuples << ") !";
                        throw INTERP_KERNEL::Exception(oss.str().c_str());
                      }
                  }
                newNb++;
              }
          }
      }
    newNbOfTuples=newNb;
    return ret.retn();
  }

  /*!
   * Returns a new DataArrayInt which is a minimal partition of elements of \a groups.
   * The i-th item of the result array is an ID of a set of elements belonging to a
   * unique set of groups, which the i-th element is a part of. This set of elements
   * belonging to a unique set of groups is called \a family, so the result array contains
   * IDs of families each element belongs to.
   *
   * \b Example: if we have two groups of elements: \a group1 [0,4] and \a group2 [ 0,1,2 ],
   * then there are 3 families:
   * - \a family1 (with ID 1) contains element [0] belonging to ( \a group1 + \a group2 ),
   * - \a family2 (with ID 2) contains elements [4] belonging to ( \a group1 ),
   * - \a family3 (with ID 3) contains element [1,2] belonging to ( \a group2 ), <br>
   * and the result array contains IDs of families [ 1,3,3,0,2 ]. <br> Note a family ID 0 which
   * stands for the element #3 which is in none of groups.
   *
   *  \param [in] groups - sequence of groups of element IDs.
   *  \param [in] newNb - total number of elements; it must be more than max ID of element
   *         in \a groups.
   *  \param [out] fidsOfGroups - IDs of families the elements of each group belong to.
   *  \return DataArrayInt * - a new instance of DataArrayInt containing IDs of families
   *         each element with ID from range [0, \a newNb ) belongs to. The caller is to
   *         delete this array using decrRef() as it is no more needed.
   *  \throw If any element ID in \a groups violates condition ( 0 <= ID < \a newNb ).
   */
  template <class T>
  DataArrayIdType *DataArrayDiscrete<T>::MakePartition(const std::vector<const DataArrayType *>& groups, mcIdType newNb, std::vector< std::vector<mcIdType> >& fidsOfGroups)
  {
    std::vector<const DataArrayType *> groups2;
    for(typename std::vector<const DataArrayType *>::const_iterator it4=groups.begin();it4!=groups.end();it4++)
      if(*it4)
        groups2.push_back(*it4);
    MCAuto<DataArrayIdType> ret=DataArrayIdType::New();
    ret->alloc(newNb,1);
    mcIdType *retPtr=ret->getPointer();
    std::fill(retPtr,retPtr+newNb,0);
    mcIdType fid=1;
    for(typename std::vector<const DataArrayType *>::const_iterator iter=groups2.begin();iter!=groups2.end();iter++)
      {
        const T *ptr=(*iter)->getConstPointer();
        std::size_t nbOfElem=(*iter)->getNbOfElems();
        mcIdType sfid=fid;
        for(mcIdType j=0;j<sfid;j++)
          {
            bool found=false;
            for(std::size_t i=0;i<nbOfElem;i++)
              {
                if(ptr[i]>=0 && ptr[i]<newNb)
                  {
                    if(retPtr[ptr[i]]==j)
                      {
                        retPtr[ptr[i]]=fid;
                        found=true;
                      }
                  }
                else
                  {
                    std::ostringstream oss; oss << "DataArrayInt::MakePartition : In group \"" << (*iter)->getName() << "\" in tuple #" << i << " value = " << ptr[i] << " ! Should be in [0," << newNb;
                    oss << ") !";
                    throw INTERP_KERNEL::Exception(oss.str().c_str());
                  }
              }
            if(found)
              fid++;
          }
      }
    fidsOfGroups.clear();
    fidsOfGroups.resize(groups2.size());
    mcIdType grId=0;
    for(typename std::vector<const DataArrayType *>::const_iterator iter=groups2.begin();iter!=groups2.end();iter++,grId++)
      {
        std::set<mcIdType> tmp;
        const T *ptr=(*iter)->getConstPointer();
        std::size_t nbOfElem=(*iter)->getNbOfElems();
        for(const T *p=ptr;p!=ptr+nbOfElem;p++)
          tmp.insert(retPtr[*p]);
        fidsOfGroups[grId].insert(fidsOfGroups[grId].end(),tmp.begin(),tmp.end());
      }
    return ret.retn();
  }

}

/// @cond INTERNAL
namespace MEDCouplingImpl
{
  template <class T>
  class OpSwitchedOn
  {
  public:
    OpSwitchedOn(T *pt):_pt(pt),_cnt(0) { }
    void operator()(const bool& b) { if(b) *_pt++=FromIdType<T>(_cnt); _cnt++; }
  private:
    T *_pt;
    MEDCoupling::mcIdType _cnt;
  };

  template <class T>
  class OpSwitchedOff
  {
  public:
    OpSwitchedOff(T *pt):_pt(pt),_cnt(0) { }
    void operator()(const bool& b) { if(!b) *_pt++=FromIdType<T>(_cnt); _cnt++; }
  private:
    T *_pt;
    MEDCoupling::mcIdType _cnt;
  };
}
/// @endcond

namespace MEDCoupling
{
  /*!
   * This method returns the list of ids in ascending mode so that v[id]==true.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::BuildListOfSwitchedOn(const std::vector<bool>& v)
  {
    std::size_t sz(std::count(v.begin(),v.end(),true));
    MCAuto<DataArrayType> ret(DataArrayType::New()); ret->alloc(sz,1);
    std::for_each(v.begin(),v.end(),MEDCouplingImpl::OpSwitchedOn<T>(ret->getPointer()));
    return ret.retn();
  }

  /*!
   * This method returns the list of ids in ascending mode so that v[id]==false.
   */
  template <class T>
  typename Traits<T>::ArrayType *DataArrayDiscrete<T>::BuildListOfSwitchedOff(const std::vector<bool>& v)
  {
    std::size_t sz(std::count(v.begin(),v.end(),false));
    MCAuto<DataArrayType> ret(DataArrayType::New()); ret->alloc(sz,1);
    std::for_each(v.begin(),v.end(),MEDCouplingImpl::OpSwitchedOff<T>(ret->getPointer()));
    return ret.retn();
  }
}

namespace MEDCoupling
{
  /*!
   * This method compares content of input vector \a v and \a this.
   * If for each id in \a this v[id]==True and for all other ids id2 not in \a this v[id2]==False, true is returned.
   * For performance reasons \a this is expected to be sorted ascendingly. If not an exception will be thrown.
   *
   * \param [in] v - the vector of 'flags' to be compared with \a this.
   *
   * \throw If \a this is not sorted ascendingly.
   * \throw If \a this has not exactly one component.
   * \throw If \a this is not allocated.
   */
  template<class T>
  bool DataArrayDiscreteSigned<T>::isFittingWith(const std::vector<bool>& v) const
  {
    this->checkAllocated();
    if(this->getNumberOfComponents()!=1)
      throw INTERP_KERNEL::Exception("DataArrayInt::isFittingWith : number of components of this should be equal to one !");
    const T *w(this->begin()),*end2(this->end());
    T refVal=-std::numeric_limits<T>::max();
    T i=0;
    std::vector<bool>::const_iterator it(v.begin());
    for(;it!=v.end();it++,i++)
      {
        if(*it)
          {
            if(w!=end2)
              {
                if(*w++==i)
                  {
                    if(i>refVal)
                      refVal=i;
                    else
                      {
                        std::ostringstream oss; oss << "DataArrayInt::isFittingWith : At pos #" << std::distance(this->begin(),w-1) << " this is not sorted ascendingly !";
                        throw INTERP_KERNEL::Exception(oss.str().c_str());
                      }
                  }
                else
                  return false;
              }
            else
              return false;
          }
      }
    return w==end2;
  }
}

#endif
