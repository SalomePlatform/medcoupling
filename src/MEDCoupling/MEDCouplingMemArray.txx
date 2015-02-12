// Copyright (C) 2007-2015  CEA/DEN, EDF R&D
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

#include <sstream>
#include <cstdlib>
#include <algorithm>

namespace ParaMEDMEM
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
}

#endif
