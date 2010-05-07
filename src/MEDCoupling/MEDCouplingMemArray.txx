//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#ifndef __PARAMEDMEM_MEDCOUPLINGMEMARRAY_TXX__
#define __PARAMEDMEM_MEDCOUPLINGMEMARRAY_TXX__

#include "MEDCouplingMemArray.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"

#include <sstream>
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
  MemArray<T>::MemArray(const MemArray<T>& other):_nb_of_elem(-1),_ownership(false),_dealloc(CPP_DEALLOC)
  {
    if(!other._pointer.isNull())
      {
        T *pointer=new T[other._nb_of_elem];
        std::copy(other._pointer.getConstPointer(),other._pointer.getConstPointer()+other._nb_of_elem,pointer);
        useArray(pointer,true,CPP_DEALLOC,other._nb_of_elem);
      }
  }

  template<class T>
  void MemArray<T>::useArray(const T *array, bool ownership, DeallocType type, int nbOfElem)
  {
    _nb_of_elem=nbOfElem;
    destroy();
    if(ownership)
      _pointer.setInternal((T *)array);
    else
      _pointer.setExternal(array);
    _ownership=ownership;
    _dealloc=type;
  }

  template<class T>
  void MemArray<T>::writeOnPlace(int id, T element0, const T *others, int sizeOfOthers)
  {
    if(id+sizeOfOthers>=_nb_of_elem)
      reAlloc(2*_nb_of_elem+sizeOfOthers+1);
    T *pointer=_pointer.getPointer();
    pointer[id]=element0;
    std::copy(others,others+sizeOfOthers,pointer+id+1);
  }

  template<class T>
  bool MemArray<T>::isEqual(const MemArray<T>& other, T prec) const
  {
    if(_nb_of_elem!=other._nb_of_elem)
      return false;
    const T *pt1=_pointer.getConstPointer();
    const T *pt2=other._pointer.getConstPointer();
    if(pt1==0 && pt2==0)
      return true;
    if(pt1==0 || pt2==0)
      return false;
    for(int i=0;i<_nb_of_elem;i++)
      if(pt1[i]-pt2[i]<-prec || (pt1[i]-pt2[i])>prec)
        return false;
    return true;
  }

  template<class T>
  void MemArray<T>::alloc(int nbOfElements)
  {
    destroy();
    _nb_of_elem=nbOfElements;
    _pointer.setInternal(new T[_nb_of_elem]);
    _ownership=true;
    _dealloc=CPP_DEALLOC;
  }
  
  template<class T>
  void MemArray<T>::reAlloc(int newNbOfElements)
  {
    T *pointer=new T[newNbOfElements];
    std::copy(_pointer.getConstPointer(),_pointer.getConstPointer()+std::min<int>(_nb_of_elem,newNbOfElements),pointer);
    if(_ownership)
      destroyPointer((T *)_pointer.getConstPointer(),_dealloc);
    _pointer.setInternal(pointer);
    _nb_of_elem=newNbOfElements;
    _ownership=true;
    _dealloc=CPP_DEALLOC;
  }

  template<class T>
  void MemArray<T>::destroyPointer(T *pt, DeallocType type)
  {
    switch(type)
      {
      case CPP_DEALLOC:
        {
          delete [] pt;
          return ;
        }
      case C_DEALLOC:
        {
          free(pt);
          return ;
        }
      default:
        std::stringstream stream;
        stream << "Invalid deallocation requested for pointer " << pt;
        throw INTERP_KERNEL::Exception(stream.str().c_str());
      }
  }

  template<class T>
  void MemArray<T>::destroy()
  {
    if(_ownership)
      destroyPointer((T *)_pointer.getConstPointer(),_dealloc);
    _pointer.null();
    _ownership=false;
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
