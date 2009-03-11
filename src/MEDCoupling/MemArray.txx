//  Copyright (C) 2007-2008  CEA/DEN, EDF R&D
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
#ifndef __PARAMEDMEM_MEMARRAY_TXX__
#define __PARAMEDMEM_MEMARRAY_TXX__

#include "MemArray.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"

#include <sstream>
#include <algorithm>

namespace ParaMEDMEM
{
  template<class T>
  MemArray<T>::MemArray(const MemArray<T>& other):_nb_of_elem(-1),_ownership(false),_pointer(0),_dealloc(CPP_DEALLOC)
  {
    if(other._pointer)
      {
        T *pointer=new T[other._nb_of_elem];
        std::copy(other._pointer,other._pointer+other._nb_of_elem,pointer);
        useArray(pointer,true,CPP_DEALLOC,other._nb_of_elem);
      }
  }

  template<class T>
  void MemArray<T>::useArray(void *array, bool ownership, DeallocType type, int nbOfElem)
  {
    _nb_of_elem=nbOfElem;
    destroy();
    _pointer=(T *)array;
    _ownership=ownership;
    _dealloc=type;
  }

  template<class T>
  void MemArray<T>::writeOnPlace(int id, T element0, const T *others, int sizeOfOthers)
  {
    if(id+sizeOfOthers>=_nb_of_elem)
      reAlloc(2*_nb_of_elem+sizeOfOthers+1);
    _pointer[id]=element0;
    memcpy(_pointer+id+1,others,sizeOfOthers*sizeof(T));
  }

  template<class T>
  void MemArray<T>::alloc(int nbOfElements)
  {
    destroy();
    _nb_of_elem=nbOfElements;
    _pointer=new T[_nb_of_elem];
    _ownership=true;
    _dealloc=CPP_DEALLOC;
  }
  
  template<class T>
  void MemArray<T>::reAlloc(int newNbOfElements)
  {
    T *pointer=new T[newNbOfElements];
    memcpy(pointer,_pointer,std::min<int>(_nb_of_elem,newNbOfElements)*sizeof(int));
    destroyPointer(_pointer,_dealloc);
    _pointer=pointer;
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
      destroyPointer(_pointer,_dealloc);
    _pointer=0;
    _ownership=false;
  }
  
  template<class T>
  MemArray<T> &MemArray<T>::operator=(const MemArray<T>& other)
  {
    alloc(other._nb_of_elem);
    memcpy(_pointer,other._pointer,_nb_of_elem*sizeof(T));
    return *this;
  }
}

#endif
