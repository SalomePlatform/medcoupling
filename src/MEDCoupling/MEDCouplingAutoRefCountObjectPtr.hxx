// Copyright (C) 2007-2012  CEA/DEN, EDF R&D
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

#ifndef __PARAMEDMEM_MEDCOUPLINGAUTOREFCOUNTOBJECTPTR_HXX__
#define __PARAMEDMEM_MEDCOUPLINGAUTOREFCOUNTOBJECTPTR_HXX__

#include "MEDCouplingRefCountObject.hxx"

namespace ParaMEDMEM
{
  template<class T>
  class MEDCouplingAutoRefCountObjectPtr
  {
  public:
    MEDCouplingAutoRefCountObjectPtr(const MEDCouplingAutoRefCountObjectPtr& other):_ptr(0) { referPtr(other._ptr); }
    MEDCouplingAutoRefCountObjectPtr(T *ptr=0):_ptr(ptr) { }
    ~MEDCouplingAutoRefCountObjectPtr() { destroyPtr(); }
    bool operator==(const MEDCouplingAutoRefCountObjectPtr& other) { return _ptr==other._ptr; }
    MEDCouplingAutoRefCountObjectPtr &operator=(const MEDCouplingAutoRefCountObjectPtr& other) { if(_ptr!=other._ptr) { destroyPtr(); referPtr(other._ptr); } return *this; }
    MEDCouplingAutoRefCountObjectPtr &operator=(T *ptr) { if(_ptr!=ptr) { destroyPtr(); _ptr=ptr; } return *this; }
    T *operator->() { return _ptr ; }
    const T *operator->() const { return _ptr; }
    T& operator*() { return *_ptr; }
    const T& operator*() const { return *_ptr; }
    operator T *() { return _ptr; }
    operator const T *() const { return _ptr; }
  private:
    void referPtr(T *ptr) { _ptr=ptr; if(_ptr) _ptr->incrRef(); }
    void destroyPtr() { if(_ptr) _ptr->decrRef(); }
  private:
    T *_ptr;
  };
}

#endif
