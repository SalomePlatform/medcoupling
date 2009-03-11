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
#ifndef __PARAMEDMEM_REFCOUNTOBJECT_HXX__
#define __PARAMEDMEM_REFCOUNTOBJECT_HXX__

#include "TimeLabel.hxx"

namespace ParaMEDMEM
{
  typedef enum
    {
      C_DEALLOC = 2,
      CPP_DEALLOC = 3
    } DeallocType;

  typedef enum
    {
      ON_CELLS = 0,
      ON_NODES = 1
    } TypeOfField;

  class RefCountObject : public TimeLabel
  {
  protected:
    RefCountObject():_cnt(1) { }
    RefCountObject(const RefCountObject& other):_cnt(1) { }
  public:
    bool decrRef() { bool ret=((--_cnt)==0); if(ret)delete this; return ret; }
    void incrRef() const { _cnt++; }
  protected:
    virtual ~RefCountObject() { }
  private:
    mutable int _cnt;
  };
}

#endif
