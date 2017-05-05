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
// Author : Anthony Geay (EDF R&D)

#ifndef __MEDCOUPLING_MEDCOUPLINGMAP_TXX__
#define __MEDCOUPLING_MEDCOUPLINGMAP_TXX__

#include "MEDCouplingMap.hxx"

namespace MEDCoupling
{
  template<class T>
  MCAuto< MapKeyVal<T> > MapKeyVal<T>::New()
  {
    MCAuto< MapKeyVal<T> > ret(new MapKeyVal<T>);
    return ret;
  }
  
  template<class T>
  std::size_t MapKeyVal<T>::getHeapMemorySizeWithoutChildren() const
  {
    return _m.size()*sizeof(std::pair<T,T>);
  }
  
  template<class T>
  std::vector<const BigMemoryObject*> MapKeyVal<T>::getDirectChildrenWithNull() const
  {
    return std::vector<const BigMemoryObject*>();//not a bug no child. Leaf object !
  }
}

#endif
