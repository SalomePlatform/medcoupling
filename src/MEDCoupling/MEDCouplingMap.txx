// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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
  template<class ID, class T>
  MCAuto< MapKeyVal<ID, T> > MapKeyVal<ID, T>::New()
  {
    MCAuto< MapKeyVal<ID, T> > ret(new MapKeyVal<ID, T>);
    return ret;
  }
  
  template<class ID, class T>
  std::size_t MapKeyVal<ID, T>::getHeapMemorySizeWithoutChildren() const
  {
    return _m.size()*sizeof(std::pair<ID, T>);
  }
  
  template<class ID, class T>
  std::vector<const BigMemoryObject*> MapKeyVal<ID, T>::getDirectChildrenWithNull() const
  {
    return std::vector<const BigMemoryObject*>();//not a bug no child. Leaf object !
  }
}

#endif
