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

#ifndef __MEDCOUPLINGMCAUTO_TXX__
#define __MEDCOUPLINGMCAUTO_TXX__

#include "MCAuto.hxx"

#include <vector>

namespace MEDCoupling
{
  template<class T>
  std::vector< const T * > VecAutoToVecOfCstPt(const std::vector< typename MEDCoupling::MCAuto<T> >& vect)
  {
    std::size_t sz(vect.size());
    std::vector< const T * > ret(sz);
    for(std::size_t i=0;i<sz;i++)
      ret[i]=(const T *)vect[i];
    return ret;
  }

  template<class T>
  std::vector< T * > VecAutoToVecOfPt(const std::vector< typename MEDCoupling::MCAuto<T> >& vect)
  {
    std::size_t sz(vect.size());
    std::vector< T * > ret(sz);
    for(std::size_t i=0;i<sz;i++)
      ret[i]=(T *)vect[i];
    return ret;
  }
}

#endif
