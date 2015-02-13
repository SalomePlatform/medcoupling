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

#ifndef __MEDPARTITIONER_MEDSKYLINEARRAY_H__
#define __MEDPARTITIONER_MEDSKYLINEARRAY_H__

#include "MEDPARTITIONER.hxx"

#include <vector>

namespace MEDPARTITIONER 
{
  class MEDPARTITIONER_EXPORT SkyLineArray
  {
  private:
    std::vector<int> _index;
    std::vector<int> _value;
  public:
    /*! if used SkyLineArray will keep empty */
    SkyLineArray();
    SkyLineArray( const SkyLineArray &myArray );
    SkyLineArray( const std::vector<int>& index, const std::vector<int>& value );
    ~SkyLineArray();
  
    int getNumberOf() const { return _index.size()-1; }
    int getLength() const { return _value.size() ; }
    const int* getIndex() const { return (const int*)(&_index[0]); }
    const int* getValue() const { return (const int*)(&_value[0]); }
  };
}
# endif
