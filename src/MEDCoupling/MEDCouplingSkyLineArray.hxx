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

#ifndef __PARAMEDMEM_MEDCOUPLINGSKYLINEARRAY_HXX__
#define __PARAMEDMEM_MEDCOUPLINGSKYLINEARRAY_HXX__

#include "MEDCoupling.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MCAuto.hxx"

#include <vector>

namespace MEDCoupling
{
  class MEDCOUPLING_EXPORT MEDCouplingSkyLineArray
  {
  private:
    MCAuto<DataArrayInt> _index;
    MCAuto<DataArrayInt> _value;
  public:
    MEDCouplingSkyLineArray();
    MEDCouplingSkyLineArray( const MEDCouplingSkyLineArray &myArray );
    MEDCouplingSkyLineArray( const std::vector<int>& index, const std::vector<int>& value );
    MEDCouplingSkyLineArray( DataArrayInt* index, DataArrayInt* value );
    ~MEDCouplingSkyLineArray();

    void set( DataArrayInt* index, DataArrayInt* value );

    int getNumberOf() const { return _index->getNbOfElems()-1; }
    int getLength()   const { return _value->getNbOfElems(); }
    const int* getIndex() const { return _index->begin(); }
    const int* getValue() const { return _value->begin(); }

    DataArrayInt* getIndexArray() const;
    DataArrayInt* getValueArray() const;

    std::string simpleRepr() const;
  };
}
# endif
