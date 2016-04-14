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

#include "MEDCouplingSkyLineArray.hxx"

using namespace MEDCoupling;

MEDCouplingSkyLineArray::MEDCouplingSkyLineArray():
  _index( DataArrayInt::New() ), _value( DataArrayInt::New() )
{
}

MEDCouplingSkyLineArray::MEDCouplingSkyLineArray(const MEDCouplingSkyLineArray &myArray)
{
  _index=myArray._index;
  _value=myArray._value;
}

MEDCouplingSkyLineArray::~MEDCouplingSkyLineArray()
{
}

MEDCouplingSkyLineArray::MEDCouplingSkyLineArray(DataArrayInt* index, DataArrayInt* value)
{
  set( index, value );
}

MEDCouplingSkyLineArray::MEDCouplingSkyLineArray( const std::vector<int>& index,
                                                  const std::vector<int>& value ):
  _index( DataArrayInt::New() ), _value( DataArrayInt::New() )
{
  _index->reserve( index.size() );
  _index->insertAtTheEnd( index.begin(), index.end() );
  _value->reserve( value.size() );
  _value->insertAtTheEnd( value.begin(), value.end() );
}

void MEDCouplingSkyLineArray::set( DataArrayInt* index, DataArrayInt* value )
{
  _index=index;
  _value=value;
  if ( (DataArrayInt*)_index ) _index->incrRef();
  else                         _index = DataArrayInt::New();
  if ( (DataArrayInt*)_value ) _value->incrRef();
  else                         _value = DataArrayInt::New();
}

DataArrayInt* MEDCouplingSkyLineArray::getIndexArray() const
{
  return ((MEDCouplingSkyLineArray*)this)->_index;
}

DataArrayInt* MEDCouplingSkyLineArray::getValueArray() const
{
  return ((MEDCouplingSkyLineArray*)this)->_value;
}

std::string MEDCouplingSkyLineArray::simpleRepr() const
{
  std::ostringstream oss;
  oss << "MEDCouplingSkyLineArray" << std::endl;
  oss << "   Nb of items: " << getNumberOf() << std::endl;
  oss << "   Nb of values: " << getLength() << std::endl;
  oss << "   Index:" << std::endl;
  oss << "   ";
  const int * i = _index->begin();
  for ( ; i != _index->end(); ++i )
    oss << *i << " ";
  oss << std::endl;
  oss << "   Value:" << std::endl;
  oss << "   ";
  const int * v = _value->begin();
  int cnt = 0;
  for ( i = _index->begin(); v != _value->end(); ++v, ++cnt )
    {
      if ( cnt == *i )
        {
          oss << "| ";
          ++i;
        }
      oss << *v << " ";
    }
  oss << std::endl;

  return oss.str();
}
