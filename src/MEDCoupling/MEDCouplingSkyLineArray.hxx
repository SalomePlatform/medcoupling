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
#include "NormalizedGeometricTypes"

#include <vector>

namespace MEDCoupling
{
  /**!
   * Class allowing the easy manipulation of the indexed array format, where the first array is a set of offsets to
   * be used in the second array, to extract packs of values.
   *
   * This class allows to pursuie this logic up to 3 levels, i.e. the first array points to packs in the second, which
   * itself points to identifiers in the third array.
   *
   * This particularly useful for connectivity of pure polygonal/polyhedral meshes.
   */
  class MEDCOUPLING_EXPORT MEDCouplingSkyLineArray : public RefCountObject
  {
  public:
    static MEDCouplingSkyLineArray * New();
    static MEDCouplingSkyLineArray * New( const std::vector<int>& index, const std::vector<int>& value);
    static MEDCouplingSkyLineArray * New( DataArrayInt* index, DataArrayInt* value );

    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;

    void set( DataArrayInt* index, DataArrayInt* value );

    int getNumberOf() const { return _index->getNbOfElems()-1; }
    int getLength()   const { return _values->getNbOfElems(); }
    const int* getIndex() const { return _index->begin(); }
    const int* getValues() const { return _values->begin(); }

    DataArrayInt* getIndexArray() const;
    DataArrayInt* getValuesArray() const;

    std::string simpleRepr() const;

//    replaceWithPackFromOther()

  private:
    MEDCouplingSkyLineArray();
    ~MEDCouplingSkyLineArray();

    MCAuto<DataArrayInt> _index;
    MCAuto<DataArrayInt> _values;
    MCAuto<DataArrayInt> _sub_values;
  };
}
# endif
