// Copyright (C) 2007-2020  CEA/DEN, EDF R&D
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
   * Class allowing the easy manipulation of the indexed array format, where the first array (the "indices") is a set of offsets to
   * be used in the second array (the "values"), to extract packs of values.
   *
   * Example:
   *    index = [0,2,7,10]
   *    values = [1,24,2,33,6,7,10,11,9,28]
   * which describes 3 packs of (integer) values : [1,24] and [2,33,6,7,10] and [11,9,28]
   *
   * Thus the index array is always monotic ascendant.
   *
   * This class allows to pursue this logic up to 3 levels, i.e. the first array (the "indices") points to packs in the
   * second (the "values"), which itself points to identifiers in a third array (the "sub-values").
   *
   * This particularly useful for connectivity of pure polygonal/polyhedral meshes.
   *
   * Example:
   *     super-index = [0,1,3]
   *     index = [0,3,6,10]
   *     values = [28,1,4,2,35,8,9,10,1,12]
   * which represent two 3 packs and two super-packs. The first super-pack is [[28,1,4]] and has only one pack [28,1,4].
   * The second super-pack is [[2,35,8], [9,10,1,12]] and has two packs [2,35,8] and [9,10,1,12].
   * Note that contrary to index, the integers in super-index are interpreted as being inclusive: the first super-pack
   * goes from offset 0 to offset 1, inclusive. This is not the same for index, where the upper bound is exclusive.
   */
  class MEDCOUPLING_EXPORT MEDCouplingSkyLineArray : public RefCountObject
  {
  public:
    static MEDCouplingSkyLineArray * New();
    static MEDCouplingSkyLineArray * New( const std::vector<mcIdType>& index, const std::vector<mcIdType>& value);
    static MEDCouplingSkyLineArray * New( DataArrayIdType* index, DataArrayIdType* value );
    static MEDCouplingSkyLineArray * New( const MEDCouplingSkyLineArray & other );

    static MEDCouplingSkyLineArray * BuildFromPolyhedronConn( const DataArrayIdType* c, const DataArrayIdType* cI );
    std::string getClassName() const override { return std::string("MEDCouplingSkyLineArray"); }
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;

    void set( DataArrayIdType* index, DataArrayIdType* value );
    void set3( DataArrayIdType* superIndex, DataArrayIdType* index, DataArrayIdType* value );

    mcIdType getSuperNumberOf()   const { return ToIdType(_super_index->getNbOfElems())-1; }
    mcIdType getNumberOf() const { return ToIdType(_index->getNbOfElems())-1; }
    mcIdType getLength()   const { return ToIdType(_values->getNbOfElems()); }

    const mcIdType* getSuperIndex() const { return _super_index->begin(); }
    const mcIdType* getIndex() const { return _index->begin(); }
    const mcIdType* getValues() const { return _values->begin(); }

    DataArrayIdType* getSuperIndexArray() const;
    DataArrayIdType* getIndexArray() const;
    DataArrayIdType* getValuesArray() const;

    std::string simpleRepr() const;

    void getSimplePackSafe(const mcIdType absolutePackId, std::vector<mcIdType> & pack) const;
    const mcIdType * getSimplePackSafePtr(const mcIdType absolutePackId, mcIdType & packSize) const;
    void findPackIds(const std::vector<mcIdType> & superPackIndices, const mcIdType *packBg, const mcIdType *packEnd,
                     std::vector<mcIdType>& out) const;

    void deletePack(const mcIdType superIdx, const mcIdType idx);
    void deleteSimplePack(const mcIdType idx);
    void pushBackPack(const mcIdType superIdx, const mcIdType * packBg, const mcIdType * packEnd);

    void replaceSimplePack(const mcIdType idx, const mcIdType * packBg, const mcIdType * packEnd);
    void replacePack(const mcIdType superIdx, const mcIdType idx, const mcIdType * packBg, const mcIdType * packEnd);

    void deleteSimplePacks(const DataArrayIdType* idx);
    void replaceSimplePacks(const DataArrayIdType* idx, const std::vector<const DataArrayIdType*>& packs);
    
    void convertToPolyhedronConn( MCAuto<DataArrayIdType>& c,  MCAuto<DataArrayIdType>& cI) const;

  private:
    MEDCouplingSkyLineArray();
    ~MEDCouplingSkyLineArray();

    void checkSuperIndex(const std::string& func) const;
    void validSuperIndex(const std::string& func, mcIdType superIndex) const;
    void validIndex(const std::string& func, mcIdType index) const;
    void validSuperIndexAndIndex(const std::string& func, mcIdType superIndex, mcIdType index) const;

    MCAuto<DataArrayIdType> _super_index;
    MCAuto<DataArrayIdType> _index;
    MCAuto<DataArrayIdType> _values;
  };

}
# endif
