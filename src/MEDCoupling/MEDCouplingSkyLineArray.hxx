// Copyright (C) 2007-2021  CEA/DEN, EDF R&D
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

#pragma once

#include "MEDCoupling.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MCAuto.hxx"
#include "NormalizedGeometricTypes"

#include <vector>
#include <functional>

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

    static std::vector< MCAuto<DataArrayIdType> > RetrieveVecIndex(const std::vector< MCAuto<MEDCouplingSkyLineArray> >& vecSka)
    {
       auto fct = [](MEDCouplingSkyLineArray *ska) { return ska->getIndexArray(); };
       return RetrieveVecOfSkyLineArrayGen(vecSka,fct);
    }
    
    static std::vector< MCAuto<DataArrayIdType> > RetrieveVecValues(const std::vector< MCAuto<MEDCouplingSkyLineArray> >& vecSka)
    {
       auto fct = [](MEDCouplingSkyLineArray *ska) { return ska->getValuesArray(); };
       return RetrieveVecOfSkyLineArrayGen(vecSka,fct);
    }
    
    static std::vector< MCAuto<DataArrayIdType> > RetrieveVecOfSkyLineArrayGen(const std::vector< MCAuto<MEDCouplingSkyLineArray> >& vecSka, std::function<DataArrayIdType *(MEDCouplingSkyLineArray *)> fct)
    {
       std::size_t sz(vecSka.size());
       std::vector< MCAuto<DataArrayIdType> > ret(sz);
       std::vector< MCAuto<DataArrayIdType> >::iterator it(ret.begin());
       std::for_each(vecSka.begin(),vecSka.end(),[&it,fct](MCAuto<MEDCouplingSkyLineArray> elt) { *it++ = MCAuto<DataArrayIdType>::TakeRef(fct(elt)); } );
       return ret;
    }

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

    MEDCouplingSkyLineArray *deepCopy() const;

    std::string simpleRepr() const;

    void thresholdPerPack(mcIdType threshold, MCAuto<MEDCouplingSkyLineArray>& left, MCAuto<MEDCouplingSkyLineArray>& right) const;

    MEDCouplingSkyLineArray *groupPacks(const DataArrayIdType *indexedPacks) const;
    MEDCouplingSkyLineArray *uniqueNotSortedByPack() const;
    static MEDCouplingSkyLineArray *AggregatePacks(const std::vector<const MEDCouplingSkyLineArray *>& sks);

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

  template<typename T>
  class SkyLineArrayGenIterator : public std::iterator< std::input_iterator_tag, const mcIdType *, mcIdType, const mcIdType **, const mcIdType *>
  {
    std::size_t _num = 0;
    std::vector<const MEDCouplingSkyLineArray *> *_data = nullptr;
  public:
    explicit SkyLineArrayGenIterator(std::size_t num , std::vector<const MEDCouplingSkyLineArray *> *data) : _num(num),_data(data) {}
    SkyLineArrayGenIterator<T>& operator++() { ++_num; return *this; }
    bool operator==(const SkyLineArrayGenIterator& other) const { return _num == other._num; }
    bool operator!=(const SkyLineArrayGenIterator& other) const { return !(*this == other); }
    reference operator*() const { T tt; return tt((*_data)[_num]); }
  };

  struct SkyLineArrayIndexPtrFunctor { const mcIdType *operator()(const MEDCouplingSkyLineArray *ska) { return ska->getIndex(); } };

  using SkyLineArrayIndexIterator = SkyLineArrayGenIterator<SkyLineArrayIndexPtrFunctor>;
  
  struct SkyLineArrayValuesPtrFunctor { const mcIdType *operator()(const MEDCouplingSkyLineArray *ska) { return ska->getValues(); } };

  using SkyLineArrayValuesIterator = SkyLineArrayGenIterator<SkyLineArrayValuesPtrFunctor>;
}
