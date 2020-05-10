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

#include "MEDCouplingSkyLineArray.hxx"

#include <sstream>
#include <deque>
#include <set>

using namespace MEDCoupling;

MEDCouplingSkyLineArray::MEDCouplingSkyLineArray():
  _super_index( DataArrayIdType::New() ), _index( DataArrayIdType::New() ), _values( DataArrayIdType::New() )
{
}

MEDCouplingSkyLineArray::~MEDCouplingSkyLineArray()
{
}

MEDCouplingSkyLineArray* MEDCouplingSkyLineArray::New()
{
  return new MEDCouplingSkyLineArray();
}

MEDCouplingSkyLineArray* MEDCouplingSkyLineArray::New( const std::vector<mcIdType>& index,
                                                       const std::vector<mcIdType>& value )
{
  MEDCouplingSkyLineArray * ret = new MEDCouplingSkyLineArray();
  ret->_index->reserve( index.size() );
  ret->_index->insertAtTheEnd( index.begin(), index.end() );
  ret->_values->reserve( value.size() );
  ret->_values->insertAtTheEnd( value.begin(), value.end() );
  return ret;
}

MEDCouplingSkyLineArray* MEDCouplingSkyLineArray::New( DataArrayIdType* index, DataArrayIdType* value )
{
  MEDCouplingSkyLineArray* ret = new MEDCouplingSkyLineArray();
  ret->set(index, value);
  return ret;
}

MEDCouplingSkyLineArray* MEDCouplingSkyLineArray::New( const MEDCouplingSkyLineArray & other )
{
  MEDCouplingSkyLineArray* ret = new MEDCouplingSkyLineArray();
  ret->_super_index = other._super_index;
  ret->_index = other._index;
  ret->_values = other._values;
  return ret;
}

/**! Build a three level SkyLine array from the dynamic connectivity of a dynamic mesh (i.e. containing only
 * polyhedrons or polygons).
 * The input arrays are deep copied, contrary to the other ctors.
 */
MEDCouplingSkyLineArray * MEDCouplingSkyLineArray::BuildFromPolyhedronConn( const DataArrayIdType* c, const DataArrayIdType* cI )
{
  using namespace std;

  MEDCouplingSkyLineArray* ret = new MEDCouplingSkyLineArray();

  const mcIdType * cP(c->begin()), * cIP(cI->begin());
  mcIdType prev = -1;
  if (c->getNbOfElems() != *(cI->end()-1))
    throw INTERP_KERNEL::Exception("MEDCouplingSkyLineArray::BuildFromDynamicConn: misformatted connectivity (wrong nb of tuples)!");
  for (mcIdType i=0; i < cI->getNbOfElems(); i++)
    {
      mcIdType j = cIP[i];
      if (cIP[i] < prev)
        throw INTERP_KERNEL::Exception("MEDCouplingSkyLineArray::BuildFromDynamicConn: misformatted connectivity (indices not monotonic ascending)!");
      prev = cIP[i];
      if (i!=cI->getNbOfElems()-1)
        if (cP[j] != INTERP_KERNEL::NORM_POLYHED)
          throw INTERP_KERNEL::Exception("MEDCouplingSkyLineArray::BuildFromDynamicConn: connectivity containing other types than POLYHED!");
    }

  vector<mcIdType> superIdx, idx, vals;
  mcIdType cnt = 0, cnt2 = 0;
  superIdx.reserve(cI->getNbOfElems());
  superIdx.push_back(0);
  idx.push_back(0);
  vals.resize(c->getNbOfElems()); // too much because of the type and the -1, but still better than push_back().
  for (mcIdType i=0; i < cI->getNbOfElems()-1; i++)
    {
      mcIdType start = cIP[i]+1, end = cIP[i+1];
      mcIdType * work = vals.data() + cnt;
      const mcIdType * w = cP+start;
      const mcIdType * w2 = find(w, cP+end, -1);
      while (w2 != cP+end)
        {
          copy(w, w2, work);
          mcIdType d = ToIdType(distance(w, w2));
          cnt += d; work +=d;
          idx.push_back(cnt); cnt2++;
          w = w2+1;  // skip the -1
          w2 = find(w, cP+end, -1);
        }
      copy(w, cP+end, work);
      cnt += ToIdType(distance(w, cP+end));
      idx.push_back(cnt); cnt2++;
      superIdx.push_back(cnt2);
    }
  ret->_super_index->alloc(superIdx.size(),1);
  copy(superIdx.begin(), superIdx.end(), ret->_super_index->getPointer());
  ret->_index->alloc(idx.size(),1);
  copy(idx.begin(), idx.end(), ret->_index->getPointer());
  ret->_values->alloc(cnt,1);
  copy(vals.begin(), vals.begin()+cnt, ret->_values->getPointer());

  return ret;
}

/**
 * Convert a three-level SkyLineArray into a polyhedral connectivity.
 * The super-packs are interpreted as cell description, and the packs represent the face connectivity.
 */
void MEDCouplingSkyLineArray::convertToPolyhedronConn( MCAuto<DataArrayIdType>& c,  MCAuto<DataArrayIdType>& cI) const
{
  // TODO: in this case an iterator would be nice
  using namespace std;

  checkSuperIndex("convertToPolyhedronConn");

  const mcIdType * siP(_super_index->begin()), * iP(_index->begin()), *vP(_values->begin());
  mcIdType cnt = 0;
  cI->alloc(_super_index->getNbOfElems(),1);  // same number of super packs as number of cells
  mcIdType * cIVecP(cI->getPointer());
  MCAuto <DataArrayIdType> dsi = _index->deltaShiftIndex();
  mcIdType sz = dsi->accumulate((std::size_t)0) + ToIdType(dsi->getNbOfElems());  // think about it: one slot for the type, -1 at the end of each face of the cell
  c->alloc(sz, 1);
  mcIdType * cVecP(c->getPointer());

  for ( mcIdType i=0; i < _super_index->getNbOfElems()-1; i++)
     {
       cIVecP[i]= cnt;
       mcIdType endId = siP[i+1];
       cVecP[cnt++] = INTERP_KERNEL::NORM_POLYHED;
       for (mcIdType j=siP[i]; j < endId; j++)
         {
           mcIdType startId2 = iP[j], endId2 = iP[j+1];
           copy(vP+startId2, vP+endId2, cVecP+cnt);
           cnt += endId2-startId2;
           if(j != endId-1)
             cVecP[cnt++] = -1;
         }
     }
  cIVecP[_super_index->getNbOfElems()-1] = cnt;
}

std::size_t MEDCouplingSkyLineArray::getHeapMemorySizeWithoutChildren() const
{
  return _index->getHeapMemorySizeWithoutChildren()+_values->getHeapMemorySizeWithoutChildren()+_super_index->getHeapMemorySizeWithoutChildren();
}

std::vector<const BigMemoryObject *> MEDCouplingSkyLineArray::getDirectChildrenWithNull() const
{
  std::vector<const BigMemoryObject *> ret;
  ret.push_back(_super_index);
  ret.push_back(_index);
  ret.push_back(_values);
  return ret;
}


void MEDCouplingSkyLineArray::set( DataArrayIdType* index, DataArrayIdType* value )
{
  _index=index;
  _values=value;
  if ( (DataArrayIdType*)_index ) _index->incrRef();
  else                            _index = DataArrayIdType::New();
  if ( (DataArrayIdType*)_values ) _values->incrRef();
  else                             _values = DataArrayIdType::New();
}

void MEDCouplingSkyLineArray::set3( DataArrayIdType* superIndex, DataArrayIdType* index, DataArrayIdType* value )
{
  _super_index=superIndex;
  if ( (DataArrayIdType*)_super_index ) _super_index->incrRef();
  else                                  _super_index = DataArrayIdType::New();
  set(index, value);
}

DataArrayIdType* MEDCouplingSkyLineArray::getSuperIndexArray() const
{
  return const_cast<MEDCouplingSkyLineArray*>(this)->_super_index;
}


DataArrayIdType* MEDCouplingSkyLineArray::getIndexArray() const
{
  return const_cast<MEDCouplingSkyLineArray*>(this)->_index;
}

DataArrayIdType* MEDCouplingSkyLineArray::getValuesArray() const
{
  return const_cast<MEDCouplingSkyLineArray*>(this)->_values;
}

MEDCouplingSkyLineArray *MEDCouplingSkyLineArray::deepCopy() const
{
  MCAuto<DataArrayIdType> indexCpy(this->_index->deepCopy());
  MCAuto<DataArrayIdType> valuesCpy(this->_values->deepCopy());
  MCAuto<MEDCouplingSkyLineArray> ret(MEDCouplingSkyLineArray::New(indexCpy,valuesCpy));
  if(_super_index.isNotNull())
  {
    MCAuto<DataArrayIdType> superIndexCpy(this->_super_index->deepCopy());
    ret->_super_index = superIndexCpy;
  }
  return ret.retn();
}

void MEDCouplingSkyLineArray::checkSuperIndex(const std::string& func) const
{
  if (!_super_index->getNbOfElems())
    {
      std::ostringstream oss;
      oss << "MEDCouplingSkyLineArray::"<< func << ": not a three level SkyLineArray! Method is not available for two-level SkyLineArray.";
      throw INTERP_KERNEL::Exception(oss.str());
    }
}

void MEDCouplingSkyLineArray::validSuperIndex(const std::string& func, mcIdType superIndex) const
{
  if(superIndex < 0 || superIndex >= _super_index->getNbOfElems())
    {
      std::ostringstream oss;
      oss << "MEDCouplingSkyLineArray::" << func <<  ": invalid super index!";
      throw INTERP_KERNEL::Exception(oss.str());
    }
}

void MEDCouplingSkyLineArray::validIndex(const std::string& func, mcIdType idx) const
{
  if(idx < 0 || idx >= _index->getNbOfElems())
    {
      std::ostringstream oss;
      oss << "MEDCouplingSkyLineArray::" << func <<  ": invalid index!";
      throw INTERP_KERNEL::Exception(oss.str());
    }
}

void MEDCouplingSkyLineArray::validSuperIndexAndIndex(const std::string& func, mcIdType superIndex, mcIdType index) const
{
  validSuperIndex(func, superIndex);
  mcIdType idx = _super_index->begin()[superIndex] + index;
  if(idx < 0 || idx >= _index->getNbOfElems())
    {
      std::ostringstream oss;
      oss << "MEDCouplingSkyLineArray::" << func <<  ": invalid index!";
      throw INTERP_KERNEL::Exception(oss.str());
    }
}

std::string MEDCouplingSkyLineArray::simpleRepr() const
{
  std::ostringstream oss;
  oss << "MEDCouplingSkyLineArray (" << this << ")" << std::endl;
  MCAuto<DataArrayIdType> super_index = _super_index->deepCopy();
  if (_super_index->getNbOfElems())
    oss << "   Nb of super-packs: " << getSuperNumberOf() << std::endl;
  else
    {
      super_index->alloc(2,1);
      super_index->setIJSilent(0,0,0);
      super_index->setIJSilent(1,0,_index->getNbOfElems()-1);
    }
  oss << "   Nb of packs: " << getNumberOf() << std::endl;
  oss << "   Nb of values: " << getLength() << std::endl;

  if (_super_index->getNbOfElems())
    {
      oss << "   Super-indices:" << std::endl;
      oss << "   ";
      const mcIdType * i = _super_index->begin();
      for ( ; i != _super_index->end(); ++i )
        oss << *i << " ";
      oss << std::endl;
    }

  oss << "   Indices:" << std::endl;
  oss << "   ";
  const mcIdType * i = _index->begin();
  for ( ; i != _index->end(); ++i )
    oss << *i << " ";
  oss << std::endl;
  oss << "   Values:" << std::endl;
  oss << "     ";
  const mcIdType * v = _values->begin();
  mcIdType cnt = 0, cntI = 0;
  i = _index->begin();
  for ( const mcIdType * si = super_index->begin()+1; v != _values->end(); ++v, ++cnt )
    {
      if ( cnt == *i )
        {
          if ( cntI == *si && cnt != 0)
            {
              oss << std::endl << "     ";
              ++si;
            }

          oss << "| ";
          ++i; ++cntI;
        }
      oss << *v << " ";
    }
  oss << std::endl;

  return oss.str();
}

/*!
 * Returns 2 SkyLineArrays with same number of packs than \a this.
 * Each pack in \a this is split in 2 parts using \a threshold parameter as cut point.
 * \a left part contains ids in \a this pack strictly lower than \a threshold
 * \a right part contains ids in \a this pack greater or equal to \a threshold
 */
void MEDCouplingSkyLineArray::thresholdPerPack(mcIdType threshold, MCAuto<MEDCouplingSkyLineArray>& left, MCAuto<MEDCouplingSkyLineArray>& right) const
{
  mcIdType nbPacks(this->getNumberOf());
  MCAuto<DataArrayIdType> lCount(DataArrayIdType::New()); lCount->alloc(nbPacks,1); lCount->fillWithZero();
  mcIdType *lCountPtr(lCount->getPointerSilent());
  const mcIdType *valuesPtr(this->_values->begin()),*indexPtr(this->_index->begin());
  for(mcIdType i = 0 ; i < nbPacks ; ++i, ++lCountPtr)
  {
    *lCountPtr = ToIdType(std::count_if(valuesPtr+indexPtr[i],valuesPtr+indexPtr[i+1],[threshold](mcIdType elt) { return elt<threshold; }));
  }
  MCAuto<DataArrayIdType> sizeOfPacks(this->_index->deltaShiftIndex());
  sizeOfPacks->substractEqual(lCount);
  mcIdType leftNbOfVal(lCount->accumulate(std::size_t(0))),rightNbOfVal(sizeOfPacks->accumulate(std::size_t(0)));
  lCount->computeOffsetsFull(); sizeOfPacks->computeOffsetsFull();
  MCAuto<DataArrayIdType> leftValues(DataArrayIdType::New()); leftValues->alloc(leftNbOfVal,1);
  MCAuto<DataArrayIdType> rightValues(DataArrayIdType::New()); rightValues->alloc(rightNbOfVal,1);
  mcIdType *rvPtr(rightValues->getPointerSilent()),*lvPtr(leftValues->getPointerSilent());
  for(mcIdType i = 0 ; i < nbPacks ; ++i)
  {
    std::for_each(valuesPtr+indexPtr[i],valuesPtr+indexPtr[i+1],[threshold,&rvPtr,&lvPtr](mcIdType elt) { if(elt<threshold) { *lvPtr++ = elt; } else { *rvPtr++ = elt; } });
  }
  left = MEDCouplingSkyLineArray::New(lCount,leftValues); right = MEDCouplingSkyLineArray::New(sizeOfPacks,rightValues);
}

MEDCouplingSkyLineArray *MEDCouplingSkyLineArray::groupPacks(const DataArrayIdType *indexedPacks) const
{
  indexedPacks->checkAllocated();
  if( indexedPacks->getNumberOfComponents() != 1 )
    throw INTERP_KERNEL::Exception("MEDCouplingSkyLineArray::groupPacks : number of components must be 1 !");
  std::size_t nbTuples(indexedPacks->getNumberOfTuples());
  if( nbTuples == 0 )
    throw INTERP_KERNEL::Exception("MEDCouplingSkyLineArray::groupPacks : number of tuples must be > 0 !");
  const DataArrayIdType *index(this->getIndexArray());
  MCAuto<DataArrayIdType> partIndex(index->selectByTupleIdSafe(indexedPacks->begin(),indexedPacks->end()));
  MCAuto<MEDCouplingSkyLineArray> ret(MEDCouplingSkyLineArray::New(partIndex,this->getValuesArray()));
  return ret.retn();
}

MEDCouplingSkyLineArray *MEDCouplingSkyLineArray::uniqueNotSortedByPack() const
{
  mcIdType nbPacks(this->getNumberOf());
  MCAuto<DataArrayIdType> retIndex(DataArrayIdType::New()); retIndex->alloc(nbPacks+1,1);
  const mcIdType *valuesPtr(this->_values->begin()),*indexPtr(this->_index->begin());
  mcIdType *retIndexPtr(retIndex->getPointer()); *retIndexPtr = 0;
  for(mcIdType i = 0 ; i < nbPacks ; ++i, ++retIndexPtr)
  {
    std::set<mcIdType> s(valuesPtr+indexPtr[i],valuesPtr+indexPtr[i+1]);
    retIndexPtr[1] = retIndexPtr[0] + ToIdType(s.size());
  }
  MCAuto<DataArrayIdType> retValues(DataArrayIdType::New()); retValues->alloc(retIndex->back(),1);
  mcIdType *retValuesPtr(retValues->getPointer());
  for(mcIdType i = 0 ; i < nbPacks ; ++i)
  {
    std::set<mcIdType> s(valuesPtr+indexPtr[i],valuesPtr+indexPtr[i+1]);
    retValuesPtr = std::copy(s.begin(),s.end(),retValuesPtr);
  }
  MCAuto<MEDCouplingSkyLineArray> ret(MEDCouplingSkyLineArray::New(retIndex,retValues));
  return ret.retn();
}
/*!
 * Take as input skylinearrays containing the same number of packs ( \a this->getNumberOf ).
 * For each packs, this method aggregates corresponding pack in \a sks.
 * 
 * \throw if either a lenght of not nullptr instances in \a sks is zero or if number of packs of not nullptr instances in \a sks is not the same.
 * \return a newly allocated skyline array that is the aggregated packs of inputs
 */
MEDCouplingSkyLineArray *MEDCouplingSkyLineArray::AggregatePacks(const std::vector<const MEDCouplingSkyLineArray *>& sks)
{
  std::vector<const MEDCouplingSkyLineArray *>sksEff;
  mcIdType nbOfPacks(std::numeric_limits<mcIdType>::max());
  constexpr char MSG[]="MEDCouplingSkyLineArray::AggregatePacks : ";
  for(auto sk : sks)
  {
    if(sk)
    {
      mcIdType curNbPacks(sk->getNumberOf());
      if(sksEff.empty())
        nbOfPacks = curNbPacks;
      if(nbOfPacks != curNbPacks)
      {
        std::ostringstream oss; oss << MSG << "first not null input ska has " << nbOfPacks << " whereas there is presence of ska with " << curNbPacks << " !";
        throw INTERP_KERNEL::Exception(oss.str());
      }
      sksEff.push_back(sk);
    }
  }
  if(sksEff.empty())
  {
    std::ostringstream oss; oss << MSG << "input vector contains no not nullptr elements !";
    throw INTERP_KERNEL::Exception(oss.str());
  }
  //
  MCAuto<DataArrayIdType> index(DataArrayIdType::New()); index->alloc(nbOfPacks+1,1);
  mcIdType *indexPtr(index->getPointer()); *indexPtr=0;
  std::vector<const mcIdType *> indicesIn(SkyLineArrayIndexIterator(0,&sksEff),SkyLineArrayIndexIterator(sksEff.size(),&sksEff));
  for( mcIdType packId = 0 ; packId < nbOfPacks ; ++packId, ++indexPtr )
  {
    mcIdType nbOfAggPacks(0);
    std::for_each(indicesIn.begin(),indicesIn.end(),[packId,&nbOfAggPacks](const mcIdType *elt) { nbOfAggPacks+=elt[packId+1]-elt[packId]; });
    indexPtr[1] = indexPtr[0] + nbOfAggPacks;
  }
  mcIdType nbOfTuplesOut(index->back());
  MCAuto<DataArrayIdType> values(DataArrayIdType::New()); values->alloc(nbOfTuplesOut,1);
  mcIdType *valuesPtr(values->getPointer());
  // let's go to populate values array
  std::vector<const mcIdType *> valuesIn(SkyLineArrayValuesIterator(0,&sksEff),SkyLineArrayValuesIterator(sksEff.size(),&sksEff));
  for( mcIdType packId = 0 ; packId < nbOfPacks ; ++packId )
  {
    std::size_t pos(0);
    std::for_each(valuesIn.begin(),valuesIn.end(),[packId,&indicesIn,&valuesPtr,&pos](const mcIdType *elt)
    { valuesPtr=std::copy(elt+indicesIn[pos][packId],elt+indicesIn[pos][packId+1],valuesPtr); ++pos; }
    );
  }
  //
  MCAuto<MEDCouplingSkyLineArray> ret(MEDCouplingSkyLineArray::New(index,values));
  return ret.retn();
}

/**
 * For a 2- or 3-level SkyLine array, return a copy of the absolute pack with given identifier.
 */
void MEDCouplingSkyLineArray::getSimplePackSafe(const mcIdType absolutePackId, std::vector<mcIdType> & pack) const
{
  if(absolutePackId < 0 || absolutePackId >= _index->getNbOfElems())
    throw INTERP_KERNEL::Exception("MEDCouplingSkyLineArray::getPackSafe: invalid index!");
  const mcIdType * iP(_index->begin()), *vP(_values->begin());
  mcIdType sz = iP[absolutePackId+1]-iP[absolutePackId];
  pack.resize(sz);
  std::copy(vP+iP[absolutePackId], vP+iP[absolutePackId+1],pack.begin());
}

/**
 * Same as getPackSafe, but directly returns a pointer to the internal data with the size of the pack.
 */
const mcIdType * MEDCouplingSkyLineArray::getSimplePackSafePtr(const mcIdType absolutePackId, mcIdType & packSize) const
{
  if(absolutePackId < 0 || absolutePackId >= _index->getNbOfElems())
    throw INTERP_KERNEL::Exception("MEDCouplingSkyLineArray::getPackSafe: invalid index!");
  const mcIdType * iP(_index->begin()), *vP(_values->begin());
  packSize = iP[absolutePackId+1]-iP[absolutePackId];
  return vP+iP[absolutePackId];
}


/**!
 * For each given super-pack ID, provide the sub-index of the first matching pack. If no matching pack is found for the
 * given super-pack -1 is returned.
 * \param[in] superPackIndices the list of super-packs that should be inspected
 * \param[in] packBg the pack that the function is looking for in each of the provided super-pack
 * \param[in] packEnd the pack that the function is looking for in each of the provided super-pack
 * \param[out] a vector of mcIdType, having the same size as superPackIndices and containing for each inspected super-pack
 * the index of the first matching pack, or -1 if none found.
 */
void MEDCouplingSkyLineArray::findPackIds(const std::vector<mcIdType> & superPackIndices,
                                          const mcIdType *packBg, const mcIdType *packEnd,
                                          std::vector<mcIdType>& out) const
{
  using namespace std;

  checkSuperIndex("findPackIds");

  mcIdType packSz = ToIdType(std::distance(packBg, packEnd));
  if (!packSz)
    throw INTERP_KERNEL::Exception("MEDCouplingSkyLineArray::findPackIds: void pack!");

  out.resize(superPackIndices.size());
  mcIdType i = 0;
  const mcIdType * siP(_super_index->begin()), * iP(_index->begin()), *vP(_values->begin());
  for(vector<mcIdType>::const_iterator it=superPackIndices.begin(); it!=superPackIndices.end(); ++it, i++)
    {
      out[i] = -1;
      const mcIdType sPackIdx = *it;
      // for each pack
      for (mcIdType idx=siP[sPackIdx], j=0; idx < siP[sPackIdx+1]; idx++, j++)
        {
          if (packSz == (iP[idx+1] - iP[idx]))
            if (equal(&vP[iP[idx]], &vP[iP[idx+1]], packBg))
              {
                out[i] = j;
                break;
              }
        }
    }
}

/**!
 * Delete pack number 'idx' in super-pack number 'superIdx'.
 * \param[in] superIdx is the super-pack number
 * \param[in] idx is the pack index inside the super-pack 'superIdx'.
 */
void MEDCouplingSkyLineArray::deletePack(const mcIdType superIdx, const mcIdType idx)
{
  checkSuperIndex("deletePack");
  validSuperIndexAndIndex("deletePack", superIdx, idx);

  mcIdType * vP = _values->getPointer();
  mcIdType * siP(_super_index->getPointer()), *iP(_index->getPointer());
  const mcIdType start = iP[siP[superIdx]+idx], end = iP[siP[superIdx]+idx+1];
  // _values
  std::copy(vP+end, vP+_values->getNbOfElems(), vP+start);
  _values->reAlloc(_values->getNbOfElems() - (end-start));

  // _index
  mcIdType nt = _index->getNbOfElems();
  std::copy(iP+siP[superIdx]+idx+1, iP+nt, iP+siP[superIdx]+idx);
  _index->reAlloc(nt-1); iP = _index->getPointer();  // better not forget this ...
  for(mcIdType ii = siP[superIdx]+idx; ii < nt-1; ii++)
    iP[ii] -= (end-start);

  // _super_index
  for(mcIdType ii = superIdx+1; ii < _super_index->getNbOfElems(); ii++)
    (siP[ii])--;
}

void MEDCouplingSkyLineArray::deleteSimplePack(const mcIdType idx)
{
  validIndex("deleteSimplePack", idx);
  
  mcIdType* iP(_index->getPointer());
  const mcIdType start(iP[idx]), end(iP[idx+1]);

  // _values
  mcIdType initValSz=_values->getNbOfElems();
  mcIdType deltaSz( start-end );  // should be negative
  mcIdType *vP(_values->getPointer());
  if (deltaSz < 0)
    {
      std::copy(vP+end, vP+initValSz, vP+start);
      _values->reAlloc(initValSz+deltaSz);
    }
  else
    throw INTERP_KERNEL::Exception("MEDCouplingSkyLineArray::deleteSimplePack");
  // _index
  mcIdType nt=_index->getNbOfElems();
  std::copy(iP+idx+1, iP+nt, iP+idx);
  for(mcIdType ii = idx; ii < nt-1; ii++)
    iP[ii] += deltaSz;
  _index->reAlloc(nt-1);
}

void MEDCouplingSkyLineArray::replaceSimplePacks(const DataArrayIdType* idx, const std::vector<const DataArrayIdType*>& packs)
{    
  if (idx->empty())
    return;
    
  for (const mcIdType * id = idx->begin(); id != idx->end(); id++)
    validIndex("deleteSimplePacks", *id);
    
  if (idx->getNbOfElems() != ToIdType( packs.size()))
    throw INTERP_KERNEL::Exception("MEDCouplingSkyLineArray::deleteSimplePacks: size of list of pack is incorrect");
    
  // copy _index, _values into a deque<set<mcIdType>>
  std::deque< std::set<mcIdType> > valuesByIdx;
  mcIdType* vP(_values->getPointer());
  mcIdType* iP(_index->getPointer());
  mcIdType nt = _index->getNbOfElems();
  for (mcIdType ii = 0; ii < nt-1; ii++)
    valuesByIdx.push_back(std::set<mcIdType>(vP+iP[ii], vP+iP[ii+1]));
    
  // modify the deque<set<mcIdType>> according to idx and packs
  mcIdType ii(0);
  for (const mcIdType *id = idx->begin(); id != idx->end(); id++)
    {
      valuesByIdx[*id] = std::set<mcIdType>(packs[ii]->begin(), packs[ii]->end());
      ii++;
    }
  // copy back the deque<set<mcIdType>> into _index, _values
  mcIdType valSz(0);
  *iP = 0;
  for (std::deque< std::set<mcIdType> >::const_iterator values=valuesByIdx.begin();values!=valuesByIdx.end();values++)
    {
      valSz += ToIdType((*values).size());
      *(++iP) = valSz;
    }
  _values->reAlloc(valSz);
  iP = _index->getPointer();
  vP = _values->getPointer();
  for (auto values : valuesByIdx)
    {
      std::copy(values.begin(), values.end(), vP+(*iP));
      iP++;
    }
}

void MEDCouplingSkyLineArray::deleteSimplePacks(const DataArrayIdType* idx)
{    
  for (auto id = idx->begin(); id != idx->end(); id++)
    validIndex("deleteSimplePacks", *id);
  
  std::set<mcIdType> packsToDelete(idx->begin(), idx->end());
    
  // _values
  mcIdType* iP(_index->getPointer());
  mcIdType initValSz = _values->getNbOfElems();
  mcIdType *vP(_values->getPointer());
  mcIdType end_prec(0),start_prec(0);
  for(std::set<mcIdType>::const_iterator ii=packsToDelete.begin();ii!=packsToDelete.end();ii++)
    {
      mcIdType start = iP[*ii];
      if (end_prec != 0)
        std::copy(vP+end_prec, vP+start, vP+start_prec);
      start_prec += start-end_prec;
      end_prec = iP[*ii+1];
    }
  if (end_prec != 0)
    std::copy(vP+end_prec, vP+initValSz, vP+start_prec);
  _values->reAlloc(initValSz-(end_prec-start_prec));
    
  // _index
  mcIdType nt = _index->getNbOfElems();
  mcIdType offset = 0;
  end_prec = 0;
  start_prec = 0;
  mcIdType deleted = 0;
  for(std::set<mcIdType>::const_iterator ii=packsToDelete.begin();ii!=packsToDelete.end();ii++)
    {
      if (end_prec != 0)
        {
          std::copy(iP+end_prec, iP+*ii, iP+start_prec);
          for (mcIdType i=start_prec; i<*ii; i++)
            iP[i] -= offset;
        }
      offset += iP[*ii+1] - iP[*ii];
      start_prec = *ii-deleted;
      end_prec = *ii+1;
      deleted += 1;
    }
  if (end_prec != 0)
    {
      std::copy(iP+end_prec, iP+nt, iP+start_prec);
      for (mcIdType i=start_prec; i<nt; i++)
        iP[i] -= offset;
    }
  _index->reAlloc(nt-deleted);
}

/**!
 * Insert a new pack in super-pack at index 'superIdx'. The pack is inserted at the end of the pack list of the chosen super-pack.
 */
void MEDCouplingSkyLineArray::pushBackPack(const mcIdType superIdx, const mcIdType * packBg, const mcIdType * packEnd)
{
  using namespace std;

  checkSuperIndex("pushBackPack");
  validSuperIndex("pushBackPack", superIdx);

  mcIdType *siP(_super_index->getPointer()), *iP(_index->getPointer());
  const mcIdType sz(ToIdType(distance(packBg, packEnd)));

  // _values
  _values->reAlloc(_values->getNbOfElems()+sz);
  mcIdType * vPE(_values->getPointer()+_values->getNbOfElems());
  mcIdType *vP(_values->getPointer());
  copy(vP+iP[siP[superIdx+1]], vPE-sz, vP+iP[siP[superIdx+1]]+sz);
  // insert pack
  copy(packBg, packEnd, vP+iP[siP[superIdx+1]]);

  // _index
  mcIdType nt = ToIdType(_index->getNbOfElems());
  _index->reAlloc(nt+1); iP = _index->getPointer();
  copy(iP+siP[superIdx+1]+1, iP+nt, iP+siP[superIdx+1]+2);
  iP[siP[superIdx+1]+1] = iP[siP[superIdx+1]] + sz;
  for(mcIdType ii = siP[superIdx+1]+2; ii < nt+1; ii++)
    iP[ii] += sz;

  // _super_index
  for(mcIdType ii = superIdx+1; ii < _super_index->getNbOfElems(); ii++)
    (siP[ii])++;
}

/**
 * Replace pack with absolute index 'idx' with the provided new pack. Function can be used either
 * for 2-level SkyLine or 3-level SkyLine.
 */
void MEDCouplingSkyLineArray::replaceSimplePack(const mcIdType idx, const mcIdType * packBg, const mcIdType * packEnd)
{
  validIndex("replaceSimplePack", idx);

  mcIdType * iP(_index->getPointer());
  mcIdType newSz = ToIdType(std::distance(packBg, packEnd));
  const mcIdType start = iP[idx], end = iP[idx+1];

  // _values
  mcIdType initValSz = _values->getNbOfElems();
  mcIdType deltaSz = newSz-(end-start);  // can be negative
  if (deltaSz)
    {
      if (deltaSz > 0)
        _values->reAlloc(initValSz+deltaSz);
      mcIdType *vP(_values->getPointer());
      std::copy(vP+end, vP+initValSz, vP+end+deltaSz);
      if (deltaSz < 0)
        _values->reAlloc(initValSz+deltaSz);
    }

  // copy new pack
  std::copy(packBg, packEnd, _values->getPointer()+start);

  // _index
  for(mcIdType ii = idx+1; ii < _index->getNbOfElems(); ii++)
    iP[ii] += deltaSz;
}

/**
 * Replace pack with super index 'superIdx' and index 'idx' with the provided new pack.
 * Function can be used only for 3-level SkyLine.
 */
void MEDCouplingSkyLineArray::replacePack(const mcIdType superIdx, const mcIdType idx, const mcIdType *packBg, const mcIdType *packEnd)
{
  checkSuperIndex("replacePack");
  validSuperIndexAndIndex("replacePack", superIdx, idx);

  mcIdType * siP(_super_index->getPointer()), *iP(_index->getPointer());
  mcIdType newSz = ToIdType(std::distance(packBg, packEnd));
  const mcIdType start = iP[siP[superIdx]+idx], end = iP[siP[superIdx]+idx+1];

  // _values
  mcIdType initValSz = _values->getNbOfElems();
  mcIdType deltaSz = newSz-(end-start);  // can be negative
  if (deltaSz)
    {
      if (deltaSz > 0)
        _values->reAlloc(initValSz+deltaSz);
      mcIdType *vP(_values->getPointer());
      std::copy(vP+end, vP+initValSz, vP+end+deltaSz);
      if (deltaSz < 0)
        _values->reAlloc(initValSz+deltaSz);
    }

  // copy new pack
  std::copy(packBg, packEnd, _values->getPointer()+start);

  // _index
  for(mcIdType ii = siP[superIdx]+idx+1; ii < _index->getNbOfElems(); ii++)
    iP[ii] += deltaSz;
}
