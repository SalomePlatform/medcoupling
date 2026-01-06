//
// Copyright (C) 2020-2026  CEA, EDF
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

#include "ParaSkyLineArray.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "CommInterface.hxx"
#include "MEDCouplingMemArray.hxx"

#include "mpi.h"

#include <fstream>
#include <sstream>
#include <numeric>
#include <memory>
#include <vector>

using namespace MEDCoupling;

ParaSkyLineArray *
ParaSkyLineArray::New(MEDCouplingSkyLineArray *ska, DataArrayIdType *globalIds)
{
    return new ParaSkyLineArray(ska, globalIds);
}

MEDCouplingSkyLineArray *
ParaSkyLineArray::getSkyLineArray() const
{
    return this->_ska.iAmATrollConstCast();
}

DataArrayIdType *
ParaSkyLineArray::getGlobalIdsArray() const
{
    return this->_global_ids.iAmATrollConstCast();
}

ParaSkyLineArray::ParaSkyLineArray(MEDCouplingSkyLineArray *ska, DataArrayIdType *globalIds)
{
    _ska.takeRef(ska);
    _global_ids.takeRef(globalIds);
    _ska.checkNotNull();
    _global_ids.checkNotNull();
    if (_ska->getNumberOf() != _global_ids->getNumberOfTuples())
    {
        std::ostringstream oss;
        oss << "ParaSkyLineArray constructor : mismatch between # globalIds (" << _global_ids->getNumberOfTuples()
            << ") and len of indices in SkyLineArray (" << _ska->getNumberOf() << ").";
        throw INTERP_KERNEL::Exception(oss.str());
    }
}

std::size_t
ParaSkyLineArray::getHeapMemorySizeWithoutChildren() const
{
    return 0;
}

std::vector<const BigMemoryObject *>
ParaSkyLineArray::getDirectChildrenWithNull() const
{
    return {_ska, _global_ids};
}

MCAuto<ParaSkyLineArray>
ParaSkyLineArray::equiRedistribute(mcIdType nbOfEntities) const
{
    MPI_Comm comm(MPI_COMM_WORLD);
    CommInterface ci;
    int size;
    ci.commSize(comm, &size);
    std::vector<MCAuto<MEDCouplingSkyLineArray> > skToBeSent(size);
    std::vector<MCAuto<DataArrayIdType> > idsCaptured(size);
    for (int curRk = 0; curRk < size; ++curRk)
    {
        mcIdType curStart(0), curEnd(0);
        DataArrayIdType::GetSlice(0, nbOfEntities, 1, curRk, size, curStart, curEnd);
        MCAuto<DataArrayIdType> idsInGlobalIds(_global_ids->findIdsInRange(curStart, curEnd));
        idsCaptured[curRk] = _global_ids->selectByTupleIdSafe(idsInGlobalIds->begin(), idsInGlobalIds->end());
        {
            DataArrayIdType *tmpValues(nullptr), *tmpIndex(nullptr);
            DataArrayIdType::ExtractFromIndexedArrays(
                idsInGlobalIds->begin(),
                idsInGlobalIds->end(),
                this->_ska->getValuesArray(),
                this->_ska->getIndexArray(),
                tmpValues,
                tmpIndex
            );
            MCAuto<DataArrayIdType> tmpValues2(tmpValues), tmpIndex2(tmpIndex);
            skToBeSent[curRk] = MEDCouplingSkyLineArray::New(tmpIndex, tmpValues);
        }
    }
    // communication : 3 arrays are going to be all2allized : ids, SkyLine index and SyLine values
    MCAuto<DataArrayIdType> aggregatedIds, indices, values;
    {
        std::vector<MCAuto<DataArrayIdType> > myRkIdsCaptured;
        ci.allToAllArrays(comm, idsCaptured, myRkIdsCaptured);
        aggregatedIds = DataArrayIdType::Aggregate(FromVecAutoToVecOfConst<DataArrayIdType>(myRkIdsCaptured));
    }
    {
        std::vector<MCAuto<DataArrayIdType> > myRkSkIndex;
        std::vector<MCAuto<DataArrayIdType> > indexToBeSent(MEDCouplingSkyLineArray::RetrieveVecIndex(skToBeSent));
        ci.allToAllArrays(comm, indexToBeSent, myRkSkIndex);
        indices = DataArrayIdType::AggregateIndexes(FromVecAutoToVecOfConst<DataArrayIdType>(myRkSkIndex));
    }
    {
        std::vector<MCAuto<DataArrayIdType> > myRkSkValues;
        std::vector<MCAuto<DataArrayIdType> > valuesToBeSent(MEDCouplingSkyLineArray::RetrieveVecValues(skToBeSent));
        ci.allToAllArrays(comm, valuesToBeSent, myRkSkValues);
        values = DataArrayIdType::Aggregate(FromVecAutoToVecOfConst<DataArrayIdType>(myRkSkValues));
    }
    // Reorder results coming from other procs
    MCAuto<DataArrayIdType> aggregatedIdsSort(aggregatedIds->deepCopy());
    aggregatedIdsSort->sort();
    MCAuto<DataArrayIdType> idsIntoAggregatedIds(
        DataArrayIdType::FindPermutationFromFirstToSecondDuplicate(aggregatedIdsSort, aggregatedIds)
    );
    MCAuto<DataArrayIdType> indicesSorted, valuesSorted;
    {
        DataArrayIdType *indicesSortedTmp(nullptr), *valuesSortedTmp(nullptr);
        DataArrayIdType::ExtractFromIndexedArrays(
            idsIntoAggregatedIds->begin(),
            idsIntoAggregatedIds->end(),
            values,
            indices,
            valuesSortedTmp,
            indicesSortedTmp
        );
        indicesSorted = indicesSortedTmp;
        valuesSorted = valuesSortedTmp;
    }
    MCAuto<DataArrayIdType> idxOfSameIds(aggregatedIdsSort->indexOfSameConsecutiveValueGroups());
    //
    MCAuto<DataArrayIdType> globalIdsOut(aggregatedIdsSort->buildUnique());
    MCAuto<MEDCouplingSkyLineArray> skOut(MEDCouplingSkyLineArray::New(indicesSorted, valuesSorted));
    skOut = skOut->groupPacks(idxOfSameIds);  // group partial packs coming from different procs
    skOut = skOut->uniqueNotSortedByPack();   // remove duplicates
    MCAuto<ParaSkyLineArray> ret(ParaSkyLineArray::New(skOut, globalIdsOut));
    return ret.retn();
}
