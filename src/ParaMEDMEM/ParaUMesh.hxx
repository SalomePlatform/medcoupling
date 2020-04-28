// Copyright (C) 2020  CEA/DEN, EDF R&D
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

#pragma once

#include "MEDCouplingUMesh.hxx"
#include "ProcessorGroup.hxx"
#include "MEDCouplingMemArray.hxx"

#include <string>
#include <vector>

namespace MEDCoupling
{
  /*!
   * Parallel representation of an unstructured mesh.
   *
   * This class is very specific to the requirement of parallel code computations.
   */
  class ParaUMesh : public RefCountObject
  {
  public:
    static ParaUMesh *New(MEDCouplingUMesh *mesh, DataArrayIdType *globalCellIds, DataArrayIdType *globalNodeIds);
    MCAuto<DataArrayIdType> getCellIdsLyingOnNodes(const DataArrayIdType *globalNodeIds, bool fullyIn) const;
    ParaUMesh *redistributeCells(const DataArrayIdType *globalCellIds) const;
    DataArrayDouble *redistributeCellField(const DataArrayIdType *globalCellIds, const DataArrayDouble *fieldValueToRed) const;
    DataArrayIdType *redistributeCellField(const DataArrayIdType *globalCellIds, const DataArrayIdType *fieldValueToRed) const;
    DataArrayDouble *redistributeNodeField(const DataArrayIdType *globalCellIds, const DataArrayDouble *fieldValueToRed) const;
    DataArrayIdType *redistributeNodeField(const DataArrayIdType *globalCellIds, const DataArrayIdType *fieldValueToRed) const;
    MEDCouplingUMesh *getMesh() { return _mesh; }
    DataArrayIdType *getGlobalCellIds() { return _cell_global; }
    DataArrayIdType *getGlobalNodeIds() { return _node_global; }
  protected:
    virtual ~ParaUMesh() { }
    ParaUMesh(MEDCouplingUMesh *mesh, DataArrayIdType *globalCellIds, DataArrayIdType *globalNodeIds);
    std::string getClassName() const override { return "ParaUMesh"; }
    std::size_t getHeapMemorySizeWithoutChildren() const override;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const override;
  private:
    MCAuto<MEDCouplingUMesh> _mesh;
    MCAuto<DataArrayIdType> _cell_global;
    MCAuto<DataArrayIdType> _node_global;
  private:
    template<class T>
    typename Traits<T>::ArrayType *redistributeCellFieldT(const DataArrayIdType *globalCellIds, const typename Traits<T>::ArrayType *fieldValueToRed) const
    {
      using DataArrayT = typename Traits<T>::ArrayType;
      MPI_Comm comm(MPI_COMM_WORLD);
      CommInterface ci;
      if( _cell_global->getNumberOfTuples() != fieldValueToRed->getNumberOfTuples() )
        throw INTERP_KERNEL::Exception("PAraUMesh::redistributeCellFieldT : invalid input length of array !");
      std::unique_ptr<mcIdType[]> allGlobalCellIds,allGlobalCellIdsIndex;
      int size(ci.allGatherArrays(comm,globalCellIds,allGlobalCellIds,allGlobalCellIdsIndex));
      // Prepare ParaUMesh parts to be sent : compute for each proc the contribution of current rank.
      std::vector< MCAuto<DataArrayIdType> > globalCellIdsToBeSent(size);
      std::vector< MCAuto<DataArrayT> > fieldToBeSent(size);
      for(int curRk = 0 ; curRk < size ; ++curRk)
      {
        mcIdType offset(allGlobalCellIdsIndex[curRk]);
        MCAuto<DataArrayIdType> globalCellIdsOfCurProc(DataArrayIdType::New());
        globalCellIdsOfCurProc->useArray(allGlobalCellIds.get()+offset,false,DeallocType::CPP_DEALLOC,allGlobalCellIdsIndex[curRk+1]-offset,1);
        // the key call is here : compute for rank curRk the cells to be sent
        MCAuto<DataArrayIdType> globalCellIdsCaptured(_cell_global->buildIntersection(globalCellIdsOfCurProc));// OK for the global cellIds
        MCAuto<DataArrayIdType> localCellIdsCaptured(_cell_global->findIdForEach(globalCellIdsCaptured->begin(),globalCellIdsCaptured->end()));
        globalCellIdsToBeSent[curRk] = globalCellIdsCaptured;
        fieldToBeSent[curRk] = fieldValueToRed->selectByTupleIdSafe(localCellIdsCaptured->begin(),localCellIdsCaptured->end());
      }
      // Receive
      std::vector< MCAuto<DataArrayIdType> > globalCellIdsReceived;
      ci.allToAllArrays(comm,globalCellIdsToBeSent,globalCellIdsReceived);
      std::vector< MCAuto<DataArrayT> > fieldValueReceived;
      ci.allToAllArrays(comm,fieldToBeSent,fieldValueReceived);
      // use globalCellIdsReceived to reorganize everything
      MCAuto<DataArrayIdType> aggregatedCellIds( DataArrayIdType::Aggregate(FromVecAutoToVecOfConst<DataArrayIdType>(globalCellIdsReceived)) );
      MCAuto<DataArrayIdType> aggregatedCellIdsSorted(aggregatedCellIds->copySorted());
      MCAuto<DataArrayIdType> idsIntoAggregatedIds(DataArrayIdType::FindPermutationFromFirstToSecondDuplicate(aggregatedCellIdsSorted,aggregatedCellIds));
      MCAuto<DataArrayIdType> cellIdsOfSameNodeIds(aggregatedCellIdsSorted->indexOfSameConsecutiveValueGroups());
      MCAuto<DataArrayIdType> n2o_cells(idsIntoAggregatedIds->selectByTupleIdSafe(cellIdsOfSameNodeIds->begin(),cellIdsOfSameNodeIds->end()-1));//new == new ordering so that global cell ids are sorted . old == coarse ordering implied by the aggregation
      //
      MCAuto<DataArrayT> fieldAggregated(DataArrayT::Aggregate(FromVecAutoToVecOfConst<DataArrayT>(fieldValueReceived)));
      MCAuto<DataArrayT> ret(fieldAggregated->selectByTupleIdSafe(n2o_cells->begin(),n2o_cells->end()));
      return ret.retn();
    }
    
    template<class T>
    typename Traits<T>::ArrayType *redistributeNodeFieldT(const DataArrayIdType *globalCellIds, const typename Traits<T>::ArrayType *fieldValueToRed) const
    {
      using DataArrayT = typename Traits<T>::ArrayType;
      MPI_Comm comm(MPI_COMM_WORLD);
      CommInterface ci;
      if( _node_global->getNumberOfTuples() != fieldValueToRed->getNumberOfTuples() )
        throw INTERP_KERNEL::Exception("PAraUMesh::redistributeNodeFieldT : invalid input length of array !");
      std::unique_ptr<mcIdType[]> allGlobalCellIds,allGlobalCellIdsIndex;
      int size(ci.allGatherArrays(comm,globalCellIds,allGlobalCellIds,allGlobalCellIdsIndex));
      // Prepare ParaUMesh parts to be sent : compute for each proc the contribution of current rank.
      std::vector< MCAuto<DataArrayIdType> > globalNodeIdsToBeSent(size);
      std::vector< MCAuto<DataArrayT> > fieldToBeSent(size);
      for(int curRk = 0 ; curRk < size ; ++curRk)
      {
        mcIdType offset(allGlobalCellIdsIndex[curRk]);
        MCAuto<DataArrayIdType> globalCellIdsOfCurProc(DataArrayIdType::New());
        globalCellIdsOfCurProc->useArray(allGlobalCellIds.get()+offset,false,DeallocType::CPP_DEALLOC,allGlobalCellIdsIndex[curRk+1]-offset,1);
        // the key call is here : compute for rank curRk the cells to be sent
        MCAuto<DataArrayIdType> globalCellIdsCaptured(_cell_global->buildIntersection(globalCellIdsOfCurProc));// OK for the global cellIds
        MCAuto<DataArrayIdType> localCellIdsCaptured(_cell_global->findIdForEach(globalCellIdsCaptured->begin(),globalCellIdsCaptured->end()));
        MCAuto<MEDCouplingUMesh> meshPart(_mesh->buildPartOfMySelf(localCellIdsCaptured->begin(),localCellIdsCaptured->end(),true));
        MCAuto<DataArrayIdType> o2n(meshPart->zipCoordsTraducer());// OK for the mesh
        MCAuto<DataArrayIdType> n2o(o2n->invertArrayO2N2N2O(meshPart->getNumberOfNodes()));
        MCAuto<DataArrayIdType> globalNodeIdsPart(_node_global->selectByTupleIdSafe(n2o->begin(),n2o->end())); // OK for the global nodeIds
        globalNodeIdsToBeSent[curRk] = globalNodeIdsPart;
        fieldToBeSent[curRk] = fieldValueToRed->selectByTupleIdSafe(n2o->begin(),n2o->end());
      }
      // Receive
      std::vector< MCAuto<DataArrayIdType> > globalNodeIdsReceived;
      ci.allToAllArrays(comm,globalNodeIdsToBeSent,globalNodeIdsReceived);
      std::vector< MCAuto<DataArrayT> > fieldValueReceived;
      ci.allToAllArrays(comm,fieldToBeSent,fieldValueReceived);
      // firstly deal with nodes.
      MCAuto<DataArrayIdType> aggregatedNodeIds( DataArrayIdType::Aggregate(FromVecAutoToVecOfConst<DataArrayIdType>(globalNodeIdsReceived)) );
      MCAuto<DataArrayIdType> aggregatedNodeIdsSorted(aggregatedNodeIds->copySorted());
      MCAuto<DataArrayIdType> nodeIdsIntoAggregatedIds(DataArrayIdType::FindPermutationFromFirstToSecondDuplicate(aggregatedNodeIdsSorted,aggregatedNodeIds));
      MCAuto<DataArrayIdType> idxOfSameNodeIds(aggregatedNodeIdsSorted->indexOfSameConsecutiveValueGroups());
      MCAuto<DataArrayIdType> n2o_nodes(nodeIdsIntoAggregatedIds->selectByTupleIdSafe(idxOfSameNodeIds->begin(),idxOfSameNodeIds->end()-1));//new == new ordering so that global node ids are sorted . old == coarse ordering implied by the aggregation
      //
      MCAuto<DataArrayT> fieldAggregated(DataArrayT::Aggregate(FromVecAutoToVecOfConst<DataArrayT>(fieldValueReceived)));
      MCAuto<DataArrayT> ret(fieldAggregated->selectByTupleIdSafe(n2o_nodes->begin(),n2o_nodes->end()));
      //
      return ret.retn();
    }
  };
}
