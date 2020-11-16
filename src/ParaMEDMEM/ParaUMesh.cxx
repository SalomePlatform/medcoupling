//
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

#include "ParaUMesh.hxx"
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

ParaUMesh *ParaUMesh::New(MEDCouplingUMesh *mesh, DataArrayIdType *globalCellIds, DataArrayIdType *globalNodeIds)
{
  return new ParaUMesh(mesh,globalCellIds,globalNodeIds);
}

ParaUMesh::ParaUMesh(MEDCouplingUMesh *mesh, DataArrayIdType *globalCellIds, DataArrayIdType *globalNodeIds)
{
  _mesh.takeRef(mesh);
  _cell_global.takeRef(globalCellIds);
  _node_global.takeRef(globalNodeIds);
  _mesh.checkNotNull();
  _cell_global.checkNotNull();
  _node_global.checkNotNull();
  _mesh->checkConsistencyLight();
  if(_mesh->getNumberOfNodes() != _node_global->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("ParaUMesh constructor : mismatch between # nodes and len of global # nodes.");
  if(_mesh->getNumberOfCells() != _cell_global->getNumberOfTuples())
    throw INTERP_KERNEL::Exception("ParaUMesh constructor : mismatch between # cells and len of global # cells.");
}

std::size_t ParaUMesh::getHeapMemorySizeWithoutChildren() const
{
  return 0;
}

std::vector<const BigMemoryObject *> ParaUMesh::getDirectChildrenWithNull() const
{
  return {_mesh,_cell_global,_node_global};
}

/*!
* This method computes the cells part of distributed mesh lying on \a globalNodeIds nodes.
* The input \a globalNodeIds are not supposed to reside on the current process.
*/
MCAuto<DataArrayIdType> ParaUMesh::getCellIdsLyingOnNodes(const DataArrayIdType *globalNodeIds, bool fullyIn) const
{
  if(fullyIn)
    return this->getCellIdsLyingOnNodesTrue(globalNodeIds);
  else
    return this->getCellIdsLyingOnNodesFalse(globalNodeIds);
}

MCAuto<DataArrayIdType> ParaUMesh::getCellIdsLyingOnNodesTrue(const DataArrayIdType *globalNodeIds) const
{
  MPI_Comm comm(MPI_COMM_WORLD);
  CommInterface ci;
  int size;
  ci.commSize(comm,&size);
  std::unique_ptr<mcIdType[]> nbOfElems(new mcIdType[size]),nbOfElems2(new mcIdType[size]),nbOfElems3(new mcIdType[size]);
  mcIdType nbOfNodeIdsLoc(globalNodeIds->getNumberOfTuples());
  ci.allGather(&nbOfNodeIdsLoc,1,MPI_ID_TYPE,nbOfElems.get(),1,MPI_ID_TYPE,comm);
  std::vector< MCAuto<DataArrayIdType> > tabs(size);
  //store for each proc the local nodeids intercepted by current proc
  int nbOfCollectiveCalls = 1;// this parameter controls the memory peak
  // loop to avoid to all procs to have all the nodes per proc
  for(int subDiv = 0 ; subDiv < nbOfCollectiveCalls ; ++subDiv)
  {
    std::unique_ptr<mcIdType[]> nbOfElemsSp(CommInterface::SplitArrayOfLength(nbOfElems,size,subDiv,nbOfCollectiveCalls));
    mcIdType nbOfNodeIdsSum(std::accumulate(nbOfElemsSp.get(),nbOfElemsSp.get()+size,0));
    std::unique_ptr<mcIdType[]> allGlobalNodeIds(new mcIdType[nbOfNodeIdsSum]);
    std::unique_ptr<int[]> nbOfElemsInt( CommInterface::ToIntArray<mcIdType>(nbOfElemsSp,size) );
    std::unique_ptr<int[]> offsetsIn( CommInterface::ComputeOffset(nbOfElemsInt,size) );
    mcIdType startGlobalNodeIds,endGlobalNodeIds;
    DataArray::GetSlice(0,globalNodeIds->getNumberOfTuples(),1,subDiv,nbOfCollectiveCalls,startGlobalNodeIds,endGlobalNodeIds);
    ci.allGatherV(globalNodeIds->begin()+startGlobalNodeIds,FromIdType<int>(endGlobalNodeIds-startGlobalNodeIds),MPI_ID_TYPE,allGlobalNodeIds.get(),nbOfElemsInt.get(),offsetsIn.get(),MPI_ID_TYPE,comm);
    mcIdType offset(0);
    for(int curRk = 0 ; curRk < size ; ++curRk)
    {
      MCAuto<DataArrayIdType> globalNodeIdsOfCurProc(DataArrayIdType::New());
      globalNodeIdsOfCurProc->useArray(allGlobalNodeIds.get()+offset,false,DeallocType::CPP_DEALLOC,nbOfElemsSp[curRk],1);
      offset += nbOfElemsSp[curRk];
      MCAuto<DataArrayIdType> globalNodeIdsCaptured(_node_global->buildIntersection(globalNodeIdsOfCurProc));
      MCAuto<DataArrayIdType> localNodeIdsToLocate(_node_global->findIdForEach(globalNodeIdsCaptured->begin(),globalNodeIdsCaptured->end()));
      if(tabs[curRk].isNull())
        tabs[curRk] = localNodeIdsToLocate;
      else
        tabs[curRk]->insertAtTheEnd(localNodeIdsToLocate->begin(),localNodeIdsToLocate->end());
    }
  }

  for(int curRk = 0 ; curRk < size ; ++curRk)
  {
    MCAuto<DataArrayIdType> localNodeIds(tabs[curRk]);
    localNodeIds->sort();
    MCAuto<DataArrayIdType> localNodeIdsUnique(localNodeIds->buildUnique());
    MCAuto<DataArrayIdType> localCellCaptured(_mesh->getCellIdsLyingOnNodes(localNodeIdsUnique->begin(),localNodeIdsUnique->end(),true));
    MCAuto<DataArrayIdType> localCellCapturedGlob(_cell_global->selectByTupleIdSafe(localCellCaptured->begin(),localCellCaptured->end()));
    tabs[curRk] = localCellCapturedGlob;
  }
  
  for(int curRk = 0 ; curRk < size ; ++curRk)
  {
    tabs[curRk] = tabs[curRk]->buildUniqueNotSorted();
    nbOfElems3[curRk] = tabs[curRk]->getNumberOfTuples();
  }
  std::vector<const DataArrayIdType *> tabss(tabs.begin(),tabs.end());
  MCAuto<DataArrayIdType> cells(DataArrayIdType::Aggregate(tabss));
  ci.allToAll(nbOfElems3.get(),1,MPI_ID_TYPE,nbOfElems2.get(),1,MPI_ID_TYPE,comm);
  mcIdType nbOfCellIdsSum(std::accumulate(nbOfElems2.get(),nbOfElems2.get()+size,0));
  MCAuto<DataArrayIdType> cellIdsFromProcs(DataArrayIdType::New());
  cellIdsFromProcs->alloc(nbOfCellIdsSum,1);
  {
    std::unique_ptr<int[]> nbOfElemsInt( CommInterface::ToIntArray<mcIdType>(nbOfElems3,size) ),nbOfElemsOutInt( CommInterface::ToIntArray<mcIdType>(nbOfElems2,size) );
    std::unique_ptr<int[]> offsetsIn( CommInterface::ComputeOffset(nbOfElemsInt,size) ), offsetsOut( CommInterface::ComputeOffset(nbOfElemsOutInt,size) );
    ci.allToAllV(cells->begin(),nbOfElemsInt.get(),offsetsIn.get(),MPI_ID_TYPE,
                 cellIdsFromProcs->getPointer(),nbOfElemsOutInt.get(),offsetsOut.get(),MPI_ID_TYPE,comm);
  }
  cellIdsFromProcs->sort();
  return cellIdsFromProcs;
}

MCAuto<DataArrayIdType> ParaUMesh::getCellIdsLyingOnNodesFalse(const DataArrayIdType *globalNodeIds) const
{
  MPI_Comm comm(MPI_COMM_WORLD);
  CommInterface ci;
  int size;
  ci.commSize(comm,&size);
  std::unique_ptr<mcIdType[]> nbOfElems(new mcIdType[size]),nbOfElems2(new mcIdType[size]),nbOfElems3(new mcIdType[size]);
  mcIdType nbOfNodeIdsLoc(globalNodeIds->getNumberOfTuples());
  ci.allGather(&nbOfNodeIdsLoc,1,MPI_ID_TYPE,nbOfElems.get(),1,MPI_ID_TYPE,comm);
  // loop to avoid to all procs to have all the nodes per proc
  int nbOfCollectiveCalls = 1;// this parameter controls the memory peak
  std::vector< MCAuto<DataArrayIdType> > tabs(size);
  for(int subDiv = 0 ; subDiv < nbOfCollectiveCalls ; ++subDiv)
  {
    std::unique_ptr<mcIdType[]> nbOfElemsSp(CommInterface::SplitArrayOfLength(nbOfElems,size,subDiv,nbOfCollectiveCalls));
    mcIdType nbOfNodeIdsSum(std::accumulate(nbOfElemsSp.get(),nbOfElemsSp.get()+size,0));
    std::unique_ptr<mcIdType[]> allGlobalNodeIds(new mcIdType[nbOfNodeIdsSum]);
    std::unique_ptr<int[]> nbOfElemsInt( CommInterface::ToIntArray<mcIdType>(nbOfElemsSp,size) );
    std::unique_ptr<int[]> offsetsIn( CommInterface::ComputeOffset(nbOfElemsInt,size) );
    mcIdType startGlobalNodeIds,endGlobalNodeIds;
    DataArray::GetSlice(0,globalNodeIds->getNumberOfTuples(),1,subDiv,nbOfCollectiveCalls,startGlobalNodeIds,endGlobalNodeIds);
    ci.allGatherV(globalNodeIds->begin()+startGlobalNodeIds,FromIdType<int>(endGlobalNodeIds-startGlobalNodeIds),MPI_ID_TYPE,allGlobalNodeIds.get(),nbOfElemsInt.get(),offsetsIn.get(),MPI_ID_TYPE,comm);
    mcIdType offset(0);
    for(int curRk = 0 ; curRk < size ; ++curRk)
    {
      MCAuto<DataArrayIdType> globalNodeIdsOfCurProc(DataArrayIdType::New());
      globalNodeIdsOfCurProc->useArray(allGlobalNodeIds.get()+offset,false,DeallocType::CPP_DEALLOC,nbOfElemsSp[curRk],1);
      offset += nbOfElemsSp[curRk];
      MCAuto<DataArrayIdType> globalNodeIdsCaptured(_node_global->buildIntersection(globalNodeIdsOfCurProc));
      MCAuto<DataArrayIdType> localNodeIdsToLocate(_node_global->findIdForEach(globalNodeIdsCaptured->begin(),globalNodeIdsCaptured->end()));
      MCAuto<DataArrayIdType> localCellCaptured(_mesh->getCellIdsLyingOnNodes(localNodeIdsToLocate->begin(),localNodeIdsToLocate->end(),false));
      MCAuto<DataArrayIdType> localCellCapturedGlob(_cell_global->selectByTupleIdSafe(localCellCaptured->begin(),localCellCaptured->end()));
      if(tabs[curRk].isNull())
        tabs[curRk] = localCellCapturedGlob;
      else
        tabs[curRk]->insertAtTheEnd(localCellCapturedGlob->begin(),localCellCapturedGlob->end());
    }
  }
  for(int curRk = 0 ; curRk < size ; ++curRk)
  {
    tabs[curRk] = tabs[curRk]->buildUniqueNotSorted();
    nbOfElems3[curRk] = tabs[curRk]->getNumberOfTuples();
  }
  std::vector<const DataArrayIdType *> tabss(tabs.begin(),tabs.end());
  MCAuto<DataArrayIdType> cells(DataArrayIdType::Aggregate(tabss));
  ci.allToAll(nbOfElems3.get(),1,MPI_ID_TYPE,nbOfElems2.get(),1,MPI_ID_TYPE,comm);
  mcIdType nbOfCellIdsSum(std::accumulate(nbOfElems2.get(),nbOfElems2.get()+size,0));
  MCAuto<DataArrayIdType> cellIdsFromProcs(DataArrayIdType::New());
  cellIdsFromProcs->alloc(nbOfCellIdsSum,1);
  {
    std::unique_ptr<int[]> nbOfElemsInt( CommInterface::ToIntArray<mcIdType>(nbOfElems3,size) ),nbOfElemsOutInt( CommInterface::ToIntArray<mcIdType>(nbOfElems2,size) );
    std::unique_ptr<int[]> offsetsIn( CommInterface::ComputeOffset(nbOfElemsInt,size) ), offsetsOut( CommInterface::ComputeOffset(nbOfElemsOutInt,size) );
    ci.allToAllV(cells->begin(),nbOfElemsInt.get(),offsetsIn.get(),MPI_ID_TYPE,
                 cellIdsFromProcs->getPointer(),nbOfElemsOutInt.get(),offsetsOut.get(),MPI_ID_TYPE,comm);
  }
  cellIdsFromProcs->sort();
  return cellIdsFromProcs;
}

DataArrayIdType *ParaUMesh::redistributeCellField(const DataArrayIdType *globalCellIds, const DataArrayIdType *fieldValueToRed) const
{
  return this->redistributeCellFieldT<mcIdType>(globalCellIds,fieldValueToRed);
}

DataArrayDouble *ParaUMesh::redistributeCellField(const DataArrayIdType *globalCellIds, const DataArrayDouble *fieldValueToRed) const
{
  return this->redistributeCellFieldT<double>(globalCellIds,fieldValueToRed);
}

DataArrayIdType *ParaUMesh::redistributeNodeField(const DataArrayIdType *globalCellIds, const DataArrayIdType *fieldValueToRed) const
{
  return this->redistributeNodeFieldT<mcIdType>(globalCellIds,fieldValueToRed);
}

DataArrayDouble *ParaUMesh::redistributeNodeField(const DataArrayIdType *globalCellIds, const DataArrayDouble *fieldValueToRed) const
{
  return this->redistributeNodeFieldT<double>(globalCellIds,fieldValueToRed);
}

/*!
 * Return part of \a this mesh split over COMM_WORLD. Part is defined by global cell ids array \a globaCellIds.
 */
ParaUMesh *ParaUMesh::redistributeCells(const DataArrayIdType *globalCellIds) const
{
  MPI_Comm comm(MPI_COMM_WORLD);
  CommInterface ci;
  std::unique_ptr<mcIdType[]> allGlobalCellIds,allGlobalCellIdsIndex;
  int size(ci.allGatherArrays(comm,globalCellIds,allGlobalCellIds,allGlobalCellIdsIndex));
  // Prepare ParaUMesh parts to be sent : compute for each proc the contribution of current rank.
  std::vector< MCAuto<DataArrayIdType> > globalCellIdsToBeSent(size),globalNodeIdsToBeSent(size);
  std::vector< MCAuto<MEDCouplingUMesh> > meshPartsToBeSent(size);
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
    meshPartsToBeSent[curRk] = meshPart;
    globalCellIdsToBeSent[curRk] = globalCellIdsCaptured;
    globalNodeIdsToBeSent[curRk] = globalNodeIdsPart;
  }
  // Receive
  std::vector< MCAuto<DataArrayIdType> > globalCellIdsReceived,globalNodeIdsReceived;
  ci.allToAllArrays(comm,globalCellIdsToBeSent,globalCellIdsReceived);
  ci.allToAllArrays(comm,globalNodeIdsToBeSent,globalNodeIdsReceived);
  //now exchange the 3 arrays for the umesh : connectivity, connectivityindex and coordinates
  std::vector<const MEDCouplingUMesh *> meshPartsToBeSent2(FromVecAutoToVecOfConst<MEDCouplingUMesh>(meshPartsToBeSent));
  //connectivityindex
  std::vector< MCAuto<DataArrayIdType> > connectivityIndexReceived,connectivityReceived;
  {
    std::vector<const DataArrayIdType *> connectivityIndexToBeSent(UMeshConnectivityIndexIterator(0,&meshPartsToBeSent2),UMeshConnectivityIndexIterator(meshPartsToBeSent2.size(),&meshPartsToBeSent2));
    ci.allToAllArrays(comm,FromVecConstToVecAuto<DataArrayIdType>(connectivityIndexToBeSent),connectivityIndexReceived);
  }
  //connectivity
  {
    std::vector<const DataArrayIdType *> connectivityToBeSent(UMeshConnectivityIterator(0,&meshPartsToBeSent2),UMeshConnectivityIterator(meshPartsToBeSent2.size(),&meshPartsToBeSent2));
    ci.allToAllArrays(comm,FromVecConstToVecAuto<DataArrayIdType>(connectivityToBeSent),connectivityReceived);
  }
  //coordinates
  MCAuto<DataArrayDouble> coords;
  {
    std::vector<const DataArrayDouble *> coordsToBeSent(UMeshCoordsIterator(0,&meshPartsToBeSent2),UMeshCoordsIterator(meshPartsToBeSent2.size(),&meshPartsToBeSent2));
    ci.allToAllArrays(comm,FromVecConstToVecAuto<DataArrayDouble>(coordsToBeSent),coords);
  }
  /////// Sort it all !
  // firstly deal with nodes.
  MCAuto<DataArrayIdType> aggregatedNodeIds( DataArrayIdType::Aggregate(FromVecAutoToVecOfConst<DataArrayIdType>(globalNodeIdsReceived)) );
  MCAuto<DataArrayIdType> aggregatedNodeIdsSorted(aggregatedNodeIds->copySorted());
  MCAuto<DataArrayIdType> nodeIdsIntoAggregatedIds(DataArrayIdType::FindPermutationFromFirstToSecondDuplicate(aggregatedNodeIdsSorted,aggregatedNodeIds));
  MCAuto<DataArrayIdType> idxOfSameNodeIds(aggregatedNodeIdsSorted->indexOfSameConsecutiveValueGroups());
  MCAuto<DataArrayIdType> n2o_nodes(nodeIdsIntoAggregatedIds->selectByTupleIdSafe(idxOfSameNodeIds->begin(),idxOfSameNodeIds->end()-1));//new == new ordering so that global node ids are sorted . old == coarse ordering implied by the aggregation
  MCAuto<DataArrayIdType> finalGlobalNodeIds(aggregatedNodeIdsSorted->selectByTupleIdSafe(idxOfSameNodeIds->begin(),idxOfSameNodeIds->end()-1));
  MCAuto<DataArrayDouble> finalCoords(coords->selectByTupleIdSafe(n2o_nodes->begin(),n2o_nodes->end()));
  finalCoords->copyStringInfoFrom(*_mesh->getCoords());
  // secondly renumbering of node ids in connectivityReceived
  for(int curRk = 0 ; curRk < size ; ++curRk)
  {
    auto current(globalNodeIdsReceived[curRk]);
    MCAuto<DataArrayIdType> aa(finalGlobalNodeIds->findIdForEach(current->begin(),current->end()));
    // work on connectivityReceived[curRk] with transformWithIndArr but do not forget type of cells that should be excluded !
    auto connectivityToModify(connectivityReceived[curRk]);
    auto connectivityIndex(connectivityIndexReceived[curRk]);
    MCAuto<DataArrayIdType> types(connectivityToModify->selectByTupleIdSafe(connectivityIndex->begin(),connectivityIndex->end()-1));
    connectivityToModify->setPartOfValuesSimple3(0,connectivityIndex->begin(),connectivityIndex->end()-1,0,1,1);
    connectivityToModify->transformWithIndArr(aa->begin(),aa->end());
    connectivityToModify->setPartOfValues3(types,connectivityIndex->begin(),connectivityIndex->end()-1,0,1,1,true);
  }
  // thirdly renumber cells
  MCAuto<DataArrayIdType> aggregatedCellIds( DataArrayIdType::Aggregate(FromVecAutoToVecOfConst<DataArrayIdType>(globalCellIdsReceived)) );
  MCAuto<DataArrayIdType> aggregatedCellIdsSorted(aggregatedCellIds->copySorted());
  MCAuto<DataArrayIdType> idsIntoAggregatedIds(DataArrayIdType::FindPermutationFromFirstToSecondDuplicate(aggregatedCellIdsSorted,aggregatedCellIds));
  MCAuto<DataArrayIdType> cellIdsOfSameNodeIds(aggregatedCellIdsSorted->indexOfSameConsecutiveValueGroups());
  MCAuto<DataArrayIdType> n2o_cells(idsIntoAggregatedIds->selectByTupleIdSafe(cellIdsOfSameNodeIds->begin(),cellIdsOfSameNodeIds->end()-1));//new == new ordering so that global cell ids are sorted . old == coarse ordering implied by the aggregation
  // TODO : check coordsReceived==globalCellIds
  MCAuto<DataArrayIdType> connSorted,indicesSorted;
  {
    MCAuto<DataArrayIdType> conn(DataArrayIdType::Aggregate(FromVecAutoToVecOfConst<DataArrayIdType>(connectivityReceived)));
    MCAuto<DataArrayIdType> connIndex(DataArrayIdType::AggregateIndexes(FromVecAutoToVecOfConst<DataArrayIdType>(connectivityIndexReceived))); 
    {
      DataArrayIdType *indicesSortedTmp(nullptr),*valuesSortedTmp(nullptr);
      DataArrayIdType::ExtractFromIndexedArrays(n2o_cells->begin(),n2o_cells->end(),conn,connIndex,valuesSortedTmp,indicesSortedTmp);
      indicesSorted = indicesSortedTmp; connSorted=valuesSortedTmp;
    }
  }
  // finalize all
  MCAuto<MEDCouplingUMesh> mesh(MEDCouplingUMesh::New(_mesh->getName(),_mesh->getMeshDimension()));
  mesh->setConnectivity(connSorted,indicesSorted,true);
  mesh->setCoords(finalCoords);
  mesh->setDescription(_mesh->getDescription());
  MCAuto<ParaUMesh> ret(ParaUMesh::New(mesh,aggregatedCellIdsSorted,finalGlobalNodeIds));
  return ret.retn();
}
