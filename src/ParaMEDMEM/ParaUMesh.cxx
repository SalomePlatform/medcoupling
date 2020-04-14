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

/*!
* This method computes the cells part of distributed mesh lying on \a globalNodeIds nodes.
* The input \a globalNodeIds are not supposed to reside on the current process.
*/
MCAuto<DataArrayIdType> ParaUMesh::getCellIdsLyingOnNodes(const DataArrayIdType *globalNodeIds, bool fullyIn) const
{
  if(fullyIn)
    throw INTERP_KERNEL::Exception("ParaUMesh::getCellIdsLyingOnNodes : not implemented yet for fullyIn == True !");
  MPI_Comm comm(MPI_COMM_WORLD);
  CommInterface ci;
  int size;
  ci.commSize(comm,&size);
  std::unique_ptr<mcIdType[]> nbOfElems(new mcIdType[size]),nbOfElems2(new mcIdType[size]),nbOfElems3(new mcIdType[size]);
  mcIdType nbOfNodeIdsLoc(globalNodeIds->getNumberOfTuples());
  ci.allGather(&nbOfNodeIdsLoc,1,MPI_ID_TYPE,nbOfElems.get(),1,MPI_ID_TYPE,comm);
  std::vector< MCAuto<DataArrayIdType> > tabs(size);
  // loop to avoid to all procs to have all the nodes per proc
  for(int subDiv = 0 ; subDiv < size ; ++subDiv)
  {
    std::unique_ptr<mcIdType[]> nbOfElemsSp(CommInterface::SplitArrayOfLength(nbOfElems,size,subDiv,size));
    mcIdType nbOfNodeIdsSum(std::accumulate(nbOfElemsSp.get(),nbOfElemsSp.get()+size,0));
    std::unique_ptr<mcIdType[]> allGlobalNodeIds(new mcIdType[nbOfNodeIdsSum]);
    std::unique_ptr<int[]> nbOfElemsInt( CommInterface::ToIntArray<mcIdType>(nbOfElemsSp,size) );
    std::unique_ptr<int[]> offsetsIn( CommInterface::ComputeOffset(nbOfElemsInt,size) );
    mcIdType startGlobalNodeIds,endGlobalNodeIds;
    DataArray::GetSlice(0,globalNodeIds->getNumberOfTuples(),1,subDiv,size,startGlobalNodeIds,endGlobalNodeIds);
    ci.allGatherV(globalNodeIds->begin()+startGlobalNodeIds,endGlobalNodeIds-startGlobalNodeIds,MPI_ID_TYPE,allGlobalNodeIds.get(),nbOfElemsInt.get(),offsetsIn.get(),MPI_ID_TYPE,comm);
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

/*!
 */
MCAuto<ParaUMesh> ParaUMesh::redistributeCells(const DataArrayIdType *globalCellIds) const
{
  MPI_Comm comm(MPI_COMM_WORLD);
  CommInterface ci;
  int size;
  ci.commSize(comm,&size);
  std::unique_ptr<mcIdType[]> nbOfElems(new mcIdType[size]);
  mcIdType nbOfCellsRequested(globalCellIds->getNumberOfTuples());
  ci.allGather(&nbOfCellsRequested,1,MPI_ID_TYPE,nbOfElems.get(),1,MPI_ID_TYPE,comm);
  mcIdType nbOfCellIdsSum(std::accumulate(nbOfElems.get(),nbOfElems.get()+size,0));
  std::unique_ptr<mcIdType[]> allGlobalCellIds(new mcIdType[nbOfCellIdsSum]);
  std::unique_ptr<int[]> nbOfElemsInt( CommInterface::ToIntArray<mcIdType>(nbOfElems,size) );
  std::unique_ptr<int[]> offsetsIn( CommInterface::ComputeOffset(nbOfElemsInt,size) );
  ci.allGatherV(globalCellIds->begin(),nbOfCellsRequested,MPI_ID_TYPE,allGlobalCellIds.get(),nbOfElemsInt.get(),offsetsIn.get(),MPI_ID_TYPE,comm);
  mcIdType offset(0);
  // Prepare ParaUMesh parts to be sent : compute for each proc the contribution of current rank.
  std::vector< MCAuto<DataArrayIdType> > globalCellIdsToBeSent(size),globalNodeIdsToBeSent(size);
  std::vector< MCAuto<MEDCouplingUMesh> > meshPartsToBeSent(size);
  for(int curRk = 0 ; curRk < size ; ++curRk)
  {
    MCAuto<DataArrayIdType> globalCellIdsOfCurProc(DataArrayIdType::New());
    globalCellIdsOfCurProc->useArray(allGlobalCellIds.get()+offset,false,DeallocType::CPP_DEALLOC,nbOfElems[curRk],1);
    offset += nbOfElems[curRk];
    // the key call is here : compute for rank curRk the cells to be sent
    MCAuto<DataArrayIdType> globalCellIdsCaptured(_cell_global->buildIntersection(globalCellIdsOfCurProc));// OK for the global cellIds
    MCAuto<DataArrayIdType> localCellIdsCaptured(_node_global->findIdForEach(globalCellIdsCaptured->begin(),globalCellIdsCaptured->end()));
    MCAuto<MEDCouplingUMesh> meshPart(_mesh->buildPartOfMySelf(localCellIdsCaptured->begin(),localCellIdsCaptured->end(),true));
    MCAuto<DataArrayIdType> o2n(meshPart->zipCoordsTraducer());// OK for the mesh
    MCAuto<DataArrayIdType> n2o(o2n->invertArrayO2N2N2O(meshPart->getNumberOfNodes()));
    MCAuto<DataArrayIdType> globalNodeIdsPart(_node_global->selectByTupleIdSafe(n2o->begin(),n2o->end())); // OK for the global nodeIds
    meshPartsToBeSent[curRk] = meshPart;
    globalCellIdsToBeSent[curRk] = globalCellIdsCaptured;
    globalNodeIdsToBeSent[curRk] = globalNodeIdsPart;
  }
  // Receive 
}
