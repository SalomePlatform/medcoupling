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

#include "MEDPARTITIONER_MEDPartitioner.hxx"
#include "MEDPARTITIONER_MeshCollection.hxx"
#include "MEDPARTITIONER_Topology.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_ParallelTopology.hxx"
#include "MEDPARTITIONER_Utils.hxx"
#include "MEDPARTITIONER_Graph.hxx"
#ifdef MED_ENABLE_METIS
#  include "MEDPARTITIONER_MetisGraph.hxx"
#endif
#ifdef MED_ENABLE_SCOTCH
#  include "MEDPARTITIONER_ScotchGraph.hxx"
#endif
#include "MEDPARTITIONER_MeshCollectionDriver.hxx"

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingSkyLineArray.hxx"

#include <iostream>
#include <vector>

MEDPARTITIONER::MEDPartitioner::MEDPartitioner(const std::string& filename, int ndomains, const std::string& library,bool create_boundary_faces, bool create_joints, bool mesure_memory):
  _input_collection( 0 ), _output_collection( 0 ), _new_topology( 0 )
{
  MyGlobals::_World_Size = 1;
  MyGlobals::_Rank = 0;
  MyGlobals::_Create_Boundary_Faces = create_boundary_faces;
  MyGlobals::_Create_Joints = create_joints;

  ParaDomainSelector parallelizer(mesure_memory);
  _input_collection=new MeshCollection(filename,parallelizer);
  _input_collection->setParaDomainSelector( &parallelizer );

  MEDPARTITIONER::ParallelTopology* aPT =
    (MEDPARTITIONER::ParallelTopology*) _input_collection->getTopology();
  aPT->setGlobalNumerotationDefault( _input_collection->getParaDomainSelector() );
  _input_collection->prepareFieldDescriptions();
  createPartitionCollection(ndomains, library, create_boundary_faces, create_joints, mesure_memory);

  parallelizer.evaluateMemory();
}

MEDPARTITIONER::MEDPartitioner::MEDPartitioner(const MEDCoupling::MEDFileData* filedata, int ndomains, const std::string& library,bool create_boundary_faces, bool create_joints, bool mesure_memory):
  _input_collection( 0 ), _output_collection( 0 ), _new_topology( 0 )
{
  MyGlobals::_World_Size = 1;
  MyGlobals::_Rank = 0;
  MyGlobals::_Create_Boundary_Faces = create_boundary_faces;
  MyGlobals::_Create_Joints = create_joints;

  ParaDomainSelector parallelizer(mesure_memory);
  _input_collection=new MeshCollection();
  _input_collection->setParaDomainSelector( &parallelizer );
  _input_collection->retrieveDriver()->readMEDFileData(filedata);

  MEDPARTITIONER::ParallelTopology* aPT =
    (MEDPARTITIONER::ParallelTopology*) _input_collection->getTopology();
  aPT->setGlobalNumerotationDefault( _input_collection->getParaDomainSelector() );
  _input_collection->prepareFieldDescriptions();
  createPartitionCollection(ndomains, library, create_boundary_faces, create_joints, mesure_memory);

  parallelizer.evaluateMemory();
}

MEDPARTITIONER::MEDPartitioner::MEDPartitioner(const MEDCoupling::MEDFileData* filedata, MEDPARTITIONER ::Graph* graph, bool create_boundary_faces, bool create_joints, bool mesure_memory):
  _input_collection( 0 ), _output_collection( 0 ), _new_topology( 0 )
{
  MyGlobals::_World_Size = 1;
  MyGlobals::_Rank = 0;
  MyGlobals::_Create_Boundary_Faces = create_boundary_faces;
  MyGlobals::_Create_Joints = create_joints;

  ParaDomainSelector parallelizer(mesure_memory);
  _input_collection=new MeshCollection();
  _input_collection->setParaDomainSelector( &parallelizer );
  _input_collection->retrieveDriver()->readMEDFileData(filedata);

  MEDPARTITIONER::ParallelTopology* aPT =
    (MEDPARTITIONER::ParallelTopology*) _input_collection->getTopology();
  aPT->setGlobalNumerotationDefault( _input_collection->getParaDomainSelector() );
  _input_collection->prepareFieldDescriptions();

  _new_topology = new MEDPARTITIONER::ParallelTopology( graph, aPT, graph->nbDomains(), _input_collection->getMeshDimension() );
  _output_collection=new MeshCollection(*_input_collection,_new_topology,false,false);
  _output_collection->filterFaceOnCell();

  parallelizer.evaluateMemory();
}

MEDPARTITIONER::MEDPartitioner::~MEDPartitioner()
{
  delete _input_collection; _input_collection = 0;
  delete _output_collection; _output_collection = 0;
  delete _new_topology; _new_topology = 0;
}

void MEDPARTITIONER::MEDPartitioner::createPartitionCollection(int ndomains, const std::string& library,bool create_boundary_faces, bool create_joints, bool mesure_memory)
{
  //ParallelTopology* aPT = (ParallelTopology*) _input_collection->getTopology();
  if (library == "metis")
    _new_topology = _input_collection->createPartition(ndomains,MEDPARTITIONER::Graph::METIS);
  else
    _new_topology = _input_collection->createPartition(ndomains,MEDPARTITIONER::Graph::SCOTCH);
  _output_collection=new MeshCollection(*_input_collection,_new_topology,false,false);
  _output_collection->filterFaceOnCell();
}

void MEDPARTITIONER::MEDPartitioner::write(const std::string& filename)
{
  ParaDomainSelector parallelizer(false);
  _output_collection->setParaDomainSelector( &parallelizer );
  _output_collection->write(filename);
  parallelizer.evaluateMemory();
}

MEDCoupling::MEDFileData* MEDPARTITIONER::MEDPartitioner::getMEDFileData()
{
  return _output_collection->retrieveDriver()->getMEDFileData();
}

MEDPARTITIONER::Graph* MEDPARTITIONER::MEDPartitioner::Graph(MEDCoupling::MEDCouplingSkyLineArray* graph, Graph::splitter_type split, int* edgeweight)
{
  MEDPARTITIONER::Graph* cellGraph=0;
  // will be destroyed by XXXGraph class:
  MEDCoupling::MEDCouplingSkyLineArray* arr = MEDCoupling::MEDCouplingSkyLineArray::New(graph->getIndexArray(), graph->getValuesArray());
  switch (split)
    {
    case Graph::METIS:
      if ( !cellGraph )
        {
#ifdef MED_ENABLE_METIS
          cellGraph=new METISGraph(arr,edgeweight);
#endif
        }
      if ( !cellGraph )
        throw INTERP_KERNEL::Exception("MEDPartitioner::Graph : PARMETIS/METIS is not available. Check your products, please.");
      break;
    case Graph::SCOTCH:
#ifdef MED_ENABLE_SCOTCH
      cellGraph=new SCOTCHGraph(arr,edgeweight);
#else
      throw INTERP_KERNEL::Exception("MEDPartitioner::Graph : SCOTCH is not available. Check your products, please.");
#endif
      break;
    }
  return cellGraph;
}
