// Copyright (C) 2017-2019  CEA/DEN, EDF R&D
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

#include "MEDPARTITIONER_PTScotchGraph.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "MEDCouplingSkyLineArray.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MCType.hxx"

#include <cstdio>
#include <mpi.h>

#ifdef MED_ENABLE_PTSCOTCH
extern "C"
{
#define restrict
#include "ptscotch.h"
}
#endif

using namespace MEDPARTITIONER;


PTSCOTCHGraph::PTSCOTCHGraph(MEDCoupling::MEDCouplingSkyLineArray *graph, int *edgeweight, DataArrayInt *vlbloctab):Graph(graph,edgeweight),_vlbloctab(vlbloctab)
{
}

PTSCOTCHGraph::~PTSCOTCHGraph()
{
}

void PTSCOTCHGraph::partGraph(int ndomain, const std::string& options_string, ParaDomainSelector* sel)
{
  if (MyGlobals::_Verbose>10)
    std::cout << "proc " << MyGlobals::_Rank << " : PTSCOTCHGraph::partGraph" << std::endl;
  
  //number of graph vertices
  int n = _graph->getNumberOf();
  //graph
  int * xadj=const_cast<int*>(_graph->getIndex());
  int * adjncy=const_cast<int*>(_graph->getValues());
  //ndomain
  int nparts=ndomain;

#if !defined(MED_ENABLE_PTSCOTCH)
  throw INTERP_KERNEL::Exception("PTSCOTCHGraph::partGraph : PTSCOTCH is not available. Check your products, please.");
#else
  //output parameters
  int* partition = new int[n+1];
  
  int* vlbloctab = _vlbloctab?const_cast<int*>(_vlbloctab->begin()):0;
  
  SCOTCH_Dgraph scotch_graph;
  SCOTCH_dgraphInit(&scotch_graph, MPI_COMM_WORLD);
  SCOTCH_dgraphBuild(&scotch_graph,
                     0,             // baseval               , base first indice 0
                     n,             // vertlocnbr            , nb of local graph nodes
                     n,             // vertlocmax            , should be set to vertlocnbr for graphs without holes
                     xadj,          // vertloctab[vertnbr+1] , index vertex table
                     0,             // vendloctab            , index end vertex table if disjoint, set to zero
                     _cell_weight,  // veloloctab            , graph vertices loads, set to zero
                     vlbloctab,     // vlblocltab            , vertex label array : global vertex index
                     xadj[n],       // edgelocnbr            , number of edges
                     xadj[n],       // edgelocsiz            , same as edgelocnbr if edgeloctab is compact
                     adjncy,        // edgeloctab[edgelocnbr], global indexes of edges
                     0,             // edgegsttab            , optional, should be computed internally, set to zero
                     _edge_weight); // edloloctab            , graph edges loads, set to zero
  
  SCOTCH_Strat scotch_strategy;
  SCOTCH_stratInit(&scotch_strategy);
  
  //!user-defined options for the strategy
  if (options_string!="")
    SCOTCH_stratGraphMap(&scotch_strategy,options_string.c_str());

  if (nparts>1)
    {
      if (MyGlobals::_Verbose>10) std::cout << "SCOTCHGraph::graphPart SCOTCH_graphPart" << std::endl;
      SCOTCH_dgraphPart(&scotch_graph,nparts,&scotch_strategy,partition);
    }
  else  //partition for 1 subdomain
    {
    for (int i=0; i<n+1; i++)
      partition[i]=0;
    }
  
  SCOTCH_stratExit(&scotch_strategy);
  SCOTCH_dgraphExit(&scotch_graph);

  std::vector<int> index(n+1);
  std::vector<int> value(n);
  index[0]=0;
  for (int i=0; i<n; i++)
    {
      index[i+1]=index[i]+1;
      value[i]=partition[i];
    }
  delete [] partition;
  
  //creating a skylinearray with no copy of the index and partition array
  //the fifth argument true specifies that only the pointers are passed 
  //to the object
  _partition = MEDCoupling::MEDCouplingSkyLineArray::New(index,value);
#endif
}

