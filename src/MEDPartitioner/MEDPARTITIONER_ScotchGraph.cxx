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

#include "MEDPARTITIONER_Graph.hxx"
#include "MEDPARTITIONER_ScotchGraph.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "MEDCouplingSkyLineArray.hxx"

#include <cstdio>

#ifdef MED_ENABLE_SCOTCH
extern "C"
{
#define restrict
#include "scotch.h"
}
#endif

using namespace MEDPARTITIONER;
  
SCOTCHGraph::SCOTCHGraph():Graph()
{
}

SCOTCHGraph::SCOTCHGraph(MEDCoupling::MEDCouplingSkyLineArray* graph, int* edgeweight):Graph(graph,edgeweight)
{
}

SCOTCHGraph::~SCOTCHGraph()
{
}

void SCOTCHGraph::partGraph(int ndomain, const std::string& options_string, ParaDomainSelector* sel)
{
  if (MyGlobals::_Verbose>10)
    std::cout << "proc " << MyGlobals::_Rank << " : SCOTCHGraph::partGraph" << std::endl;
  
  //number of graph vertices
  int n = _graph->getNumberOf();
  //graph
  int * xadj=const_cast<int*>(_graph->getIndex());
  int * adjncy=const_cast<int*>(_graph->getValues());
  //ndomain
  int nparts=ndomain;

#if !defined(MED_ENABLE_SCOTCH)
  throw INTERP_KERNEL::Exception("SCOTCHGraph::partGraph : SCOTCH is not available. Check your products, please.");
#else
  //output parameters
  int* partition = new int[n+1];
  
  SCOTCH_Graph scotch_graph;
  SCOTCH_graphInit(&scotch_graph);
  SCOTCH_graphBuild(&scotch_graph,
                    0, //base first indice 0
                    n, //nb of graph nodes
                    xadj,
                    0,
                    _cell_weight, //graph vertices loads
                    0,
                    xadj[n], //number of edges
                    adjncy,
                    _edge_weight);
  SCOTCH_Strat scotch_strategy;
  SCOTCH_stratInit(&scotch_strategy);
  
  //!user-defined options for the strategy
  if (options_string!="")
    SCOTCH_stratGraphMap(&scotch_strategy,options_string.c_str());

  if (nparts>1)
    {
      if (MyGlobals::_Verbose>10) std::cout << "SCOTCHGraph::graphPart SCOTCH_graphPart" << std::endl;
      SCOTCH_graphPart(&scotch_graph,nparts,&scotch_strategy,partition);
    }
  else  //partition for 1 subdomain
    {
    for (int i=0; i<n+1; i++)
      partition[i]=0;
    }
  
  SCOTCH_stratExit(&scotch_strategy);
  SCOTCH_graphExit(&scotch_graph);

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
