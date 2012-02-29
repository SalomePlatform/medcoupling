//  Copyright (C) 2007-2010  CEA/DEN, EDF R&D
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307 USA
//
//  See http://www.salome-platform.org/ or email : webmaster.salome@opencascade.com
//
#include <cstdio>
extern "C" {
#define restrict
#include "scotch.h"
}
#include "MEDPARTITIONER_Graph.hxx"
#include "MEDPARTITIONER_SCOTCHGraph.hxx"

using namespace MEDPARTITIONER;
  
SCOTCHGraph::SCOTCHGraph():Graph()
{
}

SCOTCHGraph::SCOTCHGraph(MEDPARTITIONER::MEDSKYLINEARRAY* graph, int* edgeweight):Graph(graph,edgeweight)
{
}

SCOTCHGraph::~SCOTCHGraph()
{
}

void SCOTCHGraph::partGraph(int ndomain, const std::string& options_string, ParaDomainSelector* sel)
{
  // number of graph vertices
  int n = _graph->getNumberOf();

  //graph
  int * xadj=const_cast<int*>(_graph->getIndex());
  int * adjncy = const_cast<int*>(_graph->getValue());

  //ndomain
  int nparts = ndomain;

  // output parameters
  int* partition = new int[n+1];

  SCOTCH_Graph scotch_graph;

  SCOTCH_graphInit(&scotch_graph);


  SCOTCH_graphBuild(&scotch_graph,
                    1, //premier indice 1
                    n, // nb of graph nodes
                    xadj,
                    0,
                    _cellweight, //graph vertices loads
                    0,
                    xadj[n], // number of edges
                    adjncy,
                    _edgeweight);

  SCOTCH_Strat scotch_strategy;           
  SCOTCH_stratInit(&scotch_strategy);

  //!user-defined options for the strategy
  if (options_string!="")
    SCOTCH_stratGraphMap(&scotch_strategy,options_string.c_str());


  if (nparts>1)           
    SCOTCH_graphPart(&scotch_graph,nparts,&scotch_strategy,partition);
  else
    // partition for 1 subdomain
    for (int i=0; i<n+1; i++)
      partition[i]=0;

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

  //creating a skylinearray with no copy of the index and partition array
  // the fifth argument true specifies that only the pointers are passed 
  //to the object
  
  _partition = new MEDPARTITIONER::MEDSKYLINEARRAY(index,value);

}
