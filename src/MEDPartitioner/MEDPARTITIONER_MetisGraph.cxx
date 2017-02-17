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

#include "MEDPARTITIONER_MetisGraph.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "MEDCouplingSkyLineArray.hxx"
#include "InterpKernelException.hxx"

#include <iostream>

extern "C"
{
#include "MEDPARTITIONER_metis.h"
}

using namespace MEDPARTITIONER;

METISGraph::METISGraph():Graph()
{
}

METISGraph::METISGraph(MEDCoupling::MEDCouplingSkyLineArray* graph, int* edgeweight)
  :Graph(graph,edgeweight)
{
}

METISGraph::~METISGraph()
{
}

void METISGraph::partGraph(int ndomain,
                           const std::string& options_string,
                           ParaDomainSelector *parallelizer)
{
  using std::vector;
  if (MyGlobals::_Verbose>10)
    std::cout << "proc " << MyGlobals::_Rank << " : METISGraph::partGraph" << std::endl;
  
  //number of graph vertices
  int n=_graph->getNumberOf();
  //graph
  int * xadj=const_cast<int*>(_graph->getIndex());
  int * adjncy=const_cast<int*>(_graph->getValues());
  //constraints
  int * vwgt=_cell_weight;
  int * adjwgt=_edge_weight;
  int wgtflag=(_edge_weight!=0)?1:0+(_cell_weight!=0)?2:0;
  //base 0 or 1
  int base=0;
  //ndomain
  int nparts=ndomain;
  //options
  /*
    (0=default_option,option,random_seed) see defs.h
    #define PMV3_OPTION_DBGLVL 1
    #define PMV3_OPTION_SEED 2
    #define PMV3_OPTION_IPART 3
    #define PMV3_OPTION_PSR 3
    seems no changes int options[4]={1,0,33,0}; //test for a random seed of 33
  */
  int options[20]={0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0};
#if !defined(MED_ENABLE_METIS)
  throw INTERP_KERNEL::Exception("METISGraph::partGraph : METIS is not available. Check your products, please.");
#else
  //output parameters
  int edgecut;
  int* partition=new int[n];

  if(nparts >1)
    {
      if (MyGlobals::_Verbose>10) 
        std::cout << "METISGraph::partGraph METIS_PartGraph METIS_PartGraph(RecursiveOrKway)" << std::endl;
      if (options_string != "k")
        MEDPARTITIONER_METIS_PartGraphRecursive(&n, xadj, adjncy, vwgt, adjwgt, &wgtflag,
                                                &base, &nparts, options, &edgecut, partition);
      else
        MEDPARTITIONER_METIS_PartGraphKway(&n, xadj, adjncy, vwgt, adjwgt, &wgtflag,
                                           &base, &nparts, options, &edgecut, partition);
    }
  else  //force this case because METIS send all 1 in value
    {
      for (int i=0; i<n; i++)
        partition[i]=0;
    }
  vector<int> index(n+1);
  vector<int> value(n);
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

