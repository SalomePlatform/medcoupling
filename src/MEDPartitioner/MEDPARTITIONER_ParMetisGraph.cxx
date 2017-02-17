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

#include "MEDPARTITIONER_ParMetisGraph.hxx"
#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "MEDCouplingSkyLineArray.hxx"
#include "InterpKernelException.hxx"

#include <iostream>

#ifdef MED_ENABLE_PARMETIS
#include <parmetis.h>
// #if PARMETIS_MAJOR_VERSION == 4
//    #define ParMETIS_PartKway ParMETIS_V3_PartKway
// #endif
#endif

using namespace MEDPARTITIONER;

ParMETISGraph::ParMETISGraph():Graph()
{
}

ParMETISGraph::ParMETISGraph(MEDCoupling::MEDCouplingSkyLineArray* graph, int* edgeweight)
  :Graph(graph,edgeweight)
{
}

ParMETISGraph::~ParMETISGraph()
{
}

void ParMETISGraph::partGraph(int ndomain,
                           const std::string& options_string,
                           ParaDomainSelector *parallelizer)
{
  using std::vector;
  vector<int> ran,vx,va; //for randomize
  
  if (MyGlobals::_Verbose>10)
    std::cout << "proc " << MyGlobals::_Rank << " : ParMETISGraph::partGraph" << std::endl;
  
  // number of graph vertices
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
  int options[4]={0,0,0,0};
  // output parameters
  int edgecut;
#if !defined(MED_ENABLE_PARMETIS)
  throw INTERP_KERNEL::Exception("ParMETISGraph::partGraph : PARMETIS is not available. Check your products, please.");
#else
  int* partition=new int[n];
  
  if (MyGlobals::_Verbose>10) 
    std::cout << "proc " << MyGlobals::_Rank << " : ParMETISGraph::partGraph ParMETIS_PartKway new" << std::endl;
  int * vtxdist=parallelizer->getProcVtxdist();
  MPI_Comm comm=MPI_COMM_WORLD;
  ParMETIS_PartKway(vtxdist, xadj, adjncy, vwgt, 
                                    adjwgt, &wgtflag, &base, &nparts, options, 
                                    &edgecut, partition, &comm );


  /*doc from parmetis.h
    void __cdecl ParMETIS_PartKway(
    idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
    idxtype *adjwgt, int *wgtflag, int *numflag, int *nparts, int *options, 
    int *edgecut, idxtype *part, MPI_Comm *comm);

    void __cdecl ParMETIS_V3_PartKway(
    idxtype *vtxdist, idxtype *xadj, idxtype *adjncy, idxtype *vwgt, 
    idxtype *adjwgt, int *wgtflag, int *numflag, int *ncon, int *nparts, 
    float *tpwgts, float *ubvec, int *options, int *edgecut, idxtype *part, 
    MPI_Comm *comm);
  */

  vector<int> index(n+1);
  vector<int> value(n);
  index[0]=0;
  if (ran.size()>0 && MyGlobals::_Atomize==0) //there is randomize
    {
      if (MyGlobals::_Is0verbose>100)
        std::cout << "randomize" << std::endl;
      for (int i=0; i<n; i++)
        {
          index[i+1]=index[i]+1;
          value[ran[i]]=partition[i];
        }
    }
  else
    {
      for (int i=0; i<n; i++)
        {
          index[i+1]=index[i]+1;
          value[i]=partition[i];
        }
    }
  delete [] partition;

  //creating a skylinearray with no copy of the index and partition array
  //the fifth argument true specifies that only the pointers are passed 
  //to the object
  
  _partition = MEDCoupling::MEDCouplingSkyLineArray::New(index,value);
#endif
}

