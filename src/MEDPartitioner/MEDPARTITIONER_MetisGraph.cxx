// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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

#include "InterpKernelException.hxx"

#include <iostream>

#ifdef ENABLE_PARMETIS
#include <parmetis.h>
#endif
extern "C" {
#include <metis.h>
}

using namespace std;
using namespace ParaMEDMEM;
using namespace MEDPARTITIONER;

METISGraph::METISGraph():Graph()
{
}

METISGraph::METISGraph(MEDPARTITIONER::SkyLineArray* graph, int* edgeweight)
  :Graph(graph,edgeweight)
{
}

METISGraph::~METISGraph()
{
}

void METISGraph::partGraph(int ndomain,
                           const std::string& options_string,
                           ParaDomainSelector* parallelizer)
{
  using std::vector;
  vector<int> ran,vx,va; //for randomize
  
  if (MyGlobals::_Verbose>10) cout<<"proc "<<MyGlobals::_Rank<<" : METISGraph::partGraph"<<endl;
  
  // number of graph vertices
  int n=_graph->getNumberOf();
  //graph
  int * xadj=const_cast<int*>(_graph->getIndex());
  int * adjncy=const_cast<int*>(_graph->getValue());
  //constraints
  int * vwgt=_cellweight;
  int * adjwgt=_edgeweight;
  int wgtflag=(_edgeweight!=0)?1:0+(_cellweight!=0)?2:0;
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
  int* partition=new int[n];

  //if (MyGlobals::_Verbose>10) cout<<"proc "<<MyGlobals::_Rank<<" : METISGraph::partGraph n="<<n<<endl;
  if (nparts >1)
    {
      if ( parallelizer )
        {
#ifdef ENABLE_PARMETIS
          // distribution of vertices of the graph among the processors
          if (MyGlobals::_Verbose>100) 
            cout<<"proc "<<MyGlobals::_Rank
                <<" : METISGraph::partGraph ParMETIS_PartKway"<<endl;
          int * vtxdist=parallelizer->getProcVtxdist();
          MPI_Comm comm=MPI_COMM_WORLD;
          try
            {
              if (MyGlobals::_Verbose>200) 
                {
                  cout<<"proc "<<MyGlobals::_Rank<<" : vtxdist :";
                  for (int i=0; i<MyGlobals::_World_Size+1; ++i) cout<<vtxdist[i]<<" ";
                  cout<<endl;
          
                  int lgxadj=vtxdist[MyGlobals::_Rank+1]-vtxdist[MyGlobals::_Rank];
                  //cout<<"lgxadj "<<lgxadj<<" "<<n<<endl;
          
                  if (lgxadj>0)
                    {
                      cout<<"\nproc "<<MyGlobals::_Rank<<" : lgxadj "<<lgxadj<<" lgadj "<<xadj[lgxadj+1]<<endl;
                      for (int i=0; i<10; ++i) cout<<xadj[i]<<" ";
                      cout<<"... "<<xadj[lgxadj]<<endl;
                      for (int i=0; i<15; ++i) cout<<adjncy[i]<<" ";
                      int ll=xadj[lgxadj]-1;
                      cout<<"... ["<<ll<<"] "<<adjncy[ll-1]<<" "<<adjncy[ll]<<endl;
                      /*for (int i=0; i<=ll; ++i) {
                        if (adjncy[i]<0) cout<<"***cvw00 error: adjncy[i]<0 "<<i<<endl;
                        }*/
                      int imaxx=0;
                      //for (int ilgxadj=0; ilgxadj<lgxadj; ilgxadj++)
                      for (int ilgxadj=0; ilgxadj<lgxadj; ilgxadj++)
                        {
                          int ilg=xadj[ilgxadj+1]-xadj[ilgxadj];
                          /*if (ilg<0) cout<<"***cvw01 error: ilg<0 in xadj "<<ilgxadj<<endl;
                            if (MyGlobals::_Is0verbose>1000) 
                            {
                            cout<<"\n -cell "<<ilgxadj<<" "<<ilg<<" :";
                            for (int i=0; i<ilg; i++) cout<<" "<<adjncy[xadj[ilgxadj]+i];
                            }*/
                          if (ilg>imaxx) imaxx=ilg;
                        }
                      cout<<"\nproc "<<MyGlobals::_Rank
                          <<" : on "<<lgxadj<<" cells, max neighbourg number (...for one cell) is "<<imaxx<<endl;
                    }
          
                }
              if ((MyGlobals::_Randomize!=0 || MyGlobals::_Atomize!=0) && MyGlobals::_World_Size==1)
                {
                  //randomize initially was for test on ParMETIS error (sometimes)
                  //due to : seems no changes int options[4]={1,0,33,0}; //test for a random seed of 33
                  //it was keeped
                  ran=CreateRandomSize(n);
                  RandomizeAdj(&xadj[0],&adjncy[0],ran,vx,va);
                  ParMETIS_PartKway( //cvwat11
                                    vtxdist, &vx[0], &va[0], vwgt, 
                                    adjwgt, &wgtflag, &base, &nparts, options, 
                                    &edgecut, partition, &comm );
                }
              else
                {
                  //MPI_Barrier(MPI_COMM_WORLD);
                  //cout<<"proc "<<MyGlobals::_Rank<<" : barrier ParMETIS_PartKway done"<<endl;
                  ParMETIS_PartKway( //cvwat11
                                    vtxdist, xadj, adjncy, vwgt, 
                                    adjwgt, &wgtflag, &base, &nparts, options, 
                                    &edgecut, partition, &comm );
                }

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

            }
          catch(...)
            {
              //shit ParMETIS "Error! Key -2 not found!" not catched...
              throw INTERP_KERNEL::Exception(LOCALIZED("Problem in ParMETIS_PartKway"));
            }
          if (n<8 && nparts==3)
            {
              for (int i=0; i<n; i++) partition[i]=i%3;
            }
#else
          throw INTERP_KERNEL::Exception(LOCALIZED("ParMETIS is not available. Check your products, please."));
#endif
          //throw INTERP_KERNEL::Exception(LOCALIZED("ParMETIS is not available. Check your products, please."));
        }
      else
        {
          if (MyGlobals::_Verbose>10) 
            cout<<"proc "<<MyGlobals::_Rank
                <<" : METISGraph::partGraph METIS_PartGraph Recursive or Kway"<<endl;
          if (options_string != "k")
            METIS_PartGraphRecursive(&n, xadj, adjncy, vwgt, adjwgt, &wgtflag,
                                     &base, &nparts, options, &edgecut, partition);
          else
            METIS_PartGraphKway(&n, xadj, adjncy, vwgt, adjwgt, &wgtflag,
                                &base, &nparts, options, &edgecut, partition);
        }
    }
  else
    {
      for (int i=0; i<n; i++) partition[i]=0;
    }
  
  vector<int> index(n+1);
  vector<int> value(n);
  index[0]=0;
  if (ran.size()>0 && MyGlobals::_Atomize==0) //there is randomize
    {
      if (MyGlobals::_Is0verbose>100) cout<<"randomize"<<endl;
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
  delete[]partition;

  //creating a skylinearray with no copy of the index and partition array
  //the fifth argument true specifies that only the pointers are passed 
  //to the object
  
  _partition = new MEDPARTITIONER::SkyLineArray(index,value);
}

