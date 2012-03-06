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

#include "MEDPARTITIONER_MeshCollection.hxx"
#include "MEDPARTITIONER_Topology.hxx"
#include "MEDPARTITIONER_Graph.hxx"
#include "MEDPARTITIONER_ParallelTopology.hxx"
#include "MEDPARTITIONER_ConnectZone.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "MEDCouplingUMesh.hxx"
#include "InterpKernelHashMap.hxx"

#include <set>
#include <map>
#include <vector>
#include <iostream>

#ifdef HAVE_MPI2
#include <mpi.h>
#endif

using namespace MEDPARTITIONER;
using namespace std;

ParallelTopology::ParallelTopology():_nb_domain(0),_mesh_dimension(0)
{
}

//constructing topology according to mesh collection without global numerotation (use setGlobalNumerotation later)
ParallelTopology::ParallelTopology(const vector<ParaMEDMEM::MEDCouplingUMesh*>& meshes)
{
  _nb_domain=meshes.size();
  _nb_cells.resize(_nb_domain);
  _nb_nodes.resize(_nb_domain);
  //  _nb_faces.resize(_nb_domain);
  
  if (MyGlobals::_Is0verbose>100) cout<<"new ParallelTopology\n";
  _loc_to_glob.resize(0);      //precaution, need gatherNbOf() setGlobalNumerotation()
  _node_loc_to_glob.resize(0); //precaution, need gatherNbOf() setGlobalNumerotation()
  //_face_loc_to_glob.resize(_nb_domain);
  _mesh_dimension = -1;
  bool parallel_mode = false;
  for (int idomain=0; !parallel_mode && idomain<_nb_domain; idomain++)
    parallel_mode = (!meshes[idomain]);

  if (MyGlobals::_Is0verbose>20 && !parallel_mode) cout<<"WARNING : ParallelTopology contructor without parallel_mode"<<endl;
  for (int idomain=0; idomain<_nb_domain; idomain++)
    {
      if ( !meshes[idomain] ) continue;
      if (_mesh_dimension==-1)
        {
          _mesh_dimension = meshes[idomain]->getMeshDimension();
        }
      else
        {
          if (_mesh_dimension!=meshes[idomain]->getMeshDimension())
            throw INTERP_KERNEL::Exception(LOCALIZED("meshes dimensions incompatible"));
        }
      _nb_cells[idomain]=meshes[idomain]->getNumberOfCells();
      _nb_nodes[idomain]=meshes[idomain]->getNumberOfNodes();
      //note: in parallel mode _nb_cells and _nb_nodes are not complete now, needs gatherNbOf()
    }
}

//constructing _loc_to_glob etc by default, needs gatherNbOf() done
void ParallelTopology::setGlobalNumerotationDefault(ParaDomainSelector* domainSelector)
{
  if (MyGlobals::_Is0verbose>100) cout<<"setGlobalNumerotationDefault on "<<_nb_domain<<" domains\n";
  if (_loc_to_glob.size()!=0) throw INTERP_KERNEL::Exception(LOCALIZED("a global numerotation is done yet"));
  _loc_to_glob.resize(_nb_domain);
  _node_loc_to_glob.resize(_nb_domain);
  
  //warning because _nb_cells[idomain] is 0 if not my domain(s)
  //we set loc_to_glob etc.. only for my domain(s)
  if (MyGlobals::_Is0verbose>500) cout<<"(c)idomain|ilocalCell|iglobalCell"<<endl;
  for (int idomain=0; idomain<_nb_domain; idomain++)
    {
      _loc_to_glob[idomain].resize(_nb_cells[idomain]);
      int domainCellShift=domainSelector->getDomainCellShift(idomain);
      for (int i=0; i<_nb_cells[idomain]; i++)
        {
          int global=domainCellShift+i ;
          _glob_to_loc.insert(make_pair(global,make_pair(idomain,i)));
          _loc_to_glob[idomain][i]=global;
          if (MyGlobals::_Verbose>500) cout<<"c"<<idomain<<"|"<<i<<"|"<<global<<" ";
        }
    }
  if (MyGlobals::_Verbose>500) MPI_Barrier(MPI_COMM_WORLD);
  if (MyGlobals::_Is0verbose>500) cout<<endl;
  
  if (MyGlobals::_Is0verbose>500) cout<<"(n)idomain|ilocalNode|iglobalNode"<<endl;
  for (int idomain=0; idomain<_nb_domain; idomain++)
    {
      _node_loc_to_glob[idomain].resize(_nb_nodes[idomain]);
      int domainNodeShift=domainSelector->getDomainNodeShift(idomain);
      for (int i=0; i<_nb_nodes[idomain]; i++)
        {
          int global=domainNodeShift+i ;
          _node_glob_to_loc.insert(make_pair(global,make_pair(idomain,i)));
          _node_loc_to_glob[idomain][i]=global;
          if (MyGlobals::_Verbose>500) cout<<"n"<<idomain<<"|"<<i<<"|"<<global<<" ";
        }
    }
  if (MyGlobals::_Verbose>500) MPI_Barrier(MPI_COMM_WORLD);
  if (MyGlobals::_Is0verbose>500) cout<<endl;
  
  _nb_total_cells=domainSelector->getNbTotalCells();
  _nb_total_nodes=domainSelector->getNbTotalNodes();
  _nb_total_faces=domainSelector->getNbTotalFaces();
  if (MyGlobals::_Is0verbose>200) cout<<"globalNumerotation default done meshDimension "<<_mesh_dimension<<" nbTotalCells "<<_nb_total_cells<<" nbTotalNodes "<<_nb_total_nodes<<endl;
}

//constructing topology according to mesh collection
ParallelTopology::ParallelTopology(const vector<ParaMEDMEM::MEDCouplingUMesh*>& meshes, 
                                   const vector<MEDPARTITIONER::ConnectZone*>& cz,
                                   vector<int*>& cellglobal,
                                   vector<int*>& nodeglobal,
                                   vector<int*>& faceglobal)
{
  _nb_domain=meshes.size();
  int index_global=0;
  int index_node_global=0;
  int index_face_global=0;

  _nb_cells.resize(_nb_domain);
  _nb_nodes.resize(_nb_domain);
  //  _nb_faces.resize(_nb_domain);
  
  _loc_to_glob.resize(_nb_domain);
  _node_loc_to_glob.resize(_nb_domain);
  //  _face_loc_to_glob.resize(_nb_domain);

  bool parallel_mode = false;
  for (int idomain=0; !parallel_mode && idomain<_nb_domain; idomain++)
    parallel_mode = (!meshes[idomain]);

  for (int idomain=0; idomain<_nb_domain; idomain++)
    {
      if ( !meshes[idomain] ) continue;
      _mesh_dimension = meshes[idomain]->getMeshDimension();
    
      //creating cell maps
      _nb_cells[idomain]=meshes[idomain]->getNumberOfCells();
      //    cout << "Nb cells (domain "<<idomain<<") = "<<_nb_cells[idomain];
      _loc_to_glob[idomain].resize(_nb_cells[idomain]);

      if (cellglobal[idomain]==0 || parallel_mode)
        {
          //int cellDomainShift=_cell_shift_by_domain[idomain];
          //creating global numbering from scratch
          for (int i=0; i<_nb_cells[idomain]; i++)
            {
              int global=i ;//cellDomainShift+i;
              _glob_to_loc.insert(make_pair(global,make_pair(idomain,i)));
              _loc_to_glob[idomain][i]=global;
              //cvw cout<<idomain<<"|"<<i<<"|"<<global<<" ";
              index_global++;
            }
        }
      //using global numbering coming from a previous numbering
      else
        {
          for (int i=0; i<_nb_cells[idomain]; i++)
            {
              int global=cellglobal[idomain][i];
              _glob_to_loc.insert(make_pair(global,make_pair(idomain,i)));
              //_loc_to_glob[make_pair(idomain,i+1)]=global;
              _loc_to_glob[idomain][i]=global;
              index_global++;
              //        cout<<"glob:"<<global<<" --> ("<<idomain<<","<<i+1<<")"<<endl;
            }
        }

      //cas sequentiel
      if (_nb_domain==1)
        {
          _nb_total_cells=index_global;
          _nb_cells[0]=index_global;
          _node_loc_to_glob[idomain].resize(meshes[idomain]->getNumberOfNodes());
          for (int i=0; i<meshes[idomain]->getNumberOfNodes(); i++)
            {
              _node_glob_to_loc.insert(make_pair(i,make_pair(0,i)));
              _node_loc_to_glob[0][i]=i;
            }
          _nb_total_nodes=meshes[idomain]->getNumberOfNodes();   
          _nb_nodes[0]=_nb_total_nodes; 
          return;
        }

      //creating node maps
      _nb_nodes[idomain]=meshes[idomain]->getNumberOfNodes();
      INTERP_KERNEL::HashMap <int,pair<int,int> > local2distant;
      _node_loc_to_glob[idomain].resize(_nb_nodes[idomain]);
      for (int icz=0; icz<cz.size(); icz++)
        {
          if (cz[icz]->getLocalDomainNumber() == idomain && 
              cz[icz]->getLocalDomainNumber()>cz[icz]->getDistantDomainNumber())
            {
              int nb_node= cz[icz]->getNodeNumber();
              const int* node_corresp=cz[icz]->getNodeCorrespValue();
              int distant_ip = cz[icz]->getDistantDomainNumber();
              for (int i=0; i< nb_node; i++)
                {
                  int local= node_corresp[i*2];
                  int distant = node_corresp[i*2+1];
                  local2distant.insert(make_pair(local, make_pair(distant_ip,distant)));    
                }
            }
        }
      // setting mappings for all nodes
      if (nodeglobal[idomain]==0)
        {
          for (int inode=0; inode<_nb_nodes[idomain]; inode++)
            {
              if (local2distant.find(inode)==local2distant.end())
                {
                  index_node_global++;
                  _node_glob_to_loc.insert(make_pair(index_node_global,make_pair(idomain,inode)));
                  //_node_loc_to_glob[make_pair(idomain,inode+1)]=index_node_global;
                  _node_loc_to_glob[idomain][inode]=index_node_global;
                }   
              else
                {
                  int ip = (local2distant.find(inode)->second).first;
                  int distant = (local2distant.find(inode)->second).second;
                  int global_number=_loc_to_glob[ip][distant];
                  _node_glob_to_loc.insert(make_pair(global_number,make_pair(idomain,inode)));
                  _node_loc_to_glob[idomain][inode]=global_number;
                } 
            }
        }
      //using former node numbering
      else
        {
          for (int inode=0; inode<_nb_nodes[idomain]; inode++)
            {
              int global_number=nodeglobal[idomain][inode];
              _node_glob_to_loc.insert(make_pair(global_number,make_pair(idomain,inode)));
              _node_loc_to_glob[idomain][inode]=global_number;
            }
        }
    }

  _nb_total_cells=index_global;
  _nb_total_nodes=index_node_global;   
  _nb_total_faces=index_face_global;
}


//constructing ParallelTopology from an old topology and a graph
ParallelTopology::ParallelTopology(Graph* graph, Topology* oldTopology, int nb_domain, int mesh_dimension)
{

  _nb_domain=nb_domain;
  //cvw !!whatisit! _nb_cells=graph->nbVertices();
  _mesh_dimension=mesh_dimension;
  
  if (MyGlobals::_Verbose>200)
    cout<<"proc "<<MyGlobals::_Rank<<" : new topology oldNbDomain "<<
      oldTopology->nbDomain()<<" newNbDomain "<<_nb_domain<<endl;
  _nb_cells.resize(_nb_domain,0);
  _nb_nodes.resize(_nb_domain,0);
  _nb_faces.resize(_nb_domain,0);

  _loc_to_glob.resize(_nb_domain);
  _node_loc_to_glob.resize(_nb_domain);
  _face_loc_to_glob.resize(_nb_domain);
  
  const int* part=graph->getPart(); //all cells for this proc (may be more domains)
  _nb_total_cells=graph->nbVertices(); //all cells for this proc (may be more domains)
  if (MyGlobals::_Verbose>300)
    cout<<"proc "<<MyGlobals::_Rank<<" : topology from partition, nbTotalCells "<<_nb_total_cells<<endl;
  
  int icellProc=0; //all cells of my domains are concatenated in part
  for (int iold=0; iold<oldTopology->nbDomain(); iold++)
    {
      int ioldNbCell=oldTopology->getCellNumber(iold);
      //cout<<"proc "<<MyGlobals::_Rank<<" : cell number old domain "<<iold<<" : "<<ioldNbCell<<endl;
      //if not my old domains getCellNumber is 0
      std::vector<int> globalids(ioldNbCell);
      oldTopology->getCellList(iold, &globalids[0]); //unique global numerotation
      for (int icell=0; icell<ioldNbCell; icell++)
        {
          int idomain=part[icellProc];
          _nb_cells[idomain]++;
          icellProc++;
          int iGlobalCell=globalids[icell];
          _loc_to_glob[idomain].push_back(iGlobalCell);
          _glob_to_loc.insert(make_pair(iGlobalCell, make_pair(idomain, _nb_cells[idomain])));
        }
    }

  if (MyGlobals::_Verbose>300)
    for (int idomain=0; idomain<_nb_domain; idomain++)
      cout<<"proc "<<MyGlobals::_Rank<<" : nbCells in new domain "<<idomain<<" : "<<_nb_cells[idomain]<<endl; 
}

ParallelTopology::~ParallelTopology()
{
} 

/*!Converts a list of global node numbers
 * to a distributed array with local cell numbers.
 * 
 * If a node in the list is represented on several domains,
 * only the first value is returned
 * */
void ParallelTopology::convertGlobalNodeList(const int* node_list, int nbnode, int* local, int* ip)
{
  if (_node_glob_to_loc.empty()) 
    throw INTERP_KERNEL::Exception(LOCALIZED("Node mapping has not yet been built"));
  for (int i=0; i< nbnode; i++)
    {
      pair<int,int> local_node = _node_glob_to_loc.find(node_list[i])->second;
      ip[i]=local_node.first;
      local[i]=local_node.second;
    }
}

/*!Converts a list of global node numbers on domain ip
 * to a distributed array with local cell numbers.
 * 
 * If a node in the list is represented on several domains,
 * only the value with domain ip is returned
 * 
 * */
void ParallelTopology::convertGlobalNodeList(const int* node_list, int nbnode, int* local, int ip)
{
  if (_node_glob_to_loc.empty()) 
    throw INTERP_KERNEL::Exception(LOCALIZED("Node mapping has not yet been built"));

  for (int i=0; i< nbnode; i++)
    {
      typedef INTERP_KERNEL::HashMultiMap<int,pair<int,int> >::iterator mmiter;
      pair<mmiter,mmiter> range=_node_glob_to_loc.equal_range(node_list[i]);
      for (mmiter it=range.first; it !=range.second; it++)
        { 
          int ipfound=(it->second).first;
          if (ipfound==ip)
            local[i]=(it->second).second;
        }
    }
} 

/*!Converts a list of global node numbers
 * to a distributed array with local cell numbers.
 * 
 * If a node in the list is represented on several domains,
 * all the values are put in the array
 * */
void ParallelTopology::convertGlobalNodeListWithTwins(const int* node_list, int nbnode, int*& local, int*& ip,int*& full_array, int& size)
{
  if (_node_glob_to_loc.empty()) 
    throw INTERP_KERNEL::Exception(LOCALIZED("Node mapping has not yet been built"));

  size=0;
  for (int i=0; i< nbnode; i++)
    {
      int count= _node_glob_to_loc.count(node_list[i]);
      //      if (count > 1) 
      //        cout << "noeud " << node_list[i]<< " doublon d'ordre " << count<<endl;
      size+=count;
    }
  int index=0;
  ip=new int[size];
  local=new int[size];
  full_array=new int[size];
  for (int i=0; i< nbnode; i++)
    {
      typedef INTERP_KERNEL::HashMultiMap<int,pair<int,int> >::iterator mmiter;
      pair<mmiter,mmiter> range=_node_glob_to_loc.equal_range(node_list[i]);
      for (mmiter it=range.first; it !=range.second; it++)
        { 
          ip[index]=(it->second).first;
          local[index]=(it->second).second;
          full_array [index]=node_list[i];
          index++;
        }

    }
}

/*!Converts a list of global face numbers
 * to a distributed array with local face numbers.
 * 
 * If a face in the list is represented on several domains,
 * all the values are put in the array
 * */
void ParallelTopology::convertGlobalFaceListWithTwins(const int* face_list, int nbface, int*& local, int*& ip, int*& full_array,int& size)
{
  size=0;
  for (int i=0; i< nbface; i++)
    {
      //int count = _face_glob_to_loc.count(face_list[i]);
      //if (count >1) MESSAGE_MED("face en doublon "<<face_list[i]);
      size+= _face_glob_to_loc.count(face_list[i]);
    }
  int index=0;
  ip=new int[size];
  local=new int[size];
  full_array=new int[size];
  for (int i=0; i< nbface; i++)
    {
      typedef INTERP_KERNEL::HashMultiMap<int,pair<int,int> >::iterator mmiter;
      pair<mmiter,mmiter> range=_face_glob_to_loc.equal_range(face_list[i]);
      for (mmiter it=range.first; it !=range.second; it++)
        { 
          ip[index]=(it->second).first;
          local[index]=(it->second).second;
          full_array[index]=face_list[i];
          index++;
        }

    }
}

//!converts a list of global cell numbers
//!to a distributed array with local cell numbers 
void ParallelTopology::convertGlobalCellList(const int* cell_list, int nbcell, int* local, int* ip)
{
  for (int i=0; i<nbcell; i++)
    {
      //cvw INTERP_KERNEL::HashMap<int, pair<int,int> >::const_iterator iter = _glob_to_loc.find(cell_list[i]);
      INTERP_KERNEL::HashMap<int, pair<int,int> >::const_iterator iter = _glob_to_loc.find(cell_list[i]);
      if (iter == _glob_to_loc.end())
        {
          cerr<<"proc "<<MyGlobals::_Rank<<" : KO cell_list["<<i<<"] : "<<cell_list[i]<<endl;
          throw INTERP_KERNEL::Exception(LOCALIZED("ParallelTopology::convertGlobalCellList : Cell not found"));
        }
      else
        {
          ip[i]=(iter->second).first;     //no domain
          local[i]=(iter->second).second; //no local cell
          //cout<<"proc "<<MyGlobals::_Rank<<" : OK cell_list["<<i<<"] : "<<cell_list[i]<<" "<<ip[i]<<" "<<local[i]<<endl;
        }
    }
}

/*!Converts a list of global face numbers
 * to a distributed array with local face numbers
 */ 
void ParallelTopology::convertGlobalFaceList(const int* face_list, int nbface, int* local, int* ip)
{
  for (int i=0; i< nbface; i++)
    {
      INTERP_KERNEL::HashMap<int, pair<int,int> >::const_iterator iter = _face_glob_to_loc.find(face_list[i]);
      if (iter == _face_glob_to_loc.end())
        {
          throw INTERP_KERNEL::Exception(LOCALIZED("ParallelTopology::convertGlobalFaceList : Face not found"));
        }
      ip[i]=(iter->second).first;
      local[i]=(iter->second).second;
      //    cout << " in convertGlobalFAceList face global "<<face_list[i]<<" -> ("<<ip[i]<<","<<local[i]<<")"<<endl;
    }
}

/*!Converts a list of global node numbers on domain ip
 * to a distributed array with local cell numbers.
 * 
 * If a node in the list is represented on several domains,
 * only the value with domain ip is returned
 * 
 */
void ParallelTopology::convertGlobalFaceList(const int* face_list, int nbface, int* local, int ip)
{
  for (int i=0; i< nbface; i++)
    {
      typedef INTERP_KERNEL::HashMultiMap<int,pair<int,int> >::iterator mmiter;
      pair<mmiter,mmiter> range=_face_glob_to_loc.equal_range(face_list[i]);
      for (mmiter it=range.first; it !=range.second; it++)
        { 
          int ipfound=(it->second).first;
          if (ipfound==ip)
            local[i]=(it->second).second; 

        }
    }
} 

//replacing a table of global numbering with a table with local numberings
// type_connectivity contains global connectivity for each type in input
// type_connectivity contains local connectivity for each type in output
void ParallelTopology::convertToLocal2ndVersion(int* nodes, int nbnodes, int idomain)
{
  for (int inode=0; inode<nbnodes; inode++)
    {
      //      cout <<" inode :"<<inode<< " global = "<<type_connectivity[type][inode];
      int global = nodes[inode];
      typedef INTERP_KERNEL::HashMultiMap<int,pair<int,int> >::iterator mmiter;
      pair<mmiter,mmiter> range=_node_glob_to_loc.equal_range(global);
      for (mmiter it=range.first; it !=range.second; it++)
        {
          if ((it->second).first==idomain)
            nodes[inode]=(it->second).second;
        }
    }
}

/*!
 * \brief Return max global face number
 */
int ParallelTopology::getMaxGlobalFace() const
{
  int max = 0;
  TGlob2LocsMap::const_iterator g_l_l = _face_glob_to_loc.begin();
  for ( ; g_l_l != _face_glob_to_loc.end(); ++g_l_l )
    if ( g_l_l->first > max )
      max = g_l_l->first;
  return max;
}

