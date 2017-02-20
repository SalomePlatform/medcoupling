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

#include "MEDPARTITIONER_ParaDomainSelector.hxx"
#include "MEDPARTITIONER_UserGraph.hxx"
#include "MEDPARTITIONER_Utils.hxx"

#include "MEDCouplingUMesh.hxx"
#include "MEDCouplingSkyLineArray.hxx"

#include <iostream>
#include <numeric>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*!
 * \brief Constructor. Find out my rank and world size
 */
MEDPARTITIONER::ParaDomainSelector::ParaDomainSelector(bool mesure_memory)
  :_rank(0),_world_size(1), _nb_result_domains(-1), _init_time(0.0),
   _mesure_memory(mesure_memory), _init_memory(0), _max_memory(0)
{
#ifdef HAVE_MPI
  if (MyGlobals::_Rank==-1)
    {
      MPI_Init(0,0);  //do once only
      MPI_Comm_size(MPI_COMM_WORLD,&_world_size) ;
      MPI_Comm_rank(MPI_COMM_WORLD,&_rank) ;
    }
  else
    {
      _world_size=MyGlobals::_World_Size;
      _rank=MyGlobals::_Rank;
    }
  _init_time = MPI_Wtime();
#else
  //sequential : no MPI
  _world_size=1;
  _rank=0;
  if (MyGlobals::_Verbose>10)
    std::cout << "WARNING : ParaDomainSelector contructor without parallel_mode World_Size=1 by default" << std::endl;
#endif
  MyGlobals::_World_Size=_world_size;
  MyGlobals::_Rank=_rank;
  
  if (MyGlobals::_Verbose>200) std::cout << "proc " << MyGlobals::_Rank << " of " << MyGlobals::_World_Size << std::endl;
  evaluateMemory();
}

MEDPARTITIONER::ParaDomainSelector::~ParaDomainSelector()
{
}

/*!
 * \brief Return true if is running on different hosts
 */
bool MEDPARTITIONER::ParaDomainSelector::isOnDifferentHosts() const
{
  evaluateMemory();
  if ( _world_size < 2 ) return false;

#ifdef HAVE_MPI
  char name_here[ MPI_MAX_PROCESSOR_NAME+1 ], name_there[ MPI_MAX_PROCESSOR_NAME+1 ];
  int size;
  MPI_Get_processor_name( name_here, &size);

  int next_proc = (rank() + 1) % nbProcs();
  int prev_proc = (rank() - 1 + nbProcs()) % nbProcs();
  int tag  = 1111111;

  MPI_Status status;
  MPI_Sendrecv((void*)&name_here[0],  MPI_MAX_PROCESSOR_NAME, MPI_CHAR, next_proc, tag,
               (void*)&name_there[0], MPI_MAX_PROCESSOR_NAME, MPI_CHAR, prev_proc, tag,
               MPI_COMM_WORLD, &status);
               
  //bug: (isOnDifferentHosts here and there) is not (isOnDifferentHosts somewhere)
  //return string(name_here) != string(name_there);
  
  int sum_same = -1;
  int same = 1;
  if (std::string(name_here) != std::string(name_there))
    same=0;
  MPI_Allreduce( &same, &sum_same, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
  return (sum_same != nbProcs());
#endif
  return false;
}

/*!
 * \brief Return true if the domain with domainIndex is to be loaded on this proc
 *  \param domainIndex - index of mesh domain
 *  \retval bool - to load or not
 */
bool MEDPARTITIONER::ParaDomainSelector::isMyDomain(int domainIndex) const
{
  evaluateMemory();
  return (_rank == getProcessorID( domainIndex ));
}

/*!
 * \brief Return processor id where the domain with domainIndex resides
 *  \param domainIndex - index of mesh domain
 *  \retval int - processor id
 */
int MEDPARTITIONER::ParaDomainSelector::getProcessorID(int domainIndex) const
{
  evaluateMemory();
  return ( domainIndex % _world_size );
}

/*!
 * \brief Gather info on nb of cell entities on each processor and return total nb.
 *
 * Is called
 * 1) for MED_CELL to know global id shift for domains at graph construction;
 * 2) for sub-entity to know total nb of sub-entities before creating those of joints
 */
void MEDPARTITIONER::ParaDomainSelector::gatherNbOf(const std::vector<MEDCoupling::MEDCouplingUMesh*>& domain_meshes)
{
  evaluateMemory();
  // get nb of elems of each domain mesh
  int nb_domains=domain_meshes.size();
  std::vector<int> nb_elems(nb_domains*2, 0); //NumberOfCells & NumberOfNodes
  for (int i=0; i<nb_domains; ++i)
    if ( domain_meshes[i] )
      {
        nb_elems[i*2] = domain_meshes[i]->getNumberOfCells();
        nb_elems[i*2+1] = domain_meshes[i]->getNumberOfNodes();
      }
  // receive nb of elems from other procs
  std::vector<int> all_nb_elems;
  if (MyGlobals::_World_Size==1)
    {
      all_nb_elems=nb_elems;
    }
  else
    {
#ifdef HAVE_MPI
      all_nb_elems.resize( nb_domains*2 );
      MPI_Allreduce((void*)&nb_elems[0], (void*)&all_nb_elems[0], nb_domains*2, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#else
      throw INTERP_KERNEL::Exception("not(HAVE_MPI) incompatible with MPI_World_Size>1");
#endif
   }
  int total_nb_cells=0, total_nb_nodes=0;
  for (int i=0; i<nb_domains; ++i)
    {
      total_nb_cells+=all_nb_elems[i*2];
      total_nb_nodes+=all_nb_elems[i*2+1];
    }
  
  if (MyGlobals::_Is0verbose>10)
    std::cout << "totalNbCells " << total_nb_cells << " totalNbNodes " << total_nb_nodes << std::endl;
  
  std::vector<int>& cell_shift_by_domain=_cell_shift_by_domain;
  std::vector<int>& node_shift_by_domain=_node_shift_by_domain;
  std::vector<int>& face_shift_by_domain=_face_shift_by_domain;
 
  std::vector< int > ordered_nbs_cell, ordered_nbs_node, domain_order(nb_domains);
  ordered_nbs_cell.push_back(0);
  ordered_nbs_node.push_back(0);
  for (int iproc=0; iproc<nbProcs(); ++iproc)
    for (int idomain=0; idomain<nb_domains; ++idomain)
      if (getProcessorID( idomain )==iproc)
        {
          domain_order[idomain] = ordered_nbs_cell.size() - 1;
          ordered_nbs_cell.push_back( ordered_nbs_cell.back() + all_nb_elems[idomain*2] );
          ordered_nbs_node.push_back( ordered_nbs_node.back() + all_nb_elems[idomain*2+1] );
        }
  cell_shift_by_domain.resize( nb_domains+1, 0 );
  node_shift_by_domain.resize( nb_domains+1, 0 );
  face_shift_by_domain.resize( nb_domains+1, 0 );
  for (int idomain=0; idomain<nb_domains; ++idomain)
    {
      cell_shift_by_domain[ idomain ] = ordered_nbs_cell[ domain_order[ idomain ]];
      node_shift_by_domain[ idomain ] = ordered_nbs_node[ domain_order[ idomain ]];
    }
  cell_shift_by_domain.back() = ordered_nbs_cell.back(); // to know total nb of elements
  node_shift_by_domain.back() = ordered_nbs_node.back(); // to know total nb of elements
  
  if (MyGlobals::_Is0verbose>300)
    {
      std::cout << "proc " << MyGlobals::_Rank << " : cellShiftByDomain ";
      for (int i=0; i<=nb_domains; ++i)
        std::cout << cell_shift_by_domain[i] << "|";
      std::cout << std::endl;
      std::cout << "proc " << MyGlobals::_Rank << " : nodeShiftBy_domain ";
      for (int i=0; i<=nb_domains; ++i)
        std::cout << node_shift_by_domain[i] << "|";
      std::cout << std::endl;
    }
  // fill _nb_vert_of_procs (is Vtxdist)
  _nb_vert_of_procs.resize(_world_size+1);
  _nb_vert_of_procs[0] = 0; // base = 0
  for (int i=0; i<nb_domains; ++i)
    {
      int rankk = getProcessorID(i);
      _nb_vert_of_procs[rankk+1] += all_nb_elems[i*2];
    }
  for (std::size_t i=1; i<_nb_vert_of_procs.size(); ++i)
    _nb_vert_of_procs[i] += _nb_vert_of_procs[i-1]; // to CSR format : cumulated
  
  if (MyGlobals::_Is0verbose>200)
    {
      std::cout << "proc " << MyGlobals::_Rank << " : gatherNbOf : vtxdist is ";
      for (int i = 0; i <= _world_size; ++i)
        std::cout << _nb_vert_of_procs[i] << " ";
      std::cout << std::endl;
    }
  
  evaluateMemory();
  return;
}

/*!
 * \brief Return distribution of the graph vertices among the processors
 *  \retval int* - array containing nb of vertices (=cells) on all processors
 *
 * gatherNbOf() must be called before.
 * The result array is to be used as the first arg of ParMETIS_V3_PartKway() and
 * is freed by ParaDomainSelector.
 */
int *MEDPARTITIONER::ParaDomainSelector::getProcVtxdist() const
{
  evaluateMemory();
  if (_nb_vert_of_procs.empty())
    throw INTERP_KERNEL::Exception("_nb_vert_of_procs not set");
  return const_cast<int*>(& _nb_vert_of_procs[0]);
}

/*!
 * \brief Return nb of cells in domains with lower index.
 *
 * gatherNbOf() must be called before.
 * Result added to local id on given domain gives id in the whole distributed mesh
 */
int MEDPARTITIONER::ParaDomainSelector::getDomainCellShift(int domainIndex) const
{
  evaluateMemory();
  if (_cell_shift_by_domain.empty())
    throw INTERP_KERNEL::Exception("_cell_shift_by_domain not set");
  return _cell_shift_by_domain[domainIndex];
}

int MEDPARTITIONER::ParaDomainSelector::getDomainNodeShift(int domainIndex) const
{
  evaluateMemory();
  if (_node_shift_by_domain.empty())
    throw INTERP_KERNEL::Exception("_node_shift_by_domain not set");
  return _node_shift_by_domain[domainIndex];
}

/*!
 * \brief Return nb of nodes on processors with lower rank.
 *
 * gatherNbOf() must be called before.
 * Result added to global id on this processor gives id in the whole distributed mesh
 */
int MEDPARTITIONER::ParaDomainSelector::getProcNodeShift() const
{
  evaluateMemory();
  if (_nb_vert_of_procs.empty())
    throw INTERP_KERNEL::Exception("_nb_vert_of_procs not set");
  return _nb_vert_of_procs[_rank];
}

/*!
 * \brief Gather graphs from all processors into one
 */
std::auto_ptr<MEDPARTITIONER::Graph> MEDPARTITIONER::ParaDomainSelector::gatherGraph(const Graph* graph) const
{
  Graph* glob_graph = 0;

  evaluateMemory();
#ifdef HAVE_MPI

  // ---------------
  // Gather indices
  // ---------------

  std::vector<int> index_size_of_proc( nbProcs() ); // index sizes - 1
  for ( std::size_t i = 1; i < _nb_vert_of_procs.size(); ++i )
    index_size_of_proc[i-1] = _nb_vert_of_procs[ i ] - _nb_vert_of_procs[ i-1 ];

  int index_size = 1 + _cell_shift_by_domain.back();
  int *graph_index = new int[ index_size ];
  const int *index = graph->getGraph()->getIndex();
  int *proc_index_displacement = const_cast<int*>( & _nb_vert_of_procs[0] );

  MPI_Allgatherv((void*) (index+1),         // send local index except first 0 (or 1)
                 index_size_of_proc[_rank], // index size on this proc
                 MPI_INT,
                 (void*) graph_index,       // receive indices
                 & index_size_of_proc[0],   // index size on each proc
                 proc_index_displacement,   // displacement of each proc index
                 MPI_INT,
                 MPI_COMM_WORLD);
  graph_index[0] = index[0]; // it is not overwritten thanks to proc_index_displacement[0]==1

  // get sizes of graph values on each proc by the got indices of graphs
  std::vector< int > value_size_of_proc( nbProcs() ), proc_value_displacement(1,0);
  for ( int i = 0; i < nbProcs(); ++i )
    {
      if ( index_size_of_proc[i] > 0 )
        value_size_of_proc[i] = graph_index[ proc_index_displacement[ i+1 ]-1 ] - graph_index[0];
      else
        value_size_of_proc[i] = 0;
      proc_value_displacement.push_back( proc_value_displacement.back() + value_size_of_proc[i] );
    }
  
  // update graph_index
  for ( int i = 1; i < nbProcs(); ++i )
    {
      int shift = graph_index[ proc_index_displacement[i]-1 ]-graph_index[0];
      for ( int j = proc_index_displacement[i]; j < proc_index_displacement[i+1]; ++j )
        graph_index[ j ] += shift;
    }
  
  // --------------
  // Gather values
  // --------------

  int value_size = graph_index[ index_size-1 ] - graph_index[ 0 ];
  int *graph_value = new int[ value_size ];
  const int *value = graph->getGraph()->getValues();

  MPI_Allgatherv((void*) value,                // send local value
                 value_size_of_proc[_rank],    // value size on this proc
                 MPI_INT,
                 (void*) graph_value,          // receive values
                 & value_size_of_proc[0],      // value size on each proc
                 & proc_value_displacement[0], // displacement of each proc value
                 MPI_INT,
                 MPI_COMM_WORLD);

  // -----------------
  // Gather partition
  // -----------------

  int * partition = new int[ _cell_shift_by_domain.back() ];
  const int* part = graph->getPart();
  
  MPI_Allgatherv((void*) part,              // send local partition
                 index_size_of_proc[_rank], // index size on this proc
                 MPI_INT,
                 (void*)(partition-1),      // -1 compensates proc_index_displacement[0]==1
                 & index_size_of_proc[0],   // index size on each proc
                 proc_index_displacement,   // displacement of each proc index
                 MPI_INT,
                 MPI_COMM_WORLD);

  // -----------
  // Make graph
  // -----------

  //   MEDCouplingSkyLineArray* array =
  //     new MEDCouplingSkyLineArray( index_size-1, value_size, graph_index, graph_value, true );

  //   glob_graph = new UserGraph( array, partition, index_size-1 );

  evaluateMemory();

  delete [] partition;

#endif // HAVE_MPI

  return std::auto_ptr<Graph>( glob_graph );
}


/*!
 * \brief Set nb of cell/cell pairs in a joint between domains
 */
void MEDPARTITIONER::ParaDomainSelector::setNbCellPairs( int nb_cell_pairs, int dist_domain, int loc_domain )
{
  // This method is needed for further computing global numbers of faces in joint.
  // Store if both domains are on this proc else on one of procs only
  if ( isMyDomain( dist_domain ) || dist_domain < loc_domain )
    {
      if ( _nb_cell_pairs_by_joint.empty() )
        _nb_cell_pairs_by_joint.resize( _nb_result_domains*(_nb_result_domains+1), 0);

      int joint_id = jointId( loc_domain, dist_domain );
      _nb_cell_pairs_by_joint[ joint_id ] = nb_cell_pairs;
    }
  evaluateMemory();
}

//================================================================================
/*!
 * \brief Return nb of cell/cell pairs in a joint between domains on different procs
 */
//================================================================================

int MEDPARTITIONER::ParaDomainSelector::getNbCellPairs( int dist_domain, int loc_domain ) const
{
  evaluateMemory();

  int joint_id = jointId( loc_domain, dist_domain );
  return _nb_cell_pairs_by_joint[ joint_id ];
}

//================================================================================
/*!
 * \brief Gather size of each joint
 */
//================================================================================

void MEDPARTITIONER::ParaDomainSelector::gatherNbCellPairs()
{
  if ( _nb_cell_pairs_by_joint.empty() )
    _nb_cell_pairs_by_joint.resize( _nb_result_domains*(_nb_result_domains+1), 0);
  evaluateMemory();

  std::vector< int > send_buf = _nb_cell_pairs_by_joint;
#ifdef HAVE_MPI
  MPI_Allreduce((void*)&send_buf[0],
                (void*)&_nb_cell_pairs_by_joint[0],
                _nb_cell_pairs_by_joint.size(),
                MPI_INT, MPI_SUM, MPI_COMM_WORLD);
#endif
  // check that the set nbs of cell pairs are correct,
  // namely that each joint is treated on one proc only
  for ( std::size_t j = 0; j < _nb_cell_pairs_by_joint.size(); ++j )
    if ( _nb_cell_pairs_by_joint[j] != send_buf[j] && send_buf[j]>0 )
      throw INTERP_KERNEL::Exception("invalid nb of cell pairs");
}

//================================================================================
/*!
 * \brief Return the first global id of sub-entity for the joint
 */
//================================================================================

int MEDPARTITIONER::ParaDomainSelector::getFisrtGlobalIdOfSubentity( int loc_domain, int dist_domain ) const
{
  // total_nb_faces includes faces existing before creation of joint faces
  // (got in gatherNbOf( MED_FACE )).
  evaluateMemory();

  int total_nb_faces = _face_shift_by_domain.empty() ? 0 : _face_shift_by_domain.back();
  int id = total_nb_faces + 1;

  if ( _nb_cell_pairs_by_joint.empty() )
    throw INTERP_KERNEL::Exception("gatherNbCellPairs() must be called before");
  int joint_id = jointId( loc_domain, dist_domain );
  for ( int j = 0; j < joint_id; ++j )
    id += _nb_cell_pairs_by_joint[ j ];

  return id;
}

//================================================================================
/*!
 * \brief Send-receive local ids of joint faces
 */
//================================================================================

int *MEDPARTITIONER::ParaDomainSelector::exchangeSubentityIds( int loc_domain, int dist_domain,
                                               const std::vector<int>& loc_ids_here ) const
{
  int* loc_ids_dist = new int[ loc_ids_here.size()];
#ifdef HAVE_MPI
  int dest = getProcessorID( dist_domain );
  int tag  = 2002 + jointId( loc_domain, dist_domain );
  MPI_Status status;
  MPI_Sendrecv((void*)&loc_ids_here[0], loc_ids_here.size(), MPI_INT, dest, tag,
               (void*) loc_ids_dist,    loc_ids_here.size(), MPI_INT, dest, tag,
               MPI_COMM_WORLD, &status);  
#endif
  evaluateMemory();

  return loc_ids_dist;
}

//================================================================================
/*!
 * \brief Return identifier for a joint
 */
//================================================================================

int MEDPARTITIONER::ParaDomainSelector::jointId( int local_domain, int distant_domain ) const
{
  evaluateMemory();
  if (_nb_result_domains < 0)
    throw INTERP_KERNEL::Exception("setNbDomains() must be called before");

  if ( local_domain < distant_domain )
    std::swap( local_domain, distant_domain );
  return local_domain * _nb_result_domains + distant_domain;
}


//================================================================================
/*!
 * \brief Return time passed from construction in seconds
 */
//================================================================================

double MEDPARTITIONER::ParaDomainSelector::getPassedTime() const
{
#ifdef HAVE_MPI
  return MPI_Wtime() - _init_time;
#else
  return 0.0;
#endif
}

/*!
  Sends content of \a mesh to processor \a target. To be used with \a recvMesh method.
  \param mesh mesh to be sent
  \param target processor id of the target
*/

void MEDPARTITIONER::ParaDomainSelector::sendMesh(const MEDCoupling::MEDCouplingUMesh& mesh, int target) const
{
#ifndef HAVE_MPI
  throw INTERP_KERNEL::Exception("ParaDomainSelector::sendMesh : incoherent call in non_MPI mode");
#else
  if (MyGlobals::_Verbose>600)
    std::cout << "proc " << _rank << " : sendMesh '" << mesh.getName() << "' size " << mesh.getNumberOfCells() << " to " << target << std::endl;
  // First stage : sending sizes
  // ------------------------------
  std::vector<int> tinyInfoLocal;
  std::vector<std::string> tinyInfoLocalS;
  std::vector<double> tinyInfoLocalD;
  //Getting tiny info of local mesh to allow the distant proc to initialize and allocate
  //the transmitted mesh.
  mesh.getTinySerializationInformation(tinyInfoLocalD,tinyInfoLocal,tinyInfoLocalS);
  tinyInfoLocal.push_back(mesh.getNumberOfCells());
  int tinySize=tinyInfoLocal.size();
  MPI_Send(&tinySize, 1, MPI_INT, target, 1113, MPI_COMM_WORLD);
  MPI_Send(&tinyInfoLocal[0], tinyInfoLocal.size(), MPI_INT, target, 1112, MPI_COMM_WORLD);

  if (mesh.getNumberOfCells()>0) //no sends if empty
    {
      MEDCoupling::DataArrayInt *v1Local=0;
      MEDCoupling::DataArrayDouble *v2Local=0;
      //serialization of local mesh to send data to distant proc.
      mesh.serialize(v1Local,v2Local);
      int nbLocalElems=0;
      int* ptLocal=0;
      if(v1Local) //if empty getNbOfElems() is 1!
        {
          nbLocalElems=v1Local->getNbOfElems(); // if empty be 1!
          ptLocal=v1Local->getPointer();
        }
      MPI_Send(ptLocal, nbLocalElems, MPI_INT, target, 1111, MPI_COMM_WORLD);
      int nbLocalElems2=0;
      double *ptLocal2=0;
      if(v2Local) //if empty be 0!
        {
          nbLocalElems2=v2Local->getNbOfElems();
          ptLocal2=v2Local->getPointer();
        }
      MPI_Send(ptLocal2, nbLocalElems2, MPI_DOUBLE, target, 1110, MPI_COMM_WORLD);
      if(v1Local) v1Local->decrRef();
      if(v2Local) v2Local->decrRef();
    }
#endif
}

/*! Receives messages from proc \a source to fill mesh \a mesh.
  To be used with \a sendMesh method.
  \param mesh  pointer to mesh that is filled
  \param source processor id of the incoming messages
*/
void MEDPARTITIONER::ParaDomainSelector::recvMesh(MEDCoupling::MEDCouplingUMesh*& mesh, int source)const
{
#ifndef HAVE_MPI
  throw INTERP_KERNEL::Exception("ParaDomainSelector::recvMesh : incoherent call in non_MPI mode");
#else
  // First stage : exchanging sizes
  // ------------------------------
  std::vector<int> tinyInfoDistant;
  std::vector<std::string> tinyInfoLocalS;
  std::vector<double> tinyInfoDistantD(1);
  //Getting tiny info of local mesh to allow the distant proc to initialize and allocate
  //the transmitted mesh.
  MPI_Status status; 
  int tinyVecSize;
  MPI_Recv(&tinyVecSize, 1, MPI_INT, source, 1113, MPI_COMM_WORLD, &status);
  tinyInfoDistant.resize(tinyVecSize);
  std::fill(tinyInfoDistant.begin(),tinyInfoDistant.end(),0);

  MPI_Recv(&tinyInfoDistant[0], tinyVecSize, MPI_INT,source,1112,MPI_COMM_WORLD, &status);
  //there was tinyInfoLocal.push_back(mesh.getNumberOfCells());
  int NumberOfCells=tinyInfoDistant[tinyVecSize-1];
  if (NumberOfCells>0)
    {
      MEDCoupling::DataArrayInt *v1Distant=MEDCoupling::DataArrayInt::New();
      MEDCoupling::DataArrayDouble *v2Distant=MEDCoupling::DataArrayDouble::New();
      //Building the right instance of copy of distant mesh.
      MEDCoupling::MEDCouplingPointSet *distant_mesh_tmp=
        MEDCoupling::MEDCouplingPointSet::BuildInstanceFromMeshType(
                                                                   (MEDCoupling::MEDCouplingMeshType) tinyInfoDistant[0]);
      std::vector<std::string> unusedTinyDistantSts;
      mesh=dynamic_cast<MEDCoupling::MEDCouplingUMesh*> (distant_mesh_tmp);
 
      mesh->resizeForUnserialization(tinyInfoDistant,v1Distant,v2Distant,unusedTinyDistantSts);
      int nbDistElem=0;
      int *ptDist=0;
      if(v1Distant)
        {
          nbDistElem=v1Distant->getNbOfElems();
          ptDist=v1Distant->getPointer();
        }
      MPI_Recv(ptDist, nbDistElem, MPI_INT, source,1111, MPI_COMM_WORLD, &status);
      double *ptDist2=0;
      nbDistElem=0;
      if(v2Distant)
        {
          nbDistElem=v2Distant->getNbOfElems();
          ptDist2=v2Distant->getPointer();
        }
      MPI_Recv(ptDist2, nbDistElem, MPI_DOUBLE,source, 1110, MPI_COMM_WORLD, &status);
      //finish unserialization
      mesh->unserialization(tinyInfoDistantD,tinyInfoDistant,v1Distant,v2Distant,unusedTinyDistantSts);
      if(v1Distant) v1Distant->decrRef();
      if(v2Distant) v2Distant->decrRef();
    }
  else
    {
      mesh=CreateEmptyMEDCouplingUMesh();
    }
  if (MyGlobals::_Verbose>600)
    std::cout << "proc " << _rank << " : recvMesh '" << mesh->getName() << "' size " << mesh->getNumberOfCells() << " from " << source << std::endl;
#endif
}

#if !defined WIN32 && !defined __APPLE__
#include <sys/sysinfo.h>
#endif

/*!
 * \brief Evaluate current memory usage and return the maximal one in KB
 */
int MEDPARTITIONER::ParaDomainSelector::evaluateMemory() const
{
  if ( _mesure_memory )
    {
      int used_memory = 0;
#if !defined WIN32 && !defined __APPLE__
      struct sysinfo si;
      int err = sysinfo( &si );
      if ( !err )
        used_memory = (( si.totalram - si.freeram + si.totalswap - si.freeswap ) * si.mem_unit ) / 1024;
#endif
      if ( used_memory > _max_memory )
        _max_memory = used_memory;

      if ( !_init_memory )
        _init_memory = used_memory;
    }
  return _max_memory - _init_memory;
}
