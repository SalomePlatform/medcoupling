// Copyright (C) 2007-2022  CEA/DEN, EDF R&D
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

#include "BlockTopology.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingCMesh.hxx"
#include "CommInterface.hxx"
#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "ComponentTopology.hxx"
#include "InterpKernelUtilities.hxx"

#include <vector>
#include <algorithm>
#include <utility>
#include <iostream>

using namespace std;

namespace MEDCoupling
{
  /*!
   * Default ctor.
   */
  BlockTopology::BlockTopology() :
    _dimension(0), _nb_procs_per_dim(0),
    _local_array_indices(0), _cycle_type(0),
    _proc_group(nullptr),_nb_elems(0),
    _owns_processor_group(false)
  {}

  /*!
   * Constructor of a block topology from a grid. 
   * This preliminary version simply splits along the first axis
   * instead of making the best choice with respect to the 
   * values of the different axes. 
   */
  BlockTopology::BlockTopology(const ProcessorGroup& group, MEDCouplingCMesh *grid):
    _dimension(grid->getSpaceDimension()), _proc_group(&group), _owns_processor_group(false)
  {
    vector <mcIdType> axis_length(_dimension);
    _nb_elems=1;
    for (int idim=0; idim <_dimension; idim++)
      {
        DataArrayDouble *arr=grid->getCoordsAt(idim);
        axis_length[idim]=arr->getNbOfElems();
        _nb_elems*=axis_length[idim];
      }  
    //default splitting along 1st dimension
    _local_array_indices.resize(_dimension);
    _nb_procs_per_dim.resize(_dimension);
  
    _local_array_indices[0].resize(_proc_group->size()+1);
    _local_array_indices[0][0]=0;
    _nb_procs_per_dim[0]=_proc_group->size();
  
    for (int i=1; i<=_proc_group->size(); i++)
      {
        _local_array_indices[0][i]=_local_array_indices[0][i-1]+
          axis_length[0]/_proc_group->size();
        if (i<= axis_length[0]%_proc_group->size())
          _local_array_indices[0][i]+=1;
      }
    for (int i=1; i<_dimension; i++)
      {
        _local_array_indices[i].resize(2);
        _local_array_indices[i][0]=0;
        _local_array_indices[i][1]=axis_length[i];
        _nb_procs_per_dim[i]=1;
      }
    _cycle_type.resize(_dimension);
    for (int i=0; i<_dimension; i++)
      _cycle_type[i]=MEDCoupling::Block;  
  }

  /*!
   * Creation of a block topology by composing 
   * a geometrical topology and a component topology.
   * This constructor is intended for creating fields 
   * for which the parallel distribution is made on the
   * components of the field rather than on the geometrical 
   * partitioning of the underlying mesh.
   * 
   */ 
  BlockTopology::BlockTopology(const BlockTopology& geom_topo, const ComponentTopology& comp_topo):_owns_processor_group(false)
  {
    // so far, the block topology can only be created if the proc group 
    // is either on geom_topo or on comp_topo
    if (geom_topo.getProcGroup()->size()>1 && comp_topo.nbBlocks()>1)
      throw INTERP_KERNEL::Exception(LOCALIZED("BlockTopology cannot yet be constructed with both complex geo and components topology"));

    if (comp_topo.nbComponents()==1)
      {
        *this=geom_topo;
        return;
      }
    else
      {
        _dimension = geom_topo.getDimension()+1;
        if (comp_topo.nbBlocks()>1)
          _proc_group=comp_topo.getProcGroup();
        else
          _proc_group=geom_topo.getProcGroup();
        _local_array_indices=geom_topo._local_array_indices;
        vector<int> comp_indices = *(comp_topo.getBlockIndices());
        _local_array_indices.emplace_back( comp_indices.begin(), comp_indices.end() );
        _nb_procs_per_dim=geom_topo._nb_procs_per_dim;
        _nb_procs_per_dim.push_back(comp_topo.nbBlocks());
        _cycle_type=geom_topo._cycle_type;
        _cycle_type.push_back(Block);
        _nb_elems=geom_topo.getNbElements()*comp_topo.nbComponents();
      }  
  }

  /*! Constructor for creating a one-dimensional
   * topology from a processor group and a local 
   * number of elements on each processor
   * 
   * The function must be called only by the processors belonging
   * to group \a group. Calling it from a processor not belonging
   * to \a group will cause an MPI error, while calling from a subset
   * of \a group will result in a deadlock. 
   */
  BlockTopology::BlockTopology(const ProcessorGroup& group, mcIdType nb_elem):_dimension(1),_proc_group(&group),_owns_processor_group(false)
  {
    mcIdType* nbelems_per_proc = new mcIdType[group.size()];
    const MPIProcessorGroup* mpi_group=dynamic_cast<const MPIProcessorGroup*>(_proc_group);
    const MPI_Comm* comm=mpi_group->getComm();
    mcIdType nbtemp=nb_elem;
    mpi_group->getCommInterface().allGather(&nbtemp, 1, MPI_ID_TYPE, 
                                            nbelems_per_proc, 1, MPI_ID_TYPE, 
                                            *comm);
    _nb_elems=0;
  
    //splitting along only dimension
    _local_array_indices.resize(1);
    _nb_procs_per_dim.resize(1);  
          
    _local_array_indices[0].resize(_proc_group->size()+1);
    _local_array_indices[0][0]=0;
    _nb_procs_per_dim[0]=_proc_group->size();
  
    for (int i=1; i<=_proc_group->size(); i++)
      {
        _local_array_indices[0][i]=_local_array_indices[0][i-1]+
          nbelems_per_proc[i-1];
        _nb_elems+=nbelems_per_proc[i-1];
      }
    _cycle_type.resize(1);
    _cycle_type[0]=MEDCoupling::Block;
    delete[] nbelems_per_proc;
  }

  BlockTopology::~BlockTopology()
  {
    release();
  }

  /** Destructor involves MPI operations: make sure this is accessible from a proper
   * method for Python wrapping.
   */
  void BlockTopology::release()
  {
    if (_owns_processor_group)
      delete _proc_group;
    _proc_group = nullptr;
  }

  //!converts a pair <subdomainid,local> to a global number
  std::pair<int,mcIdType> BlockTopology::globalToLocal(const mcIdType global) const
  {
    int subdomain_id=0;
    mcIdType position=global;
    mcIdType size=_nb_elems;
    std::size_t size_procs=_proc_group->size();
    mcIdType increment=size;
    vector<mcIdType>axis_position(_dimension);
    vector<mcIdType>axis_offset(_dimension);
    for (int idim=0; idim<_dimension; idim++)
      {
        std::size_t axis_size=_local_array_indices[idim].size()-1;
        mcIdType axis_nb_elem=_local_array_indices[idim][axis_size];
        increment=increment/axis_nb_elem;
        int proc_increment = (int)(size_procs/axis_size);
        mcIdType axis_pos=position/increment;
        position=position%increment;
        int iaxis=1;
        while (_local_array_indices[idim][iaxis]<=axis_pos)
          {
            subdomain_id+=proc_increment;
            iaxis++;
          }
        axis_position[idim]=axis_pos-_local_array_indices[idim][iaxis-1];
        axis_offset[idim]=iaxis;
      }
    mcIdType local=0;
    mcIdType local_increment=1;
    for (int idim=_dimension-1; idim>=0; idim--)
      {
        local+=axis_position[idim]*local_increment;
        local_increment*=_local_array_indices[idim][axis_offset[idim]]-_local_array_indices[idim][axis_offset[idim]-1];
      }
    return make_pair(subdomain_id,local);
  }

  //!converts local number to a global number
  mcIdType BlockTopology::localToGlobal(const pair<int,mcIdType> local) const
  {

    std::size_t subdomain_id=local.first;
    mcIdType global=0;
    mcIdType loc=local.second;
    mcIdType increment=_nb_elems;
    std::size_t proc_increment=_proc_group->size();
    mcIdType local_increment=getNbLocalElements();
    for (int idim=0; idim < _dimension; idim++)
      {
        std::size_t axis_size=_local_array_indices[idim].size()-1;
        mcIdType axis_nb_elem=_local_array_indices[idim][axis_size];
        increment=axis_nb_elem==0?0:increment/axis_nb_elem;
        proc_increment = proc_increment/axis_size;
        std::size_t proc_axis=subdomain_id/proc_increment;
        subdomain_id=subdomain_id%proc_increment;
        mcIdType local_axis_nb_elem=_local_array_indices[idim][proc_axis+1]-_local_array_indices[idim][proc_axis];
        local_increment = (local_axis_nb_elem==0)?0:(local_increment/local_axis_nb_elem);
        mcIdType iaxis=((local_increment==0)?0:(loc/local_increment))+_local_array_indices[idim][proc_axis];
        global+=increment*iaxis;
        loc = (local_increment==0)?0:(loc%local_increment);
      }
    return global;
  }

  //Retrieves the local number of elements
  mcIdType BlockTopology::getNbLocalElements()const
  {
    int position=_proc_group->myRank();
    mcIdType nb_elem = 1;
    int increment=1;
    for (int i=_dimension-1; i>=0; i--)
      {
        increment *=_nb_procs_per_dim[i];
        int idim=position%increment;
        position=position/increment;
        mcIdType imin=_local_array_indices[i][idim];
        mcIdType imax=_local_array_indices[i][idim+1];
        nb_elem*=(imax-imin);
      }
    return nb_elem;
  }

  /*! Retrieves the min and max indices of the domain stored locally
   * for each dimension. The output vector has the topology dimension
   * as a size and each pair <int,int> contains min and max. Indices 
   * range from min to max-1.
   */
  std::vector<std::pair<int,mcIdType> > BlockTopology::getLocalArrayMinMax() const
  {
    vector<pair<int,mcIdType> > local_indices (_dimension);
    int myrank=_proc_group->myRank();
    int increment=1;
    for (int i=_dimension-1; i>=0; i--)
      {  
        increment *=_nb_procs_per_dim[i];
        int idim=myrank%increment;
        local_indices[i].first=(int)_local_array_indices[i][idim];
        local_indices[i].second=_local_array_indices[i][idim+1];
        cout << local_indices[i].first << " "<< local_indices[i].second<<endl;
      }
    return local_indices;
  }

  /*! Serializes the data contained in the Block Topology
   * for communication purposes*/
  void BlockTopology::serialize(mcIdType* & serializer, mcIdType& size) const 
  {
    vector<mcIdType> buffer;
  
    buffer.push_back(_dimension);
    buffer.push_back(_nb_elems);
    for (int i=0; i<_dimension; i++)
      {
        buffer.push_back(_nb_procs_per_dim[i]);
        buffer.push_back(_cycle_type[i]);
        buffer.push_back(ToIdType(_local_array_indices[i].size()));
        for (std::size_t j=0; j<_local_array_indices[i].size(); j++)
          buffer.push_back(_local_array_indices[i][j]);
      }
  
    //serializing the comm group
    mcIdType size_comm=_proc_group->size();
    buffer.push_back(size_comm);
    MPIProcessorGroup world_group(_proc_group->getCommInterface());
    for (int i=0; i<size_comm;i++)
      {
        int world_rank=world_group.translateRank(_proc_group, i);
        buffer.push_back(world_rank);
      }
  
    serializer=new mcIdType[buffer.size()];
    size=ToIdType(buffer.size());
    copy(buffer.begin(), buffer.end(), serializer);
  }

  /*!
   *
   * Unserializes the data contained in the Block Topology
   * after communication. Uses the same structure as the one used for serialize() 
   *
   */
  void BlockTopology::unserialize(const mcIdType* serializer,const CommInterface& comm_interface)
  {
    const mcIdType* ptr_serializer=serializer;
    _dimension=(int)*(ptr_serializer++);
    _nb_elems=*(ptr_serializer++);
    _nb_procs_per_dim.resize(_dimension);
    _cycle_type.resize(_dimension);
    _local_array_indices.resize(_dimension);
    for (int i=0; i<_dimension; i++)
      {
        _nb_procs_per_dim[i]=(int)*(ptr_serializer++);
        _cycle_type[i]=(CYCLE_TYPE)*(ptr_serializer++);
        _local_array_indices[i].resize(*(ptr_serializer++));
        for (std::size_t j=0; j<_local_array_indices[i].size(); j++)
          _local_array_indices[i][j]=*(ptr_serializer++);
      }
    set<int> procs;
    mcIdType size_comm=*(ptr_serializer++);
    for (int i=0; i<size_comm; i++)
      procs.insert((int)*(ptr_serializer++));

    if (_owns_processor_group)
      delete _proc_group;
    _proc_group=new MPIProcessorGroup(comm_interface,procs);
    _owns_processor_group=true;
  }
}
