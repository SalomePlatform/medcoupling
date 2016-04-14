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
    _proc_group(NULL),_nb_elems(0),
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
    vector <int> axis_length(_dimension);
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
        _local_array_indices.push_back(comp_indices);
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
  BlockTopology::BlockTopology(const ProcessorGroup& group, int nb_elem):_dimension(1),_proc_group(&group),_owns_processor_group(false)
  {
    int* nbelems_per_proc = new int[group.size()];
    const MPIProcessorGroup* mpi_group=dynamic_cast<const MPIProcessorGroup*>(_proc_group);
    const MPI_Comm* comm=mpi_group->getComm();
    int nbtemp=nb_elem;
    mpi_group->getCommInterface().allGather(&nbtemp, 1, MPI_INT, 
                                            nbelems_per_proc, 1, MPI_INT, 
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
    if (_owns_processor_group)
      delete _proc_group;
  }

  //!converts a pair <subdomainid,local> to a global number
  std::pair<int,int> BlockTopology::globalToLocal(const int global) const
  {
    int subdomain_id=0;
    int position=global;
    int size=_nb_elems;
    int size_procs=_proc_group->size();
    int increment=size;
    vector<int>axis_position(_dimension);
    vector<int>axis_offset(_dimension);
    for (int idim=0; idim<_dimension; idim++)
      {
        int axis_size=_local_array_indices[idim].size()-1;
        int axis_nb_elem=_local_array_indices[idim][axis_size];
        increment=increment/axis_nb_elem;
        int proc_increment = size_procs/(axis_size);
        int axis_pos=position/increment;
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
    int local=0;
    int local_increment=1;
    for (int idim=_dimension-1; idim>=0; idim--)
      {
        local+=axis_position[idim]*local_increment;
        local_increment*=_local_array_indices[idim][axis_offset[idim]]-_local_array_indices[idim][axis_offset[idim]-1];
      }
    return make_pair(subdomain_id,local);
  }

  //!converts local number to a global number
  int BlockTopology::localToGlobal(const pair<int,int> local) const
  {

    int subdomain_id=local.first;
    int global=0;
    int loc=local.second;
    int increment=_nb_elems;
    int proc_increment=_proc_group->size();
    int local_increment=getNbLocalElements();
    for (int idim=0; idim < _dimension; idim++)
      {
        int axis_size=_local_array_indices[idim].size()-1;
        int axis_nb_elem=_local_array_indices[idim][axis_size];
        increment=axis_nb_elem==0?0:increment/axis_nb_elem;
        proc_increment = proc_increment/(axis_size);
        int proc_axis=subdomain_id/proc_increment;
        subdomain_id=subdomain_id%proc_increment;
        int local_axis_nb_elem=_local_array_indices[idim][proc_axis+1]-_local_array_indices[idim][proc_axis];
        local_increment = (local_axis_nb_elem==0)?0:(local_increment/local_axis_nb_elem);
        int iaxis=((local_increment==0)?0:(loc/local_increment))+_local_array_indices[idim][proc_axis];
        global+=increment*iaxis;
        loc = (local_increment==0)?0:(loc%local_increment);
      }
    return global;
  }

  //Retrieves the local number of elements
  int BlockTopology::getNbLocalElements()const
  {
    int position=_proc_group->myRank();
    int nb_elem = 1;
    int increment=1;
    for (int i=_dimension-1; i>=0; i--)
      {
        increment *=_nb_procs_per_dim[i];
        int idim=position%increment;
        position=position/increment;
        int imin=_local_array_indices[i][idim];
        int imax=_local_array_indices[i][idim+1];
        nb_elem*=(imax-imin);
      }
    return nb_elem;
  }

  /*! Retrieves the min and max indices of the domain stored locally
   * for each dimension. The output vector has the topology dimension
   * as a size and each pair <int,int> contains min and max. Indices 
   * range from min to max-1.
   */
  std::vector<std::pair<int,int> > BlockTopology::getLocalArrayMinMax() const
  {
    vector<pair<int,int> > local_indices (_dimension);
    int myrank=_proc_group->myRank();
    int increment=1;
    for (int i=_dimension-1; i>=0; i--)
      {  
        increment *=_nb_procs_per_dim[i];
        int idim=myrank%increment;
        local_indices[i].first=_local_array_indices[i][idim];
        local_indices[i].second=_local_array_indices[i][idim+1];
        cout << local_indices[i].first << " "<< local_indices[i].second<<endl;
      }
    return local_indices;
  }

  /*! Serializes the data contained in the Block Topology
   * for communication purposes*/
  void BlockTopology::serialize(int* & serializer, int& size) const 
  {
    vector<int> buffer;
  
    buffer.push_back(_dimension);
    buffer.push_back(_nb_elems);
    for (int i=0; i<_dimension; i++)
      {
        buffer.push_back(_nb_procs_per_dim[i]);
        buffer.push_back(_cycle_type[i]);
        buffer.push_back(_local_array_indices[i].size());
        for (int j=0; j<(int)_local_array_indices[i].size(); j++)
          buffer.push_back(_local_array_indices[i][j]);
      }
  
    //serializing the comm group
    int size_comm=_proc_group->size();
    buffer.push_back(size_comm);
    MPIProcessorGroup world_group(_proc_group->getCommInterface());
    for (int i=0; i<size_comm;i++)
      {
        int world_rank=world_group.translateRank(_proc_group, i);
        buffer.push_back(world_rank);
      }
  
    serializer=new int[buffer.size()];
    size=buffer.size();
    copy(buffer.begin(), buffer.end(), serializer);
  }

  /*!
   *
   * Unserializes the data contained in the Block Topology
   * after communication. Uses the same structure as the one used for serialize() 
   *
   */
  void BlockTopology::unserialize(const int* serializer,const CommInterface& comm_interface)
  {
    const int* ptr_serializer=serializer;
    cout << "unserialize..."<<endl;
    _dimension=*(ptr_serializer++);
    cout << "dimension "<<_dimension<<endl;
    _nb_elems=*(ptr_serializer++);
    cout << "nbelems "<<_nb_elems<<endl;
    _nb_procs_per_dim.resize(_dimension);
    _cycle_type.resize(_dimension);
    _local_array_indices.resize(_dimension);
    for (int i=0; i<_dimension; i++)
      {
        _nb_procs_per_dim[i]=*(ptr_serializer++);
        _cycle_type[i]=(CYCLE_TYPE)*(ptr_serializer++);
        _local_array_indices[i].resize(*(ptr_serializer++));
        for (int j=0; j<(int)_local_array_indices[i].size(); j++)
          _local_array_indices[i][j]=*(ptr_serializer++);
      }
    set<int> procs;
    int size_comm=*(ptr_serializer++);
    for (int i=0; i<size_comm; i++)
      procs.insert(*(ptr_serializer++));
    cout << "unserialize..."<<procs.size()<<endl;
    _proc_group=new MPIProcessorGroup(comm_interface,procs);
    _owns_processor_group=true;
    //TODO manage memory ownership of _proc_group  
  }
}
