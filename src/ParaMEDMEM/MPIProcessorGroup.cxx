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

#include "ProcessorGroup.hxx"
#include "MPIProcessorGroup.hxx"
#include "CommInterface.hxx"
#include "InterpolationUtils.hxx"

#include <iostream>
#include <set>
#include <algorithm>
#include "mpi.h"

using namespace std;


namespace MEDCoupling
{
  /*!
   \anchor MPIProcessorGroup-det
   \class MPIProcessorGroup

   The MPIProcessorGroup class represents a set of distinct "processors" (computation nodes)
   in a MPI code. It is used to define the MPI topology of code couplings.

   Groups can be set up in various ways, the most common being
   the use of the \c MPIProcessorGroup(Comminterface, int pfirst, int plast)
   constructor.

   The following code excerpt creates two processor groups on respectively 3 and 2 processors.
   \verbatim
   int main()
   {
   MPI_Init(&argc,&argv);
   CommInterface comm_interface;
   MPIProcessorGroup codeA_group(comm_interface, 0, 2);  // groups processors 0, 1 and 2
   MPIProcessorGroup codeB_group(comm_interface, 3, 4);  // groups processors 3 and 4

   ...
   }
   \endverbatim
  */


  /*! 
   * Creates a processor group that is based on all the
   processors of MPI_COMM_WORLD .This routine must be called by all processors in MPI_COMM_WORLD.
   \param interface CommInterface object giving access to the MPI
   communication layer
  */
  MPIProcessorGroup::MPIProcessorGroup(const CommInterface& interface):
    ProcessorGroup(interface),_world_comm(MPI_COMM_WORLD)
  {
    _comm=_world_comm;
    _comm_interface.commGroup(_world_comm, &_group);
    int size;
    _comm_interface.commSize(_world_comm,&size);
    for (int i=0; i<size; i++)
      _proc_ids.insert(i);

  }

  /*! Creates a processor group that is based on the processors included in \a proc_ids.
    This routine must be called by all processors in MPI_COMM_WORLD.

    \param interface CommInterface object giving access to the MPI
    communication layer
    \param proc_ids set of ids that are to be integrated in the group. The ids number are 
    to be understood in terms of MPI_COMM_WORLD ranks.
  */

  MPIProcessorGroup::MPIProcessorGroup(const CommInterface& interface, set<int> proc_ids, const MPI_Comm& world_comm):
    ProcessorGroup(interface, proc_ids), _world_comm(world_comm)
  {
    updateMPISpecificAttributes();
  }


  void MPIProcessorGroup::updateMPISpecificAttributes()
  {
    //Creation of a communicator 
    MPI_Group group_world;
  
    int size_world;
    _comm_interface.commSize(_world_comm,&size_world);
    int rank_world;
    _comm_interface.commRank(_world_comm,&rank_world);
    _comm_interface.commGroup(_world_comm, &group_world);

    int* ranks=new int[_proc_ids.size()];
   
    // copying proc_ids in ranks
    copy<set<int>::const_iterator,int*> (_proc_ids.begin(), _proc_ids.end(), ranks);
    for (int i=0; i< (int)_proc_ids.size();i++)
      if (ranks[i]>size_world-1)
        {
          delete[] ranks;
          _comm_interface.groupFree(&group_world);  // MPI_Group is a C structure and won't get de-allocated automatically?
          throw INTERP_KERNEL::Exception("invalid rank in set<int> argument of MPIProcessorGroup constructor");
        }
      
    _comm_interface.groupIncl(group_world, _proc_ids.size(), ranks, &_group);
  
    _comm_interface.commCreate(_world_comm, _group, &_comm);

    // clean-up
    delete[] ranks;
    _comm_interface.groupFree(&group_world);  // MPI_Group is a C structure and won't get de-allocated automatically?
  }

  /*! Creates a processor group that is based on the processors between \a pstart and \a pend.
    This routine must be called by all processors in MPI_COMM_WORLD.

    \param comm_interface CommInterface object giving access to the MPI
    communication layer
    \param pstart id in MPI_COMM_WORLD of the first processor in the group
    \param pend id in MPI_COMM_WORLD of the last processor in the group
  */
  MPIProcessorGroup::MPIProcessorGroup (const CommInterface& comm_interface, int pstart, int pend, const MPI_Comm& world_comm): ProcessorGroup(comm_interface,pstart,pend),_world_comm(world_comm)
  {
    //Creation of a communicator 
    MPI_Group group_world;
  
    int size_world;
    _comm_interface.commSize(_world_comm,&size_world);
    int rank_world;
    _comm_interface.commRank(_world_comm,&rank_world);
    _comm_interface.commGroup(_world_comm, &group_world);

    if (pend>size_world-1 || pend <pstart || pstart<0)
      {
        _comm_interface.groupFree(&group_world);
        throw INTERP_KERNEL::Exception("invalid argument in MPIProcessorGroup constructor (comm,pfirst,plast)");
      }
    int nprocs=pend-pstart+1;
    int* ranks=new int[nprocs];
    for (int i=pstart; i<=pend;i++)
      {
        ranks[i-pstart]=i;
      }

    _comm_interface.groupIncl(group_world, nprocs, ranks, &_group);
  
    _comm_interface.commCreate(_world_comm, _group, &_comm);

    // clean-up
    delete[] ranks;
    _comm_interface.groupFree(&group_world);  // MPI_Group is a C structured and won't get de-allocated automatically?
  }

  MPIProcessorGroup::MPIProcessorGroup (const ProcessorGroup& proc_group, set<int> proc_ids) :
    ProcessorGroup(proc_group.getCommInterface()),
    _world_comm(MPI_COMM_WORLD), _group(MPI_GROUP_NULL), _comm(MPI_COMM_NULL)
  {
    cout << "MPIProcessorGroup (const ProcessorGroup& proc_group, set<int> proc_ids)" <<endl;
    cout << "Not implemented yet !"<<endl;
    exit(1);
  }

  MPIProcessorGroup::MPIProcessorGroup(const MPIProcessorGroup& other):
      ProcessorGroup(other),_world_comm(other._world_comm)
  {
    updateMPISpecificAttributes();
  }

  MPIProcessorGroup::~MPIProcessorGroup()
  {
    _comm_interface.groupFree(&_group);
    if (_comm!=_world_comm && _comm !=MPI_COMM_NULL)
      _comm_interface.commFree(&_comm);
  
  }

  /*! Translation of the rank id between two processor groups. This method translates rank \a rank
    on the current processor group to the rank on group pointed by \a group.
    \param group group from which the rank is expected
    \param rank rank on group \a group of the processor which is to be translated
    \return rank on local group
  */
  int MPIProcessorGroup::translateRank(const ProcessorGroup* group, int rank) const
  {
    const MPIProcessorGroup* targetgroup=dynamic_cast<const MPIProcessorGroup*>(group);
    int local_rank;
    MPI_Group_translate_ranks(targetgroup->_group, 1, &rank, _group, &local_rank);
    return local_rank;
  }
  
  /*!Creates a processor group that is the complement of the current group 
    inside MPI_COMM_WORLD
    \return pointer to the new ProcessorGroup structure.
  */
  ProcessorGroup* MPIProcessorGroup::createComplementProcGroup() const
  {
    set <int> procs;
    int world_size=_comm_interface.worldSize();
    for (int i=0; i<world_size; i++)
      procs.insert(i);
    for (set<int>::const_iterator iter=_proc_ids.begin(); iter!= _proc_ids.end(); iter++)
      procs.erase(*iter);
    
    return new MPIProcessorGroup(_comm_interface, procs, _world_comm);
    
  }

  MPIProcessorGroup *MPIProcessorGroup::deepCopy() const
  {
    return new MPIProcessorGroup(*this);
  }

  /*!Adding processors of group \a group to local group.
    \param group group that is to be fused with current group
    \return new group formed by the fusion of local group and \a group.
  */
  ProcessorGroup*  MPIProcessorGroup::fuse (const ProcessorGroup& group) const
  {
    set <int> procs = _proc_ids;
    const set<int>& distant_proc_ids = group.getProcIDs();
    for (set<int>::const_iterator iter=distant_proc_ids.begin(); iter!=distant_proc_ids.end(); iter++)
      {
        procs.insert(*iter);
      }
    return new MPIProcessorGroup(_comm_interface, procs, _world_comm);
  }

  int MPIProcessorGroup::myRank() const
  { 
    int rank;
    MPI_Comm_rank(_comm,&rank);
    return rank;
  }
  
  ProcessorGroup* MPIProcessorGroup::createProcGroup() const
  {
    set <int> procs;
    for (set<int>::const_iterator iter=_proc_ids.begin(); iter!= _proc_ids.end(); iter++)
      procs.insert(*iter);
  
    return new MPIProcessorGroup(_comm_interface, procs, _world_comm);

  }
}
