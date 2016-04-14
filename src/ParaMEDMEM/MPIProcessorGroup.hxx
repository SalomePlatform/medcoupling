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

#ifndef __MPIPROCESSORGROUP_HXX__
#define __MPIPROCESSORGROUP_HXX__

#include "ProcessorGroup.hxx"

#include <set>
#include <mpi.h>

namespace MEDCoupling
{
  class CommInterface;

  class MPIProcessorGroup : public ProcessorGroup
  {
  public:
    MPIProcessorGroup(const CommInterface& interface);
    MPIProcessorGroup(const CommInterface& interface, std::set<int> proc_ids, const MPI_Comm& world_comm=MPI_COMM_WORLD);
    MPIProcessorGroup (const ProcessorGroup& proc_group, std::set<int> proc_ids);
    MPIProcessorGroup(const CommInterface& interface,int pstart, int pend, const MPI_Comm& world_comm=MPI_COMM_WORLD);
    MPIProcessorGroup(const MPIProcessorGroup& other);
    virtual ~MPIProcessorGroup();
    virtual MPIProcessorGroup *deepCopy() const;
    virtual ProcessorGroup* fuse (const ProcessorGroup&) const;
    void intersect (ProcessorGroup&) { }
    int myRank() const;
    bool containsMyRank() const { int rank; MPI_Group_rank(_group, &rank); return (rank!=MPI_UNDEFINED); }
    int translateRank(const ProcessorGroup* group, int rank) const;
    const MPI_Comm* getComm() const { return &_comm; }
    ProcessorGroup* createComplementProcGroup() const;
    ProcessorGroup* createProcGroup() const;
    MPI_Comm getWorldComm() { return _world_comm; }
  private:
    void updateMPISpecificAttributes();
  private:
    const MPI_Comm _world_comm;  // just an observer - current instance is not responsible for the management of this comm
    MPI_Group _group;
    MPI_Comm _comm;
  };
}

#endif
