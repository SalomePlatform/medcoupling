// Copyright (C) 2025  CEA, EDF
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

#include "DisjointDECAbstract.hxx"
#include "CommInterface.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "ComponentTopology.hxx"
#include "MPIProcessorGroup.hxx"

#include <cmath>
#include <iostream>

namespace MEDCoupling
{
DisjointDECAbstract::DisjointDECAbstract(ProcessorGroup &source_group, ProcessorGroup &target_group)
    : _union_comm(MPI_COMM_NULL),
      _source_group(&source_group),
      _target_group(&target_group),
      _owns_groups(false),
      _comm_interface(nullptr)
{
    checkPartitionGroup();
    _union_group = source_group.fuse(target_group);
}

DisjointDECAbstract::DisjointDECAbstract(const DisjointDECAbstract &s)
    : DEC(s),
      _owns_groups(false),
      _union_group(nullptr),
      _source_group(nullptr),
      _target_group(nullptr),
      _comm_interface(nullptr),
      _union_comm(MPI_COMM_NULL)
{
    copyInstance(s);
}

DisjointDECAbstract &
DisjointDECAbstract::operator=(const DisjointDECAbstract &s)
{
    cleanInstance();
    copyInstance(s);
    return *this;
}

void
DisjointDECAbstract::copyInstance(const DisjointDECAbstract &other)
{
    DEC::copyFrom(other);
    if (other._union_comm != MPI_COMM_NULL)
    {
        // Tricky: the DEC is responsible for the management of _union_comm. And this comm is referenced by
        // the MPIProcGroups (source/targets). In the case where _union_comm is not NULL we must take care of rebuilding
        // the MPIProcGroups with a communicator that will survive the destruction of 'other'.
        _owns_groups = true;
        MPI_Comm_dup(other._union_comm, &_union_comm);
        //        std::cout << "DUP union comm - new is "<< _union_comm << "\n";
        _target_group = new MPIProcessorGroup(*_comm_interface, other._target_group->getProcIDs(), _union_comm);
        _source_group = new MPIProcessorGroup(*_comm_interface, other._source_group->getProcIDs(), _union_comm);
    }
    else
    {
        if (other._target_group)
        {
            _target_group = other._target_group->deepCopy();
            _owns_groups = true;
        }
        if (other._source_group)
        {
            _source_group = other._source_group->deepCopy();
            _owns_groups = true;
        }
    }
    if (_source_group && _target_group)
        _union_group = _source_group->fuse(*_target_group);
}

DisjointDECAbstract::DisjointDECAbstract(
    const std::set<int> &source_ids, const std::set<int> &target_ids, const MPI_Comm &world_comm
)
    : _comm_interface(nullptr), _owns_groups(true), _union_comm(MPI_COMM_NULL)
{
    MEDCoupling::CommInterface comm;
    // Create the list of procs including source and target
    std::set<int> union_ids;  // source and target ids in world_comm
    union_ids.insert(source_ids.begin(), source_ids.end());
    union_ids.insert(target_ids.begin(), target_ids.end());
    if (union_ids.size() != (source_ids.size() + target_ids.size()))
        throw INTERP_KERNEL::Exception(
            "DisjointDECAbstract constructor : source_ids and target_ids overlap partially or fully. This type of DEC "
            "does not "
            "support it! OverlapDEC class could be the solution!"
        );
    int *union_ranks_world = new int[union_ids.size()];  // ranks of sources and targets in world_comm
    std::copy(union_ids.begin(), union_ids.end(), union_ranks_world);

    // Create a communicator on these procs
    MPI_Group union_group, world_group;
    comm.commGroup(world_comm, &world_group);
    comm.groupIncl(world_group, (int)union_ids.size(), union_ranks_world, &union_group);
    comm.commCreate(world_comm, union_group, &_union_comm);
    delete[] union_ranks_world;
    if (_union_comm == MPI_COMM_NULL)
    {  // This process is not in union
        _source_group = 0;
        _target_group = 0;
        _union_group = 0;
        comm.groupFree(&union_group);
        comm.groupFree(&world_group);
        return;
    }

    // Translate source_ids and target_ids from world_comm to union_comm
    int *source_ranks_world = new int[source_ids.size()];  // ranks of sources in world_comm
    std::copy(source_ids.begin(), source_ids.end(), source_ranks_world);
    int *source_ranks_union = new int[source_ids.size()];  // ranks of sources in union_comm
    int *target_ranks_world = new int[target_ids.size()];  // ranks of targets in world_comm
    std::copy(target_ids.begin(), target_ids.end(), target_ranks_world);
    int *target_ranks_union = new int[target_ids.size()];  // ranks of targets in union_comm
    MPI_Group_translate_ranks(world_group, (int)source_ids.size(), source_ranks_world, union_group, source_ranks_union);
    MPI_Group_translate_ranks(world_group, (int)target_ids.size(), target_ranks_world, union_group, target_ranks_union);
    std::set<int> source_ids_union;
    for (int i = 0; i < (int)source_ids.size(); i++) source_ids_union.insert(source_ranks_union[i]);
    std::set<int> target_ids_union;
    for (int i = 0; i < (int)target_ids.size(); i++) target_ids_union.insert(target_ranks_union[i]);
    delete[] source_ranks_world;
    delete[] source_ranks_union;
    delete[] target_ranks_world;
    delete[] target_ranks_union;

    // Create the MPIProcessorGroups
    _source_group = new MPIProcessorGroup(comm, source_ids_union, _union_comm);
    _target_group = new MPIProcessorGroup(comm, target_ids_union, _union_comm);
    _union_group = _source_group->fuse(*_target_group);
    comm.groupFree(&union_group);
    comm.groupFree(&world_group);
}

DisjointDECAbstract::~DisjointDECAbstract() { cleanInstance(); }

void
DisjointDECAbstract::cleanInstance()
{
    if (_owns_groups)
    {
        delete _source_group;
        delete _target_group;
    }
    _owns_groups = false;
    _source_group = nullptr;
    _target_group = nullptr;
    delete _union_group;
    _union_group = nullptr;
    if (_union_comm != MPI_COMM_NULL)
        _comm_interface->commFree(&_union_comm);
    _union_comm = MPI_COMM_NULL;
}

/**
 * Check that the sources and targets procs form a partition of the world communicator referenced in the groups.
 * This world communicator is not necessarily MPI_WORLD_COMM, but it has to be covered completely for the DECs to work.
 */
void
DisjointDECAbstract::checkPartitionGroup() const
{
    int size = -1;
    MPIProcessorGroup *tgt = static_cast<MPIProcessorGroup *>(_target_group);
    MPIProcessorGroup *src = static_cast<MPIProcessorGroup *>(_source_group);
    MPI_Comm comm_t = tgt->getWorldComm();
    MPI_Comm comm_s = src->getWorldComm();
    if (comm_t != comm_s)
        throw INTERP_KERNEL::Exception(
            "DisjointDECAbstract constructor: Inconsistent world communicator when building DisjointDECAbstract"
        );
    MPI_Comm_size(comm_t, &size);

    std::set<int> union_ids;  // source and target ids in world_comm
    union_ids.insert(src->getProcIDs().begin(), src->getProcIDs().end());
    union_ids.insert(tgt->getProcIDs().begin(), tgt->getProcIDs().end());
    if ((int)union_ids.size() != size)
        throw INTERP_KERNEL::Exception(
            "DisjointDECAbstract constructor: source_ids and target_ids do not form a partition of the communicator! "
            "Restrain "
            "the world communicator passed to MPIProcessorGroup ctor."
        );
}

bool
DisjointDECAbstract::isInSourceSide() const
{
    if (!_source_group)
        return false;
    return _source_group->containsMyRank();
}

bool
DisjointDECAbstract::isInTargetSide() const
{
    if (!_target_group)
        return false;
    return _target_group->containsMyRank();
}

bool
DisjointDECAbstract::isInUnion() const
{
    if (!_union_group)
        return false;
    return _union_group->containsMyRank();
}

}  // namespace MEDCoupling
