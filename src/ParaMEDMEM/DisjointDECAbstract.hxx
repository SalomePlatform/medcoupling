// Copyright (C) 2007-2026  CEA, EDF
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

#pragma once

#include "MEDCouplingFieldDouble.hxx"
#include "NormalizedUnstructuredMesh.hxx"
#include "DEC.hxx"

#include <mpi.h>
#include <set>

namespace ICoCo
{
class MEDDoubleField;
}

namespace MEDCoupling
{
class ProcessorGroup;
class ParaFIELD;

class DisjointDECAbstract : public DEC
{
   public:
    DisjointDECAbstract()
        : _union_group(nullptr),
          _source_group(nullptr),
          _target_group(nullptr),
          _comm_interface(nullptr),
          _owns_groups(false),
          _union_comm(MPI_COMM_NULL)
    {
    }
    DisjointDECAbstract(ProcessorGroup &source_group, ProcessorGroup &target_group);
    DisjointDECAbstract(const DisjointDECAbstract &);
    DisjointDECAbstract &operator=(const DisjointDECAbstract &s);
    DisjointDECAbstract(
        const std::set<int> &src_ids, const std::set<int> &trg_ids, const MPI_Comm &world_comm = MPI_COMM_WORLD
    );
    virtual ~DisjointDECAbstract();

    virtual void computeProcGroup() {}
    //
    ProcessorGroup *getSourceGrp() const { return _source_group; }
    ProcessorGroup *getTargetGrp() const { return _target_group; }
    ProcessorGroup *getUnionGrp() const { return _union_group; }
    bool isInSourceSide() const;
    bool isInTargetSide() const;
    bool isInUnion() const;

   protected:
    void compareFieldAndMethod() const;
    void cleanInstance();
    void copyInstance(const DisjointDECAbstract &other);
    void checkPartitionGroup() const;

   protected:
    //! Processor group representing the union of target and source processors
    ProcessorGroup *_union_group;
    ProcessorGroup *_source_group;
    ProcessorGroup *_target_group;
    bool _owns_groups;

    const CommInterface *_comm_interface;
    MPI_Comm _union_comm;
};

}  // namespace MEDCoupling
