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
#include "MEDCouplingUMesh.hxx"
#include "ParaFIELD.hxx"
#include "ParaMESH.hxx"

#include <map>
#include <vector>

namespace MEDCoupling
{
class ProcessorGroup;
class MPIProcessorGroup;

class OverlapManagerDEC
{
   public:
    OverlapManagerDEC();

    void attachLocalField(const ParaFIELD *field, const MPI_Comm &comm);
    const ParaFIELD *getRestrictedField(bool &originalOrComputed) const;
    void doNotOwnRestrictedFieldAnymore();
    void ownOriginalField();

    void synchronizeGhosts();
    void updateLocalArray();

    ~OverlapManagerDEC();
    void release();

   private:
    void computeOwnership();
    bool recomputeOwnership() const;

    ParaFIELD *restrict(const ParaFIELD *pfield) const;
    ParaMESH *restrict(const MEDCouplingMesh *mesh) const;

    const ProcessorGroup *_local_group;
    MPIProcessorGroup *_mpi_group;

    MCAuto<DataArrayIdType> _global_ids;
    std::vector<int> _owner;
    std::vector<mcIdType> _kept, _removed, _shared, _o2r;
    std::vector<std::vector<mcIdType>> _joints_send, _joints_recv;

    bool _own_originalField;
    const ParaFIELD *_originalField;

    ParaFIELD *_restrictedField;

    bool _hasLocalOverlap, _hasGlobalOverlap;
    bool _own_restrictedField;
    void releaseRestrictedField();
    void releaseOriginalField();
    void releaseMPIGroup();
};
}  // namespace MEDCoupling
