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

#include "DisjointDEC.hxx"
#include "CommInterface.hxx"
#include "Topology.hxx"
#include "BlockTopology.hxx"
#include "ComponentTopology.hxx"
#include "ParaFIELD.hxx"
#include "ParaMESH.hxx"
#include "ICoCoField.hxx"
#include "ICoCoMEDDoubleField.hxx"
#include "MPIProcessorGroup.hxx"

#include <cmath>
#include <iostream>

namespace MEDCoupling
{

DisjointDEC::DisjointDEC(ProcessorGroup &source_group, ProcessorGroup &target_group)
    : DisjointDECAbstract(source_group, target_group), _local_field(nullptr), _owns_field(false)
{
    checkPartitionGroup();
    _union_group = source_group.fuse(target_group);
}

DisjointDEC::DisjointDEC(const DisjointDEC &s) : DisjointDECAbstract(s), _local_field(nullptr), _owns_field(false)
{
    copyInstance(s);
}

DisjointDEC &
DisjointDEC::operator=(const DisjointDEC &s)
{
    cleanInstance();
    copyInstance(s);
    return *this;
}

DisjointDEC::DisjointDEC(const std::set<int> &source_ids, const std::set<int> &target_ids, const MPI_Comm &world_comm)
    : DisjointDECAbstract(source_ids, target_ids, world_comm), _local_field(nullptr), _owns_field(false)
{
}

DisjointDEC::~DisjointDEC() {}

void
DisjointDEC::setNature(NatureOfField nature)
{
    if (_local_field)
        _local_field->getField()->setNature(nature);
}

/*! Attaches a local field to a DEC.
  If the processor is on the receiving end of the DEC, the field
  will be updated by a recvData() call.
  Reversely, if the processor is on the sending end, the field will be read, possibly transformed, and sent
  appropriately to the other side.
*/
void
DisjointDEC::attachLocalField(const ParaFIELD *field, bool ownPt)
{
    if (!isInUnion())
        return;
    if (_owns_field)
        delete _local_field;
    _local_field = field;
    _owns_field = ownPt;
    _comm_interface = &(field->getTopology()->getProcGroup()->getCommInterface());
    compareFieldAndMethod();
}

/*! Attaches a local field to a DEC. The method will test whether the processor
  is on the source or the target side and will associate the mesh underlying the
  field to the local side.

  If the processor is on the receiving end of the DEC, the field
  will be updated by a recvData() call.
  Reversely, if the processor is on the sending end, the field will be read, possibly transformed,
  and sent appropriately to the other side.
*/

void
DisjointDEC::attachLocalField(MEDCouplingFieldDouble *field)
{
    if (!isInUnion())
        return;
    ProcessorGroup *local_group;
    if (_source_group->containsMyRank())
        local_group = _source_group;
    else if (_target_group->containsMyRank())
        local_group = _target_group;
    else
        throw INTERP_KERNEL::Exception("Invalid procgroup for field attachment to DEC");
    ParaMESH *paramesh = new ParaMESH(
        static_cast<MEDCouplingPointSet *>(const_cast<MEDCouplingMesh *>(field->getMesh())),
        *local_group,
        field->getMesh()->getName()
    );
    ParaFIELD *tmp = new ParaFIELD(field, paramesh, *local_group);
    tmp->setOwnSupport(true);
    attachLocalField(tmp, true);
    //_comm_interface=&(local_group->getCommInterface());
}

/*!
  Attaches a local field to a DEC.
  If the processor is on the receiving end of the DEC, the field
  will be updated by a recvData() call.
  Reversely, if the processor is on the sending end, the field will be read, possibly transformed, and sent
  appropriately to the other side. The field type is a generic ICoCo Field, so that the DEC can couple a number of
  different fields :
  - a ICoCo::MEDDoubleField, that is created from a MEDCoupling structure

*/
void
DisjointDEC::attachLocalField(const ICoCo::MEDDoubleField *field)
{
    if (!isInUnion())
        return;
    if (!field)
        throw INTERP_KERNEL::Exception("DisjointDEC::attachLocalField : ICoCo::MEDDoubleField pointer is NULL !");
    attachLocalField(field->getMCField());
}

/*!
  Computes the field norm over its support
  on the source side and renormalizes the field on the target side
  so that the norms match.

  \f[
  I_{source}=\sum_{i=1}^{n_{source}}V_{i}.|\Phi^{source}_{i}|^2,
  \f]

  \f[
  I_{target}=\sum_{i=1}^{n_{target}}V_{i}.|\Phi^{target}_{i}|^2,
  \f]

  \f[
  \Phi^{target}:=\Phi^{target}.\sqrt{I_{source}/I_{target}}.
  \f]

*/
void
DisjointDEC::renormalizeTargetField(bool isWAbs)
{
    if (_source_group->containsMyRank())
        for (int icomp = 0; icomp < (int)_local_field->getField()->getArray()->getNumberOfComponents(); icomp++)
        {
            double total_norm = _local_field->getVolumeIntegral(icomp + 1, isWAbs);
            double source_norm = total_norm;
            _comm_interface->broadcast(
                &source_norm, 1, MPI_DOUBLE, 0, *dynamic_cast<MPIProcessorGroup *>(_union_group)->getComm()
            );
        }
    if (_target_group->containsMyRank())
    {
        for (int icomp = 0; icomp < (int)_local_field->getField()->getArray()->getNumberOfComponents(); icomp++)
        {
            double total_norm = _local_field->getVolumeIntegral(icomp + 1, isWAbs);
            double source_norm = total_norm;
            _comm_interface->broadcast(
                &source_norm, 1, MPI_DOUBLE, 0, *dynamic_cast<MPIProcessorGroup *>(_union_group)->getComm()
            );

            if (fabs(total_norm) > 1e-100)
                _local_field->getField()->applyLin(source_norm / total_norm, 0.0, icomp + 1);
        }
    }
}

void
DisjointDEC::cleanInstance()
{
    if (_owns_field)
    {
        delete _local_field;
    }
    _local_field = nullptr;
    _owns_field = false;
    DisjointDECAbstract::cleanInstance();
}

void
DisjointDEC::compareFieldAndMethod() const
{
    if (_local_field)
    {
        TypeOfField entity = _local_field->getField()->getTypeOfField();
        if (getMethod() == "P0")
        {
            if (entity != ON_CELLS)
                throw INTERP_KERNEL::Exception(
                    "Field support and interpolation method mismatch."
                    " For P0 interpolation, field must be on MED_CELL's"
                );
        }
        else if (getMethod() == "P1")
        {
            if (entity != ON_NODES)
                throw INTERP_KERNEL::Exception(
                    "Field support and interpolation method mismatch."
                    " For P1 interpolation, field must be on MED_NODE's"
                );
        }
        else if (getMethod() == "P1d")
        {
            if (entity != ON_CELLS)
                throw INTERP_KERNEL::Exception(
                    "Field support and interpolation method mismatch."
                    " For P1d interpolation, field must be on MED_CELL's"
                );
            if (_target_group->containsMyRank())
                throw INTERP_KERNEL::Exception("Projection to P1d field not supported");
        }
        else
        {
            throw INTERP_KERNEL::Exception("Unknown interpolation method. Possible methods: P0, P1, P1d");
        }
    }
}

/*!
  If way==true, source procs call sendData() and target procs call recvData().
  if way==false, it's the other way round.
*/
void
DisjointDEC::sendRecvData(bool way)
{
    if (!isInUnion())
        return;
    if (isInSourceSide())
    {
        if (way)
            sendData();
        else
            recvData();
    }
    else if (isInTargetSide())
    {
        if (way)
            recvData();
        else
            sendData();
    }
}

}  // namespace MEDCoupling
