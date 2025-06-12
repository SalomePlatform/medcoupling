// Copyright (C) 2007-2025  CEA, EDF
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
// Author : Anthony Geay (CEA/DEN)

#include "OverlapDEC.hxx"
#include "CommInterface.hxx"
#include "ParaMESH.hxx"
#include "ParaFIELD.hxx"
#include "MPIProcessorGroup.hxx"
#include "OverlapElementLocator.hxx"
#include "OverlapInterpolationMatrix.hxx"
#include "ICoCoMEDDoubleField.hxx"

namespace MEDCoupling
{
OverlapDEC::OverlapDEC(const std::set<int> &procIds, const MPI_Comm &world_comm)
    : _load_balancing_algo(1),
      _own_group(true),
      _interpolation_matrix(0),
      _locator(0),
      _default_field_value(0.0),
      _source_field(0),
      _own_source_field(false),
      _target_field(0),
      _own_target_field(false),
      _comm(MPI_COMM_NULL)
{
    MEDCoupling::CommInterface comm;
    int *ranks_world = new int[procIds.size()];  // ranks of sources and targets in world_comm
    std::copy(procIds.begin(), procIds.end(), ranks_world);
    MPI_Group group, world_group;
    comm.commGroup(world_comm, &world_group);
    comm.groupIncl(world_group, (int)procIds.size(), ranks_world, &group);
    delete[] ranks_world;
    comm.commCreate(world_comm, group, &_comm);
    comm.groupFree(&group);
    comm.groupFree(&world_group);
    if (_comm == MPI_COMM_NULL)
    {
        _group = 0;
        return;
    }
    std::set<int> idsUnion;
    for (unsigned int i = 0; i < procIds.size(); i++) idsUnion.insert(i);
    _group = new MPIProcessorGroup(comm, idsUnion, _comm);
}

OverlapDEC::~OverlapDEC() { release(); }

/** Destructor involves MPI operations: make sure this is accessible from a proper
 * method for Python wrapping.
 */
void
OverlapDEC::release()
{
    if (_own_group)
    {
        delete _group;
        _group = nullptr;
    }
    if (_own_source_field)
    {
        delete _source_field;
        _source_field = nullptr;
    }
    if (_own_target_field)
    {
        delete _target_field;
        _target_field = nullptr;
    }
    delete _interpolation_matrix;
    _interpolation_matrix = nullptr;
    delete _locator;
    _locator = nullptr;
    if (_comm != MPI_COMM_NULL)
    {
        MEDCoupling::CommInterface comm;
        comm.commFree(&_comm);
    }
    _comm = MPI_COMM_NULL;
}

void
OverlapDEC::sendRecvData(bool way)
{
    if (way)
        sendData();
    else
        recvData();
}

void
OverlapDEC::sendData()
{
    _interpolation_matrix->multiply(_default_field_value);
}

void
OverlapDEC::recvData()
{
    throw INTERP_KERNEL::Exception("Not implemented yet !!!!");
    //_interpolation_matrix->transposeMultiply();
}

void
OverlapDEC::synchronize()
{
    if (!isInGroup())
        return;
    // Check number of components of field on both side (for now allowing void field/mesh on one proc is not allowed)
    if (!_source_field || !_source_field->getField())
        throw INTERP_KERNEL::Exception(
            "OverlapDEC::synchronize(): currently, having a void source field on a proc is not allowed!"
        );
    if (!_target_field || !_target_field->getField())
        throw INTERP_KERNEL::Exception(
            "OverlapDEC::synchronize(): currently, having a void target field on a proc is not allowed!"
        );
    if (_target_field->getField()->getNumberOfComponents() != _source_field->getField()->getNumberOfComponents())
        throw INTERP_KERNEL::Exception(
            "OverlapDEC::synchronize(): source and target field have different number of components!"
        );
    delete _interpolation_matrix;
    _locator = new OverlapElementLocator(
        _source_field, _target_field, *_group, getBoundingBoxAdjustmentAbs(), _load_balancing_algo
    );
    _interpolation_matrix =
        new OverlapInterpolationMatrix(_source_field, _target_field, *_group, *this, *this, *_locator);
    _locator->copyOptions(*this);
    _locator->exchangeMeshes(*_interpolation_matrix);
    std::vector<std::pair<int, int> > jobs = _locator->getToDoList();
    std::string srcMeth = _locator->getSourceMethod();
    std::string trgMeth = _locator->getTargetMethod();
    for (std::vector<std::pair<int, int> >::const_iterator it = jobs.begin(); it != jobs.end(); it++)
    {
        const MEDCouplingPointSet *src = _locator->getSourceMesh((*it).first);
        const DataArrayIdType *srcIds = _locator->getSourceIds((*it).first);
        const MEDCouplingPointSet *trg = _locator->getTargetMesh((*it).second);
        const DataArrayIdType *trgIds = _locator->getTargetIds((*it).second);
        _interpolation_matrix->computeLocalIntersection(
            src, srcIds, srcMeth, (*it).first, trg, trgIds, trgMeth, (*it).second
        );
    }
    _interpolation_matrix->prepare(_locator->getProcsToSendFieldData());
    _interpolation_matrix->computeSurfacesAndDeno();
}

void
OverlapDEC::attachSourceLocalField(ParaFIELD *field, bool ownPt)
{
    if (!isInGroup())
        return;
    if (_own_source_field)
        delete _source_field;
    _source_field = field;
    _own_source_field = ownPt;
}

void
OverlapDEC::attachTargetLocalField(ParaFIELD *field, bool ownPt)
{
    if (!isInGroup())
        return;
    if (_own_target_field)
        delete _target_field;
    _target_field = field;
    _own_target_field = ownPt;
}

void
OverlapDEC::attachSourceLocalField(MEDCouplingFieldDouble *field)
{
    if (!isInGroup())
        return;

    ParaMESH *paramesh = new ParaMESH(
        static_cast<MEDCouplingPointSet *>(const_cast<MEDCouplingMesh *>(field->getMesh())),
        *_group,
        field->getMesh()->getName()
    );
    ParaFIELD *tmpField = new ParaFIELD(field, paramesh, *_group);
    tmpField->setOwnSupport(true);
    attachSourceLocalField(tmpField, true);
}

void
OverlapDEC::attachTargetLocalField(MEDCouplingFieldDouble *field)
{
    if (!isInGroup())
        return;

    ParaMESH *paramesh = new ParaMESH(
        static_cast<MEDCouplingPointSet *>(const_cast<MEDCouplingMesh *>(field->getMesh())),
        *_group,
        field->getMesh()->getName()
    );
    ParaFIELD *tmpField = new ParaFIELD(field, paramesh, *_group);
    tmpField->setOwnSupport(true);
    attachTargetLocalField(tmpField, true);
}

void
OverlapDEC::attachSourceLocalField(ICoCo::MEDDoubleField *field)
{
    attachSourceLocalField(field->getMCField());
}

void
OverlapDEC::attachTargetLocalField(ICoCo::MEDDoubleField *field)
{
    attachTargetLocalField(field->getMCField());
}

bool
OverlapDEC::isInGroup() const
{
    if (!_group)
        return false;
    return _group->containsMyRank();
}

void
OverlapDEC::debugPrintWorkSharing(std::ostream &ostr) const
{
    _locator->debugPrintWorkSharing(ostr);
}
}  // namespace MEDCoupling
