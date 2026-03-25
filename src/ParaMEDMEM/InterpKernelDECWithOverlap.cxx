// Copyright (C) 2026  CEA, EDF
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

#include "InterpKernelDECWithOverlap.hxx"

using namespace MEDCoupling;

InterpKernelDECWithOverlap::~InterpKernelDECWithOverlap() { this->release(); }

void
InterpKernelDECWithOverlap::release()
{
    _overlap.release();
}

void
InterpKernelDECWithOverlap::attachLocalField(const ParaFIELD *field, bool ownPt)
{
    const MPIProcessorGroup *grp(dynamic_cast<const MPIProcessorGroup *>(this->getUnionGrp()));
    _overlap.attachLocalField(field, *grp->getComm());
    bool originalOrComputed;
    const ParaFIELD *field2(_overlap.getRestrictedField(originalOrComputed));
    bool myOwnPt(ownPt);
    if (originalOrComputed)  // if Computed tell _overlap to stop owning computed ParaField instance -> transfert
                             // ownership to this
    {
        if (ownPt)
        {
            _overlap.ownOriginalField();
        }
        myOwnPt = true;
        _overlap.doNotOwnRestrictedFieldAnymore();
    }
    InterpKernelDEC::attachLocalField(field2, myOwnPt);
};

void
InterpKernelDECWithOverlap::attachLocalField(MEDCouplingFieldDouble *field, DataArrayIdType *globalIds)
{
    if (!isInUnion())
        return;
    MPIProcessorGroup *local_group(this->getMyProcGroup());

    ParaMESH *paramesh = new ParaMESH(
        static_cast<MEDCouplingPointSet *>(const_cast<MEDCouplingMesh *>(field->getMesh())),
        *local_group,
        field->getMesh()->getName()
    );
    ParaFIELD *tmp = new ParaFIELD(field, paramesh, *local_group);
    tmp->setOwnSupport(true);
    tmp->setGlobalNumbering(globalIds);
    this->attachLocalField(tmp, true);
};

void
InterpKernelDECWithOverlap::recvData()
{
    InterpKernelDEC::recvData();
    /* Synchonize ghost if needed */
    _overlap.synchronizeGhosts();
}

void
InterpKernelDECWithOverlap::sendData()
{
    /* Since array can be modified without call to attachLocalField,
     *  we need to update with new values before to send it
     *  even if no modification
     *  Called only if overlapping is detected
     */
    _overlap.updateLocalArray();
    InterpKernelDEC::sendData();
}
