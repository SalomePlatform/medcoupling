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

#include "MEDFileUtilities.hxx"
#include "MEDFileSafeCaller.txx"
#include "MEDLoaderBase.hxx"
#include "MEDLoader.hxx"
#include "MEDFileBasis.hxx"

using namespace MEDCoupling;

void
MEDFileUtilities::WrapperOf_MEDfieldQuantityKindWr(
    med_idt fid,
    const std::string &fieldName,
    const MEDCoupling::QuantityKindAbstract *qk,
    const MEDCoupling::MEDFileWritable &opts
)
{
}

void
MEDFileUtilities::WrapperOf_MEDfieldQuantityKindRd(
    med_idt fid, const std::string &fieldName, MEDCoupling::MCAuto<MEDCoupling::QuantityKindAbstract> &qk
)
{
    qk = QuantityKindUnDef::New().retn();
}
