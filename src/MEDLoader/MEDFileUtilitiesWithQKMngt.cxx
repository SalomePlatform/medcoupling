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

#include "InterpKernelAutoPtr.hxx"

#include <sstream>
#include <fstream>

namespace
{
/*!
 * convertion from medcoupling ids in QuantityKindEnum::MCIdOfValue to MEDFile med_quantity_kind
 * manupulate with care !
 */
med_quantity_kind
QKToMedFileEnum_ForQuantityKindEnum(int idOfMedCouplingQKEnum)
{
    return med_quantity_kind(idOfMedCouplingQKEnum + 1);
}

/*!
 * convertion to medcoupling ids in QuantityKindEnum::MCIdOfValue from MEDFile med_quantity_kind
 * manupulate with care !
 */
int
QKFromMedFileEnum_ToMCId(med_quantity_kind mqk)
{
    if (mqk == MED_PQK_UNDEF || mqk == MED_PQK_USER)
    {
        THROW_IK_EXCEPTION("Invalid input med_quantity_kind enum for conversion");
    }
    return (int)mqk - 1;
}

med_quantity_kind
QKToMedFileEnum(const MEDCoupling::QuantityKindAbstract *qk, std::string &mqkDesc)
{
    if (dynamic_cast<const MEDCoupling::QuantityKindUnDef *>(qk))
    {
        return MED_PQK_UNDEF;
    }
    const MEDCoupling::QuantityKindUser *qkus(dynamic_cast<const MEDCoupling::QuantityKindUser *>(qk));
    if (qkus)
    {
        mqkDesc = qkus->getValue();
        return MED_PQK_USER;
    }
    const MEDCoupling::QuantityKindEnum *qke(dynamic_cast<const MEDCoupling::QuantityKindEnum *>(qk));
    if (qke)
    {
        mqkDesc = qke->getDescription();
        med_quantity_kind ret(
            QKToMedFileEnum_ForQuantityKindEnum(MEDCoupling::QuantityKindEnum::MCIdOfValue(qke->getValue()))
        );
        med_quantity_kind ret2(MEDphysicalQuantityKindFromStr(qke->getValue().c_str()));
        if (ret == ret2)
        {
            return ret;
        }
        const char *medfileStr(MEDphysicalQuantityKindStr(ret));
        if (!medfileStr)
        {
            return MED_PQK_UNDEF;
        }
        if (qke->getValue() != std::string(medfileStr))
        {  // smells bad. See QKToMedFileEnum_ForQuantityKindEnum implementation
            return MED_PQK_UNDEF;
        }
        return ret;
    }
    else
    {
        return MED_PQK_UNDEF;
    }
}

MEDCoupling::MCAuto<MEDCoupling::QuantityKindAbstract>
QKFromMedFile(med_quantity_kind mqk, const std::string &value)
{
    switch (mqk)
    {
        case MED_PQK_UNDEF:
        {
            MEDCoupling::MCAuto<MEDCoupling::QuantityKindUnDef> ret(MEDCoupling::QuantityKindUnDef::New());
            return MEDCoupling::StaticCast<MEDCoupling::QuantityKindUnDef, MEDCoupling::QuantityKindAbstract>(
                std::move(ret)
            );
        }
        case MED_PQK_USER:
        {
            return MEDCoupling::StaticCast<MEDCoupling::QuantityKindUser, MEDCoupling::QuantityKindAbstract>(
                MEDCoupling::QuantityKindUser::New(value)
            );
        }
        default:
        {
            int id(QKFromMedFileEnum_ToMCId(mqk));
            MEDCoupling::MCAuto<MEDCoupling::QuantityKindEnum> ret(
                MEDCoupling::QuantityKindEnum::New(MEDCoupling::QuantityKindEnum::AllowedValuesAt(id))
            );
            ret->setDescription(value);
            return MEDCoupling::StaticCast<MEDCoupling::QuantityKindEnum, MEDCoupling::QuantityKindAbstract>(
                std::move(ret)
            );
        }
    }
}

std::string
FromMCQKEnumStrToFrenchDescr(const std::string &valueEnum)
{
    med_quantity_kind ret(QKToMedFileEnum_ForQuantityKindEnum(MEDCoupling::QuantityKindEnum::MCIdOfValue(valueEnum)));
    return MEDphysicalQuantityKindDescription(ret);
}

std::string
FromMCQKEnumStrToQUDTUri(const std::string &valueEnum)
{
    med_quantity_kind ret(QKToMedFileEnum_ForQuantityKindEnum(MEDCoupling::QuantityKindEnum::MCIdOfValue(valueEnum)));
    return MEDphysicalQuantityKindQUDTUri(ret);
}

}  // namespace

void
MEDFileUtilities::WrapperOf_MEDfieldQuantityKindWr(
    med_idt fid,
    const std::string &fieldName,
    const MEDCoupling::QuantityKindAbstract *qk,
    const MEDCoupling::MEDFileWritable &opts
)
{
    if (!qk)
        return;
    INTERP_KERNEL::AutoPtr<char> nomcha(MEDLoaderBase::buildEmptyString(MED_NAME_SIZE));
    MEDLoaderBase::safeStrCpy(fieldName.c_str(), MED_NAME_SIZE, nomcha, opts.getTooLongStrPolicy());
    std::string mqkDesc;
    med_quantity_kind mqk(QKToMedFileEnum(qk, mqkDesc));
    if (mqk != MED_PQK_UNDEF)
    {  // see src/ci/MEDfieldQuantityKindWr.c:67. If mqk==MED_PQK_UNDEF do not call MEDfieldQuantityKindWr
        INTERP_KERNEL::AutoPtr<char> nomdesc(MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE));
        MEDLoaderBase::safeStrCpy(mqkDesc.c_str(), MED_COMMENT_SIZE, nomdesc, opts.getTooLongStrPolicy());
        MEDFILESAFECALLERWR0(MEDfieldQuantityKindWr, (fid, nomcha, mqk, nomdesc));
    }
}

void
MEDFileUtilities::WrapperOf_MEDfieldQuantityKindRd(
    med_idt fid, const std::string &fieldName, MEDCoupling::MCAuto<MEDCoupling::QuantityKindAbstract> &qk
)
{
    med_quantity_kind mqk;
    INTERP_KERNEL::AutoPtr<char> nomdesc(MEDLoaderBase::buildEmptyString(MED_COMMENT_SIZE));
    MEDFILESAFECALLERRD0(MEDfieldQuantityKindRd, (fid, fieldName.c_str(), &mqk, nomdesc));
    std::string qkValue(MEDLoaderBase::buildStringFromFortran(nomdesc, MED_COMMENT_SIZE));
    qk = QKFromMedFile(mqk, qkValue);
}

std::string
MEDFileUtilities::WrapperOf_ToFrenchDescr(MEDCoupling::MCAuto<MEDCoupling::QuantityKindAbstract> obj)
{
    if (obj.isNull())
    {
        THROW_IK_EXCEPTION("Nullptr QuantityKind : no description");
    }
    MEDCoupling::QuantityKindUnDef *obj1(
        dynamic_cast<MEDCoupling::QuantityKindUnDef *>((MEDCoupling::QuantityKindAbstract *)obj)
    );
    if (obj1)
    {
        return std::string();
    }
    MEDCoupling::QuantityKindUser *obj2(
        dynamic_cast<MEDCoupling::QuantityKindUser *>((MEDCoupling::QuantityKindAbstract *)obj)
    );
    if (obj2)
    {
        return std::string();
    }
    MEDCoupling::QuantityKindEnum *obj3(
        dynamic_cast<MEDCoupling::QuantityKindEnum *>((MEDCoupling::QuantityKindAbstract *)obj)
    );
    if (obj3)
    {
        return FromMCQKEnumStrToFrenchDescr(obj3->getValue());
    }
    THROW_IK_EXCEPTION("Not managed QuantityKind class for FrenchDescription.");
}

std::string
MEDFileUtilities::WrapperOf_ToQUDTUri(MEDCoupling::MCAuto<MEDCoupling::QuantityKindAbstract> obj)
{
    if (obj.isNull())
    {
        THROW_IK_EXCEPTION("Nullptr QuantityKind : no description");
    }
    MEDCoupling::QuantityKindUnDef *obj1(
        dynamic_cast<MEDCoupling::QuantityKindUnDef *>((MEDCoupling::QuantityKindAbstract *)obj)
    );
    if (obj1)
    {
        THROW_IK_EXCEPTION("No QUDTUri for UndefKind.");
    }
    MEDCoupling::QuantityKindUser *obj2(
        dynamic_cast<MEDCoupling::QuantityKindUser *>((MEDCoupling::QuantityKindAbstract *)obj)
    );
    if (obj2)
    {
        THROW_IK_EXCEPTION("No QUDTUri for UserKind.");
    }
    MEDCoupling::QuantityKindEnum *obj3(
        dynamic_cast<MEDCoupling::QuantityKindEnum *>((MEDCoupling::QuantityKindAbstract *)obj)
    );
    if (obj3)
    {
        return FromMCQKEnumStrToQUDTUri(obj3->getValue());
    }
    THROW_IK_EXCEPTION("Not managed QuantityKind class for QUDTUri.");
}
