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
// Author : Anthony Geay (CEA/DEN)

#ifndef __MEDFILEMESHELT_HXX__
#define __MEDFILEMESHELT_HXX__

#include "MEDCouplingMemArray.hxx"
#include "MEDCoupling1GTUMesh.hxx"
#include "MEDCouplingPartDefinition.hxx"
#include "MCAuto.hxx"

#include "NormalizedUnstructuredMesh.hxx"

#include "med.h"

namespace MEDCoupling
{
class MEDCouplingUMesh;
class MEDFileMeshReadSelector;

class MEDFileUMeshPerTypeCommon : public RefCountObject
{
   public:
    static MEDFileUMeshPerTypeCommon *New();
    std::string getClassName() const override { return std::string("MEDFileUMeshPerTypeCommon"); }
    void loadCommonPart(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        mcIdType curNbOfElem,
        med_geometry_type geoElt,
        med_entity_type entity,
        MEDFileMeshReadSelector *mrs
    );
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    const DataArrayIdType *getFam() const { return _fam; }
    const DataArrayIdType *getNum() const { return _num; }
    const DataArrayAsciiChar *getNames() const { return _names; }

   protected:
    MCAuto<DataArrayIdType> _num;
    MCAuto<DataArrayIdType> _fam;
    MCAuto<DataArrayAsciiChar> _names;
};

class MEDFileUMeshPerType : public MEDFileUMeshPerTypeCommon
{
   public:
    static MEDFileUMeshPerType *New(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        int mdim,
        med_geometry_type geoElt,
        INTERP_KERNEL::NormalizedCellType geoElt2,
        MEDFileMeshReadSelector *mrs
    );
    static MEDFileUMeshPerType *NewPart(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        int mdim,
        INTERP_KERNEL::NormalizedCellType geoElt2,
        mcIdType strt,
        mcIdType stp,
        mcIdType step,
        MEDFileMeshReadSelector *mrs
    );
    static MEDFileUMeshPerType *NewPart(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        int mdim,
        INTERP_KERNEL::NormalizedCellType geoElt2,
        med_geometry_type geoElt,
        med_entity_type whichEntity,
        const std::vector<mcIdType> &distrib,
        MEDFileMeshReadSelector *mrs
    );

    std::string getClassName() const override { return std::string("MEDFileUMeshPerType"); }
    static bool isExisting(
        med_idt fid, const char *mName, int dt, int it, med_geometry_type geoElt, med_entity_type &whichEntity
    );
    static void getMedTypes(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        INTERP_KERNEL::NormalizedCellType MCgeoElt,
        med_geometry_type &geoElt,
        med_entity_type &whichEntity
    );
    std::size_t getHeapMemorySizeWithoutChildren() const;
    std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    int getDim() const;
    MEDCoupling1GTUMesh *getMesh() const { return const_cast<MEDCoupling1GTUMesh *>((const MEDCoupling1GTUMesh *)_m); }
    const PartDefinition *getPartDef() const { return _pd; }
    static void Write(
        med_idt fid,
        const std::string &mname,
        int mdim,
        const MEDCoupling1GTUMesh *m,
        const DataArrayIdType *fam,
        const DataArrayIdType *num,
        const DataArrayAsciiChar *names
    );

   private:
    MEDFileUMeshPerType();
    MEDFileUMeshPerType(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        int mdim,
        med_geometry_type geoElt,
        INTERP_KERNEL::NormalizedCellType type,
        med_entity_type entity,
        MEDFileMeshReadSelector *mrs
    );
    void loadPart(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        int mdim,
        med_geometry_type geoElt,
        INTERP_KERNEL::NormalizedCellType type,
        med_entity_type entity,
        mcIdType strt,
        mcIdType end,
        mcIdType step,
        MEDFileMeshReadSelector *mrs
    );
    void loadPart(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        int mdim,
        med_geometry_type geoElt,
        INTERP_KERNEL::NormalizedCellType type,
        med_entity_type entity,
        const std::vector<mcIdType> &distrib,
        MEDFileMeshReadSelector *mrs
    );
    void loadFromStaticType(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        int mdim,
        mcIdType curNbOfElem,
        med_geometry_type geoElt,
        INTERP_KERNEL::NormalizedCellType type,
        med_entity_type entity,
        MEDFileMeshReadSelector *mrs
    );
    void loadPartStaticType(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        int mdim,
        mcIdType curNbOfElem,
        med_geometry_type geoElt,
        INTERP_KERNEL::NormalizedCellType type,
        med_entity_type entity,
        MEDFileMeshReadSelector *mrs
    );
    void loadPolyg(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        int mdim,
        mcIdType arraySize,
        med_geometry_type geoElt,
        med_entity_type entity,
        MEDFileMeshReadSelector *mrs
    );
    void loadPolyh(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        int mdim,
        mcIdType connFaceLgth,
        med_geometry_type geoElt,
        med_entity_type entity,
        MEDFileMeshReadSelector *mrs
    );
    void loadPartOfCellCommonPart(
        med_idt fid,
        const char *mName,
        int dt,
        int it,
        int mdim,
        mcIdType curNbOfElem,
        med_geometry_type geoElt,
        med_entity_type entity,
        MEDFileMeshReadSelector *mrs
    );

   private:
    MCAuto<MEDCoupling1GTUMesh> _m;
    MCAuto<PartDefinition> _pd;
};
}  // namespace MEDCoupling
#endif
