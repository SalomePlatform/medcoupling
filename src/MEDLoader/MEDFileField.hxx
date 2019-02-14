// Copyright (C) 2007-2019  CEA/DEN, EDF R&D
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
// Author : Anthony Geay (EDF R&D)

#ifndef __MEDFILEFIELD_HXX__
#define __MEDFILEFIELD_HXX__

#include "MEDLoaderDefines.hxx"

#include "MEDFileFieldInternal.hxx"
#include "MEDFileFieldGlobs.hxx"
#include "MEDFileField1TS.hxx"
#include "MEDFileFieldMultiTS.hxx"

#include "MEDFileFieldOverView.hxx"
#include "MEDFileUtilities.txx"
#include "MEDFileEntities.hxx"

#include "MCAuto.hxx"
#include "MEDLoaderTraits.hxx"
#include "MEDCouplingTraits.hxx"
#include "MEDCouplingRefCountObject.hxx"
#include "MEDCouplingFieldInt.hxx"
#include "MEDCouplingMemArray.hxx"
#include "MEDCouplingPartDefinition.hxx"

#include "NormalizedUnstructuredMesh.hxx"
#include "InterpKernelException.hxx"

#include <vector>
#include <string>
#include <list>
#include <set>

#include "med.h"

namespace MEDCoupling
{
  class MEDFileFieldGlobs;
  class MEDCouplingMesh;
  class MEDCouplingFieldDouble;
  class MEDFileMesh;
  class MEDFileFieldVisitor;

  class MEDFileMeshes;

  class MEDFileFieldsIterator;
  class MEDFileStructureElements;
  
  /*!
   * Use class.
   */
  class MEDFileFields : public RefCountObject, public MEDFileFieldGlobsReal, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT static MEDFileFields *New();
    MEDLOADER_EXPORT static MEDFileFields *New(const std::string& fileName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileFields *New(med_idt fid, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileFields *NewAdv(const std::string& fileName, bool loadAll, const MEDFileEntities *entities);
    MEDLOADER_EXPORT static MEDFileFields *NewAdv(med_idt fid, bool loadAll, const MEDFileEntities *entities);
    MEDLOADER_EXPORT static MEDFileFields *NewWithDynGT(const std::string& fileName, const MEDFileStructureElements *se, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileFields *NewWithDynGT(med_idt fid, const MEDFileStructureElements *se, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileFields *New(DataArrayByte *db) { return BuildFromMemoryChunk<MEDFileFields>(db); }
    MEDLOADER_EXPORT static MEDFileFields *LoadPartOf(const std::string& fileName, bool loadAll=true, const MEDFileMeshes *ms=0);
    MEDLOADER_EXPORT static MEDFileFields *LoadSpecificEntities(const std::string& fileName, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >& entities, bool loadAll=true);
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT MEDFileFields *deepCopy() const;
    MEDLOADER_EXPORT MEDFileFields *shallowCpy() const;
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    MEDLOADER_EXPORT void loadArrays();
    MEDLOADER_EXPORT void loadArraysIfNecessary();
    MEDLOADER_EXPORT void unloadArrays();
    MEDLOADER_EXPORT void unloadArraysWithoutDataLoss();
    MEDLOADER_EXPORT int getNumberOfFields() const;
    MEDLOADER_EXPORT std::vector< std::pair<int,int> > getCommonIterations(bool& areThereSomeForgottenTS) const;
    MEDLOADER_EXPORT std::vector<std::string> getFieldsNames() const;
    MEDLOADER_EXPORT std::vector<std::string> getMeshesNames() const;
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT void simpleRepr(int bkOffset, std::ostream& oss) const;
    //
    MEDLOADER_EXPORT void resize(int newSize);
    MEDLOADER_EXPORT void pushField(MEDFileAnyTypeFieldMultiTS *field);
    MEDLOADER_EXPORT void pushFields(const std::vector<MEDFileAnyTypeFieldMultiTS *>& fields);
    MEDLOADER_EXPORT void setFieldAtPos(int i, MEDFileAnyTypeFieldMultiTS *field);
    MEDLOADER_EXPORT int getPosFromFieldName(const std::string& fieldName) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *getFieldAtPos(int i) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *getFieldWithName(const std::string& fieldName) const;
    MEDLOADER_EXPORT MEDFileFields *buildSubPart(const int *startIds, const int *endIds) const;
    MEDLOADER_EXPORT bool removeFieldsWithoutAnyTimeStep();
    MEDLOADER_EXPORT MEDFileFields *partOfThisLyingOnSpecifiedMeshName(const std::string& meshName) const;
    MEDLOADER_EXPORT MEDFileFields *partOfThisLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const;
    MEDLOADER_EXPORT MEDFileFields *partOfThisNotLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const;
    MEDLOADER_EXPORT bool presenceOfStructureElements() const;
    MEDLOADER_EXPORT void killStructureElements();
    MEDLOADER_EXPORT void keepOnlyStructureElements();
    MEDLOADER_EXPORT void keepOnlyOnMeshSE(const std::string& meshName, const std::string& seName);
    MEDLOADER_EXPORT void getMeshSENames(std::vector< std::pair<std::string,std::string> >& ps) const;
    MEDLOADER_EXPORT void blowUpSE(MEDFileMeshes *ms, const MEDFileStructureElements *ses);
    MEDLOADER_EXPORT MCAuto<MEDFileFields> partOfThisOnStructureElements() const;
    MEDLOADER_EXPORT MCAuto<MEDFileFields> partOfThisLyingOnSpecifiedMeshSEName(const std::string& meshName, const std::string& seName) const;
    MEDLOADER_EXPORT void aggregate(const MEDFileFields& other);
    MEDLOADER_EXPORT MEDFileFieldsIterator *iterator();
    MEDLOADER_EXPORT void destroyFieldAtPos(int i);
    MEDLOADER_EXPORT void destroyFieldsAtPos(const int *startIds, const int *endIds);
    MEDLOADER_EXPORT void destroyFieldsAtPos2(int bg, int end, int step);
    MEDLOADER_EXPORT bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT bool renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N);
    MEDLOADER_EXPORT void accept(MEDFileFieldVisitor& visitor) const;
    MEDLOADER_EXPORT MCAuto<MEDFileFields> linearToQuadratic(const MEDFileMeshes *oldLin, const MEDFileMeshes *newQuad) const;
  public:
    MEDLOADER_EXPORT MEDFileFields *extractPart(const std::map<int, MCAuto<DataArrayInt> >& extractDef, MEDFileMesh *mm) const;
  public:
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsed() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsed() const;
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsedMulti() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsedMulti() const;
    MEDLOADER_EXPORT void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
  private:
    ~MEDFileFields() { }
    MEDFileFields();
    MEDFileFields(med_idt fid, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities);
  private:
    std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > _fields;
  };

  class MEDFileFieldsIterator
  {
  public:
    MEDLOADER_EXPORT MEDFileFieldsIterator(MEDFileFields *fs);
    MEDLOADER_EXPORT ~MEDFileFieldsIterator();
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *nextt();
  private:
    MCAuto<MEDFileFields> _fs;
    int _iter_id;
    int _nb_iter;
  };
}

#endif
