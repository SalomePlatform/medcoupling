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

#ifndef __MEDFILEFIELDMULTITS_HXX__
#define __MEDFILEFIELDMULTITS_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDFileField1TS.hxx"
#include "MEDFileFieldGlobs.hxx"
#include "MEDLoaderTraits.hxx"
#include "MEDFileUtilities.hxx"

namespace MEDCoupling
{
  class MEDFileMesh;
  class MEDFileMeshes;
  class MEDCouplingMesh;
  class MEDFileFieldVisitor;
  class MEDFileAnyTypeField1TS;
  class MEDFileAnyTypeField1TSWithoutSDA;
  
  class MEDFileAnyTypeFieldMultiTSWithoutSDA : public RefCountObject, public MEDFileFieldNameScope
  {
  protected:
    MEDFileAnyTypeFieldMultiTSWithoutSDA();
    MEDFileAnyTypeFieldMultiTSWithoutSDA(const std::string& fieldName, const std::string& meshName);
    MEDFileAnyTypeFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities);
    MEDFileAnyTypeFieldMultiTSWithoutSDA(med_idt fid, const std::string& fieldName, const std::string& meshName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities);
  public:
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeFieldMultiTSWithoutSDA *deepCopy() const;
    MEDLOADER_EXPORT virtual std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > splitComponents() const;
    MEDLOADER_EXPORT virtual std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > splitDiscretizations() const;
    MEDLOADER_EXPORT virtual std::vector< MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> > splitMultiDiscrPerGeoTypes() const;
    MEDLOADER_EXPORT virtual const char *getTypeStr() const = 0;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeFieldMultiTSWithoutSDA *shallowCpy() const = 0;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeFieldMultiTSWithoutSDA *createNew() const = 0;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeField1TSWithoutSDA *createNew1TSWithoutSDAEmptyInstance() const = 0;
    MEDLOADER_EXPORT virtual void checkCoherencyOfType(const MEDFileAnyTypeField1TSWithoutSDA *f1ts) const = 0;
    MEDLOADER_EXPORT const std::vector<std::string>& getInfo() const;
    MEDLOADER_EXPORT bool presenceOfMultiDiscPerGeoType() const;
    MEDLOADER_EXPORT void setInfo(const std::vector<std::string>& info);
    MEDLOADER_EXPORT int getTimeStepPos(int iteration, int order) const;
    MEDLOADER_EXPORT const MEDFileAnyTypeField1TSWithoutSDA& getTimeStepEntry(int iteration, int order) const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA& getTimeStepEntry(int iteration, int order);
    MEDLOADER_EXPORT bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT int getNumberOfTS() const;
    MEDLOADER_EXPORT void eraseEmptyTS();
    MEDLOADER_EXPORT void eraseTimeStepIds(const int *startIds, const int *endIds);
    MEDLOADER_EXPORT void eraseTimeStepIds2(int bg, int end, int step);
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *buildFromTimeStepIds(const int *startIds, const int *endIds) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *buildFromTimeStepIds2(int bg, int end, int step) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *partOfThisLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *partOfThisNotLyingOnSpecifiedTimeSteps(const std::vector< std::pair<int,int> >& timeSteps) const;
    MEDLOADER_EXPORT bool presenceOfStructureElements() const;
    MEDLOADER_EXPORT bool onlyStructureElements() const;
    MEDLOADER_EXPORT void killStructureElements();
    MEDLOADER_EXPORT void keepOnlyStructureElements();
    MEDLOADER_EXPORT void keepOnlyOnSE(const std::string& seName);
    MEDLOADER_EXPORT void getMeshSENames(std::vector< std::pair<std::string,std::string> >& ps) const;
    MEDLOADER_EXPORT int getPosOfTimeStep(int iteration, int order) const;
    MEDLOADER_EXPORT int getPosGivenTime(double time, double eps=1e-8) const;
    MEDLOADER_EXPORT std::vector< std::pair<int,int> > getIterations() const;
    MEDLOADER_EXPORT std::vector< std::pair<int,int> > getTimeSteps(std::vector<double>& ret1) const;
    MEDLOADER_EXPORT void pushBackTimeStep(MCAuto<MEDFileAnyTypeField1TSWithoutSDA>& tse);
    MEDLOADER_EXPORT void synchronizeNameScope();
    MEDLOADER_EXPORT void simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const;
    MEDLOADER_EXPORT int getNonEmptyLevels(int iteration, int order, const std::string& mname, std::vector<int>& levs) const;
    MEDLOADER_EXPORT void appendFieldNoProfileSBT(const MEDCouplingFieldDouble *field, const DataArray *arr, MEDFileFieldGlobsReal& glob);
    MEDLOADER_EXPORT void appendFieldProfile(const MEDCouplingFieldDouble *field, const DataArray *arr, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob, bool smartPflKiller);
    MEDLOADER_EXPORT std::vector< std::vector< std::pair<int,int> > > getFieldSplitedByType(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    MEDLOADER_EXPORT std::vector< std::vector<TypeOfField> > getTypesOfFieldAvailable() const;
    MEDLOADER_EXPORT DataArray *getUndergroundDataArray(int iteration, int order) const;
    MEDLOADER_EXPORT DataArray *getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDLOADER_EXPORT bool renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N, MEDFileFieldGlobsReal& glob);
    MEDLOADER_EXPORT void accept(MEDFileFieldVisitor& visitor) const;
    MEDLOADER_EXPORT void loadStructureOrStructureAndBigArraysRecursively(med_idt fid, int nbPdt, med_field_type fieldTyp, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities);
    MEDLOADER_EXPORT void writeLL(med_idt fid, const MEDFileWritable& opts) const;
    MEDLOADER_EXPORT void loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc);
    MEDLOADER_EXPORT void loadBigArraysRecursivelyIfNecessary(med_idt fid, const MEDFileFieldNameScope& nasc);
    MEDLOADER_EXPORT void unloadArrays();
  public:
    MEDLOADER_EXPORT const MEDFileAnyTypeField1TSWithoutSDA *getTimeStepAtPos2(int pos) const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA *getTimeStepAtPos2(int pos);
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsed2() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsed2() const;
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsedMulti2() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsedMulti2() const;
    MEDLOADER_EXPORT void changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void setIteration(int i, MCAuto<MEDFileAnyTypeField1TSWithoutSDA> ts);
  protected:
    virtual med_field_type getMEDFileFieldType() const = 0;
    void copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr);
    void checkCoherencyOfTinyInfo(const MEDCouplingFieldDouble *field, const DataArray *arr) const;
    void checkThatComponentsMatch(const std::vector<std::string>& compos) const;
    void checkThatNbOfCompoOfTSMatchThis() const;
  protected:
    std::vector<std::string> _infos;
    std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > _time_steps;
  };

  class MEDFileIntFieldMultiTSWithoutSDA;

  template<class T>
  class MEDFileTemplateFieldMultiTSWithoutSDA : public MEDFileAnyTypeFieldMultiTSWithoutSDA
  {
  public:
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::FMTSWSDAType *New(med_idt fid, const std::string& fieldName, const std::string& meshName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities);
    MEDLOADER_EXPORT const char *getTypeStr() const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *createNew() const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA *createNew1TSWithoutSDAEmptyInstance() const;
  protected:
    MEDFileTemplateFieldMultiTSWithoutSDA() { }
    MEDFileTemplateFieldMultiTSWithoutSDA(const std::string& fieldName, const std::string& meshName):MEDFileAnyTypeFieldMultiTSWithoutSDA(fieldName,meshName) { }
    /** \param [in] fieldId field id in C mode */
    MEDFileTemplateFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities):MEDFileAnyTypeFieldMultiTSWithoutSDA(fid,fieldId,loadAll,ms,entities) { }
    MEDFileTemplateFieldMultiTSWithoutSDA(med_idt fid, const std::string& fieldName, const std::string& meshName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities):MEDFileAnyTypeFieldMultiTSWithoutSDA(fid,fieldName,meshName,fieldTyp,infos,nbOfStep,dtunit,loadAll,ms,entities) { }
    void checkCoherencyOfType(const MEDFileAnyTypeField1TSWithoutSDA *f1ts) const;
  };
  
  class MEDFileFieldMultiTSWithoutSDA : public MEDFileTemplateFieldMultiTSWithoutSDA<double>
  {
    friend class MEDFileTemplateFieldMultiTSWithoutSDA<double>;
  public:
    MEDLOADER_EXPORT MEDFileFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities):MEDFileTemplateFieldMultiTSWithoutSDA<double>(fid,fieldId,loadAll,ms,entities) { }
    MEDLOADER_EXPORT std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    MEDLOADER_EXPORT MEDFileIntFieldMultiTSWithoutSDA *convertToInt() const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *shallowCpy() const { return new MEDFileFieldMultiTSWithoutSDA(*this); }
  protected:
    MEDFileFieldMultiTSWithoutSDA(const std::string& fieldName, const std::string& meshName):MEDFileTemplateFieldMultiTSWithoutSDA<double>(fieldName,meshName) { }
    MEDFileFieldMultiTSWithoutSDA(med_idt fid, const std::string& fieldName, const std::string& meshName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities):MEDFileTemplateFieldMultiTSWithoutSDA<double>(fid,fieldName,meshName,fieldTyp,infos,nbOfStep,dtunit,loadAll,ms,entities) { }
    med_field_type getMEDFileFieldType() const { return MED_FLOAT64; }
  public:
    MEDLOADER_EXPORT MEDFileFieldMultiTSWithoutSDA() { }
  };

  template<class T>
  class MEDFileNDTemplateFieldMultiTSWithoutSDA : public MEDFileTemplateFieldMultiTSWithoutSDA<T>
  {
  public:
    MEDLOADER_EXPORT MEDFileFieldMultiTSWithoutSDA *convertToDouble() const;
  protected:
    MEDFileNDTemplateFieldMultiTSWithoutSDA() { }
    MEDFileNDTemplateFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities):MEDFileTemplateFieldMultiTSWithoutSDA<T>(fid,fieldId,loadAll,ms,entities) { }
    MEDFileNDTemplateFieldMultiTSWithoutSDA(const std::string& fieldName, const std::string& meshName):MEDFileTemplateFieldMultiTSWithoutSDA<T>(fieldName,meshName) { }
    MEDFileNDTemplateFieldMultiTSWithoutSDA(med_idt fid, const std::string& fieldName, const std::string& meshName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities):MEDFileTemplateFieldMultiTSWithoutSDA<T>(fid,fieldName,meshName,fieldTyp,infos,nbOfStep,dtunit,loadAll,ms,entities) { }
  };

  class MEDFileIntFieldMultiTSWithoutSDA : public MEDFileNDTemplateFieldMultiTSWithoutSDA<int>
  {
    friend class MEDFileTemplateFieldMultiTSWithoutSDA<int>;
  public:
    MEDLOADER_EXPORT MEDFileIntFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities):MEDFileNDTemplateFieldMultiTSWithoutSDA<int>(fid,fieldId,loadAll,ms,entities) { }
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *shallowCpy() const { return new MEDFileIntFieldMultiTSWithoutSDA(*this); }
  protected:
    MEDFileIntFieldMultiTSWithoutSDA(const std::string& fieldName, const std::string& meshName):MEDFileNDTemplateFieldMultiTSWithoutSDA<int>(fieldName,meshName) { }
    MEDFileIntFieldMultiTSWithoutSDA(med_idt fid, const std::string& fieldName, const std::string& meshName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities):MEDFileNDTemplateFieldMultiTSWithoutSDA<int>(fid,fieldName,meshName,fieldTyp,infos,nbOfStep,dtunit,loadAll,ms,entities) { }
    med_field_type getMEDFileFieldType() const { return MED_INT32; }
  public:
    MEDLOADER_EXPORT MEDFileIntFieldMultiTSWithoutSDA() { }
  };

  class MEDFileFloatFieldMultiTSWithoutSDA : public MEDFileNDTemplateFieldMultiTSWithoutSDA<float>
  {
    friend class MEDFileTemplateFieldMultiTSWithoutSDA<float>;
  public:
    MEDLOADER_EXPORT MEDFileFloatFieldMultiTSWithoutSDA(med_idt fid, int fieldId, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities):MEDFileNDTemplateFieldMultiTSWithoutSDA<float>(fid,fieldId,loadAll,ms,entities) { }
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSWithoutSDA *shallowCpy() const { return new MEDFileFloatFieldMultiTSWithoutSDA(*this); }
  protected:
    MEDFileFloatFieldMultiTSWithoutSDA(const std::string& fieldName, const std::string& meshName):MEDFileNDTemplateFieldMultiTSWithoutSDA<float>(fieldName,meshName) { }
    MEDFileFloatFieldMultiTSWithoutSDA(med_idt fid, const std::string& fieldName, const std::string& meshName, med_field_type fieldTyp, const std::vector<std::string>& infos, int nbOfStep, const std::string& dtunit, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities):MEDFileNDTemplateFieldMultiTSWithoutSDA<float>(fid,fieldName,meshName,fieldTyp,infos,nbOfStep,dtunit,loadAll,ms,entities) { }
    med_field_type getMEDFileFieldType() const { return MED_FLOAT32; }
  public:
    MEDLOADER_EXPORT MEDFileFloatFieldMultiTSWithoutSDA() { }
  };

  class MEDFileAnyTypeFieldMultiTSIterator;
  class MEDFileFastCellSupportComparator;
  /*!
   * User class.
   */
  class MEDFileAnyTypeFieldMultiTS : public RefCountObject, public MEDFileWritableStandAlone, public MEDFileFieldGlobsReal
  {
  protected:
    MEDFileAnyTypeFieldMultiTS();
    MEDFileAnyTypeFieldMultiTS(med_idt fid, bool loadAll, const MEDFileMeshes *ms);
    MEDFileAnyTypeFieldMultiTS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities=0);
    MEDFileAnyTypeFieldMultiTS(const MEDFileAnyTypeFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent);
    static MEDFileAnyTypeFieldMultiTS *BuildNewInstanceFromContent(MEDFileAnyTypeFieldMultiTSWithoutSDA *c, med_idt fid);
    static MEDFileAnyTypeFieldMultiTSWithoutSDA *BuildContentFrom(med_idt fid, bool loadAll, const MEDFileMeshes *ms);
    static MEDFileAnyTypeFieldMultiTSWithoutSDA *BuildContentFrom(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities);
  public:
    MEDLOADER_EXPORT static MEDFileAnyTypeFieldMultiTS *New(const std::string& fileName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeFieldMultiTS *New(med_idt fid, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeFieldMultiTS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeFieldMultiTS *New(med_idt fid, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeFieldMultiTS *BuildNewInstanceFromContent(MEDFileAnyTypeFieldMultiTSWithoutSDA *c);
    MEDLOADER_EXPORT void loadArrays();
    MEDLOADER_EXPORT void loadArraysIfNecessary();
    MEDLOADER_EXPORT void unloadArrays();
    MEDLOADER_EXPORT void unloadArraysWithoutDataLoss();
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeFieldMultiTS *deepCopy() const;
    MEDLOADER_EXPORT std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > splitComponents() const;
    MEDLOADER_EXPORT std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > splitDiscretizations() const;
    MEDLOADER_EXPORT std::vector< MCAuto< MEDFileAnyTypeFieldMultiTS > > splitMultiDiscrPerGeoTypes() const;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeFieldMultiTS *shallowCpy() const = 0;
    MEDLOADER_EXPORT virtual void checkCoherencyOfType(const MEDFileAnyTypeField1TS *f1ts) const = 0;
    //
    MEDLOADER_EXPORT virtual MEDFileAnyTypeField1TS *getTimeStepAtPos(int pos) const = 0;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *getTimeStep(int iteration, int order) const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *getTimeStepGivenTime(double time, double eps=1e-8) const;
    MEDLOADER_EXPORT static std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > SplitIntoCommonTimeSeries(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS);
    MEDLOADER_EXPORT static std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > SplitPerCommonSupport(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS, const MEDFileMesh *mesh, std::vector< MCAuto<MEDFileFastCellSupportComparator> >& fsc);
    MEDLOADER_EXPORT static int CheckSupportAcrossTime(MEDFileAnyTypeFieldMultiTS *f0, MEDFileAnyTypeFieldMultiTS *f1, const MEDFileMesh *mesh, TypeOfField& tof0, TypeOfField& tof1);
  public:// direct forwarding to MEDFileField1TSWithoutSDA instance _content
    MEDLOADER_EXPORT std::string getName() const;
    MEDLOADER_EXPORT void setName(const std::string& name);
    MEDLOADER_EXPORT std::string getDtUnit() const;
    MEDLOADER_EXPORT void setDtUnit(const std::string& dtUnit);
    MEDLOADER_EXPORT std::string getMeshName() const;
    MEDLOADER_EXPORT void setMeshName(const std::string& newMeshName);
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT void simpleRepr(int bkOffset, std::ostream& oss, int fmtsId) const;
    MEDLOADER_EXPORT int getNumberOfTS() const;
    MEDLOADER_EXPORT void eraseEmptyTS();
    MEDLOADER_EXPORT void eraseTimeStepIds(const int *startIds, const int *endIds);
    MEDLOADER_EXPORT void eraseTimeStepIds2(int bg, int end, int step);
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *buildSubPart(const int *startIds, const int *endIds) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *buildSubPartSlice(int bg, int end, int step) const;
    MEDLOADER_EXPORT std::vector< std::pair<int,int> > getTimeSteps(std::vector<double>& ret1) const;
    MEDLOADER_EXPORT std::vector< std::pair<int,int> > getIterations() const;
    MEDLOADER_EXPORT void pushBackTimeSteps(const std::vector<MEDFileAnyTypeField1TS *>& f1ts);
    MEDLOADER_EXPORT void pushBackTimeSteps(MEDFileAnyTypeFieldMultiTS *fmts);
    MEDLOADER_EXPORT void pushBackTimeStep(MEDFileAnyTypeField1TS *f1ts);
    MEDLOADER_EXPORT void synchronizeNameScope();
    MEDLOADER_EXPORT int getPosOfTimeStep(int iteration, int order) const;
    MEDLOADER_EXPORT int getPosGivenTime(double time, double eps=1e-8) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSIterator *iterator();
    MEDLOADER_EXPORT bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT const std::vector<std::string>& getInfo() const;
    MEDLOADER_EXPORT bool presenceOfMultiDiscPerGeoType() const;
    MEDLOADER_EXPORT void setInfo(const std::vector<std::string>& info);
    MEDLOADER_EXPORT int getNumberOfComponents() const;
    MEDLOADER_EXPORT int getNonEmptyLevels(int iteration, int order, const std::string& mname, std::vector<int>& levs) const;
    MEDLOADER_EXPORT std::vector< std::vector<TypeOfField> > getTypesOfFieldAvailable() const;
    MEDLOADER_EXPORT std::vector< std::vector< std::pair<int,int> > > getFieldSplitedByType(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    MEDLOADER_EXPORT MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> getContent();
  public:
    MEDLOADER_EXPORT virtual MEDFileAnyTypeFieldMultiTS *buildNewEmpty() const = 0;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeFieldMultiTS *extractPart(const std::map<int, MCAuto<DataArrayInt> >& extractDef, MEDFileMesh *mm) const = 0;
    MEDLOADER_EXPORT static MCAuto<MEDFileAnyTypeFieldMultiTS> Aggregate(const std::vector<const MEDFileAnyTypeFieldMultiTS *>& fmtss, const std::vector< std::vector< std::pair<int,int> > >& dts);
  public:
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsed() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsed() const;
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsedMulti() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsedMulti() const;
    MEDLOADER_EXPORT void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
  protected:
    MEDFileAnyTypeFieldMultiTSWithoutSDA *contentNotNullBase();
    const MEDFileAnyTypeFieldMultiTSWithoutSDA *contentNotNullBase() const;
  private:
    static std::vector< std::vector<MEDFileAnyTypeFieldMultiTS *> > SplitPerCommonSupportNotNodesAlg(const std::vector<MEDFileAnyTypeFieldMultiTS *>& vectFMTS, const MEDFileMesh *mesh, std::vector< MCAuto<MEDFileFastCellSupportComparator> >& cmps);
  protected:
    MCAuto<MEDFileAnyTypeFieldMultiTSWithoutSDA> _content;
  };

  template<class T>
  class MEDFileTemplateFieldMultiTS : public MEDFileAnyTypeFieldMultiTS
  {
  public:
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::FMTSType *New();
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::FMTSType *New(const std::string& fileName, bool loadAll=true);
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::FMTSType *New(med_idt fid, bool loadAll=true);
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::FMTSType *New(DataArrayByte *db) { return BuildFromMemoryChunk<typename MLFieldTraits<T>::FMTSType>(db); }
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::FMTSType *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::FMTSType *New(med_idt fid, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::FMTSType *New(const typename MLFieldTraits<T>::FMTSWSDAType& other, bool shallowCopyOfContent);
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::FMTSType *LoadSpecificEntities(const std::string& fileName, const std::string& fieldName, const std::vector< std::pair<TypeOfField,INTERP_KERNEL::NormalizedCellType> >& entities, bool loadAll=true);
    MEDLOADER_EXPORT typename MLFieldTraits<T>::FMTSType *extractPartImpl(const std::map<int, MCAuto<DataArrayInt> >& extractDef, MEDFileMesh *mm) const;
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *extractPart(const std::map<int, MCAuto<DataArrayInt> >& extractDef, MEDFileMesh *mm) const { return this->extractPartImpl(extractDef,mm); }
    //
    MEDLOADER_EXPORT typename Traits<T>::FieldType *field(int iteration, int order, const MEDFileMesh *mesh) const;
    MEDLOADER_EXPORT typename Traits<T>::FieldType *getFieldAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, int renumPol=0) const;
    MEDLOADER_EXPORT typename Traits<T>::FieldType *getFieldAtTopLevel(TypeOfField type, int iteration, int order, int renumPol=0) const;
    MEDLOADER_EXPORT typename Traits<T>::FieldType *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const;
    MEDLOADER_EXPORT typename Traits<T>::FieldType *getFieldOnMeshAtLevel(TypeOfField type, int iteration, int order, const MEDCouplingMesh *mesh, int renumPol=0) const;
    MEDLOADER_EXPORT typename Traits<T>::FieldType *getFieldAtLevelOld(TypeOfField type, int iteration, int order, const std::string& mname, int meshDimRelToMax, int renumPol=0) const;
    MEDLOADER_EXPORT typename Traits<T>::ArrayType *getFieldWithProfile(TypeOfField type, int iteration, int order, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const;
    //
    MEDLOADER_EXPORT void appendFieldNoProfileSBT(const typename Traits<T>::FieldType *field);
    MEDLOADER_EXPORT void appendFieldProfile(const typename Traits<T>::FieldType *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    MEDLOADER_EXPORT void appendFieldProfileFlatly(const typename Traits<T>::FieldType *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    //
    MEDLOADER_EXPORT typename MLFieldTraits<T>::F1TSType *getTimeStepAtPos(int pos) const;
    MEDLOADER_EXPORT typename Traits<T>::ArrayType *getUndergroundDataArray(int iteration, int order) const;
    MEDLOADER_EXPORT typename Traits<T>::ArrayType *getUndergroundDataArrayExt(int iteration, int order, std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDLOADER_EXPORT typename MLFieldTraits<T>::FMTSType *buildNewEmptyImpl() const;
    MEDLOADER_EXPORT void checkCoherencyOfType(const MEDFileAnyTypeField1TS *f1ts) const;
  protected:
    const typename MLFieldTraits<T>::FMTSWSDAType *contentNotNull() const;
    typename MLFieldTraits<T>::FMTSWSDAType *contentNotNull();
    void appendFieldProfileGeneral(const typename Traits<T>::FieldType *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, bool smartPflKiller);
  protected:
    ~MEDFileTemplateFieldMultiTS() { }
    MEDFileTemplateFieldMultiTS();
    MEDFileTemplateFieldMultiTS(const typename MLFieldTraits<T>::FMTSWSDAType& other, bool shallowCopyOfContent);
    MEDFileTemplateFieldMultiTS(med_idt fid, bool loadAll, const MEDFileMeshes *ms);
    MEDFileTemplateFieldMultiTS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities=0);
  };
  
  class MEDFileIntFieldMultiTS;

  /*!
   * User class.
   */
  class MEDFileFieldMultiTS : public MEDFileTemplateFieldMultiTS<double>
  {
    friend class MEDFileTemplateFieldMultiTS<double>;
  public:
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *shallowCpy() const;
    MEDLOADER_EXPORT MEDFileIntFieldMultiTS *convertToInt(bool isDeepCpyGlobs=true) const;
    //
    MEDLOADER_EXPORT std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(int iteration, int order, const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    MEDLOADER_EXPORT MEDFileFieldMultiTS *buildNewEmpty() const { return buildNewEmptyImpl(); }
  public:
  private:
    ~MEDFileFieldMultiTS() { }
    MEDFileFieldMultiTS() { }
    MEDFileFieldMultiTS(const MEDFileFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent);
    MEDFileFieldMultiTS(med_idt fid, bool loadAll, const MEDFileMeshes *ms);
    MEDFileFieldMultiTS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities=0);
  };

  template<class T>
  class MEDFileNDTemplateFieldMultiTS : public MEDFileTemplateFieldMultiTS<T>
  {
  public:
    MEDLOADER_EXPORT MEDFileFieldMultiTS *convertToDouble(bool isDeepCpyGlobs=true) const;
  protected:
    ~MEDFileNDTemplateFieldMultiTS() { }
    MEDFileNDTemplateFieldMultiTS() { }
    MEDFileNDTemplateFieldMultiTS(const typename MLFieldTraits<T>::FMTSWSDAType& other, bool shallowCopyOfContent):MEDFileTemplateFieldMultiTS<T>(other,shallowCopyOfContent) { }
    MEDFileNDTemplateFieldMultiTS(med_idt fid, bool loadAll, const MEDFileMeshes *ms):MEDFileTemplateFieldMultiTS<T>(fid,loadAll,ms) { }
    MEDFileNDTemplateFieldMultiTS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities):MEDFileTemplateFieldMultiTS<T>(fid,fieldName,loadAll,ms,entities) { }
  };

  /*!
   * User class.
   */
  class MEDFileIntFieldMultiTS : public MEDFileNDTemplateFieldMultiTS<int>
  {
    friend class MEDFileTemplateFieldMultiTS<int>;
  public:
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *shallowCpy() const { return new MEDFileIntFieldMultiTS(*this); }
    MEDLOADER_EXPORT MEDFileIntFieldMultiTS *buildNewEmpty() const { return buildNewEmptyImpl(); }
  private:
    ~MEDFileIntFieldMultiTS() { }
    MEDFileIntFieldMultiTS() { }
    MEDFileIntFieldMultiTS(const MEDFileIntFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent):MEDFileNDTemplateFieldMultiTS<int>(other,shallowCopyOfContent) { }
    MEDFileIntFieldMultiTS(med_idt fid, bool loadAll, const MEDFileMeshes *ms):MEDFileNDTemplateFieldMultiTS<int>(fid,loadAll,ms) { }
    MEDFileIntFieldMultiTS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities=0):MEDFileNDTemplateFieldMultiTS<int>(fid,fieldName,loadAll,ms,entities) { }
  };

  /*!
   * User class.
   */
  class MEDFileFloatFieldMultiTS : public MEDFileNDTemplateFieldMultiTS<float>
  {
    friend class MEDFileTemplateFieldMultiTS<float>;
  public:
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTS *shallowCpy() const { return new MEDFileFloatFieldMultiTS(*this); }
    MEDLOADER_EXPORT MEDFileFloatFieldMultiTS *buildNewEmpty() const { return buildNewEmptyImpl(); }
  private:
    ~MEDFileFloatFieldMultiTS() { }
    MEDFileFloatFieldMultiTS() { }
    MEDFileFloatFieldMultiTS(const MEDFileFloatFieldMultiTSWithoutSDA& other, bool shallowCopyOfContent):MEDFileNDTemplateFieldMultiTS<float>(other,shallowCopyOfContent) { }
    MEDFileFloatFieldMultiTS(med_idt fid, bool loadAll, const MEDFileMeshes *ms):MEDFileNDTemplateFieldMultiTS<float>(fid,loadAll,ms) { }
    MEDFileFloatFieldMultiTS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities=0):MEDFileNDTemplateFieldMultiTS<float>(fid,fieldName,loadAll,ms,entities) { }
  };

  class MEDFileAnyTypeFieldMultiTSIterator
  {
  public:
    MEDLOADER_EXPORT MEDFileAnyTypeFieldMultiTSIterator(MEDFileAnyTypeFieldMultiTS *fmts);
    MEDLOADER_EXPORT ~MEDFileAnyTypeFieldMultiTSIterator();
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *nextt();
  private:
    MCAuto<MEDFileAnyTypeFieldMultiTS> _fmts;
    int _iter_id;
    int _nb_iter;
  };
}

#endif
