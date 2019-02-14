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

#ifndef __MEDFILEFIELD1TS_HXX__
#define __MEDFILEFIELD1TS_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDFileFieldGlobs.hxx"
#include "MEDFileFieldInternal.hxx"
#include "MEDLoaderTraits.hxx"

#include "med.h"

namespace MEDCoupling
{
  class TimeHolder;
  class MEDFileMeshes;
  
  /*!
   * SDA is for Shared Data Arrays such as profiles.
   */
  class MEDFileAnyTypeField1TSWithoutSDA : public RefCountObject, public MEDFileFieldNameScope
  {
  public:
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA();
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order);
    MEDLOADER_EXPORT int getIteration() const { return _iteration; }
    MEDLOADER_EXPORT int getOrder() const { return _order; }
    MEDLOADER_EXPORT double getTime(int& iteration, int& order) const { iteration=_iteration; order=_order; return _dt; }
    MEDLOADER_EXPORT void setTime(int iteration, int order, double val) { _dt=val; _iteration=iteration; _order=order; }
    MEDLOADER_EXPORT int getDimension() const;
    MEDLOADER_EXPORT bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT int getMeshIteration() const;
    MEDLOADER_EXPORT int getMeshOrder() const;
    MEDLOADER_EXPORT bool isDealingTS(int iteration, int order) const;
    MEDLOADER_EXPORT std::pair<int,int> getDtIt() const;
    MEDLOADER_EXPORT void fillIteration(std::pair<int,int>& p) const;
    MEDLOADER_EXPORT void fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const;
    MEDLOADER_EXPORT std::vector<TypeOfField> getTypesOfFieldAvailable() const;
    //
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsed2() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsed2() const;
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsedMulti2() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsedMulti2() const;
    MEDLOADER_EXPORT void changePflsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changeLocsRefsNamesGen2(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    //
    MEDLOADER_EXPORT int getNonEmptyLevels(const std::string& mname, std::vector<int>& levs) const;
    MEDLOADER_EXPORT void convertMedBallIntoClassic();
    MEDLOADER_EXPORT void makeReduction(INTERP_KERNEL::NormalizedCellType ct, TypeOfField tof, const DataArrayInt *pfl);
    MEDLOADER_EXPORT std::vector< std::vector<std::pair<int,int> > > getFieldSplitedByType(const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    //
    MEDLOADER_EXPORT MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId);
    MEDLOADER_EXPORT const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const;
    MEDLOADER_EXPORT void deepCpyLeavesFrom(const MEDFileAnyTypeField1TSWithoutSDA& other);
    MEDLOADER_EXPORT void accept(MEDFileFieldVisitor& visitor) const;
  public:
    MEDLOADER_EXPORT int getNumberOfComponents() const;
    MEDLOADER_EXPORT const std::vector<std::string>& getInfo() const;
    MEDLOADER_EXPORT std::vector<std::string>& getInfo();
    MEDLOADER_EXPORT bool presenceOfMultiDiscPerGeoType() const;
    MEDLOADER_EXPORT bool presenceOfStructureElements() const;
    MEDLOADER_EXPORT bool onlyStructureElements() const;
    MEDLOADER_EXPORT void killStructureElements();
    MEDLOADER_EXPORT void keepOnlyStructureElements();
    MEDLOADER_EXPORT void keepOnlyOnSE(const std::string& seName);
    MEDLOADER_EXPORT void getMeshSENames(std::vector< std::pair<std::string,std::string> >& ps) const;
    MEDLOADER_EXPORT void setInfo(const std::vector<std::string>& infos);
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT int copyTinyInfoFrom(const TimeHolder *th, const MEDCouplingFieldTemplate *field, const DataArray *arr);
    MEDLOADER_EXPORT void setFieldNoProfileSBT(const TimeHolder *th, const MEDCouplingFieldTemplate *field, const DataArray *arr, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc);
    MEDLOADER_EXPORT void setFieldProfile(const TimeHolder *th, const MEDCouplingFieldTemplate *field, const DataArray *arrOfVals, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, MEDFileFieldGlobsReal& glob, const MEDFileFieldNameScope& nasc, bool smartPflKiller=true);
    MEDLOADER_EXPORT virtual void simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeField1TSWithoutSDA *deepCopy() const = 0;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeField1TSWithoutSDA *shallowCpy() const = 0;
    MEDLOADER_EXPORT virtual std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > splitComponents() const;
    MEDLOADER_EXPORT virtual const char *getTypeStr() const = 0;
    MEDLOADER_EXPORT virtual DataArray *getUndergroundDataArray() const = 0;
    MEDLOADER_EXPORT virtual DataArray *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const = 0;
    MEDLOADER_EXPORT virtual void setArray(DataArray *arr) = 0;
    MEDLOADER_EXPORT virtual DataArray *createNewEmptyDataArrayInstance() const = 0;
    MEDLOADER_EXPORT virtual DataArray *getOrCreateAndGetArray() = 0;
    MEDLOADER_EXPORT virtual const DataArray *getOrCreateAndGetArray() const = 0;
  public:
    MEDLOADER_EXPORT MEDCouplingFieldDouble *fieldOnMesh(const MEDFileFieldGlobsReal *glob, const MEDFileMesh *mesh, MCAuto<DataArray>& arrOut, const MEDFileFieldNameScope& nasc) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, const std::string& mName, int renumPol, const MEDFileFieldGlobsReal *glob, MCAuto<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDFileMesh *mesh, MCAuto<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldAtTopLevel(TypeOfField type, const std::string& mName, int renumPol, const MEDFileFieldGlobsReal *glob, MCAuto<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    MEDLOADER_EXPORT MEDCouplingFieldDouble *getFieldOnMeshAtLevel(TypeOfField type, int renumPol, const MEDFileFieldGlobsReal *glob, const MEDCouplingMesh *mesh, const DataArrayInt *cellRenum, const DataArrayInt *nodeRenum, MCAuto<DataArray> &arrOut, const MEDFileFieldNameScope& nasc) const;
    DataArray *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl, const MEDFileFieldGlobsReal *glob, const MEDFileFieldNameScope& nasc) const;
  public:
    MEDLOADER_EXPORT bool renumberEntitiesLyingOnMesh(const std::string& meshName, const std::vector<int>& oldCode, const std::vector<int>& newCode, const DataArrayInt *renumO2N, MEDFileFieldGlobsReal& glob);
    MEDLOADER_EXPORT std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > splitDiscretizations() const;
    MEDLOADER_EXPORT std::vector< MCAuto<MEDFileAnyTypeField1TSWithoutSDA> > splitMultiDiscrPerGeoTypes() const;
    MEDLOADER_EXPORT int keepOnlySpatialDiscretization(TypeOfField tof, std::vector< std::pair<int,int> >& its);
    MEDLOADER_EXPORT int keepOnlyGaussDiscretization(std::size_t idOfDisc, std::vector< std::pair<int,int> >& its);
  public:
    MEDLOADER_EXPORT void allocNotFromFile(int newNbOfTuples);
    MEDLOADER_EXPORT bool allocIfNecessaryTheArrayToReceiveDataFromFile();
    MEDLOADER_EXPORT void loadOnlyStructureOfDataRecursively(med_idt fid, const MEDFileFieldNameScope& nasc, const MEDFileMeshes *ms, const MEDFileEntities *entities);
    MEDLOADER_EXPORT void loadBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc);
    MEDLOADER_EXPORT void loadBigArraysRecursivelyIfNecessary(med_idt fid, const MEDFileFieldNameScope& nasc);
    MEDLOADER_EXPORT void loadStructureAndBigArraysRecursively(med_idt fid, const MEDFileFieldNameScope& nasc, const MEDFileMeshes *ms, const MEDFileEntities *entities);
    MEDLOADER_EXPORT void unloadArrays();
    MEDLOADER_EXPORT void writeLL(med_idt fid, const MEDFileWritable& opts, const MEDFileFieldNameScope& nasc) const;
    MEDLOADER_EXPORT static std::string FieldNameToMEDFileConvention(const std::string& nonCorrectFieldName);
  protected:
    int getMeshIdFromMeshName(const std::string& mName) const;
    int addNewEntryIfNecessary(const MEDCouplingMesh *mesh);
    void updateData(int newLgth, const std::vector< std::pair<int,int> >& oldStartStops);
  protected:
    std::vector< MCAuto< MEDFileFieldPerMesh > > _field_per_mesh;
    int _iteration;
    int _order;
    double _dt;
  public:
    //! only useable on reading
    mutable int _csit;
    // -3 means allocated and build from scratch
    // -2 means allocated and read from a file
    // -1 means not allocated and build from scratch
    // >=0 means not allocated and read from a file
    mutable int _nb_of_tuples_to_be_allocated;
  };

  class MEDFileIntField1TSWithoutSDA;

  template<class T>
  class MEDFileField1TSTemplateWithoutSDA : public MEDFileAnyTypeField1TSWithoutSDA
  {
  protected:
    MEDFileField1TSTemplateWithoutSDA(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order):MEDFileAnyTypeField1TSWithoutSDA(fieldName,meshName,csit,iteration,order) { }
    MEDFileField1TSTemplateWithoutSDA() { }
  public:
    MEDLOADER_EXPORT void setArray(DataArray *arr);
    MEDLOADER_EXPORT DataArray *createNewEmptyDataArrayInstance() const;
    MEDLOADER_EXPORT typename Traits<T>::ArrayType *getOrCreateAndGetArrayTemplate();
    MEDLOADER_EXPORT typename Traits<T>::ArrayType const *getOrCreateAndGetArrayTemplate() const;
    MEDLOADER_EXPORT typename Traits<T>::ArrayType *getUndergroundDataArrayTemplate() const;
    MEDLOADER_EXPORT typename Traits<T>::ArrayType *getUndergroundDataArrayTemplateExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDLOADER_EXPORT DataArray *getOrCreateAndGetArray();
    MEDLOADER_EXPORT const DataArray *getOrCreateAndGetArray() const;
    MEDLOADER_EXPORT DataArray *getUndergroundDataArray() const;
    MEDLOADER_EXPORT void aggregate(const typename std::vector< typename MLFieldTraits<T>::F1TSWSDAType const * >& f1tss, const std::vector< std::vector< std::pair<int,int> > >& dts);
    MEDLOADER_EXPORT void copyTimeInfoFrom(const typename Traits<T>::FieldType *mcf);
  protected:
    MCAuto< typename Traits<T>::ArrayType > _arr;
  };

  /*!
   * SDA is for Shared Data Arrays such as profiles.
   */
  class MEDFileField1TSWithoutSDA : public MEDFileField1TSTemplateWithoutSDA<double>
  {
  public:
    MEDLOADER_EXPORT const char *getTypeStr() const;
    MEDLOADER_EXPORT DataArray *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDLOADER_EXPORT std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF, std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    MEDLOADER_EXPORT static void CheckMeshDimRel(int meshDimRelToMax);
    MEDLOADER_EXPORT static std::vector<int> CheckSBTMesh(const MEDCouplingMesh *mesh);
    MEDLOADER_EXPORT static MEDFileField1TSWithoutSDA *New(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order, const std::vector<std::string>& infos);
  public:
    MEDLOADER_EXPORT MEDFileField1TSWithoutSDA();
    MEDLOADER_EXPORT MEDFileField1TSWithoutSDA(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order, const std::vector<std::string>& infos);
    MEDLOADER_EXPORT MEDFileField1TSWithoutSDA *shallowCpy() const;
    MEDLOADER_EXPORT MEDFileField1TSWithoutSDA *deepCopy() const;
    MEDLOADER_EXPORT MEDFileIntField1TSWithoutSDA *convertToInt() const;
  public:
    static const char TYPE_STR[];
  };

  template<class T>
  class MEDFileField1TSNDTemplateWithoutSDA : public MEDFileField1TSTemplateWithoutSDA<T>
  {
  public:
    MEDLOADER_EXPORT MEDFileField1TSWithoutSDA *convertToDouble() const;
  protected:
    MEDFileField1TSNDTemplateWithoutSDA() { }
    MEDFileField1TSNDTemplateWithoutSDA(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order, const std::vector<std::string>& infos):MEDFileField1TSTemplateWithoutSDA<T>(fieldName,meshName,csit,iteration,order) { }
  };
  
  /*!
   * SDA is for Shared Data Arrays such as profiles.
   */
  class MEDFileIntField1TSWithoutSDA : public MEDFileField1TSNDTemplateWithoutSDA<int>
  {
  public:
    MEDLOADER_EXPORT MEDFileIntField1TSWithoutSDA();
    MEDLOADER_EXPORT static MEDFileIntField1TSWithoutSDA *New(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order, const std::vector<std::string>& infos);
    MEDLOADER_EXPORT MEDFileIntField1TSWithoutSDA *deepCopy() const;
    MEDLOADER_EXPORT MEDFileIntField1TSWithoutSDA *shallowCpy() const;
    MEDLOADER_EXPORT const char *getTypeStr() const;
    MEDLOADER_EXPORT DataArray *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDLOADER_EXPORT DataArrayInt *getUndergroundDataArrayIntExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
  protected:
    MEDFileIntField1TSWithoutSDA(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order, const std::vector<std::string>& infos);
  public:
    MEDLOADER_EXPORT static const char TYPE_STR[];
  };

  /*!
   * SDA is for Shared Data Arrays such as profiles.
   */
  class MEDFileFloatField1TSWithoutSDA : public MEDFileField1TSNDTemplateWithoutSDA<float>
  {
  public:
    MEDLOADER_EXPORT MEDFileFloatField1TSWithoutSDA();
    MEDLOADER_EXPORT static MEDFileFloatField1TSWithoutSDA *New(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order, const std::vector<std::string>& infos);
    MEDLOADER_EXPORT MEDFileFloatField1TSWithoutSDA *deepCopy() const;
    MEDLOADER_EXPORT MEDFileFloatField1TSWithoutSDA *shallowCpy() const;
    MEDLOADER_EXPORT const char *getTypeStr() const;
    MEDLOADER_EXPORT DataArray *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDLOADER_EXPORT DataArrayFloat *getUndergroundDataArrayFloatExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
  protected:
    MEDFileFloatField1TSWithoutSDA(const std::string& fieldName, const std::string& meshName, int csit, int iteration, int order, const std::vector<std::string>& infos);
  public:
    MEDLOADER_EXPORT static const char TYPE_STR[];
  };

  /*!
   * User class.
   */
  class MEDFileAnyTypeField1TS : public RefCountObject, public MEDFileWritableStandAlone, public MEDFileFieldGlobsReal
  {
  protected:
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS();
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS(med_idt fid, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities=0);
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities=0);
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities=0);
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS(const MEDFileAnyTypeField1TSWithoutSDA& other, bool shallowCopyOfContent);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *BuildNewInstanceFromContent(MEDFileAnyTypeField1TSWithoutSDA *c);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *BuildNewInstanceFromContent(MEDFileAnyTypeField1TSWithoutSDA *c, med_idt fid);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TSWithoutSDA *BuildContentFrom(med_idt fid, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TSWithoutSDA *BuildContentFrom(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TSWithoutSDA *BuildContentFrom(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms, const MEDFileEntities *entities);
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    // direct forwarding to MEDFileAnyTypeField1TSWithoutSDA instance _content
  public:
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *New(const std::string& fileName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *New(med_idt fid, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *New(med_idt fid, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *New(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll=true);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *NewAdv(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileEntities *entities);
    MEDLOADER_EXPORT static MEDFileAnyTypeField1TS *NewAdv(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileEntities *entities);
    MEDLOADER_EXPORT int getDimension() const;
    MEDLOADER_EXPORT int getIteration() const;
    MEDLOADER_EXPORT int getOrder() const;
    MEDLOADER_EXPORT double getTime(int& iteration, int& order) const;
    MEDLOADER_EXPORT void setTime(int iteration, int order, double val);
    MEDLOADER_EXPORT std::string getName() const;
    MEDLOADER_EXPORT void setName(const std::string& name);
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT void simpleRepr(int bkOffset, std::ostream& oss, int f1tsId) const;
    MEDLOADER_EXPORT std::string getDtUnit() const;
    MEDLOADER_EXPORT void setDtUnit(const std::string& dtUnit);
    MEDLOADER_EXPORT std::string getMeshName() const;
    MEDLOADER_EXPORT void setMeshName(const std::string& newMeshName);
    MEDLOADER_EXPORT bool changeMeshNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT int getMeshIteration() const;
    MEDLOADER_EXPORT int getMeshOrder() const;
    MEDLOADER_EXPORT int getNumberOfComponents() const;
    MEDLOADER_EXPORT bool isDealingTS(int iteration, int order) const;
    MEDLOADER_EXPORT std::pair<int,int> getDtIt() const;
    MEDLOADER_EXPORT void fillIteration(std::pair<int,int>& p) const;
    MEDLOADER_EXPORT void fillTypesOfFieldAvailable(std::vector<TypeOfField>& types) const;
    MEDLOADER_EXPORT void setInfo(const std::vector<std::string>& infos);
    MEDLOADER_EXPORT const std::vector<std::string>& getInfo() const;
    MEDLOADER_EXPORT std::vector<std::string>& getInfo();
    MEDLOADER_EXPORT bool presenceOfMultiDiscPerGeoType() const;
    MEDLOADER_EXPORT std::vector<TypeOfField> getTypesOfFieldAvailable() const;
    MEDLOADER_EXPORT std::vector< std::vector<std::pair<int,int> > > getFieldSplitedByType(const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
        std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
    MEDLOADER_EXPORT MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId);
    MEDLOADER_EXPORT const MEDFileFieldPerMeshPerTypePerDisc *getLeafGivenMeshAndTypeAndLocId(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId) const;
    MEDLOADER_EXPORT int getNonEmptyLevels(const std::string& mname, std::vector<int>& levs) const;
    MEDLOADER_EXPORT void convertMedBallIntoClassic();
    MEDLOADER_EXPORT void makeReduction(INTERP_KERNEL::NormalizedCellType ct, TypeOfField tof, const DataArrayInt *pfl);
    MEDLOADER_EXPORT virtual MEDFileAnyTypeField1TS *extractPart(const std::map<int, MCAuto<DataArrayInt> >& extractDef, MEDFileMesh *mm) const = 0;
    MEDLOADER_EXPORT virtual MEDFileAnyTypeField1TS *shallowCpy() const = 0;
  public:
    MEDLOADER_EXPORT void loadArrays();
    MEDLOADER_EXPORT void loadArraysIfNecessary();
    MEDLOADER_EXPORT void unloadArrays();
    MEDLOADER_EXPORT void unloadArraysWithoutDataLoss();
    MEDLOADER_EXPORT std::vector< MCAuto< MEDFileAnyTypeField1TS > > splitComponents() const;
    MEDLOADER_EXPORT std::vector< MCAuto< MEDFileAnyTypeField1TS > > splitDiscretizations() const;
    MEDLOADER_EXPORT std::vector< MCAuto< MEDFileAnyTypeField1TS > > splitMultiDiscrPerGeoTypes() const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *deepCopy() const;
    MEDLOADER_EXPORT int copyTinyInfoFrom(const MEDCouplingFieldDouble *field, const DataArray *arr);
    MEDLOADER_EXPORT int copyTinyInfoFrom(const TimeHolder *th, const MEDCouplingFieldTemplate *field, const DataArray *arr);
  public:
    //! underground method see MEDFileField1TSWithoutSDA::setProfileNameOnLeaf
    MEDLOADER_EXPORT void setProfileNameOnLeaf(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newPflName, bool forceRenameOnGlob=false);
    //! underground method see MEDFileField1TSWithoutSDA::setLocNameOnLeaf
    MEDLOADER_EXPORT void setLocNameOnLeaf(const std::string& mName, INTERP_KERNEL::NormalizedCellType typ, int locId, const std::string& newLocName, bool forceRenameOnGlob=false);
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsed() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsed() const;
    MEDLOADER_EXPORT std::vector<std::string> getPflsReallyUsedMulti() const;
    MEDLOADER_EXPORT std::vector<std::string> getLocsReallyUsedMulti() const;
    MEDLOADER_EXPORT void changePflsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
    MEDLOADER_EXPORT void changeLocsRefsNamesGen(const std::vector< std::pair<std::vector<std::string>, std::string > >& mapOfModif);
  public:
    MEDLOADER_EXPORT static int LocateField2(med_idt fid, int fieldIdCFormat, bool checkFieldId, std::string& fieldName, med_field_type& typcha, std::vector<std::string>& infos, std::string& dtunitOut, std::string& meshName);
    MEDLOADER_EXPORT static int LocateField(med_idt fid, const std::string& fieldName, int& posCFormat, med_field_type& typcha, std::vector<std::string>& infos, std::string& dtunitOut, std::string& meshName);
  public:
    MEDLOADER_EXPORT virtual med_field_type getMEDFileFieldType() const = 0;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TSWithoutSDA *contentNotNullBase();
    MEDLOADER_EXPORT const MEDFileAnyTypeField1TSWithoutSDA *contentNotNullBase() const;
  protected:
    MCAuto<MEDFileAnyTypeField1TSWithoutSDA> _content;
  };

  class MEDFileIntField1TS;

  template<class T>
  class MEDFileTemplateField1TS : public MEDFileAnyTypeField1TS
  {
  public:
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::F1TSType *New();
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::F1TSType *New(const std::string& fileName, bool loadAll=true);
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::F1TSType *New(med_idt fid, bool loadAll=true);
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::F1TSType *New(DataArrayByte *db) { return BuildFromMemoryChunk<typename MLFieldTraits<T>::F1TSType>(db); }
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::F1TSType *New(const std::string& fileName, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::F1TSType *New(med_idt fid, const std::string& fieldName, bool loadAll=true);
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::F1TSType *New(const std::string& fileName, const std::string& fieldName, int iteration, int order, bool loadAll=true);
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::F1TSType *New(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll=true);
    MEDLOADER_EXPORT static typename MLFieldTraits<T>::F1TSType *New(const typename MLFieldTraits<T>::F1TSWSDAType& other, bool shallowCopyOfContent);
  public:
    MEDLOADER_EXPORT static typename Traits<T>::ArrayType *ReturnSafelyTypedDataArray(MCAuto<DataArray>& arr);
    MEDLOADER_EXPORT typename Traits<T>::ArrayType *getFieldWithProfile(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, DataArrayInt *&pfl) const;
    MEDLOADER_EXPORT void setArray(DataArray *arr);
    MEDLOADER_EXPORT typename Traits<T>::ArrayType *getUndergroundDataArray() const;
    MEDLOADER_EXPORT typename Traits<T>::ArrayType *getUndergroundDataArrayExt(std::vector< std::pair<std::pair<INTERP_KERNEL::NormalizedCellType,int>,std::pair<int,int> > >& entries) const;
    MEDLOADER_EXPORT static MCAuto<typename Traits<T>::FieldType> SetDataArrayInField(MEDCouplingFieldDouble *f, MCAuto<DataArray>& arr);
    MEDLOADER_EXPORT static MCAuto<MEDCouplingFieldDouble> ToFieldTemplateWithTime(const typename Traits<T>::FieldType *f);
  public:
    MEDLOADER_EXPORT void copyTimeInfoFrom(const typename Traits<T>::FieldType *mcf);
    MEDLOADER_EXPORT typename Traits<T>::FieldType *field(const MEDFileMesh *mesh) const;
    MEDLOADER_EXPORT typename Traits<T>::FieldType *getFieldAtLevel(TypeOfField type, int meshDimRelToMax, int renumPol=0) const;
    MEDLOADER_EXPORT typename Traits<T>::FieldType *getFieldAtTopLevel(TypeOfField type, int renumPol=0) const;
    MEDLOADER_EXPORT typename Traits<T>::FieldType *getFieldOnMeshAtLevel(TypeOfField type, int meshDimRelToMax, const MEDFileMesh *mesh, int renumPol=0) const;
    MEDLOADER_EXPORT typename Traits<T>::FieldType *getFieldOnMeshAtLevel(TypeOfField type, const MEDCouplingMesh *mesh, int renumPol=0) const;
    MEDLOADER_EXPORT typename Traits<T>::FieldType *getFieldAtLevelOld(TypeOfField type, const std::string& mname, int meshDimRelToMax, int renumPol=0) const;
    MEDLOADER_EXPORT void setFieldNoProfileSBT(const typename Traits<T>::FieldType *field);
    MEDLOADER_EXPORT void setFieldProfile(const typename Traits<T>::FieldType *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    MEDLOADER_EXPORT void setFieldProfileFlatly(const typename Traits<T>::FieldType *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile);
    MEDLOADER_EXPORT typename MLFieldTraits<T>::F1TSType *extractPartImpl(const std::map<int, MCAuto<DataArrayInt> >& extractDef, MEDFileMesh *mm) const;
    MEDLOADER_EXPORT MEDFileAnyTypeField1TS *extractPart(const std::map<int, MCAuto<DataArrayInt> >& extractDef, MEDFileMesh *mm) const { return this->extractPartImpl(extractDef,mm); }
  protected:
    ~MEDFileTemplateField1TS() { }
    MEDFileTemplateField1TS();
    MEDFileTemplateField1TS(med_idt fid, bool loadAll, const MEDFileMeshes *ms):MEDFileAnyTypeField1TS(fid,loadAll,ms) { }
    MEDFileTemplateField1TS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms):MEDFileAnyTypeField1TS(fid,fieldName,loadAll,ms) { }
    MEDFileTemplateField1TS(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms):MEDFileAnyTypeField1TS(fid,fieldName,iteration,order,loadAll,ms) { }
    MEDFileTemplateField1TS(const typename MLFieldTraits<T>::F1TSWSDAType& other, bool shallowCopyOfContent):MEDFileAnyTypeField1TS(other,shallowCopyOfContent) { }
    const typename MLFieldTraits<T>::F1TSWSDAType *contentNotNull() const;
    typename MLFieldTraits<T>::F1TSWSDAType *contentNotNull();
    void setFieldProfileGeneral(const typename Traits<T>::FieldType *field, const MEDFileMesh *mesh, int meshDimRelToMax, const DataArrayInt *profile, bool smartPflKiller);
  };

  /*!
   * User class.
   */
  class MEDFileField1TS : public MEDFileTemplateField1TS<double>
  {
    friend class MEDFileTemplateField1TS<double>;
  public:
    MEDLOADER_EXPORT MEDFileIntField1TS *convertToInt(bool isDeepCpyGlobs=true) const;
  public:
    MEDLOADER_EXPORT MEDFileField1TS *shallowCpy() const;
    MEDLOADER_EXPORT std::vector< std::vector<DataArrayDouble *> > getFieldSplitedByType2(const std::string& mname, std::vector<INTERP_KERNEL::NormalizedCellType>& types, std::vector< std::vector<TypeOfField> >& typesF,
                                                                                          std::vector< std::vector<std::string> >& pfls, std::vector< std::vector<std::string> >& locs) const;
  public:
  private:
    med_field_type getMEDFileFieldType() const { return MED_FLOAT64; }
  private:
    ~MEDFileField1TS() { }
    MEDFileField1TS(med_idt fid, bool loadAll, const MEDFileMeshes *ms);
    MEDFileField1TS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms);
    MEDFileField1TS(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms);
    MEDFileField1TS(const MEDFileField1TSWithoutSDA& other, bool shallowCopyOfContent);
    MEDFileField1TS() { }
  };

  template<class T>
  class MEDFileNDTemplateField1TS : public MEDFileTemplateField1TS<T>
  {
  public:
    MEDLOADER_EXPORT MEDFileField1TS *convertToDouble(bool isDeepCpyGlobs=true) const;
  protected:
    ~MEDFileNDTemplateField1TS() { }
    MEDFileNDTemplateField1TS() { }
    MEDFileNDTemplateField1TS(med_idt fid, bool loadAll, const MEDFileMeshes *ms):MEDFileTemplateField1TS<T>(fid,loadAll,ms) { }
    MEDFileNDTemplateField1TS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms):MEDFileTemplateField1TS<T>(fid,fieldName,loadAll,ms) { }
    MEDFileNDTemplateField1TS(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms):MEDFileTemplateField1TS<T>(fid,fieldName,iteration,order,loadAll,ms) { }
    MEDFileNDTemplateField1TS(const typename MLFieldTraits<T>::F1TSWSDAType& other, bool shallowCopyOfContent):MEDFileTemplateField1TS<T>(other,shallowCopyOfContent) { }
  };

  class MEDFileIntField1TS : public MEDFileNDTemplateField1TS<int>
  {
    friend class MEDFileTemplateField1TS<int>;
  public:
    MEDLOADER_EXPORT MEDFileIntField1TS *shallowCpy() const { return new MEDFileIntField1TS(*this); }
  public:
    MEDLOADER_EXPORT static MCAuto<MEDCouplingFieldDouble> ConvertFieldIntToFieldDouble(const MEDCouplingFieldInt *f);
  private:
    med_field_type getMEDFileFieldType() const { return MED_INT32; }
  private:
    ~MEDFileIntField1TS() { }
    MEDFileIntField1TS() { }
    MEDFileIntField1TS(med_idt fid, bool loadAll, const MEDFileMeshes *ms):MEDFileNDTemplateField1TS<int>(fid,loadAll,ms) { }
    MEDFileIntField1TS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms):MEDFileNDTemplateField1TS<int>(fid,fieldName,loadAll,ms) { }
    MEDFileIntField1TS(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms):MEDFileNDTemplateField1TS<int>(fid,fieldName,iteration,order,loadAll,ms) { }
    /*!
     * This constructor is a shallow copy constructor. If \a shallowCopyOfContent is true the content of \a other is shallow copied.
     * If \a shallowCopyOfContent is false, \a other is taken to be the content of \a this.
     *
     * \warning this is a shallow copy constructor
     */
    MEDFileIntField1TS(const MEDFileIntField1TSWithoutSDA& other, bool shallowCopyOfContent):MEDFileNDTemplateField1TS<int>(other,shallowCopyOfContent) { }
  };

  class MEDFileFloatField1TS : public MEDFileNDTemplateField1TS<float>
  {
    friend class MEDFileTemplateField1TS<float>;
  private:
    med_field_type getMEDFileFieldType() const { return MED_FLOAT32; }
    MEDLOADER_EXPORT MEDFileFloatField1TS *shallowCpy() const { return new MEDFileFloatField1TS(*this); }
  private:
    ~MEDFileFloatField1TS() { }
    MEDFileFloatField1TS() { }
    MEDFileFloatField1TS(med_idt fid, bool loadAll, const MEDFileMeshes *ms):MEDFileNDTemplateField1TS<float>(fid,loadAll,ms) { }
    MEDFileFloatField1TS(med_idt fid, const std::string& fieldName, bool loadAll, const MEDFileMeshes *ms):MEDFileNDTemplateField1TS<float>(fid,fieldName,loadAll,ms) { }
    MEDFileFloatField1TS(med_idt fid, const std::string& fieldName, int iteration, int order, bool loadAll, const MEDFileMeshes *ms):MEDFileNDTemplateField1TS<float>(fid,fieldName,iteration,order,loadAll,ms) { }
    /*!
     * This constructor is a shallow copy constructor. If \a shallowCopyOfContent is true the content of \a other is shallow copied.
     * If \a shallowCopyOfContent is false, \a other is taken to be the content of \a this.
     *
     * \warning this is a shallow copy constructor
     */
    MEDFileFloatField1TS(const MEDFileFloatField1TSWithoutSDA& other, bool shallowCopyOfContent):MEDFileNDTemplateField1TS<float>(other,shallowCopyOfContent) { }
  };
}

#endif
