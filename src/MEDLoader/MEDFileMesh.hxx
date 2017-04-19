// Copyright (C) 2007-2016  CEA/DEN, EDF R&D
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

#ifndef __MEDFILEMESH_HXX__
#define __MEDFILEMESH_HXX__

#include "MEDLoaderDefines.hxx"
#include "MEDFileMeshLL.hxx"
#include "MEDFileUtilities.txx"
#include "MEDCouplingPartDefinition.hxx"
#include "MEDFileMeshReadSelector.hxx"
#include "MEDFileJoint.hxx"
#include "MEDFileEquivalence.hxx"

#include <map>
#include <list>

namespace MEDCoupling
{
  class MEDFileFieldGlobsReal;
  class MEDFileField1TSStructItem;

  class MEDFileMesh : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT static MEDFileMesh *New(const std::string& fileName, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileMesh *New(med_idt fid, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileMesh *New(DataArrayByte *db) { return BuildFromMemoryChunk<MEDFileMesh>(db); }
    MEDLOADER_EXPORT static MEDFileMesh *New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0, MEDFileJoints* joints=0);
    MEDLOADER_EXPORT static MEDFileMesh *New(med_idt fid, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0, MEDFileJoints* joints=0);
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT virtual MEDFileMesh *createNewEmpty() const = 0;
    MEDLOADER_EXPORT virtual MEDFileMesh *deepCopy() const = 0;
    MEDLOADER_EXPORT virtual MEDFileMesh *shallowCpy() const = 0;
    MEDLOADER_EXPORT virtual bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    MEDLOADER_EXPORT virtual void clearNonDiscrAttributes() const;
    MEDLOADER_EXPORT virtual void setName(const std::string& name);
    MEDLOADER_EXPORT bool changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT std::string getName() const { return _name; }
    MEDLOADER_EXPORT std::string getUnivName() const { return _univ_name; }
    MEDLOADER_EXPORT bool getUnivNameWrStatus() const { return _univ_wr_status; }
    MEDLOADER_EXPORT void setUnivNameWrStatus(bool newStatus) { _univ_wr_status=newStatus; }
    MEDLOADER_EXPORT void setDescription(const std::string& name) { _desc_name=name; }
    MEDLOADER_EXPORT std::string getDescription() const { return _desc_name; }
    MEDLOADER_EXPORT void setOrder(int order) { _order=order; }
    MEDLOADER_EXPORT int getOrder() const { return _order; }
    MEDLOADER_EXPORT void setIteration(int it) { _iteration=it; }
    MEDLOADER_EXPORT int getIteration() const { return _iteration; }
    MEDLOADER_EXPORT void setTimeValue(double time) { _time=time; }
    MEDLOADER_EXPORT void setTime(int dt, int it, double time) { _time=time; _iteration=dt; _order=it; }
    MEDLOADER_EXPORT double getTime(int& dt, int& it) const { dt=_iteration; it=_order; return _time; }
    MEDLOADER_EXPORT double getTimeValue() const { return _time; }
    MEDLOADER_EXPORT void setTimeUnit(const std::string& unit) { _dt_unit=unit; }
    MEDLOADER_EXPORT std::string getTimeUnit() const { return _dt_unit; }
    MEDLOADER_EXPORT void setAxisType(MEDCouplingAxisType at) { _axis_type=at; }
    MEDLOADER_EXPORT MEDCouplingAxisType getAxisType() const { return _axis_type; }
    MEDLOADER_EXPORT std::vector<INTERP_KERNEL::NormalizedCellType> getAllGeoTypes() const;
    MEDLOADER_EXPORT virtual int getNumberOfNodes() const = 0;
    MEDLOADER_EXPORT virtual int getNumberOfCellsAtLevel(int meshDimRelToMaxExt) const = 0;
    MEDLOADER_EXPORT virtual bool hasImplicitPart() const = 0;
    MEDLOADER_EXPORT virtual int buildImplicitPartIfAny(INTERP_KERNEL::NormalizedCellType gt) const = 0;
    MEDLOADER_EXPORT virtual void releaseImplicitPartIfAny() const = 0;
    MEDLOADER_EXPORT virtual std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypesAtLevel(int meshDimRelToMax) const = 0;
    MEDLOADER_EXPORT virtual int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType ct) const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getNonEmptyLevels() const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getNonEmptyLevelsExt() const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getFamArrNonEmptyLevelsExt() const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getNumArrNonEmptyLevelsExt() const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getNameArrNonEmptyLevelsExt() const = 0;
    MEDLOADER_EXPORT virtual int getSizeAtLevel(int meshDimRelToMaxExt) const = 0;
    MEDLOADER_EXPORT virtual MEDCouplingMesh *getMeshAtLevel(int meshDimRelToMax, bool renum=false) const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getDistributionOfTypes(int meshDimRelToMax) const;
    MEDLOADER_EXPORT virtual void whichAreNodesFetched(const MEDFileField1TSStructItem& st, const MEDFileFieldGlobsReal *globs, std::vector<bool>& nodesFetched) const = 0;
    MEDLOADER_EXPORT virtual MEDFileMesh *cartesianize() const = 0;
    MEDLOADER_EXPORT virtual bool presenceOfStructureElements() const = 0;
    MEDLOADER_EXPORT virtual void killStructureElements() { }
    //
    MEDLOADER_EXPORT bool areFamsEqual(const MEDFileMesh *other, std::string& what) const;
    MEDLOADER_EXPORT bool areGrpsEqual(const MEDFileMesh *other, std::string& what) const;
    MEDLOADER_EXPORT bool existsGroup(const std::string& groupName) const;
    MEDLOADER_EXPORT bool existsFamily(int famId) const;
    MEDLOADER_EXPORT bool existsFamily(const std::string& familyName) const;
    MEDLOADER_EXPORT void setFamilyId(const std::string& familyName, int id);
    MEDLOADER_EXPORT void setFamilyIdUnique(const std::string& familyName, int id);
    MEDLOADER_EXPORT virtual void addFamily(const std::string& familyName, int id);
    MEDLOADER_EXPORT virtual void createGroupOnAll(int meshDimRelToMaxExt, const std::string& groupName);
    MEDLOADER_EXPORT virtual bool keepFamIdsOnlyOnLevs(const std::vector<int>& famIds, const std::vector<int>& levs);
    MEDLOADER_EXPORT void addFamilyOnGrp(const std::string& grpName, const std::string& famName);
    MEDLOADER_EXPORT std::string findOrCreateAndGiveFamilyWithId(int id, bool& created);
    MEDLOADER_EXPORT void setFamilyInfo(const std::map<std::string,int>& info);
    MEDLOADER_EXPORT void setGroupInfo(const std::map<std::string, std::vector<std::string> >&info);
    MEDLOADER_EXPORT void copyFamGrpMapsFrom(const MEDFileMesh& other);
    MEDLOADER_EXPORT void clearGrpMap();
    MEDLOADER_EXPORT void clearFamMap();
    MEDLOADER_EXPORT void clearFamGrpMaps();
    MEDLOADER_EXPORT const std::map<std::string,int>& getFamilyInfo() const { return _families; }
    MEDLOADER_EXPORT const std::map<std::string, std::vector<std::string> >& getGroupInfo() const { return _groups; }
    MEDLOADER_EXPORT std::vector<std::string> getFamiliesOnGroup(const std::string& name) const;
    MEDLOADER_EXPORT std::vector<std::string> getFamiliesOnGroups(const std::vector<std::string>& grps) const;
    MEDLOADER_EXPORT std::vector<int> getFamiliesIdsOnGroup(const std::string& name) const;
    MEDLOADER_EXPORT void setFamiliesOnGroup(const std::string& name, const std::vector<std::string>& fams);
    MEDLOADER_EXPORT void setFamiliesIdsOnGroup(const std::string& name, const std::vector<int>& famIds);
    MEDLOADER_EXPORT std::vector<std::string> getGroupsOnFamily(const std::string& name) const;
    MEDLOADER_EXPORT void setGroupsOnFamily(const std::string& famName, const std::vector<std::string>& grps);
    MEDLOADER_EXPORT std::vector<std::string> getGroupsNames() const;
    MEDLOADER_EXPORT std::vector<std::string> getFamiliesNames() const;
    MEDLOADER_EXPORT std::vector<std::string> getGroupsOnSpecifiedLev(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT std::vector<int> getGrpNonEmptyLevelsExt(const std::string& grp) const;
    MEDLOADER_EXPORT std::vector<int> getGrpNonEmptyLevels(const std::string& grp) const;
    MEDLOADER_EXPORT std::vector<int> getGrpsNonEmptyLevels(const std::vector<std::string>& grps) const;
    MEDLOADER_EXPORT std::vector<int> getGrpsNonEmptyLevelsExt(const std::vector<std::string>& grps) const;
    MEDLOADER_EXPORT virtual std::vector<int> getFamsNonEmptyLevels(const std::vector<std::string>& fams) const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getFamsNonEmptyLevelsExt(const std::vector<std::string>& fams) const = 0;
    MEDLOADER_EXPORT std::vector<int> getFamNonEmptyLevels(const std::string& fam) const;
    MEDLOADER_EXPORT std::vector<int> getFamNonEmptyLevelsExt(const std::string& fam) const;
    MEDLOADER_EXPORT std::vector<std::string> getFamiliesNamesWithFilePointOfView() const;
    MEDLOADER_EXPORT static std::string GetMagicFamilyStr();
    MEDLOADER_EXPORT void assignFamilyNameWithGroupName();
    MEDLOADER_EXPORT std::vector<std::string> removeEmptyGroups();
    MEDLOADER_EXPORT void removeGroup(const std::string& name);
    MEDLOADER_EXPORT void removeFamily(const std::string& name);
    MEDLOADER_EXPORT std::vector<std::string> removeOrphanGroups();
    MEDLOADER_EXPORT std::vector<std::string> removeOrphanFamilies();
    MEDLOADER_EXPORT void removeFamiliesReferedByNoGroups();
    MEDLOADER_EXPORT void rearrangeFamilies();
    MEDLOADER_EXPORT void checkOrphanFamilyZero() const;
    MEDLOADER_EXPORT void changeGroupName(const std::string& oldName, const std::string& newName);
    MEDLOADER_EXPORT void changeFamilyName(const std::string& oldName, const std::string& newName);
    MEDLOADER_EXPORT void changeFamilyId(int oldId, int newId);
    MEDLOADER_EXPORT void changeAllGroupsContainingFamily(const std::string& familyNameToChange, const std::vector<std::string>& newFamiliesNames);
    MEDLOADER_EXPORT int getFamilyId(const std::string& name) const;
    MEDLOADER_EXPORT int getMaxAbsFamilyId() const;
    MEDLOADER_EXPORT int getMaxFamilyId() const;
    MEDLOADER_EXPORT int getMinFamilyId() const;
    MEDLOADER_EXPORT int getTheMaxAbsFamilyId() const;
    MEDLOADER_EXPORT int getTheMaxFamilyId() const;
    MEDLOADER_EXPORT int getTheMinFamilyId() const;
    MEDLOADER_EXPORT virtual int getMaxAbsFamilyIdInArrays() const = 0;
    MEDLOADER_EXPORT virtual int getMaxFamilyIdInArrays() const = 0;
    MEDLOADER_EXPORT virtual int getMinFamilyIdInArrays() const = 0;
    MEDLOADER_EXPORT DataArrayInt *getAllFamiliesIdsReferenced() const;
    MEDLOADER_EXPORT DataArrayInt *computeAllFamilyIdsInUse() const;
    MEDLOADER_EXPORT std::vector<int> getFamiliesIds(const std::vector<std::string>& famNames) const;
    MEDLOADER_EXPORT std::string getFamilyNameGivenId(int id) const;
    MEDLOADER_EXPORT bool ensureDifferentFamIdsPerLevel();
    MEDLOADER_EXPORT void normalizeFamIdsTrio();
    MEDLOADER_EXPORT void normalizeFamIdsMEDFile();
    MEDLOADER_EXPORT virtual int getMeshDimension() const = 0;
    MEDLOADER_EXPORT virtual std::string simpleRepr() const;
    MEDLOADER_EXPORT virtual std::string advancedRepr() const = 0;
    //
    MEDLOADER_EXPORT virtual void setGroupsAtLevel(int meshDimRelToMaxExt, const std::vector<const DataArrayInt *>& grps, bool renum=false);
    MEDLOADER_EXPORT virtual void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr) = 0;
    MEDLOADER_EXPORT virtual void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr) = 0;
    MEDLOADER_EXPORT virtual void setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr) = 0;
    MEDLOADER_EXPORT virtual void setGlobalNumFieldAtLevel(int meshDimRelToMaxExt, DataArrayInt *globalNumArr) = 0;
    MEDLOADER_EXPORT virtual void addNodeGroup(const DataArrayInt *ids) = 0;
    MEDLOADER_EXPORT virtual void addGroup(int meshDimRelToMaxExt, const DataArrayInt *ids) = 0;
    MEDLOADER_EXPORT virtual const DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const = 0;
    MEDLOADER_EXPORT virtual DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt) = 0;
    MEDLOADER_EXPORT DataArrayInt *getOrCreateAndGetFamilyFieldAtLevel(int meshDimRelToMaxExt);
    MEDLOADER_EXPORT virtual const DataArrayInt *getNumberFieldAtLevel(int meshDimRelToMaxExt) const = 0;
    MEDLOADER_EXPORT virtual const DataArrayInt *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const = 0;
    MEDLOADER_EXPORT virtual const DataArrayAsciiChar *getNameFieldAtLevel(int meshDimRelToMaxExt) const = 0;
    MEDLOADER_EXPORT virtual MCAuto<DataArrayInt> getGlobalNumFieldAtLevel(int meshDimRelToMaxExt) const = 0;
    MEDLOADER_EXPORT virtual DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const = 0;
    MEDLOADER_EXPORT virtual DataArrayInt *getGroupsArr(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const;
    MEDLOADER_EXPORT virtual DataArrayInt *getGroupArr(int meshDimRelToMaxExt, const std::string& grp, bool renum=false) const;
    MEDLOADER_EXPORT virtual DataArrayInt *getFamilyArr(int meshDimRelToMaxExt, const std::string& fam, bool renum=false) const;
    MEDLOADER_EXPORT virtual DataArrayInt *getNodeGroupArr(const std::string& grp, bool renum=false) const;
    MEDLOADER_EXPORT virtual DataArrayInt *getNodeGroupsArr(const std::vector<std::string>& grps, bool renum=false) const;
    MEDLOADER_EXPORT virtual DataArrayInt *getNodeFamilyArr(const std::string& fam, bool renum=false) const;
    MEDLOADER_EXPORT virtual DataArrayInt *getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum=false) const;
    // tools
    MEDLOADER_EXPORT virtual bool unPolyze(std::vector<int>& oldCode, std::vector<int>& newCode, DataArrayInt *& o2nRenumCell) = 0;
    MEDLOADER_EXPORT int getNumberOfJoints() const;
    MEDLOADER_EXPORT MEDFileJoints *getJoints() const;
    MEDLOADER_EXPORT void setJoints( MEDFileJoints* joints );
    MEDFileEquivalences *getEquivalences() { return _equiv; }
    const MEDFileEquivalences *getEquivalences() const { return _equiv; }
    void killEquivalences() { _equiv=(MEDFileEquivalences *)0; }
    void initializeEquivalences() { _equiv=MEDFileEquivalences::New(this); }
    static INTERP_KERNEL::NormalizedCellType ConvertFromMEDFileGeoType(med_geometry_type geoType);
    static TypeOfField ConvertFromMEDFileEntity(med_entity_type etype);
  protected:
    MEDFileMesh();
    //! protected because no way in MED file API to specify this name
    void setUnivName(const std::string& name) { _univ_name=name; }
    void addFamilyOnAllGroupsHaving(const std::string& famName, const std::string& otherFamName);
    void dealWithTinyInfo(const MEDCouplingMesh *m);
    virtual void synchronizeTinyInfoOnLeaves() const = 0;
    void getFamilyRepr(std::ostream& oss) const;
    virtual void appendFamilyEntries(const DataArrayInt *famIds, const std::vector< std::vector<int> >& fidsOfGrps, const std::vector<std::string>& grpNames);
    virtual void changeFamilyIdArr(int oldId, int newId) = 0;
    virtual std::list< MCAuto<DataArrayInt> > getAllNonNullFamilyIds() const = 0;
    virtual void loadLL(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs) = 0;
    void loadLLWithAdditionalItems(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void addGroupUnderground(bool isNodeGroup, const DataArrayInt *ids, DataArrayInt *famArr);
    static void TranslateFamilyIds(int offset, DataArrayInt *famArr, std::vector< std::vector<int> >& famIdsPerGrp);
    static void ChangeAllGroupsContainingFamily(std::map<std::string, std::vector<std::string> >& groups, const std::string& familyNameToChange, const std::vector<std::string>& newFamiliesNames);
    static std::string FindOrCreateAndGiveFamilyWithId(std::map<std::string,int>& families, int id, bool& created);
    static std::string CreateNameNotIn(const std::string& nameTry, const std::vector<std::string>& namesToAvoid);
    static int PutInThirdComponentOfCodeOffset(std::vector<int>& code, int strt);
    void writeJoints(med_idt fid) const;
    void loadJointsFromFile(med_idt fid, MEDFileJoints *toUseInstedOfReading=0);
    void loadEquivalences(med_idt fid);
    void deepCpyEquivalences(const MEDFileMesh& other);
    bool areEquivalencesEqual(const MEDFileMesh *other, std::string& what) const;
    void getEquivalencesRepr(std::ostream& oss) const;
    void checkCartesian() const;
  private:
    virtual void writeMeshLL(med_idt fid) const = 0;
  protected:
    int _order;
    int _iteration;
    double _time;
    std::string _dt_unit;
    std::string _name;
    //! this attribute do not impact the state of instance -> mutable
    mutable std::string _univ_name;
    bool _univ_wr_status;
    std::string _desc_name;
    MEDCouplingAxisType _axis_type;
    MCAuto<MEDFileJoints> _joints;
    MCAuto<MEDFileEquivalences> _equiv;
  protected:
    std::map<std::string, std::vector<std::string> > _groups;
    std::map<std::string,int> _families;
  public:
    MEDLOADER_EXPORT static const char DFT_FAM_NAME[];
  };

  class MEDCouplingMappedExtrudedMesh;
  
  class MEDFileUMesh : public MEDFileMesh
  {
    friend class MEDFileMesh;
  public:
    MEDLOADER_EXPORT static MEDFileUMesh *New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileUMesh *New(med_idt fid, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileUMesh *New(const std::string& fileName, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileUMesh *New(med_idt fid, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileUMesh *New(DataArrayByte *db) { return BuildFromMemoryChunk<MEDFileUMesh>(db); }
    MEDLOADER_EXPORT static MEDFileUMesh *New(const MEDCouplingMappedExtrudedMesh *mem);
    MEDLOADER_EXPORT static MEDFileUMesh *New();
    MEDLOADER_EXPORT static MEDFileUMesh *LoadPartOf(const std::string& fileName, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<int>& slicPerTyp, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileUMesh *LoadPartOf(med_idt fid, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<int>& slicPerTyp, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static const char *GetSpeStr4ExtMesh() { return SPE_FAM_STR_EXTRUDED_MESH; }
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT MEDFileMesh *createNewEmpty() const;
    MEDLOADER_EXPORT MEDFileUMesh *deepCopy() const;
    MEDLOADER_EXPORT MEDFileUMesh *shallowCpy() const;
    MEDLOADER_EXPORT bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    MEDLOADER_EXPORT void checkConsistency() const;
    MEDLOADER_EXPORT void checkSMESHConsistency() const;
    MEDLOADER_EXPORT void clearNodeAndCellNumbers();
    MEDLOADER_EXPORT void clearNonDiscrAttributes() const;
    MEDLOADER_EXPORT void setName(const std::string& name);
    MEDLOADER_EXPORT const std::vector< MCAuto<MEDFileEltStruct4Mesh> >& getAccessOfUndergroundEltStrs() const { return _elt_str; }
    //
    MEDLOADER_EXPORT int getMaxAbsFamilyIdInArrays() const;
    MEDLOADER_EXPORT int getMaxFamilyIdInArrays() const;
    MEDLOADER_EXPORT int getMinFamilyIdInArrays() const;
    MEDLOADER_EXPORT int getMeshDimension() const;
    MEDLOADER_EXPORT int getSpaceDimension() const;
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT std::string advancedRepr() const;
    MEDLOADER_EXPORT int getSizeAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT const DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt);
    MEDLOADER_EXPORT const DataArrayInt *getNumberFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT const DataArrayInt *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT const DataArrayAsciiChar *getNameFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT MCAuto<DataArrayInt> getGlobalNumFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT const PartDefinition *getPartDefAtLevel(int meshDimRelToMaxExt, INTERP_KERNEL::NormalizedCellType gt=INTERP_KERNEL::NORM_ERROR) const;
    MEDLOADER_EXPORT int getNumberOfNodes() const;
    MEDLOADER_EXPORT int getNumberOfCellsAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT bool hasImplicitPart() const;
    MEDLOADER_EXPORT int buildImplicitPartIfAny(INTERP_KERNEL::NormalizedCellType gt) const;
    MEDLOADER_EXPORT void releaseImplicitPartIfAny() const;
    MEDLOADER_EXPORT std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypesAtLevel(int meshDimRelToMax) const;
    MEDLOADER_EXPORT int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType ct) const;
    MEDLOADER_EXPORT void whichAreNodesFetched(const MEDFileField1TSStructItem& st, const MEDFileFieldGlobsReal *globs, std::vector<bool>& nodesFetched) const;
    MEDLOADER_EXPORT MEDFileMesh *cartesianize() const;
    MEDLOADER_EXPORT bool presenceOfStructureElements() const;
    MEDLOADER_EXPORT void killStructureElements();
    MEDLOADER_EXPORT std::vector<int> getNonEmptyLevels() const;
    MEDLOADER_EXPORT std::vector<int> getNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getFamArrNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getNumArrNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getNameArrNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getFamsNonEmptyLevels(const std::vector<std::string>& fams) const;
    MEDLOADER_EXPORT std::vector<int> getFamsNonEmptyLevelsExt(const std::vector<std::string>& fams) const;
    MEDLOADER_EXPORT DataArrayDouble *getCoords() const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getGroup(int meshDimRelToMaxExt, const std::string& grp, bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getGroups(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getFamily(int meshDimRelToMaxExt, const std::string& fam, bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getFamilies(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const;
    MEDLOADER_EXPORT DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getMeshAtLevel(int meshDimRelToMax, bool renum=false) const;
    MEDLOADER_EXPORT std::vector<int> getDistributionOfTypes(int meshDimRelToMax) const;
    MEDLOADER_EXPORT std::vector< std::pair<int,int> > getAllDistributionOfTypes() const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getLevel0Mesh(bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getLevelM1Mesh(bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getLevelM2Mesh(bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getLevelM3Mesh(bool renum=false) const;
    MEDLOADER_EXPORT void forceComputationOfParts() const;
    MEDLOADER_EXPORT std::vector<MEDCoupling1GTUMesh *> getDirectUndergroundSingleGeoTypeMeshes(int meshDimRelToMax) const;
    MEDLOADER_EXPORT MEDCoupling1GTUMesh *getDirectUndergroundSingleGeoTypeMesh(INTERP_KERNEL::NormalizedCellType gt) const;
    MEDLOADER_EXPORT DataArrayInt *extractFamilyFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    MEDLOADER_EXPORT DataArrayInt *extractNumberFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    MEDLOADER_EXPORT int getRelativeLevOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    //
    MEDLOADER_EXPORT void setFamilyNameAttachedOnId(int id, const std::string& newFamName);
    MEDLOADER_EXPORT void setCoords(DataArrayDouble *coords);
    MEDLOADER_EXPORT void setCoordsForced(DataArrayDouble *coords);
    MEDLOADER_EXPORT void eraseGroupsAtLevel(int meshDimRelToMaxExt);
    MEDLOADER_EXPORT void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr);
    MEDLOADER_EXPORT void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr);
    MEDLOADER_EXPORT void setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr);
    MEDLOADER_EXPORT void setGlobalNumFieldAtLevel(int meshDimRelToMaxExt, DataArrayInt *globalNumArr);
    MEDLOADER_EXPORT void addNodeGroup(const DataArrayInt *ids);
    MEDLOADER_EXPORT void addGroup(int meshDimRelToMaxExt, const DataArrayInt *ids);
    MEDLOADER_EXPORT void removeMeshAtLevel(int meshDimRelToMax);
    MEDLOADER_EXPORT void setMeshAtLevel(int meshDimRelToMax, MEDCoupling1GTUMesh *m);
    MEDLOADER_EXPORT void setMeshAtLevel(int meshDimRelToMax, MEDCouplingUMesh *m, bool newOrOld=false);
    MEDLOADER_EXPORT void setMeshes(const std::vector<const MEDCouplingUMesh *>& ms, bool renum=false);
    MEDLOADER_EXPORT void setGroupsFromScratch(int meshDimRelToMax, const std::vector<const MEDCouplingUMesh *>& ms, bool renum=false);
    MEDLOADER_EXPORT void setGroupsOnSetMesh(int meshDimRelToMax, const std::vector<const MEDCouplingUMesh *>& ms, bool renum=false);
    MEDLOADER_EXPORT void optimizeFamilies();
    // tools
    MEDLOADER_EXPORT void buildInnerBoundaryAlongM1Group(const std::string& grpNameM1, DataArrayInt *&nodesDuplicated, DataArrayInt *&cellsModified, DataArrayInt *&cellsNotModified);
    MEDLOADER_EXPORT bool unPolyze(std::vector<int>& oldCode, std::vector<int>& newCode, DataArrayInt *& o2nRenumCell);
    MEDLOADER_EXPORT DataArrayInt *zipCoords();
    MEDLOADER_EXPORT DataArrayInt *deduceNodeSubPartFromCellSubPart(const std::map<int, MCAuto<DataArrayInt> >& extractDef) const;
    MEDLOADER_EXPORT MEDFileUMesh *extractPart(const std::map<int, MCAuto<DataArrayInt> >& extractDef) const;
    MEDLOADER_EXPORT MEDFileUMesh *buildExtrudedMesh(const MEDCouplingUMesh *m1D, int policy) const;
    MEDLOADER_EXPORT MEDFileUMesh *linearToQuadratic(int conversionType=0, double eps=1e-12) const;
    MEDLOADER_EXPORT MEDFileUMesh *quadraticToLinear(double eps=1e-12) const;
    MEDLOADER_EXPORT MCAuto<MEDFileUMesh> symmetry3DPlane(const double point[3], const double normalVector[3]) const;
    MEDLOADER_EXPORT static MCAuto<MEDFileUMesh> Aggregate(const std::vector<const MEDFileUMesh *>& meshes);
    MEDLOADER_EXPORT MEDCouplingMappedExtrudedMesh *convertToExtrudedMesh() const;
    // serialization
    MEDLOADER_EXPORT void serialize(std::vector<double>& tinyDouble, std::vector<int>& tinyInt, std::vector<std::string>& tinyStr,
                                    std::vector< MCAuto<DataArrayInt> >& bigArraysI, MCAuto<DataArrayDouble>& bigArrayD);
    MEDLOADER_EXPORT void unserialize(std::vector<double>& tinyDouble, std::vector<int>& tinyInt, std::vector<std::string>& tinyStr,
                                      std::vector< MCAuto<DataArrayInt> >& bigArraysI, MCAuto<DataArrayDouble>& bigArrayD);
  private:
    MEDLOADER_EXPORT ~MEDFileUMesh();
    void writeMeshLL(med_idt fid) const;
    MEDFileUMesh();
    MEDFileUMesh(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadPartUMeshFromFile(med_idt fid, const std::string& mName, const std::vector<INTERP_KERNEL::NormalizedCellType>& types, const std::vector<int>& slicPerTyp, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    void loadLL(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void dispatchLoadedPart(med_idt fid, const MEDFileUMeshL2& loaderl2, const std::string& mName, MEDFileMeshReadSelector *mrs);
    const MEDFileUMeshSplitL1 *getMeshAtLevSafe(int meshDimRelToMaxExt) const;
    MEDFileUMeshSplitL1 *getMeshAtLevSafe(int meshDimRelToMaxExt);
    void checkMeshDimCoherency(int meshDim, int meshDimRelToMax) const;
    DataArrayDouble *checkMultiMesh(const std::vector<const MEDCouplingUMesh *>& ms) const;
    void computeRevNum() const;
    void synchronizeTinyInfoOnLeaves() const;
    void changeFamilyIdArr(int oldId, int newId);
    std::list< MCAuto<DataArrayInt> > getAllNonNullFamilyIds() const;
    MCAuto<MEDFileUMeshSplitL1>& checkAndGiveEntryInSplitL1(int meshDimRelToMax, MEDCouplingPointSet *m);
  private:
    static const char SPE_FAM_STR_EXTRUDED_MESH[];
  private:
    std::vector< MCAuto<MEDFileUMeshSplitL1> > _ms;
    MCAuto<DataArrayDouble> _coords;
    MCAuto<DataArrayInt> _fam_coords;
    MCAuto<DataArrayInt> _num_coords;
    MCAuto<DataArrayInt> _global_num_coords;
    MCAuto<DataArrayAsciiChar> _name_coords;
    mutable MCAuto<DataArrayInt> _rev_num_coords;
    MCAuto<PartDefinition> _part_coords;
    std::vector< MCAuto<MEDFileEltStruct4Mesh> > _elt_str;
  };

  class MEDFileStructuredMesh : public MEDFileMesh
  {
    friend class MEDFileMesh;
  public:
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT int getMaxAbsFamilyIdInArrays() const;
    MEDLOADER_EXPORT int getMaxFamilyIdInArrays() const;
    MEDLOADER_EXPORT int getMinFamilyIdInArrays() const;
    MEDLOADER_EXPORT bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    MEDLOADER_EXPORT void clearNonDiscrAttributes() const;
    MEDLOADER_EXPORT DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const;
    MEDLOADER_EXPORT const DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt);
    MEDLOADER_EXPORT void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr);
    MEDLOADER_EXPORT void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr);
    MEDLOADER_EXPORT void setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr);
    MEDLOADER_EXPORT void setGlobalNumFieldAtLevel(int meshDimRelToMaxExt, DataArrayInt *globalNumArr);
    MEDLOADER_EXPORT void addNodeGroup(const DataArrayInt *ids);
    MEDLOADER_EXPORT void addGroup(int meshDimRelToMaxExt, const DataArrayInt *ids);
    MEDLOADER_EXPORT const DataArrayInt *getNumberFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT const DataArrayInt *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT const DataArrayAsciiChar *getNameFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT MCAuto<DataArrayInt> getGlobalNumFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT std::vector<int> getNonEmptyLevels() const;
    MEDLOADER_EXPORT std::vector<int> getNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getFamArrNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getNumArrNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getNameArrNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT MEDCouplingMesh *getMeshAtLevel(int meshDimRelToMax, bool renum=false) const;
    MEDLOADER_EXPORT std::vector<int> getFamsNonEmptyLevels(const std::vector<std::string>& fams) const;
    MEDLOADER_EXPORT std::vector<int> getFamsNonEmptyLevelsExt(const std::vector<std::string>& fams) const;
    MEDLOADER_EXPORT int getSizeAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT int getNumberOfNodes() const;
    MEDLOADER_EXPORT int getNumberOfCellsAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT bool hasImplicitPart() const;
    MEDLOADER_EXPORT int buildImplicitPartIfAny(INTERP_KERNEL::NormalizedCellType gt) const;
    MEDLOADER_EXPORT void releaseImplicitPartIfAny() const;
    MEDLOADER_EXPORT MEDCoupling1SGTUMesh *getImplicitFaceMesh() const;
    MEDLOADER_EXPORT std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypesAtLevel(int meshDimRelToMax) const;
    MEDLOADER_EXPORT int getNumberOfCellsWithType(INTERP_KERNEL::NormalizedCellType ct) const;
    MEDLOADER_EXPORT void whichAreNodesFetched(const MEDFileField1TSStructItem& st, const MEDFileFieldGlobsReal *globs, std::vector<bool>& nodesFetched) const;
    MEDLOADER_EXPORT bool presenceOfStructureElements() const { return false; }
    MEDLOADER_EXPORT virtual const MEDCouplingStructuredMesh *getStructuredMesh() const = 0;
    // tools
    MEDLOADER_EXPORT bool unPolyze(std::vector<int>& oldCode, std::vector<int>& newCode, DataArrayInt *& o2nRenumCell);
  protected:
    ~MEDFileStructuredMesh() { }
    void changeFamilyIdArr(int oldId, int newId);
    std::list< MCAuto<DataArrayInt> > getAllNonNullFamilyIds() const;
    void deepCpyAttributes();
    void loadStrMeshFromFile(MEDFileStrMeshL2 *strm, med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void writeStructuredLL(med_idt fid, const std::string& maa) const;
    void buildImplicitPart() const;
    void buildMinusOneImplicitPartIfNeeded() const;
    static med_geometry_type GetGeoTypeFromMeshDim(int meshDim);
  private:
    static void LoadStrMeshDAFromFile(med_idt fid, int meshDim, int dt, int it, const std::string& mName, MEDFileMeshReadSelector *mrs,
                                      MCAuto<DataArrayInt>& famCells, MCAuto<DataArrayInt>& numCells, MCAuto<DataArrayAsciiChar>& namesCells);
  private:
    MCAuto<DataArrayInt> _fam_nodes;
    MCAuto<DataArrayInt> _num_nodes;
    MCAuto<DataArrayAsciiChar> _names_nodes;
    MCAuto<DataArrayInt> _fam_cells;
    MCAuto<DataArrayInt> _num_cells;
    MCAuto<DataArrayAsciiChar> _names_cells;
    MCAuto<DataArrayInt> _fam_faces;
    MCAuto<DataArrayInt> _num_faces;
    MCAuto<DataArrayAsciiChar> _names_faces;
    mutable MCAuto<DataArrayInt> _rev_num_nodes;
    mutable MCAuto<DataArrayInt> _rev_num_cells;
    mutable MCAuto<MEDCoupling1SGTUMesh> _faces_if_necessary;
  };

  class MEDFileCMesh : public MEDFileStructuredMesh
  {
    friend class MEDFileMesh;
  public:
    MEDLOADER_EXPORT static MEDFileCMesh *New();
    MEDLOADER_EXPORT static MEDFileCMesh *New(const std::string& fileName, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileCMesh *New(med_idt fid, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileCMesh *New(DataArrayByte *db) { return BuildFromMemoryChunk<MEDFileCMesh>(db); }
    MEDLOADER_EXPORT static MEDFileCMesh *New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileCMesh *New(med_idt fid, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT MEDFileMesh *createNewEmpty() const;
    MEDLOADER_EXPORT MEDFileCMesh *deepCopy() const;
    MEDLOADER_EXPORT MEDFileCMesh *shallowCpy() const;
    MEDLOADER_EXPORT bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    MEDLOADER_EXPORT int getMeshDimension() const;
    MEDLOADER_EXPORT int getSpaceDimension() const;
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT std::string advancedRepr() const;
    MEDLOADER_EXPORT void clearNonDiscrAttributes() const;
    MEDLOADER_EXPORT const MEDCouplingCMesh *getMesh() const;
    MEDLOADER_EXPORT void setMesh(MEDCouplingCMesh *m);
    MEDLOADER_EXPORT MEDFileMesh *cartesianize() const;
  private:
    ~MEDFileCMesh() { }
    const MEDCouplingStructuredMesh *getStructuredMesh() const;
    void writeMeshLL(med_idt fid) const;
    MEDFileCMesh();
    void synchronizeTinyInfoOnLeaves() const;
    MEDFileCMesh(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadLL(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
  private:
    MCAuto<MEDCouplingCMesh> _cmesh;
  };

  class MEDFileCurveLinearMesh : public MEDFileStructuredMesh
  {
    friend class MEDFileMesh;
  public:
    MEDLOADER_EXPORT static MEDFileCurveLinearMesh *New();
    MEDLOADER_EXPORT static MEDFileCurveLinearMesh *New(const std::string& fileName, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileCurveLinearMesh *New(med_idt fid, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileCurveLinearMesh *New(DataArrayByte *db) { return BuildFromMemoryChunk<MEDFileCurveLinearMesh>(db); }
    MEDLOADER_EXPORT static MEDFileCurveLinearMesh *New(const std::string& fileName, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileCurveLinearMesh *New(med_idt fid, const std::string& mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT MEDFileMesh *createNewEmpty() const;
    MEDLOADER_EXPORT MEDFileCurveLinearMesh *deepCopy() const;
    MEDLOADER_EXPORT MEDFileCurveLinearMesh *shallowCpy() const;
    MEDLOADER_EXPORT bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    MEDLOADER_EXPORT int getMeshDimension() const;
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT std::string advancedRepr() const;
    MEDLOADER_EXPORT void clearNonDiscrAttributes() const;
    MEDLOADER_EXPORT const MEDCouplingCurveLinearMesh *getMesh() const;
    MEDLOADER_EXPORT void setMesh(MEDCouplingCurveLinearMesh *m);
    MEDLOADER_EXPORT MEDFileMesh *cartesianize() const;
  private:
    ~MEDFileCurveLinearMesh() { }
    MEDFileCurveLinearMesh();
    MEDFileCurveLinearMesh(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    const MEDCouplingStructuredMesh *getStructuredMesh() const;
    void synchronizeTinyInfoOnLeaves() const;
    void writeMeshLL(med_idt fid) const;
    void loadLL(med_idt fid, const std::string& mName, int dt, int it, MEDFileMeshReadSelector *mrs);//to imp
  private:
    MCAuto<MEDCouplingCurveLinearMesh> _clmesh;
  };

  class MEDFileMeshMultiTS : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT static MEDFileMeshMultiTS *New();
    MEDLOADER_EXPORT static MEDFileMeshMultiTS *New(med_idt fid);
    MEDLOADER_EXPORT static MEDFileMeshMultiTS *New(const std::string& fileName);
    MEDLOADER_EXPORT static MEDFileMeshMultiTS *New(med_idt fid, const std::string& mName);
    MEDLOADER_EXPORT static MEDFileMeshMultiTS *New(const std::string& fileName, const std::string& mName);
    MEDLOADER_EXPORT MEDFileMeshMultiTS *deepCopy() const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT std::string getName() const;
    MEDLOADER_EXPORT void setName(const std::string& newMeshName);
    MEDLOADER_EXPORT bool changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT void cartesianizeMe();
    MEDLOADER_EXPORT MEDFileMesh *getOneTimeStep() const;
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    MEDLOADER_EXPORT void setOneTimeStep(MEDFileMesh *mesh1TimeStep);
    MEDLOADER_EXPORT MEDFileJoints *getJoints() const;
    MEDLOADER_EXPORT void setJoints(MEDFileJoints* joints);
    MEDLOADER_EXPORT bool presenceOfStructureElements() const;
    MEDLOADER_EXPORT void killStructureElements();
  private:
    ~MEDFileMeshMultiTS() { }
    void loadFromFile(med_idt fid, const std::string& mName);
    MEDFileMeshMultiTS();
    MEDFileMeshMultiTS(med_idt fid);
    MEDFileMeshMultiTS(med_idt fid, const std::string& mName);
  private:
    std::vector< MCAuto<MEDFileMesh> > _mesh_one_ts;
  };

  class MEDFileMeshesIterator;

  class MEDFileMeshes : public RefCountObject, public MEDFileWritableStandAlone
  {
  public:
    MEDLOADER_EXPORT static MEDFileMeshes *New();
    MEDLOADER_EXPORT static MEDFileMeshes *New(med_idt fid);
    MEDLOADER_EXPORT static MEDFileMeshes *New(DataArrayByte *db) { return BuildFromMemoryChunk<MEDFileMeshes>(db); }
    MEDLOADER_EXPORT static MEDFileMeshes *New(const std::string& fileName);
    MEDLOADER_EXPORT MEDFileMeshes *deepCopy() const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildrenWithNull() const;
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT void simpleReprWithoutHeader(std::ostream& oss) const;
    MEDLOADER_EXPORT void writeLL(med_idt fid) const;
    MEDLOADER_EXPORT int getNumberOfMeshes() const;
    MEDLOADER_EXPORT MEDFileMeshesIterator *iterator();
    MEDLOADER_EXPORT MEDFileMesh *getMeshAtPos(int i) const;
    MEDLOADER_EXPORT MEDFileMesh *getMeshWithName(const std::string& mname) const;
    MEDLOADER_EXPORT std::vector<std::string> getMeshesNames() const;
    MEDLOADER_EXPORT bool changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT void cartesianizeMe();
    //
    MEDLOADER_EXPORT void resize(int newSize);
    MEDLOADER_EXPORT void pushMesh(MEDFileMesh *mesh);
    MEDLOADER_EXPORT void setMeshAtPos(int i, MEDFileMesh *mesh);
    MEDLOADER_EXPORT void destroyMeshAtPos(int i);
    MEDLOADER_EXPORT bool presenceOfStructureElements() const;
    MEDLOADER_EXPORT void killStructureElements();
  private:
    ~MEDFileMeshes() { }
    void checkConsistencyLight() const;
    void loadFromFile(med_idt fid);
    MEDFileMeshes();
    MEDFileMeshes(med_idt fid);
  private:
    std::vector< MCAuto<MEDFileMeshMultiTS> > _meshes;
  };

  class MEDFileMeshesIterator
  {
  public:
    MEDLOADER_EXPORT MEDFileMeshesIterator(MEDFileMeshes *ms);
    MEDLOADER_EXPORT ~MEDFileMeshesIterator();
    MEDLOADER_EXPORT MEDFileMesh *nextt();
  private:
    MCAuto<MEDFileMeshes> _ms;
    int _iter_id;
    int _nb_iter;
  };
}

#endif
