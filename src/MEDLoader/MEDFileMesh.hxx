// Copyright (C) 2007-2013  CEA/DEN, EDF R&D
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License.
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
#include "MEDFileUtilities.hxx"
#include "MEDFileMeshReadSelector.hxx"

#include <map>
#include <list>

namespace ParaMEDMEM
{
  class MEDFileFieldGlobsReal;
  class MEDFileField1TSStructItem;
  
  class MEDFileMesh : public RefCountObject, public MEDFileWritable
  {
  public:
    MEDLOADER_EXPORT static MEDFileMesh *New(const char *fileName, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDLOADER_EXPORT virtual MEDFileMesh *createNewEmpty() const = 0;
    MEDLOADER_EXPORT virtual MEDFileMesh *deepCpy() const = 0;
    MEDLOADER_EXPORT virtual MEDFileMesh *shallowCpy() const = 0;
    MEDLOADER_EXPORT virtual bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    MEDLOADER_EXPORT virtual void clearNonDiscrAttributes() const;
    MEDLOADER_EXPORT void setName(const char *name) { _name=name; }
    MEDLOADER_EXPORT bool changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT std::string getName() const { return _name; }
    MEDLOADER_EXPORT const char *getUnivName() const { return _univ_name.c_str(); }
    MEDLOADER_EXPORT bool getUnivNameWrStatus() const { return _univ_wr_status; }
    MEDLOADER_EXPORT void setUnivNameWrStatus(bool newStatus) { _univ_wr_status=newStatus; }
    MEDLOADER_EXPORT void setDescription(const char *name) { _desc_name=name; }
    MEDLOADER_EXPORT std::string getDescription() const { return _desc_name; }
    MEDLOADER_EXPORT void setOrder(int order) { _order=order; }
    MEDLOADER_EXPORT int getOrder() const { return _order; }
    MEDLOADER_EXPORT void setIteration(int it) { _iteration=it; }
    MEDLOADER_EXPORT int getIteration() const { return _iteration; }
    MEDLOADER_EXPORT void setTimeValue(double time) { _time=time; }
    MEDLOADER_EXPORT void setTime(int dt, int it, double time) { _time=time; _iteration=dt; _order=it; }
    MEDLOADER_EXPORT double getTime(int& dt, int& it) { dt=_iteration; it=_order; return _time; }
    MEDLOADER_EXPORT double getTimeValue() const { return _time; }
    MEDLOADER_EXPORT void setTimeUnit(const char *unit) { _dt_unit=unit; }
    MEDLOADER_EXPORT const char *getTimeUnit() const { return _dt_unit.c_str(); }
    MEDLOADER_EXPORT virtual int getNumberOfNodes() const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getNonEmptyLevels() const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getNonEmptyLevelsExt() const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getFamArrNonEmptyLevelsExt() const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getNumArrNonEmptyLevelsExt() const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getNameArrNonEmptyLevelsExt() const = 0;
    MEDLOADER_EXPORT virtual void write(const char *fileName, int mode) const;
    MEDLOADER_EXPORT virtual void write(med_idt fid) const;
    MEDLOADER_EXPORT virtual int getSizeAtLevel(int meshDimRelToMaxExt) const = 0;
    MEDLOADER_EXPORT virtual MEDCouplingMesh *getGenMeshAtLevel(int meshDimRelToMax, bool renum=false) const = 0;
    MEDLOADER_EXPORT virtual std::vector<int> getDistributionOfTypes(int meshDimRelToMax) const;
    MEDLOADER_EXPORT virtual void whichAreNodesFetched(const MEDFileField1TSStructItem& st, const MEDFileFieldGlobsReal *globs, std::vector<bool>& nodesFetched) const = 0;
    //
    MEDLOADER_EXPORT bool areFamsEqual(const MEDFileMesh *other, std::string& what) const;
    MEDLOADER_EXPORT bool areGrpsEqual(const MEDFileMesh *other, std::string& what) const;
    MEDLOADER_EXPORT bool existsGroup(const char *groupName) const;
    MEDLOADER_EXPORT bool existsFamily(int famId) const;
    MEDLOADER_EXPORT bool existsFamily(const char *familyName) const;
    MEDLOADER_EXPORT void setFamilyId(const char *familyName, int id);
    MEDLOADER_EXPORT void setFamilyIdUnique(const char *familyName, int id);
    MEDLOADER_EXPORT virtual void addFamily(const char *familyName, int id);
    MEDLOADER_EXPORT virtual void createGroupOnAll(int meshDimRelToMaxExt, const char *groupName);
    MEDLOADER_EXPORT virtual bool keepFamIdsOnlyOnLevs(const std::vector<int>& famIds, const std::vector<int>& levs);
    MEDLOADER_EXPORT void addFamilyOnGrp(const char *grpName, const char *famName);
    MEDLOADER_EXPORT std::string findOrCreateAndGiveFamilyWithId(int id, bool& created);
    MEDLOADER_EXPORT void setFamilyInfo(const std::map<std::string,int>& info);
    MEDLOADER_EXPORT void setGroupInfo(const std::map<std::string, std::vector<std::string> >&info);
    MEDLOADER_EXPORT void copyFamGrpMapsFrom(const MEDFileMesh& other);
    MEDLOADER_EXPORT void clearGrpMap();
    MEDLOADER_EXPORT void clearFamMap();
    MEDLOADER_EXPORT void clearFamGrpMaps();
    MEDLOADER_EXPORT const std::map<std::string,int>& getFamilyInfo() const { return _families; }
    MEDLOADER_EXPORT const std::map<std::string, std::vector<std::string> >& getGroupInfo() const { return _groups; }
    MEDLOADER_EXPORT std::vector<std::string> getFamiliesOnGroup(const char *name) const;
    MEDLOADER_EXPORT std::vector<std::string> getFamiliesOnGroups(const std::vector<std::string>& grps) const;
    MEDLOADER_EXPORT std::vector<int> getFamiliesIdsOnGroup(const char *name) const;
    MEDLOADER_EXPORT void setFamiliesOnGroup(const char *name, const std::vector<std::string>& fams);
    MEDLOADER_EXPORT void setFamiliesIdsOnGroup(const char *name, const std::vector<int>& famIds);
    MEDLOADER_EXPORT std::vector<std::string> getGroupsOnFamily(const char *name) const;
    MEDLOADER_EXPORT void setGroupsOnFamily(const char *famName, const std::vector<std::string>& grps);
    MEDLOADER_EXPORT std::vector<std::string> getGroupsNames() const;
    MEDLOADER_EXPORT std::vector<std::string> getFamiliesNames() const;
    MEDLOADER_EXPORT void assignFamilyNameWithGroupName();
    MEDLOADER_EXPORT std::vector<std::string> removeEmptyGroups();
    MEDLOADER_EXPORT void removeGroup(const char *name);
    MEDLOADER_EXPORT void removeFamily(const char *name);
    MEDLOADER_EXPORT std::vector<std::string> removeOrphanGroups();
    MEDLOADER_EXPORT std::vector<std::string> removeOrphanFamilies();
    MEDLOADER_EXPORT void changeGroupName(const char *oldName, const char *newName);
    MEDLOADER_EXPORT void changeFamilyName(const char *oldName, const char *newName);
    MEDLOADER_EXPORT void changeFamilyId(int oldId, int newId);
    MEDLOADER_EXPORT void changeAllGroupsContainingFamily(const char *familyNameToChange, const std::vector<std::string>& newFamiliesNames);
    MEDLOADER_EXPORT int getFamilyId(const char *name) const;
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
    MEDLOADER_EXPORT virtual const DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const = 0;
    MEDLOADER_EXPORT virtual const DataArrayInt *getNumberFieldAtLevel(int meshDimRelToMaxExt) const = 0;
    MEDLOADER_EXPORT virtual const DataArrayInt *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const = 0;
    MEDLOADER_EXPORT virtual const DataArrayAsciiChar *getNameFieldAtLevel(int meshDimRelToMaxExt) const = 0;
    MEDLOADER_EXPORT virtual DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const = 0;
    MEDLOADER_EXPORT virtual DataArrayInt *getGroupsArr(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const;
    MEDLOADER_EXPORT virtual DataArrayInt *getGroupArr(int meshDimRelToMaxExt, const char *grp, bool renum=false) const;
    MEDLOADER_EXPORT virtual DataArrayInt *getFamilyArr(int meshDimRelToMaxExt, const char *fam, bool renum=false) const;
    MEDLOADER_EXPORT virtual DataArrayInt *getNodeGroupArr(const char *grp, bool renum=false) const;
    MEDLOADER_EXPORT virtual DataArrayInt *getNodeGroupsArr(const std::vector<std::string>& grps, bool renum=false) const;
    MEDLOADER_EXPORT virtual DataArrayInt *getNodeFamilyArr(const char *fam, bool renum=false) const;
    MEDLOADER_EXPORT virtual DataArrayInt *getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum=false) const;
    // tools
    MEDLOADER_EXPORT virtual bool unPolyze(std::vector<int>& oldCode, std::vector<int>& newCode, DataArrayInt *& o2nRenumCell) = 0;
  protected:
    MEDFileMesh();
    //! protected because no way in MED file API to specify this name
    void setUnivName(const char *name) { _univ_name=name; }
    void addFamilyOnAllGroupsHaving(const char *famName, const char *otherFamName);
    virtual void writeLL(med_idt fid) const = 0;
    void dealWithTinyInfo(const MEDCouplingMesh *m);
    virtual void synchronizeTinyInfoOnLeaves() const = 0;
    void getFamilyRepr(std::ostream& oss) const;
    virtual void appendFamilyEntries(const DataArrayInt *famIds, const std::vector< std::vector<int> >& fidsOfGrps, const std::vector<std::string>& grpNames);
    virtual void changeFamilyIdArr(int oldId, int newId) = 0;
    static void TranslateFamilyIds(int offset, DataArrayInt *famArr, std::vector< std::vector<int> >& famIdsPerGrp);
    static void ChangeAllGroupsContainingFamily(std::map<std::string, std::vector<std::string> >& groups, const char *familyNameToChange, const std::vector<std::string>& newFamiliesNames);
    static std::string FindOrCreateAndGiveFamilyWithId(std::map<std::string,int>& families, int id, bool& created);
    static std::string CreateNameNotIn(const std::string& nameTry, const std::vector<std::string>& namesToAvoid);
    static int PutInThirdComponentOfCodeOffset(std::vector<int>& code, int strt);
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
  protected:
    std::map<std::string, std::vector<std::string> > _groups;
    std::map<std::string,int> _families;
  public:
    MEDLOADER_EXPORT static const char DFT_FAM_NAME[];
  };

  class MEDFileUMesh : public MEDFileMesh
  {
    friend class MEDFileMesh;
  public:
    MEDLOADER_EXPORT static MEDFileUMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileUMesh *New(const char *fileName, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileUMesh *New();
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDLOADER_EXPORT MEDFileMesh *createNewEmpty() const;
    MEDLOADER_EXPORT MEDFileMesh *deepCpy() const;
    MEDLOADER_EXPORT MEDFileMesh *shallowCpy() const;
    MEDLOADER_EXPORT bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    MEDLOADER_EXPORT void clearNonDiscrAttributes() const;
    MEDLOADER_EXPORT ~MEDFileUMesh();
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
    MEDLOADER_EXPORT const DataArrayInt *getNumberFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT const DataArrayInt *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT const DataArrayAsciiChar *getNameFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT int getNumberOfNodes() const;
    MEDLOADER_EXPORT std::vector<INTERP_KERNEL::NormalizedCellType> getGeoTypesAtLevel(int meshDimRelToMax) const;
    MEDLOADER_EXPORT void whichAreNodesFetched(const MEDFileField1TSStructItem& st, const MEDFileFieldGlobsReal *globs, std::vector<bool>& nodesFetched) const;
    MEDLOADER_EXPORT std::vector<int> getNonEmptyLevels() const;
    MEDLOADER_EXPORT std::vector<int> getNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getFamArrNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getNumArrNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getNameArrNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getGrpNonEmptyLevels(const char *grp) const;
    MEDLOADER_EXPORT std::vector<int> getGrpNonEmptyLevelsExt(const char *grp) const;
    MEDLOADER_EXPORT std::vector<int> getFamNonEmptyLevels(const char *fam) const;
    MEDLOADER_EXPORT std::vector<int> getFamNonEmptyLevelsExt(const char *fam) const;
    MEDLOADER_EXPORT std::vector<int> getGrpsNonEmptyLevels(const std::vector<std::string>& grps) const;
    MEDLOADER_EXPORT std::vector<int> getGrpsNonEmptyLevelsExt(const std::vector<std::string>& grps) const;
    MEDLOADER_EXPORT std::vector<int> getFamsNonEmptyLevels(const std::vector<std::string>& fams) const;
    MEDLOADER_EXPORT std::vector<int> getFamsNonEmptyLevelsExt(const std::vector<std::string>& fams) const;
    MEDLOADER_EXPORT std::vector<std::string> getGroupsOnSpecifiedLev(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT DataArrayDouble *getCoords() const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getGroup(int meshDimRelToMaxExt, const char *grp, bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getGroups(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getFamily(int meshDimRelToMaxExt, const char *fam, bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getFamilies(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const;
    MEDLOADER_EXPORT DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getMeshAtLevel(int meshDimRelToMaxExt, bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingMesh *getGenMeshAtLevel(int meshDimRelToMax, bool renum=false) const;
    MEDLOADER_EXPORT std::vector<int> getDistributionOfTypes(int meshDimRelToMax) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getLevel0Mesh(bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getLevelM1Mesh(bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getLevelM2Mesh(bool renum=false) const;
    MEDLOADER_EXPORT MEDCouplingUMesh *getLevelM3Mesh(bool renum=false) const;
    MEDLOADER_EXPORT void forceComputationOfParts() const;
    MEDLOADER_EXPORT std::vector<MEDCoupling1GTUMesh *> getDirectUndergroundSingleGeoTypeMeshes(int meshDimRelToMax) const;
    MEDLOADER_EXPORT MEDCoupling1GTUMesh *getDirectUndergroundSingleGeoTypeMesh(INTERP_KERNEL::NormalizedCellType gt) const;
    MEDLOADER_EXPORT DataArrayInt *extractFamilyFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    MEDLOADER_EXPORT DataArrayInt *extractNumberFieldOnGeoType(INTERP_KERNEL::NormalizedCellType gt) const;
    //
    MEDLOADER_EXPORT void setFamilyNameAttachedOnId(int id, const std::string& newFamName);
    MEDLOADER_EXPORT void setCoords(DataArrayDouble *coords);
    MEDLOADER_EXPORT void eraseGroupsAtLevel(int meshDimRelToMaxExt);
    MEDLOADER_EXPORT void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr);
    MEDLOADER_EXPORT void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr);
    MEDLOADER_EXPORT void setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr);
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
    MEDLOADER_EXPORT void duplicateNodesOnM1Group(const char *grpNameM1, DataArrayInt *&nodesDuplicated, DataArrayInt *&cellsModified, DataArrayInt *&cellsNotModified);
    MEDLOADER_EXPORT bool unPolyze(std::vector<int>& oldCode, std::vector<int>& newCode, DataArrayInt *& o2nRenumCell);
    MEDLOADER_EXPORT DataArrayInt *zipCoords();
  private:
    void writeLL(med_idt fid) const;
    MEDFileUMesh();
    MEDFileUMesh(med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadUMeshFromFile(med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    const MEDFileUMeshSplitL1 *getMeshAtLevSafe(int meshDimRelToMaxExt) const;
    MEDFileUMeshSplitL1 *getMeshAtLevSafe(int meshDimRelToMaxExt);
    void checkMeshDimCoherency(int meshDim, int meshDimRelToMax) const;
    DataArrayDouble *checkMultiMesh(const std::vector<const MEDCouplingUMesh *>& ms) const;
    void computeRevNum() const;
    void synchronizeTinyInfoOnLeaves() const;
    void changeFamilyIdArr(int oldId, int newId);
    std::list< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > getAllNonNullFamilyIds() const;
    void addGroupUnderground(bool isNodeGroup, const DataArrayInt *ids, DataArrayInt *famArr);
    MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1>& checkAndGiveEntryInSplitL1(int meshDimRelToMax, MEDCouplingPointSet *m);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> > _ms;
    MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> _coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam_coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num_coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayAsciiChar> _name_coords;
    mutable MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _rev_num_coords;
  };

  class MEDFileStructuredMesh : public MEDFileMesh
  {
    friend class MEDFileMesh;
  public:
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDLOADER_EXPORT int getMaxAbsFamilyIdInArrays() const;
    MEDLOADER_EXPORT int getMaxFamilyIdInArrays() const;
    MEDLOADER_EXPORT int getMinFamilyIdInArrays() const;
    MEDLOADER_EXPORT bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    MEDLOADER_EXPORT void clearNonDiscrAttributes() const;
    MEDLOADER_EXPORT DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const;
    MEDLOADER_EXPORT const DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr);
    MEDLOADER_EXPORT void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr);
    MEDLOADER_EXPORT void setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr);
    MEDLOADER_EXPORT const DataArrayInt *getNumberFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT const DataArrayInt *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT const DataArrayAsciiChar *getNameFieldAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT std::vector<int> getNonEmptyLevels() const;
    MEDLOADER_EXPORT std::vector<int> getNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getFamArrNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getNumArrNonEmptyLevelsExt() const;
    MEDLOADER_EXPORT std::vector<int> getNameArrNonEmptyLevelsExt() const;
    MEDCouplingMesh *getGenMeshAtLevel(int meshDimRelToMax, bool renum=false) const;
    MEDLOADER_EXPORT int getSizeAtLevel(int meshDimRelToMaxExt) const;
    MEDLOADER_EXPORT int getNumberOfNodes() const;
    MEDLOADER_EXPORT void whichAreNodesFetched(const MEDFileField1TSStructItem& st, const MEDFileFieldGlobsReal *globs, std::vector<bool>& nodesFetched) const;
    // tools
    MEDLOADER_EXPORT bool unPolyze(std::vector<int>& oldCode, std::vector<int>& newCode, DataArrayInt *& o2nRenumCell);
  protected:
    void changeFamilyIdArr(int oldId, int newId);
    void deepCpyAttributes();
    void loadStrMeshFromFile(MEDFileStrMeshL2 *strm, med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void writeStructuredLL(med_idt fid, const char *maa) const;
    virtual const MEDCouplingStructuredMesh *getStructuredMesh() const = 0;
    static med_geometry_type GetGeoTypeFromMeshDim(int meshDim);
  private:
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam_nodes;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num_nodes;
    MEDCouplingAutoRefCountObjectPtr<DataArrayAsciiChar> _names_nodes;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam_cells;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num_cells;
    MEDCouplingAutoRefCountObjectPtr<DataArrayAsciiChar> _names_cells;
    mutable MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _rev_num_nodes;
    mutable MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _rev_num_cells;
  };

  class MEDFileCMesh : public MEDFileStructuredMesh
  {
    friend class MEDFileMesh;
  public:
    MEDLOADER_EXPORT static MEDFileCMesh *New();
    MEDLOADER_EXPORT static MEDFileCMesh *New(const char *fileName, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileCMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDLOADER_EXPORT MEDFileMesh *createNewEmpty() const;
    MEDLOADER_EXPORT MEDFileMesh *deepCpy() const;
    MEDLOADER_EXPORT MEDFileMesh *shallowCpy() const;
    MEDLOADER_EXPORT bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    MEDLOADER_EXPORT int getMeshDimension() const;
    MEDLOADER_EXPORT int getSpaceDimension() const;
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT std::string advancedRepr() const;
    MEDLOADER_EXPORT void clearNonDiscrAttributes() const;
    MEDLOADER_EXPORT const MEDCouplingCMesh *getMesh() const;
    MEDLOADER_EXPORT void setMesh(MEDCouplingCMesh *m);
  private:
    const MEDCouplingStructuredMesh *getStructuredMesh() const;
    void writeLL(med_idt fid) const;
    MEDFileCMesh();
    void synchronizeTinyInfoOnLeaves() const;
    MEDFileCMesh(med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    void loadCMeshFromFile(med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs);
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingCMesh> _cmesh;
  };

  class MEDFileCurveLinearMesh : public MEDFileStructuredMesh
  {
    friend class MEDFileMesh;
  public:
    MEDLOADER_EXPORT static MEDFileCurveLinearMesh *New();
    MEDLOADER_EXPORT static MEDFileCurveLinearMesh *New(const char *fileName, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT static MEDFileCurveLinearMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0);
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDLOADER_EXPORT MEDFileMesh *createNewEmpty() const;
    MEDLOADER_EXPORT MEDFileMesh *deepCpy() const;
    MEDLOADER_EXPORT MEDFileMesh *shallowCpy() const;
    MEDLOADER_EXPORT bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    MEDLOADER_EXPORT int getMeshDimension() const;
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT std::string advancedRepr() const;
    MEDLOADER_EXPORT void clearNonDiscrAttributes() const;
    MEDLOADER_EXPORT const MEDCouplingCurveLinearMesh *getMesh() const;
    MEDLOADER_EXPORT void setMesh(MEDCouplingCurveLinearMesh *m);
  private:
    MEDFileCurveLinearMesh();
    MEDFileCurveLinearMesh(med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs);
    const MEDCouplingStructuredMesh *getStructuredMesh() const;
    void synchronizeTinyInfoOnLeaves() const;
    void writeLL(med_idt fid) const;
    void loadCLMeshFromFile(med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs);//to imp
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingCurveLinearMesh> _clmesh;
  };

  class MEDFileMeshMultiTS : public RefCountObject, public MEDFileWritable
  {
  public:
    MEDLOADER_EXPORT static MEDFileMeshMultiTS *New();
    MEDLOADER_EXPORT static MEDFileMeshMultiTS *New(const char *fileName);
    MEDLOADER_EXPORT static MEDFileMeshMultiTS *New(const char *fileName, const char *mName);
    MEDLOADER_EXPORT MEDFileMeshMultiTS *deepCpy() const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDLOADER_EXPORT std::string getName() const;
    MEDLOADER_EXPORT void setName(const char *newMeshName);
    MEDLOADER_EXPORT bool changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    MEDLOADER_EXPORT MEDFileMesh *getOneTimeStep() const;
    MEDLOADER_EXPORT void write(med_idt fid) const;
    MEDLOADER_EXPORT void write(const char *fileName, int mode) const;
    MEDLOADER_EXPORT void setOneTimeStep(MEDFileMesh *mesh1TimeStep);
  private:
    void loadFromFile(const char *fileName, const char *mName);
    MEDFileMeshMultiTS();
    MEDFileMeshMultiTS(const char *fileName);
    MEDFileMeshMultiTS(const char *fileName, const char *mName);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> > _mesh_one_ts;
  };

  class MEDFileMeshesIterator;

  class MEDFileMeshes : public RefCountObject, public MEDFileWritable
  {
  public:
    MEDLOADER_EXPORT static MEDFileMeshes *New();
    MEDLOADER_EXPORT static MEDFileMeshes *New(const char *fileName);
    MEDLOADER_EXPORT MEDFileMeshes *deepCpy() const;
    MEDLOADER_EXPORT std::size_t getHeapMemorySizeWithoutChildren() const;
    MEDLOADER_EXPORT std::vector<const BigMemoryObject *> getDirectChildren() const;
    MEDLOADER_EXPORT std::string simpleRepr() const;
    MEDLOADER_EXPORT void simpleReprWithoutHeader(std::ostream& oss) const;
    MEDLOADER_EXPORT void write(const char *fileName, int mode) const;
    MEDLOADER_EXPORT void write(med_idt fid) const;
    MEDLOADER_EXPORT int getNumberOfMeshes() const;
    MEDLOADER_EXPORT MEDFileMeshesIterator *iterator();
    MEDLOADER_EXPORT MEDFileMesh *getMeshAtPos(int i) const;
    MEDLOADER_EXPORT MEDFileMesh *getMeshWithName(const char *mname) const;
    MEDLOADER_EXPORT std::vector<std::string> getMeshesNames() const;
    MEDLOADER_EXPORT bool changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab);
    //
    MEDLOADER_EXPORT void resize(int newSize);
    MEDLOADER_EXPORT void pushMesh(MEDFileMesh *mesh);
    MEDLOADER_EXPORT void setMeshAtPos(int i, MEDFileMesh *mesh);
    MEDLOADER_EXPORT void destroyMeshAtPos(int i);
  private:
    void checkCoherency() const;
    void loadFromFile(const char *fileName);
    MEDFileMeshes();
    MEDFileMeshes(const char *fileName);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMeshMultiTS> > _meshes;
  };

  class MEDFileMeshesIterator
  {
  public:
    MEDLOADER_EXPORT MEDFileMeshesIterator(MEDFileMeshes *ms);
    MEDLOADER_EXPORT ~MEDFileMeshesIterator();
    MEDLOADER_EXPORT MEDFileMesh *nextt();
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDFileMeshes> _ms;
     int _iter_id;
     int _nb_iter;
  };
}

#endif
