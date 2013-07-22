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
  class MEDFileFieldGlobs;
  class MEDFileField1TSStructItem;
  
  class MEDLOADER_EXPORT MEDFileMesh : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileMesh *New(const char *fileName, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    static MEDFileMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    virtual MEDFileMesh *createNewEmpty() const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDFileMesh *deepCpy() const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDFileMesh *shallowCpy() const throw(INTERP_KERNEL::Exception) = 0;
    virtual bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    virtual void clearNonDiscrAttributes() const;
    void setName(const char *name) { _name=name; }
    bool changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    const char *getName() const { return _name.c_str(); }
    const char *getUnivName() const { return _univ_name.c_str(); }
    bool getUnivNameWrStatus() const { return _univ_wr_status; }
    void setUnivNameWrStatus(bool newStatus) { _univ_wr_status=newStatus; }
    void setDescription(const char *name) { _desc_name=name; }
    const char *getDescription() const { return _desc_name.c_str(); }
    void setOrder(int order) { _order=order; }
    int getOrder() const { return _order; }
    void setIteration(int it) { _iteration=it; }
    int getIteration() const { return _iteration; }
    void setTimeValue(double time) { _time=time; }
    void setTime(int dt, int it, double time) { _time=time; _iteration=dt; _order=it; }
    double getTime(int& dt, int& it) { dt=_iteration; it=_order; return _time; }
    double getTimeValue() const { return _time; }
    void setTimeUnit(const char *unit) { _dt_unit=unit; }
    const char *getTimeUnit() const { return _dt_unit.c_str(); }
    virtual int getNumberOfNodes() const throw(INTERP_KERNEL::Exception) = 0;
    virtual std::vector<int> getNonEmptyLevels() const = 0;
    virtual std::vector<int> getNonEmptyLevelsExt() const = 0;
    virtual std::vector<int> getFamArrNonEmptyLevelsExt() const = 0;
    virtual std::vector<int> getNumArrNonEmptyLevelsExt() const = 0;
    virtual std::vector<int> getNameArrNonEmptyLevelsExt() const = 0;
    virtual void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    virtual void write(med_idt fid) const throw(INTERP_KERNEL::Exception);
    virtual int getSizeAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDCouplingMesh *getGenMeshAtLevel(int meshDimRelToMax, bool renum=false) const throw(INTERP_KERNEL::Exception) = 0;
    virtual void whichAreNodesFetched(const MEDFileField1TSStructItem& st, const MEDFileFieldGlobs *globs, std::vector<bool>& nodesFetched) const throw(INTERP_KERNEL::Exception) = 0;
    //
    bool areFamsEqual(const MEDFileMesh *other, std::string& what) const;
    bool areGrpsEqual(const MEDFileMesh *other, std::string& what) const;
    bool existsGroup(const char *groupName) const;
    bool existsFamily(int famId) const;
    bool existsFamily(const char *familyName) const;
    void setFamilyId(const char *familyName, int id);
    void setFamilyIdUnique(const char *familyName, int id) throw(INTERP_KERNEL::Exception);
    virtual void addFamily(const char *familyName, int id) throw(INTERP_KERNEL::Exception);
    virtual void createGroupOnAll(int meshDimRelToMaxExt, const char *groupName) throw(INTERP_KERNEL::Exception);
    virtual bool keepFamIdsOnlyOnLevs(const std::vector<int>& famIds, const std::vector<int>& levs) throw(INTERP_KERNEL::Exception);
    void addFamilyOnGrp(const char *grpName, const char *famName) throw(INTERP_KERNEL::Exception);
    std::string findOrCreateAndGiveFamilyWithId(int id, bool& created) throw(INTERP_KERNEL::Exception);
    void setFamilyInfo(const std::map<std::string,int>& info);
    void setGroupInfo(const std::map<std::string, std::vector<std::string> >&info);
    void copyFamGrpMapsFrom(const MEDFileMesh& other);
    const std::map<std::string,int>& getFamilyInfo() const { return _families; }
    const std::map<std::string, std::vector<std::string> >& getGroupInfo() const { return _groups; }
    std::vector<std::string> getFamiliesOnGroup(const char *name) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getFamiliesOnGroups(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamiliesIdsOnGroup(const char *name) const throw(INTERP_KERNEL::Exception);
    void setFamiliesOnGroup(const char *name, const std::vector<std::string>& fams) throw(INTERP_KERNEL::Exception);
    void setFamiliesIdsOnGroup(const char *name, const std::vector<int>& famIds) throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getGroupsOnFamily(const char *name) const throw(INTERP_KERNEL::Exception);
    void setGroupsOnFamily(const char *famName, const std::vector<std::string>& grps) throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getGroupsNames() const;
    std::vector<std::string> getFamiliesNames() const;
    void assignFamilyNameWithGroupName() throw(INTERP_KERNEL::Exception);
    std::vector<std::string> removeEmptyGroups() throw(INTERP_KERNEL::Exception);
    void removeGroup(const char *name) throw(INTERP_KERNEL::Exception);
    void removeFamily(const char *name) throw(INTERP_KERNEL::Exception);
    std::vector<std::string> removeOrphanGroups() throw(INTERP_KERNEL::Exception);
    std::vector<std::string> removeOrphanFamilies() throw(INTERP_KERNEL::Exception);
    void changeGroupName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception);
    void changeFamilyName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception);
    void changeFamilyId(int oldId, int newId) throw(INTERP_KERNEL::Exception);
    void changeAllGroupsContainingFamily(const char *familyNameToChange, const std::vector<std::string>& newFamiliesNames) throw(INTERP_KERNEL::Exception);
    int getFamilyId(const char *name) const throw(INTERP_KERNEL::Exception);
    int getMaxAbsFamilyId() const throw(INTERP_KERNEL::Exception);
    int getMaxFamilyId() const throw(INTERP_KERNEL::Exception);
    int getMinFamilyId() const throw(INTERP_KERNEL::Exception);
    int getTheMaxAbsFamilyId() const throw(INTERP_KERNEL::Exception);
    int getTheMaxFamilyId() const throw(INTERP_KERNEL::Exception);
    int getTheMinFamilyId() const throw(INTERP_KERNEL::Exception);
    virtual int getMaxAbsFamilyIdInArrays() const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getMaxFamilyIdInArrays() const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getMinFamilyIdInArrays() const throw(INTERP_KERNEL::Exception) = 0;
    DataArrayInt *getAllFamiliesIdsReferenced() const throw(INTERP_KERNEL::Exception);
    DataArrayInt *computeAllFamilyIdsInUse() const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamiliesIds(const std::vector<std::string>& famNames) const throw(INTERP_KERNEL::Exception);
    std::string getFamilyNameGivenId(int id) const throw(INTERP_KERNEL::Exception);
    bool ensureDifferentFamIdsPerLevel() throw(INTERP_KERNEL::Exception);
    void normalizeFamIdsTrio() throw(INTERP_KERNEL::Exception);
    void normalizeFamIdsMEDFile() throw(INTERP_KERNEL::Exception);
    virtual int getMeshDimension() const throw(INTERP_KERNEL::Exception) = 0;
    virtual std::string simpleRepr() const;
    virtual std::string advancedRepr() const = 0;
    //
    virtual void setGroupsAtLevel(int meshDimRelToMaxExt, const std::vector<const DataArrayInt *>& grps, bool renum=false) throw(INTERP_KERNEL::Exception);
    virtual void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr) throw(INTERP_KERNEL::Exception) = 0;
    virtual void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr) throw(INTERP_KERNEL::Exception) = 0;
    virtual void setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr) throw(INTERP_KERNEL::Exception) = 0;
    virtual const DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) = 0;
    virtual const DataArrayInt *getNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) = 0;
    virtual const DataArrayInt *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) = 0;
    virtual const DataArrayAsciiChar *getNameFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayInt *getGroupsArr(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getGroupArr(int meshDimRelToMaxExt, const char *grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getFamilyArr(int meshDimRelToMaxExt, const char *fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeGroupArr(const char *grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeGroupsArr(const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeFamilyArr(const char *fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    // tools
    virtual bool unPolyze(std::vector<int>& oldCode, std::vector<int>& newCode, DataArrayInt *& o2nRenumCell) throw(INTERP_KERNEL::Exception) = 0;
  protected:
    MEDFileMesh();
    //! protected because no way in MED file API to specify this name
    void setUnivName(const char *name) { _univ_name=name; }
    void addFamilyOnAllGroupsHaving(const char *famName, const char *otherFamName) throw(INTERP_KERNEL::Exception);
    virtual void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception) = 0;
    void dealWithTinyInfo(const MEDCouplingMesh *m) throw(INTERP_KERNEL::Exception);
    virtual void synchronizeTinyInfoOnLeaves() const = 0;
    void getFamilyRepr(std::ostream& oss) const;
    virtual void appendFamilyEntries(const DataArrayInt *famIds, const std::vector< std::vector<int> >& fidsOfGrps, const std::vector<std::string>& grpNames);
    virtual void changeFamilyIdArr(int oldId, int newId) throw(INTERP_KERNEL::Exception) = 0;
    static void TranslateFamilyIds(int offset, DataArrayInt *famArr, std::vector< std::vector<int> >& famIdsPerGrp);
    static void ChangeAllGroupsContainingFamily(std::map<std::string, std::vector<std::string> >& groups, const char *familyNameToChange, const std::vector<std::string>& newFamiliesNames) throw(INTERP_KERNEL::Exception);
    static std::string FindOrCreateAndGiveFamilyWithId(std::map<std::string,int>& families, int id, bool& created) throw(INTERP_KERNEL::Exception);
    static std::string CreateNameNotIn(const std::string& nameTry, const std::vector<std::string>& namesToAvoid) throw(INTERP_KERNEL::Exception);
    static int PutInThirdComponentOfCodeOffset(std::vector<int>& code, int strt) throw(INTERP_KERNEL::Exception);
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
    static const char DFT_FAM_NAME[];
  };

  class MEDLOADER_EXPORT MEDFileUMesh : public MEDFileMesh
  {
    friend class MEDFileMesh;
  public:
    static MEDFileUMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    static MEDFileUMesh *New(const char *fileName, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    static MEDFileUMesh *New();
    std::size_t getHeapMemorySize() const;
    MEDFileMesh *createNewEmpty() const throw(INTERP_KERNEL::Exception);
    MEDFileMesh *deepCpy() const throw(INTERP_KERNEL::Exception);
    MEDFileMesh *shallowCpy() const throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    void clearNonDiscrAttributes() const;
    ~MEDFileUMesh();
    //
    int getMaxAbsFamilyIdInArrays() const throw(INTERP_KERNEL::Exception);
    int getMaxFamilyIdInArrays() const throw(INTERP_KERNEL::Exception);
    int getMinFamilyIdInArrays() const throw(INTERP_KERNEL::Exception);
    int getMeshDimension() const throw(INTERP_KERNEL::Exception);
    int getSpaceDimension() const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    int getSizeAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    const DataArrayAsciiChar *getNameFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    int getNumberOfNodes() const throw(INTERP_KERNEL::Exception);
    void whichAreNodesFetched(const MEDFileField1TSStructItem& st, const MEDFileFieldGlobs *globs, std::vector<bool>& nodesFetched) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getNonEmptyLevels() const;
    std::vector<int> getNonEmptyLevelsExt() const;
    std::vector<int> getFamArrNonEmptyLevelsExt() const;
    std::vector<int> getNumArrNonEmptyLevelsExt() const;
    std::vector<int> getNameArrNonEmptyLevelsExt() const;
    std::vector<int> getGrpNonEmptyLevels(const char *grp) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getGrpNonEmptyLevelsExt(const char *grp) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamNonEmptyLevels(const char *fam) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamNonEmptyLevelsExt(const char *fam) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getGrpsNonEmptyLevels(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getGrpsNonEmptyLevelsExt(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamsNonEmptyLevels(const std::vector<std::string>& fams) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamsNonEmptyLevelsExt(const std::vector<std::string>& fams) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getGroupsOnSpecifiedLev(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *getCoords() const;
    MEDCouplingUMesh *getGroup(int meshDimRelToMaxExt, const char *grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getGroups(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getFamily(int meshDimRelToMaxExt, const char *fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getFamilies(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getMeshAtLevel(int meshDimRelToMaxExt, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingMesh *getGenMeshAtLevel(int meshDimRelToMax, bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getLevel0Mesh(bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getLevelM1Mesh(bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getLevelM2Mesh(bool renum=false) const throw(INTERP_KERNEL::Exception);
    MEDCouplingUMesh *getLevelM3Mesh(bool renum=false) const throw(INTERP_KERNEL::Exception);
    std::vector<MEDCoupling1GTUMesh *> getDirectUndergroundSingleGeoTypeMeshes(int meshDimRelToMax) const throw(INTERP_KERNEL::Exception);
    MEDCoupling1GTUMesh *getDirectUndergroundSingleGeoTypeMesh(INTERP_KERNEL::NormalizedCellType gt) const throw(INTERP_KERNEL::Exception);
    //
    void setFamilyNameAttachedOnId(int id, const std::string& newFamName) throw(INTERP_KERNEL::Exception);
    void setCoords(DataArrayDouble *coords) throw(INTERP_KERNEL::Exception);
    void eraseGroupsAtLevel(int meshDimRelToMaxExt) throw(INTERP_KERNEL::Exception);
    void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr) throw(INTERP_KERNEL::Exception);
    void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr) throw(INTERP_KERNEL::Exception);
    void setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr) throw(INTERP_KERNEL::Exception);
    void addNodeGroup(const DataArrayInt *ids) throw(INTERP_KERNEL::Exception);
    void addGroup(int meshDimRelToMaxExt, const DataArrayInt *ids) throw(INTERP_KERNEL::Exception);
    void removeMeshAtLevel(int meshDimRelToMax) throw(INTERP_KERNEL::Exception);
    void setMeshAtLevel(int meshDimRelToMax, MEDCouplingUMesh *m, bool newOrOld=false) throw(INTERP_KERNEL::Exception);
    void setMeshes(const std::vector<const MEDCouplingUMesh *>& ms, bool renum=false) throw(INTERP_KERNEL::Exception);
    void setGroupsFromScratch(int meshDimRelToMax, const std::vector<const MEDCouplingUMesh *>& ms, bool renum=false) throw(INTERP_KERNEL::Exception);
    void setGroupsOnSetMesh(int meshDimRelToMax, const std::vector<const MEDCouplingUMesh *>& ms, bool renum=false) throw(INTERP_KERNEL::Exception);
    void optimizeFamilies() throw(INTERP_KERNEL::Exception);
    // tools
    void duplicateNodesOnM1Group(const char *grpNameM1, DataArrayInt *&nodesDuplicated, DataArrayInt *&cellsModified, DataArrayInt *&cellsNotModified) throw(INTERP_KERNEL::Exception);
    bool unPolyze(std::vector<int>& oldCode, std::vector<int>& newCode, DataArrayInt *& o2nRenumCell) throw(INTERP_KERNEL::Exception);
    DataArrayInt *zipCoords() throw(INTERP_KERNEL::Exception);
  private:
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    MEDFileUMesh();
    MEDFileUMesh(med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs) throw(INTERP_KERNEL::Exception);
    void loadUMeshFromFile(med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs) throw(INTERP_KERNEL::Exception);
    const MEDFileUMeshSplitL1 *getMeshAtLevSafe(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    MEDFileUMeshSplitL1 *getMeshAtLevSafe(int meshDimRelToMaxExt) throw(INTERP_KERNEL::Exception);
    void checkMeshDimCoherency(int meshDim, int meshDimRelToMax) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *checkMultiMesh(const std::vector<const MEDCouplingUMesh *>& ms) const throw(INTERP_KERNEL::Exception);
    void computeRevNum() const;
    void synchronizeTinyInfoOnLeaves() const;
    void changeFamilyIdArr(int oldId, int newId) throw(INTERP_KERNEL::Exception);
    std::list< MEDCouplingAutoRefCountObjectPtr<DataArrayInt> > getAllNonNullFamilyIds() const;
    void addGroupUnderground(bool isNodeGroup, const DataArrayInt *ids, DataArrayInt *famArr) throw(INTERP_KERNEL::Exception);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> > _ms;
    MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> _coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam_coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num_coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayAsciiChar> _name_coords;
    mutable MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _rev_num_coords;
  };

  class MEDLOADER_EXPORT MEDFileStructuredMesh : public MEDFileMesh
  {
    friend class MEDFileMesh;
  public:
    std::size_t getHeapMemorySize() const;
    int getMaxAbsFamilyIdInArrays() const throw(INTERP_KERNEL::Exception);
    int getMaxFamilyIdInArrays() const throw(INTERP_KERNEL::Exception);
    int getMinFamilyIdInArrays() const throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    void clearNonDiscrAttributes() const;
    DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr) throw(INTERP_KERNEL::Exception);
    void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr) throw(INTERP_KERNEL::Exception);
    void setNameFieldAtLevel(int meshDimRelToMaxExt, DataArrayAsciiChar *nameArr) throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    const DataArrayAsciiChar *getNameFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getNonEmptyLevels() const;
    std::vector<int> getNonEmptyLevelsExt() const;
    std::vector<int> getFamArrNonEmptyLevelsExt() const;
    std::vector<int> getNumArrNonEmptyLevelsExt() const;
    std::vector<int> getNameArrNonEmptyLevelsExt() const;
    MEDCouplingMesh *getGenMeshAtLevel(int meshDimRelToMax, bool renum=false) const throw(INTERP_KERNEL::Exception);
    int getSizeAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    int getNumberOfNodes() const throw(INTERP_KERNEL::Exception);
    void whichAreNodesFetched(const MEDFileField1TSStructItem& st, const MEDFileFieldGlobs *globs, std::vector<bool>& nodesFetched) const throw(INTERP_KERNEL::Exception);
    // tools
    bool unPolyze(std::vector<int>& oldCode, std::vector<int>& newCode, DataArrayInt *& o2nRenumCell) throw(INTERP_KERNEL::Exception);
  protected:
    void changeFamilyIdArr(int oldId, int newId) throw(INTERP_KERNEL::Exception);
    void deepCpyAttributes() throw(INTERP_KERNEL::Exception);
    void loadStrMeshFromFile(MEDFileStrMeshL2 *strm, med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs) throw(INTERP_KERNEL::Exception);
    void writeStructuredLL(med_idt fid, const char *maa) const throw(INTERP_KERNEL::Exception);
    virtual const MEDCouplingStructuredMesh *getStructuredMesh() const = 0;
    static med_geometry_type GetGeoTypeFromMeshDim(int meshDim) throw(INTERP_KERNEL::Exception);
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

  class MEDLOADER_EXPORT MEDFileCMesh : public MEDFileStructuredMesh
  {
    friend class MEDFileMesh;
  public:
    static MEDFileCMesh *New();
    static MEDFileCMesh *New(const char *fileName, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    static MEDFileCMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    MEDFileMesh *createNewEmpty() const throw(INTERP_KERNEL::Exception);
    MEDFileMesh *deepCpy() const throw(INTERP_KERNEL::Exception);
    MEDFileMesh *shallowCpy() const throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    int getMeshDimension() const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    void clearNonDiscrAttributes() const;
    const MEDCouplingCMesh *getMesh() const;
    void setMesh(MEDCouplingCMesh *m) throw(INTERP_KERNEL::Exception);
  private:
    const MEDCouplingStructuredMesh *getStructuredMesh() const;
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    MEDFileCMesh();
    void synchronizeTinyInfoOnLeaves() const;
    MEDFileCMesh(med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs) throw(INTERP_KERNEL::Exception);
    void loadCMeshFromFile(med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs) throw(INTERP_KERNEL::Exception);
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingCMesh> _cmesh;
  };

  class MEDLOADER_EXPORT MEDFileCurveLinearMesh : public MEDFileStructuredMesh
  {
    friend class MEDFileMesh;
  public:
    static MEDFileCurveLinearMesh *New();
    static MEDFileCurveLinearMesh *New(const char *fileName, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    static MEDFileCurveLinearMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1, MEDFileMeshReadSelector *mrs=0) throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    MEDFileMesh *createNewEmpty() const throw(INTERP_KERNEL::Exception);
    MEDFileMesh *deepCpy() const throw(INTERP_KERNEL::Exception);
    MEDFileMesh *shallowCpy() const throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    int getMeshDimension() const throw(INTERP_KERNEL::Exception);
    std::string simpleRepr() const;
    std::string advancedRepr() const;
    void clearNonDiscrAttributes() const;
    const MEDCouplingCurveLinearMesh *getMesh() const;
    void setMesh(MEDCouplingCurveLinearMesh *m) throw(INTERP_KERNEL::Exception);
  private:
    MEDFileCurveLinearMesh();
    MEDFileCurveLinearMesh(med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs) throw(INTERP_KERNEL::Exception);
    const MEDCouplingStructuredMesh *getStructuredMesh() const;
    void synchronizeTinyInfoOnLeaves() const;
    void writeLL(med_idt fid) const throw(INTERP_KERNEL::Exception);
    void loadCLMeshFromFile(med_idt fid, const char *mName, int dt, int it, MEDFileMeshReadSelector *mrs) throw(INTERP_KERNEL::Exception);//to imp
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingCurveLinearMesh> _clmesh;
  };

  class MEDLOADER_EXPORT MEDFileMeshMultiTS : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileMeshMultiTS *New();
    static MEDFileMeshMultiTS *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileMeshMultiTS *New(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception);
    MEDFileMeshMultiTS *deepCpy() const throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    const char *getName() const throw(INTERP_KERNEL::Exception);
    void setName(const char *newMeshName) throw(INTERP_KERNEL::Exception);
    bool changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    MEDFileMesh *getOneTimeStep() const throw(INTERP_KERNEL::Exception);
    void write(med_idt fid) const throw(INTERP_KERNEL::Exception);
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void setOneTimeStep(MEDFileMesh *mesh1TimeStep) throw(INTERP_KERNEL::Exception);
  private:
    void loadFromFile(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception);
    MEDFileMeshMultiTS();
    MEDFileMeshMultiTS(const char *fileName) throw(INTERP_KERNEL::Exception);
    MEDFileMeshMultiTS(const char *fileName, const char *mName) throw(INTERP_KERNEL::Exception);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMesh> > _mesh_one_ts;
  };

  class MEDFileMeshesIterator;

  class MEDLOADER_EXPORT MEDFileMeshes : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileMeshes *New();
    static MEDFileMeshes *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    MEDFileMeshes *deepCpy() const throw(INTERP_KERNEL::Exception);
    std::size_t getHeapMemorySize() const;
    std::string simpleRepr() const;
    void simpleReprWithoutHeader(std::ostream& oss) const;
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    void write(med_idt fid) const throw(INTERP_KERNEL::Exception);
    int getNumberOfMeshes() const throw(INTERP_KERNEL::Exception);
    MEDFileMeshesIterator *iterator() throw(INTERP_KERNEL::Exception);
    MEDFileMesh *getMeshAtPos(int i) const throw(INTERP_KERNEL::Exception);
    MEDFileMesh *getMeshWithName(const char *mname) const throw(INTERP_KERNEL::Exception);
    std::vector<std::string> getMeshesNames() const throw(INTERP_KERNEL::Exception);
    bool changeNames(const std::vector< std::pair<std::string,std::string> >& modifTab) throw(INTERP_KERNEL::Exception);
    //
    void resize(int newSize) throw(INTERP_KERNEL::Exception);
    void pushMesh(MEDFileMesh *mesh) throw(INTERP_KERNEL::Exception);
    void setMeshAtPos(int i, MEDFileMesh *mesh) throw(INTERP_KERNEL::Exception);
    void destroyMeshAtPos(int i) throw(INTERP_KERNEL::Exception);
  private:
    void checkCoherency() const throw(INTERP_KERNEL::Exception);
    void loadFromFile(const char *fileName) throw(INTERP_KERNEL::Exception);
    MEDFileMeshes();
    MEDFileMeshes(const char *fileName) throw(INTERP_KERNEL::Exception);
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileMeshMultiTS> > _meshes;
  };

  class MEDCOUPLING_EXPORT MEDFileMeshesIterator
  {
  public:
    MEDFileMeshesIterator(MEDFileMeshes *ms);
    ~MEDFileMeshesIterator();
    MEDFileMesh *nextt();
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDFileMeshes> _ms;
     int _iter_id;
     int _nb_iter;
  };
}

#endif
