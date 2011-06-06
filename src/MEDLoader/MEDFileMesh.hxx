// Copyright (C) 2007-2011  CEA/DEN, EDF R&D
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

#ifndef __MEDFILEMESH_HXX__
#define __MEDFILEMESH_HXX__

#include "MEDFileMeshLL.hxx"
#include "MEDFileUtilities.hxx"

#include <map>

namespace ParaMEDMEM
{
  class MEDFileMesh : public RefCountObject, public MEDFileWritable
  {
  public:
    static MEDFileMesh *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1) throw(INTERP_KERNEL::Exception);
    virtual bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    virtual void clearNonDiscrAttributes() const;
    void setName(const char *name) { _name=name; }
    const char *getName() const { return _name.c_str(); }
    void setUnivName(const char *name) { _univ_name=name; }
    const char *getUnivName() const { return _univ_name.c_str(); }
    void setDescription(const char *name) { _desc_name=name; }
    const char *getDescription() const { return _desc_name.c_str(); }
    void setOrder(int order) { _order=order; }
    int getOrder() const { return _order; }
    void setIteration(int it) { _iteration=it; }
    int getIteration() const { return _iteration; }
    void setTimeValue(double time) { _time=time; }
    void setTime(double time, int dt, int it) { _time=time; _iteration=dt; _order=it; }
    double getTime(int& dt, int& it) { dt=_iteration; it=_order; return _time; }
    double getTimeValue() const { return _time; }
    void setTimeUnit(const char *unit) { _dt_unit=unit; }
    const char *getTimeUnit() const { return _dt_unit.c_str(); }
    virtual void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception) = 0;
    virtual int getSizeAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) = 0;
    virtual MEDCouplingMesh *getGenMeshAtLevel(int meshDimRelToMax, bool renum=false) const throw(INTERP_KERNEL::Exception) = 0;
    //
    bool areFamsEqual(const MEDFileMesh *other, std::string& what) const;
    bool areGrpsEqual(const MEDFileMesh *other, std::string& what) const;
    bool existsFamily(int famId) const;
    bool existsFamily(const char *familyName) const;
    void setFamilyId(const char *familyName, int id);
    void addFamily(const char *familyName, int id) throw(INTERP_KERNEL::Exception);
    void addGrpOnFamily(const char *grpName, const char *famName) throw(INTERP_KERNEL::Exception);
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
    void removeGroup(const char *name) throw(INTERP_KERNEL::Exception);
    void removeFamily(const char *name) throw(INTERP_KERNEL::Exception);
    void changeGroupName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception);
    void changeFamilyName(const char *oldName, const char *newName) throw(INTERP_KERNEL::Exception);
    int getFamilyId(const char *name) const throw(INTERP_KERNEL::Exception);
    int getMaxFamilyId() const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamiliesIds(const std::vector<std::string>& famNames) const throw(INTERP_KERNEL::Exception);
    std::string getFamilyNameGivenId(int id) const throw(INTERP_KERNEL::Exception);
    virtual int getMeshDimension() const throw(INTERP_KERNEL::Exception) = 0;
    //
    virtual void setGroupsAtLevel(int meshDimRelToMaxExt, const std::vector<const DataArrayInt *>& grps, bool renum=false) throw(INTERP_KERNEL::Exception);
    virtual void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr) throw(INTERP_KERNEL::Exception) = 0;
    virtual void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr) throw(INTERP_KERNEL::Exception) = 0;
    virtual const DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) = 0;
    virtual const DataArrayInt *getNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) = 0;
    virtual const DataArrayInt *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception) = 0;
    virtual DataArrayInt *getGroupsArr(int meshDimRelToMaxExt, const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getGroupArr(int meshDimRelToMaxExt, const char *grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getFamilyArr(int meshDimRelToMaxExt, const char *fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeGroupArr(const char *grp, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeGroupsArr(const std::vector<std::string>& grps, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeFamilyArr(const char *fam, bool renum=false) const throw(INTERP_KERNEL::Exception);
    virtual DataArrayInt *getNodeFamiliesArr(const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
  protected:
    MEDFileMesh();
    virtual void synchronizeTinyInfoOnLeaves() const = 0;
    virtual void appendFamilyEntries(const std::set<int>& famIds, const std::vector< std::vector<int> >& fidsOfGrps, const std::vector<std::string>& grpNames);
    static void TranslateFamilyIds(int offset, DataArrayInt *famArr, std::vector< std::vector<int> >& famIdsPerGrp);
  protected:
    int _order;
    int _iteration;
    double _time;
    std::string _dt_unit;
    std::string _name;
    std::string _univ_name;
    std::string _desc_name;
  protected:
    std::map<std::string, std::vector<std::string> > _groups;
    std::map<std::string,int> _families;
  public:
    static const char DFT_FAM_NAME[];
  };

  class MEDFileUMesh : public MEDFileMesh
  {
    friend class MEDFileMesh;
  public:
    static MEDFileUMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1) throw(INTERP_KERNEL::Exception);
    static MEDFileUMesh *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileUMesh *New();
    bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    void clearNonDiscrAttributes() const;
    ~MEDFileUMesh();
    //
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    int getMeshDimension() const throw(INTERP_KERNEL::Exception);
    int getSizeAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getNonEmptyLevels() const;
    std::vector<int> getNonEmptyLevelsExt() const;
    std::vector<int> getGrpNonEmptyLevels(const char *grp) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getGrpNonEmptyLevelsExt(const char *grp) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamNonEmptyLevels(const char *fam) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamNonEmptyLevelsExt(const char *fam) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getGrpsNonEmptyLevels(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getGrpsNonEmptyLevelsExt(const std::vector<std::string>& grps) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamsNonEmptyLevels(const std::vector<std::string>& fams) const throw(INTERP_KERNEL::Exception);
    std::vector<int> getFamsNonEmptyLevelsExt(const std::vector<std::string>& fams) const throw(INTERP_KERNEL::Exception);
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
    //
    void setFamilyNameAttachedOnId(int id, const std::string& newFamName) throw(INTERP_KERNEL::Exception);
    void setCoords(DataArrayDouble *coords) throw(INTERP_KERNEL::Exception);
    void eraseGroupsAtLevel(int meshDimRelToMaxExt) throw(INTERP_KERNEL::Exception);
    void setFamilyField(DataArrayInt *arr, const std::vector< std::vector< int > > &userfids, const std::vector<std::string>& grpNames, bool renum=false) throw(INTERP_KERNEL::Exception);
    void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr) throw(INTERP_KERNEL::Exception);
    void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr) throw(INTERP_KERNEL::Exception);
    void addNodeGroup(const std::string& name, const std::vector<int>& ids) throw(INTERP_KERNEL::Exception);
    void removeMeshAtLevel(int meshDimRelToMax) throw(INTERP_KERNEL::Exception);
    void setMeshAtLevel(int meshDimRelToMax, MEDCouplingUMesh *m, bool newOrOld=false) throw(INTERP_KERNEL::Exception);
    void setMeshAtLevelGen(int meshDimRelToMax, MEDCouplingUMesh *m, bool newOrOld) throw(INTERP_KERNEL::Exception);
    void setGroupsFromScratch(int meshDimRelToMax, const std::vector<const MEDCouplingUMesh *>& ms) throw(INTERP_KERNEL::Exception);
    void setGroupsOnSetMesh(int meshDimRelToMax, const std::vector<const MEDCouplingUMesh *>& ms, bool renum) throw(INTERP_KERNEL::Exception);
    void optimizeFamilies() throw(INTERP_KERNEL::Exception);
  private:
    MEDFileUMesh();
    MEDFileUMesh(med_idt fid, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception);
    void loadUMeshFromFile(med_idt fid, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception);
    const MEDFileUMeshSplitL1 *getMeshAtLevSafe(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    MEDFileUMeshSplitL1 *getMeshAtLevSafe(int meshDimRelToMaxExt) throw(INTERP_KERNEL::Exception);
    void checkMeshDimCoherency(int meshDim, int meshDimRelToMax) const throw(INTERP_KERNEL::Exception);
    DataArrayDouble *checkMultiMesh(const std::vector<const MEDCouplingUMesh *>& ms) const throw(INTERP_KERNEL::Exception);
    void computeRevNum() const;
    void synchronizeTinyInfoOnLeaves() const;
  private:
    std::vector< MEDCouplingAutoRefCountObjectPtr<MEDFileUMeshSplitL1> > _ms;
    MEDCouplingAutoRefCountObjectPtr<DataArrayDouble> _coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam_coords;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num_coords;
    mutable MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _rev_num_coords;
  };


  class MEDFileCMesh : public MEDFileMesh
  {
    friend class MEDFileMesh;
  public:
    static MEDFileCMesh *New();
    static MEDFileCMesh *New(const char *fileName) throw(INTERP_KERNEL::Exception);
    static MEDFileCMesh *New(const char *fileName, const char *mName, int dt=-1, int it=-1) throw(INTERP_KERNEL::Exception);
    bool isEqual(const MEDFileMesh *other, double eps, std::string& what) const;
    int getMeshDimension() const throw(INTERP_KERNEL::Exception);
    void clearNonDiscrAttributes() const;
    const MEDCouplingCMesh *getMesh() const;
    MEDCouplingMesh *getGenMeshAtLevel(int meshDimRelToMax, bool renum=false) const throw(INTERP_KERNEL::Exception);
    void setMesh(MEDCouplingCMesh *m);
    void write(const char *fileName, int mode) const throw(INTERP_KERNEL::Exception);
    int getSizeAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    DataArrayInt *getFamiliesArr(int meshDimRelToMaxExt, const std::vector<std::string>& fams, bool renum=false) const throw(INTERP_KERNEL::Exception);
    void setFamilyFieldArr(int meshDimRelToMaxExt, DataArrayInt *famArr) throw(INTERP_KERNEL::Exception);
    void setRenumFieldArr(int meshDimRelToMaxExt, DataArrayInt *renumArr) throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getFamilyFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
    const DataArrayInt *getRevNumberFieldAtLevel(int meshDimRelToMaxExt) const throw(INTERP_KERNEL::Exception);
  private:
    MEDFileCMesh();
    void synchronizeTinyInfoOnLeaves() const;
    MEDFileCMesh(med_idt fid, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception);
    void loadCMeshFromFile(med_idt fid, const char *mName, int dt, int it) throw(INTERP_KERNEL::Exception);
  private:
    MEDCouplingAutoRefCountObjectPtr<MEDCouplingCMesh> _cmesh;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam_nodes;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num_nodes;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _fam_cells;
    MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _num_cells;
    mutable MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _rev_num_nodes;
    mutable MEDCouplingAutoRefCountObjectPtr<DataArrayInt> _rev_num_cells;
  };
}

#endif
